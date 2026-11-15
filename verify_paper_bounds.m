function verified = verify_paper_bounds()
%VERIFY_PAPER_BOUNDS Recheck Table 2 of the paper from the saved certificates.
%   verified = VERIFY_PAPER_BOUNDS() loads, for N_u = 800 and N_u = 1300, the
%   single report results/validation_<N_u>_*.mat produced by run_proof, and
%   checks in INTLAB that
%     1. the report is a successful certificate, and the Newton-Kantorovich
%        inequalities of Theorem 1 hold with its exact binary bounds;
%     2. the decimal values printed in Table 2 are upper bounds of the
%        certified Y, Z1, Z2, r0 (they are rounded upward);
%     3. both inequalities of Theorem 1 still hold with the printed values,
%            p = Y + Z1*r0 + Z2*r0^2/2 - r0 <= 0,   q = Z1 + Z2*r0 < 1.
%   It does not rerun the full bound computation (use run_proof for that).
%   Initialize INTLAB first.
environment = check_intlab();
root = fileparts(mfilename('fullpath'));
orders = [800, 1300];
% Values printed in Table 2: Y, Z1, Z2, r0.
% Strings preserve the intended decimal numbers during interval conversion.
rows = {{'6.33171e-12','0.344185','163.177','9.655e-12'}, ...
        {'1.87641e-12','0.234178','172.690','2.451e-12'}};
verified = struct([]);
for j = 1:numel(orders)
    % Exactly one saved certificate per order is expected in results/.
    found = dir(fullfile(root,'results',sprintf('validation_%d_*.mat',orders(j))));
    assert(numel(found)==1, ...
        'Expected exactly one report for N_u=%d in results/, found %d.',orders(j),numel(found));
    file = found(1).name;
    data = load(fullfile(root,'results',file),'report');
    report = data.report;
    c = report.certificate;
    assert(strcmp(report.status,'validated') && c.success, ...
        'Saved report does not contain a successful certificate.');
    % Independently recompute the inequalities from the saved bounds.
    savedR = intval(c.r0);
    assert(sup(c.Y+c.Z1*savedR+c.Z2*savedR^2/2-savedR)<=0);
    assert(sup(c.Z1+c.Z2*savedR)<1);
    % The printed values must not be smaller than the certified ones.
    Y = intval(rows{j}{1}); Z1 = intval(rows{j}{2});
    Z2 = intval(rows{j}{3}); r = intval(rows{j}{4});
    assert(inf(Y)>=sup(c.Y) && inf(Z1)>=sup(c.Z1) ...
        && inf(Z2)>=sup(c.Z2),'Printed bound rounds downward.');
    assert(inf(r)>=c.r0,'Printed radius rounds downward.');
    % Theorem 1 with the printed values.
    p = Y+Z1*r+Z2*r^2/2-r;
    q = Z1+Z2*r;
    assert(sup(p)<=0 && sup(q)<1,'Printed inequalities do not close.');
    % Auxiliary check (see run_proof): 1 + 2*pi + eps2*u(t) > 0 for the true
    % solution u. The lag changes by at most epsilon_2*r, with epsilon_2=2.
    shiftedLag = report.shiftedLagLower-2*r;
    assert(inf(shiftedLag)>0,'Positive shifted lag not verified.');
    verified(j).Nu = report.Nu;
    verified(j).Y = Y; verified(j).Z1 = Z1; verified(j).Z2 = Z2;
    verified(j).r0 = r; verified(j).p = p; verified(j).q = q;
    verified(j).shiftedLag = shiftedLag;
    verified(j).sourceReport = file;
    verified(j).environment = environment;
    fprintf('Nu=%d: printed bounds VERIFIED (p<=%.17g, q<=%.17g).\n', ...
        report.Nu,sup(p),sup(q));
end
end
