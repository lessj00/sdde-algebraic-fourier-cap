function [Kv,Kv_mid,Rv] = Cs_grid_SDDE(a,eps2,smax,M)
%CS_GRID_SDDE Sharp bounds tildeK_s(v) of Lemma 4 (midpoint enclosure).
%   [Kv,Kv_mid,Rv] = CS_GRID_SDDE(a,eps2,smax,M)
%   For v(t) = ubar(t - 1 - eps2*ubar(t)) and s = 1..smax, returns
%       Kv(s) = tildeK_s(v) >= K_s(v) = (1/2pi) ||v^{(s)}||_{L^1(-pi,pi)},
%   where, by Lemma 4 with h = 2*pi/M and cells I_l = [t_l, t_l+h],
%       tildeK_s(v) = (1/M) sum_l |v^{(s)}(t_l + h/2)|                 (Kv_mid)
%                   + (pi/(2M)) (1/M) sum_l sup_{I_l} |v^{(s+1)}|      (Rv).
%   The derivatives are evaluated with Faa di Bruno's formula (6) (Lemma 3),
%       v^{(s)} = sum_{j=1}^s ubar^{(j)}(theta) B_{s,j}(theta',...,theta^{(s-j+1)}),
%   with theta' = 1 - eps2*ubar' and theta^{(q)} = -eps2*ubar^{(q)} (q >= 2),
%   equation (5). As described after Lemma 4, the derivatives of theta are
%   obtained by interval FFTs and ubar^{(j)}(theta(.)) by interval
%   matrix-vector products. These constants are used for the Fourier tail of
%   v in (8) (Proposition 1(b)) and for T_K in Lemma 5.
%   Values on the grid t_l = l*h, l = 0..M-1, use interval inverse FFTs:
%   verifyfft(.,-1)*M returns sum_n alpha_n e^{i n t_l}.
%   Since Ubar_{-n} = conj(Ubar_n), all functions involved are real on R and
%   |.| is applied to the real part of the enclosures. Requires INTLAB.
N   = (length(a)-1)/2;
assert(M==2^nextpow2(M) && smax>=1 && smax<=30);
if 2*N >= M, error('need M > 2N_u'), end
n   = intval((-N:N)');
a   = intval(a);
ipi = intval('pi');
h   = 2*ipi/M;
% delta = [0,h]: t_l + delta encloses the whole cell I_l.
delta = infsup(0,1)*h;
tl  = intval((0:M-1)')*h;
% Zero padding of a coefficient vector (modes -N..N) to length M, in FFT order.
pad = @(c) ifftshift([intval(zeros(M/2-N,1)); c; intval(zeros(M/2-N-1,1))]);

Wm = cell(smax+1,1);   % ubar^{(j)}(t_l + h/2), j = 0..smax      (midpoints)
Uc = cell(smax+2,1);   % ubar^{(j)}([t_l,t_l+h]), j = 0..smax+1  (whole cells)
for j = 0:smax+1
    % Fourier coefficients of ubar^{(j)}: (i n)^j Ubar_n.
    cj = (1i*n).^j .* a;
    if j <= smax
        % Shift by h/2 to evaluate at the midpoints t_l + h/2.
        Wm{j+1} = real(verifyfft(pad(cj.*exp(1i*n*(h/2))),-1)*M);
    end
    % Multiplying by exp(i n delta) encloses the values on the cell I_l.
    Uc{j+1} = real(verifyfft(pad(cj.*exp(1i*n*delta)),-1)*M);
end
Cm = intval(zeros(2*N+1,smax+1));                    % column j: (i n)^j Ubar_n
for j = 1:smax+1, Cm(:,j) = (1i*n).^j .* a; end

mid_abs = zeros(M,smax); sup_rem = zeros(M,smax);
blk = 512;                                           % grid points per block
for b = 1:blk:M
    idx = (b:min(b+blk-1,M))';
    % ---- leading term of Lemma 4: |v^{(s)}| at the midpoints t_l + h/2
    % theta(t_l + h/2), and the Bell arguments theta^{(q)} from (5).
    th = tl(idx) + h/2 - 1 - eps2*Wm{1}(idx);
    x = cell(smax,1); x{1} = 1 - eps2*Wm{2}(idx);
    for i = 2:smax, x{i} = -eps2*Wm{i+1}(idx); end
    B  = bell_table(x,smax);
    % ubar^{(j)}(theta) = sum_n (i n)^j Ubar_n e^{i n theta}, j = 1..smax.
    UM = real(exp(1i*(th*n.'))*Cm(:,1:smax));
    for s = 1:smax
        % Faa di Bruno (6): v^{(s)} = sum_j ubar^{(j)}(theta) B_{s,j}.
        vs = intval(zeros(numel(idx),1));
        for j = 1:s, vs = vs + UM(:,j).*B{s+1,j+1}; end
        mid_abs(idx,s) = sup(abs(vs));
    end
    % ---- remainder term of Lemma 4: sup of |v^{(s+1)}| over each cell I_l
    th = tl(idx) + delta - 1 - eps2*Uc{1}(idx);
    x = cell(smax+1,1); x{1} = 1 - eps2*Uc{2}(idx);
    for i = 2:smax+1, x{i} = -eps2*Uc{i+1}(idx); end
    B  = bell_table(x,smax+1);
    UM = real(exp(1i*(th*n.'))*Cm);                  % j = 1..smax+1
    for s = 1:smax
        % v^{(s+1)} = sum_{j=1}^{s+1} ubar^{(j)}(theta) B_{s+1,j}.
        vs = intval(zeros(numel(idx),1));
        for j = 1:s+1, vs = vs + UM(:,j).*B{s+2,j+1}; end
        sup_rem(idx,s) = sup(abs(vs));
    end
end
% Assemble Lemma 4 with rigorous (interval) summation of the upper bounds.
Kv_mid = sum(intval(mid_abs),1).'/M;                 % (1/M) sum_l |v^{(s)}(t_l+h/2)|
Rv     = (ipi/(2*M))*sum(intval(sup_rem),1).'/M;     % (pi/(2M))(1/M) sum_l sup|v^{(s+1)}|
Kv     = intval(sup(Kv_mid + Rv));
end
