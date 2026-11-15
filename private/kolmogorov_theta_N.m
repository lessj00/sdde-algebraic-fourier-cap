function ThN = kolmogorov_theta_N(N,opts)
%KOLMOGOROV_THETA_N   Verified enclosure of a single Kolmogorov constant.
%
%   ThN = KOLMOGOROV_THETA_N(N) returns an INTLAB interval containing
%
%       theta_N    = (1/pi) \int_0^pi |kappa_N(t)| dt ,
%       kappa_N(t) = pi - t - 2 \sum_{k=1}^N sin(kt)/k ,
%
%   the optimal constant in  ||pi^{>N} h||_{C^0} <= theta_N ||h'||_{C^0}.
%   This Kolmogorov constant [Kolmogorov 1935] is used in Sect. 2 and in
%   Lemma 5 of the paper (bounds Z_{1,1}, Z_{1,2}, Y_res, and the norm of
%   A*Ddag); its computation follows [4, Sect. 2.6.3], where [4] is the
%   companion preprint (van den Berg, Breden, Cadiot, Church, Lessard).
%
%   ThN = KOLMOGOROV_THETA_N(N,opts) passes options; see below.
%
%   The floating-point Newton iteration only proposes root centers.
%   All inclusion tests and the final integral enclosure use INTLAB.
%   Source: [4, Sect. 2.6.3], proposition "computation Kolmogorov constant".
%   No continuation in N is used in this proof driver.
%
%   Method ([4, Sect. 2.6.3]):
%     (1) kappa_N is monotone on [t_{j-1},t_j], t_j = 2*pi*j/(2N+1), with a
%         unique zero tau_j in (t_{j-1},t_j); Newton gives taubar_j.
%     (2) Newton-Kantorovich test: with
%           Y_j >= |kappa_N(taubar_j)/kappa_N'(taubar_j)| ,
%           Z_j >= |1 - kappa_N'(taubar_j+[-eps,eps])/kappa_N'(taubar_j)| ,
%         if Z_j < 1 and r_j = Y_j/(1-Z_j) <= eps then
%         tau_j \in [taubar_j-r_j, taubar_j+r_j].  kappa_N' = -D_N is
%         available in closed form, so this costs O(1) per zero.
%     (3) theta_N is assembled from
%         theta_N = (2/pi) sum_j (-1)^{j-1} Kcal_N(tau_j)
%                   - (2/pi) sum_k 1/k^2
%                   + (-1)^N ( pi/2 + (2/pi) sum_k (-1)^k/k^2 ).
%
%   Options (struct, all optional)
%     opts.eps      radius of the a priori box in the Z-bound      (1e-9)
%     opts.tol      Newton stopping tolerance                      (1e-11)
%     opts.maxit    maximum number of Newton steps                 (8)
%     opts.blkelem  block size (matrix entries) of the trig sums   (2e6)
%
%   Cost is O(N^2) in time and O(opts.blkelem) in memory.
%
%   Example
%     Th = kolmogorov_theta_N(2553);
%     fprintf('theta_2553 in [%.16e, %.16e]\n',inf(Th),sup(Th));
%
%   Requires INTLAB.

if nargin < 2 || isempty(opts), opts = struct; end
def = struct('eps',1e-9,'tol',1e-11,'maxit',8,'blkelem',2e6);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts,f{i}), opts.(f{i}) = def.(f{i}); end
end
if ~(isscalar(N) && N == round(N) && N >= 1)
    error('kolmogorov_theta_N:N','N must be a positive integer');
end

Pi = 4*atan(intval(1));             % verified enclosure of pi

% ----------------------------------------------------- grid and constants
j    = (1:N)';
kk   = N:-1:1;   jr = kk(:);        % DESCENDING k: the partial sums of the
% series are then O(1/k) instead of O(1), so the rounding accumulated in the
% length-N dot products is O(eps*log N) instead of O(eps*N).  This keeps
% rad(theta_N) growing like N rather than like N^2.
h    = 2*pi/(2*N+1);
m    = (2*j-1)*pi/(2*N+1);          % midpoints of (t_{j-1},t_j)
tlo  = (j-1)*h;   thi = j*h;        % t_{j-1} , t_j
blk  = max(1,floor(opts.blkelem/N));
q    = 51 - max(0,ceil(log2(N)));   % see QUANT: makes taubar*k exact

% (1) --------------------------------------------- starting point + Newton
tb = m;
tb = quant(tb,q);

w1 = 1./jr;   cold = inf;
for it = 1:opts.maxit
    kap  = pi - tb - 2*tsum(tb,kk,w1,blk,@sin);
    dkap = -sin(((2*N+1)/2)*tb)./sin(tb/2);
    corr = kap./dkap;
    tb   = quant(min(max(tb-corr,tlo+0.05*h),thi-0.05*h),q);
    cnew = max(abs(corr));
    if cnew <= opts.tol || cnew > 0.5*cold, break, end   % converged / stalled
    cold = cnew;
end

% (2) ------------------------------- verified enclosure of the zeros tau_j
% QUANT guarantees that every product taubar_j*k, k = 1..N, is exact in
% binary floating point; check the hypothesis rather than assume it.
if sup(intval(max(tb))*N*intval(pow2(1,q))) > 2^53
    error('kolmogorov_theta_N:exactness', ...
          'products taubar*k are not exact at N = %d',N);
end
W1 = intval(1)./intval(jr);         % 1/k   (descending order)
assert(N<2^26,'Integer squares must be exactly representable.');
W2 = intval(1)./intval(jr.^2);      % 1/k^2 (jr.^2 exact in double)

kapI = Pi - intval(tb) - 2*tsum(tb,kk,W1,blk,@sin);     % kappa_N(taubar)
dkI  = -Dn(N,intval(tb));                               % kappa_N'(taubar)
Ieps = midrad(tb,opts.eps*ones(N,1));
dkE  = -Dn(N,Ieps);                                     % kappa_N'([taubar+-eps])

if any(mig(dkI) == 0)
    error('kolmogorov_theta_N:derivative', ...
          'kappa_N''(taubar) not bounded away from 0 at N = %d',N);
end
Y = mag(kapI./dkI);                      % rigorous upper bound of |kappa_N/kappa_N'|
Z = mag(1 - dkE./dkI);
if ~all(Z < 1)
    error('kolmogorov_theta_N:Zbound','Z >= 1 at N = %d',N);
end
r = sup(intval(Y)./(1-intval(Z)));       % rigorous upper bound of Y/(1-Z)
% Identify the roots using exact-pi cell endpoints, not rounded double
% grid boundaries from the nonrigorous Newton proposal stage.
tloI = 2*Pi*intval(j-1)/(2*N+1);
thiI = 2*Pi*intval(j)/(2*N+1);
if ~( all(r <= opts.eps) && all(inf(intval(tb)-intval(r)) > sup(tloI)) ...
        && all(sup(intval(tb)+intval(r)) < inf(thiI)) )
    error('kolmogorov_theta_N:verification', ...
          'Newton-Kantorovich test failed at N = %d (max r = %.3e)',N,max(r));
end
% here:  tau_j \in [tb_j-r_j, tb_j+r_j] \subset (t_{j-1},t_j)  and, by
% monotonicity of kappa_N on each [t_{j-1},t_j], these are all its zeros.

% (3) ------------------------------------------------ enclosure of theta_N
% Kcal_N(tau_j) = Kcal_N(tb) + kappa_N(tb)*(tau-tb) + kappa_N'(xi)*(tau-tb)^2/2
Kbar = Pi.*intval(tb) - intval(tb).^2/2 + 2*tsum(tb,kk,W2,blk,@cos);
Kj   = Kbar + kapI.*infsup(-r,r) + 0.5*dkE.*infsup(0*r,sup(intval(r).^2));

% alternating sum, shifted by the constant Kcal_N(tau_j) ~ pi^2/3 so that no
% cancellation is lost in the accumulation
c   = Pi^2/3;
sg  = (-1).^(j-1);
S1  = c*mod(N,2) + sum(sg.*(Kj-c));
S2  = sum(W2);                      % sum_{k=1}^N 1/k^2
S3  = sum(((-1).^jr).*W2);          % sum_{k=1}^N (-1)^k/k^2

ThN = (2/Pi)*S1 - (2/Pi)*S2 + (-1)^N*( Pi/2 + (2/Pi)*S3 );
end
% =====================================================================
function y = Dn(N,t)
%DN   Dirichlet kernel  D_N(t) = sin((N+1/2)t)/sin(t/2)  (= -kappa_N').
M = intval(2*N+1)/2;                    % (2N+1)/2 exact, division by 2 exact
y = sin(M.*t)./sin(t/2);
end
% ---------------------------------------------------------------------
function s = tsum(tb,kk,w,blk,fun)
%TSUM   s(i) = sum_k fun(tb(i)*k)*w(k), computed blockwise.
%   tb is a double column vector quantised by QUANT, so every product
%   tb(i)*k is exact in binary floating point and intval(tb(i)*k) is a
%   point interval.  kk is a double row vector of the indices k and w the
%   matching column of weights: double (fast, non-verified) or intval
%   (verified).  NOTE: kk is passed in rather than derived from w, because
%   NUMEL is not overloaded for intval (it returns 1 for an interval array).
n = numel(tb);
rig = isa(w,'intval');
if rig, s = intval(zeros(n,1)); else, s = zeros(n,1); end
for i0 = 1:blk:n
    i1 = min(i0+blk-1,n);
    T  = tb(i0:i1)*kk;                  % exact outer product
    if rig, T = intval(T); end
    s(i0:i1) = fun(T)*w;
end
end
% ---------------------------------------------------------------------
function y = quant(x,q)
%QUANT   Round to a multiple of 2^-q (exact operations).
%   With q = 51 - ceil(log2(N)) and 0 < x < pi one has x*k exactly
%   representable in double precision for every k = 1,...,N.
s = pow2(1,q);
y = round(x*s)/s;
end
