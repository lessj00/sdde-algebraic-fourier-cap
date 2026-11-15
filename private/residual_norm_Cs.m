function [rnorm, Mg, info] = residual_norm_Cs(a, eps, Kv_g, Mp, K)
%RESIDUAL_NORM_CS Pointwise bound on the residual (second part of Lemma 5).
%   [rnorm, Mg, info] = RESIDUAL_NORM_CS(a, eps, Kv_g, Mp, K)
%   Returns rnorm >= ||r||_{C^0} for the residual of the differential equation
%       r(t) = ubar'(t) - g(ubar)(t) + M g(ubar),
%       g(u)(t) = -u(t) - eps1*u(t-1-eps2*u(t)) + cos t + sin t,
%   an enclosure Mg of the mean M g(ubar), and info.upi, an enclosure of
%   ubar(-pi). These are the quantities in Lemma 5 and in the bound (9).
%   Method (Lemma 5): q = r - eps1*pi^{>K} v is a trigonometric polynomial of
%   order K (since K >= N_u), so with t_k = 2*pi*k/Mp (Mp = M' in the paper),
%       ||q||_{C^0} <= max_k |q(t_k)| / cos(pi*K/Mp),
%   and |q(t_k)| <= |r(t_k)| + |eps1|*T_K. This gives
%       ||r||_{C^0} <= (max_k |r(t_k)| + |eps1| T_K)/cos(pi K/Mp) + |eps1| T_K,
%       T_K >= sum_{|n|>K} |V_n|   (Lemma 1 with the sharp constants of Lemma 4).
%   As described after Lemma 5, r(t_k) is evaluated in real interval
%   arithmetic, with the phases n*t_k reduced modulo 2*pi in exact integer
%   arithmetic and the sums taken from the highest modes downward. This
%   reduces, but does not remove, the accumulation of rounding errors;
%   info.rad_r records the largest interval radius of r(t_k).
%   INPUTS : a = (Ubar_n)_{n=-N_u}^{N_u} (unpadded), eps = (eps1,eps2),
%            Kv_g (Cs_grid_SDDE), Mp = number of sample points (M'),
%            K = cutoff with N_u <= K < Mp/2. Requires INTLAB.

N  = (length(a)-1)/2;
if K < N || 2*K >= Mp, error('need N_u <= K < Mp/2'), end
assert(N*Mp+Mp/2<2^53,'Grid phase products must be exact integers.');
e1 = intval(eps(1)); e2 = intval(eps(2));
% Ubar_n = A_n + i B_n for n = 1..N, and Ubar_0 = U0 (real).
A  = real(intval(a(N+2:end)));  B = imag(intval(a(N+2:end)));
U0 = real(intval(a(N+1)));
ipi = intval('pi');
kk = (-Mp/2:Mp/2-1)';                 % grid t_k = 2*pi*k/Mp in [-pi,pi)
t  = 2*ipi*intval(kk)/Mp;

% ubar(t_k)  = U0 + 2*sum_{n>=1} (A_n cos nt_k - B_n sin nt_k),
% ubar'(t_k) = -2*sum_{n>=1} n (A_n sin nt_k + B_n cos nt_k),
% summed from n = N down to 1. The phases n*t_k are reduced EXACTLY in
% integer arithmetic, n*t_k = 2*pi*m/Mp mod 2*pi with m in [-Mp/2,Mp/2), so
% that the radius of cos(n*t_k) does not grow with n.
u  = intval(zeros(Mp,1)); up = intval(zeros(Mp,1));
for n = N:-1:1
    m  = mod(n*kk + Mp/2, Mp) - Mp/2;  % exact (n*Mp < 2^53)
    ph = 2*ipi*intval(m)/Mp;
    c = cos(ph); s = sin(ph);
    u  = u  + (A(n)*c - B(n)*s);
    up = up - n*(A(n)*s + B(n)*c);
end
u = U0 + 2*u; up = 2*up;
% v(t_k) = ubar(theta(t_k)), theta(t) = t - 1 - eps2*ubar(t), equation (4).
th = t - 1 - e2*u;
v  = intval(zeros(Mp,1));
for n = N:-1:1
    v = v + (A(n)*cos(n*th) - B(n)*sin(n*th));
end
v = U0 + 2*v;

% T_K >= sum_{|n|>K} |V_n|: Lemma 1 with p = 0, minimized over s.
% Using s >= 3 (Lemma 1 allows s >= 2) is conservative.
T = inf; Kt = intval(K);
for s = 3:length(Kv_g)
    T = min(T, sup(2*Kv_g(s)/((s-1)*Kt^(s-1))));
end
% Mean of v: V_0 = (1/Mp) sum_k v(t_k) - sum_{l~=0} V_{l*Mp} (Poisson
% summation, proof of Lemma 2). The alias frequencies |l*Mp| exceed K, so
% the whole Fourier tail T also bounds this aliasing error (conservative but rigorous).
V0 = sum(v)/Mp + infsup(-1,1)*T;
% M g(ubar) = mean of (-ubar - eps1*v + cos + sin) = -Ubar_0 - eps1*V_0.
Mg = -U0 - e1*V0;

% r = ubar' - g(ubar) + M g(ubar) = ubar' + ubar + eps1*v - cos - sin + M g(ubar).
r = up + u + e1*v - cos(t) - sin(t) + Mg;
rmax  = max(sup(abs(r)));
% Sampling inequality of Lemma 5 (Ehlich-Zeller type, [4, Sect. 2.6.6]).
Cinf  = 1/cos(ipi*K/Mp);
rnorm = sup(Cinf*(rmax + abs(e1)*T) + abs(e1)*T);
info  = struct('rmax',rmax,'T',T,'Cinf',sup(Cinf),'rad_r',max(rad(r)), ...
               'upi',u(1));   % enclosure of ubar(-pi) (t_1 = -pi), used in (9)
end
