function F = F_SDDE_rig_Cs(x,eps,Kv_c,M)
%F_SDDE_RIG_CS Interval enclosure of the Galerkin part of F(xbar).
%   F = F_SDDE_RIG_CS(x,eps,Kv_c,M) encloses the finite vector Fhat of
%   Sect. 4 of the paper: the scalar component and the Fourier modes
%   |n| <= N of F(xbar), where (Sect. 2)
%       F(beta,u) = ( M g(u) , u - beta - L1 g(u) ),
%       g(u)(t)   = -u(t) - eps1*v(t) + cos t + sin t,   v = u o theta,
%   so that the Fourier coefficients of g(ubar) are G = cos + sin - Ubar - eps1*V.
%   x = (betabar, Ubar_{-N}, ..., Ubar_N) with N >= N_u + 1.
%   The coefficients V_n, |n| <= N, are enclosed by Proposition 1(a):
%       V_n in Vbar_n + calK_s(v) E_s(n,M) B,    s >= 2,
%   with Vbar_n the M-point DFT of v computed with interval arithmetic,
%   E_s(n,M) from Lemma 2 (alias_bound) and calK_s(v) = Kv_c(s+1) from (7).
%   The best s in 2..smax is used for each n.
%   The zeroth coefficient of L1 g(ubar) is sum_{j~=0} (-1)^{j+1} G_j/(ij);
%   here it is truncated to 0 < |j| <= N, which defines Fhat (Sect. 4). The
%   missing part (|j| > N) is added in computation_bounds_Cs using the tail
%   bound of Proposition 1(b).

eps = intval(eps);
beta = intval(x(1)); a = intval(x(2:end));
N = (length(a)-1)/2; n = intval((-N:N)');
if abs(log2(M)-ceil(log2(M)))>0 || N>M/4
    error('Wrong size for M')
end
% Fourier coefficients of sin t and cos t (modes -N..N).
sin_t = [intval(zeros(N-1,1));1i/intval(2);0;-1i/intval(2);intval(zeros(N-1,1))];
cos_t = [intval(zeros(N-1,1));1/intval(2);0;1/intval(2);intval(zeros(N-1,1))];

% Values ubar(t_k) on the grid t_k = 2*pi*k/M, k = -M/2..M/2-1 (inverse FFT).
ta  = [intval(zeros(M/2-N,1));a;intval(zeros(M/2-N-1,1))];
u_k = real(ifftshift(verifyfft(fftshift(ta),-1))*intval(M));    % ubar(t_k)
t   = (2*intval('pi')/M)*(-M/2:M/2-1)';
% v(t_k) = ubar(theta(t_k)), theta(t) = t - 1 - eps2*ubar(t), equation (4).
v_at_t_k = intval(zeros(M,1));
for k = 0:M-1
    v_at_t_k(k+1) = real(sum(a.*exp(1i*n*(t(k+1)-1-eps(2)*u_k(k+1)))));
end
% Vbar_n = (1/M) sum_k v(t_k) e^{-i n t_k}: DFT coefficients (Sect. 3).
Vbar = fftshift(verifyfft(ifftshift(v_at_t_k))/intval(M));

% Aliasing error of Proposition 1(a): calK_s(v) E_s(n,M), minimized over s.
smax = length(Kv_c)-1;
err = inf(2*N+1,1);
for s = 2:smax
    err = min(err, sup(Kv_c(s+1)*alias_bound((-N:N)',M,s)));
end
% V_n, |n| <= N: DFT value plus a complex disc of radius err.
V = Vbar(M/2-N+1:M/2+N+1,1) + complex_error_disc(err);

G = cos_t + sin_t - a - eps(1)*V;          % Fourier coefficients of g(ubar)
% Galerkin matrix of the periodic primitive L1 (Sect. 2):
% (L1 h)_j = H_j/(ij) for j ~= 0, and its zeroth row
% (L1 h)_0 = sum_{j~=0} (-1)^{j+1} H_j/(ij), here truncated to |j| <= N.
L1 = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*intval(j));
    L1(N+1,j+N+1)   = (-1)^(j+1)/(1i*intval(j));
end
MG   = G(N+1);                             % M g(ubar) = G_0
calB = [intval(zeros(N,1));beta;intval(zeros(N,1))];   % constant function betabar
F = [MG; a - calB - L1*G];
end
