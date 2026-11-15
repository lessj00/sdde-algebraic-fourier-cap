function DF = DF_SDDE_rig_Cs(x,eps,Kw_c,Bs,M,N0)
%DF_SDDE_RIG_CS Interval enclosure of the Galerkin matrix pi^{<=N} DF(xbar) pi^{<=N}.
%   DF = DF_SDDE_RIG_CS(x,eps,Kw_c,Bs,M,N0)
%   The Frechet derivative of F (Sect. 2) is
%       DF(beta,u) = [ 0    M dg(u)        ]
%                    [ -1   I - L1 dg(u)   ],
%   where, with theta, w and e_k as in (4),
%       dg(ubar) h = -h + eps1*eps2*w*h - eps1*(h o theta).
%   In Fourier coordinates this is -I + eps1*eps2*L2 - eps1*L3, where
%     L2(n,k) = W_{n-k}    : multiplication by w,   Proposition 1(c),
%     L3(n,k) = (L3)_{n,k} : n-th coefficient of e_k, Proposition 1(d).
%   The zeroth row of L1*dg(ubar) involves ALL Fourier coefficients H_j of
%   h = dg(ubar) e^{ik.}; it is enclosed separately with Proposition 1(e).
%   Inputs: x = (betabar, Ubar_{-N..N}), eps = (eps1,eps2),
%   Kw_c(s+1) = calK_s(w) and Bs(s+1,j+1) = B_{s,j}(alpha) from Lemma 3 (7),
%   M (DFT grid of Prop. 1(c),(d)) and N0 (DFT grid of Prop. 1(e)).
%   The beta column is exact: (0,-e_0). All other arithmetic is interval.
%   Modes are ordered -N:N.

eps = intval(eps);
a = intval(x(2:end));
N = (length(a)-1)/2; n = (-N:N)';
% Hypotheses: M >= 4N (stronger than M > 2N in Prop. 1) and N0 > 2N.
if abs(log2(M)-ceil(log2(M)))>0 || N>M/4, error('Wrong size for M'), end
if abs(log2(N0)-ceil(log2(N0)))>0 || 2*N>=N0, error('Wrong size for N0'), end
smax = size(Bs,1)-1;
ae1 = abs(intval(eps(1))); ae2 = abs(intval(eps(2)));

% Values ubar(t_k) on the grid t_k = 2*pi*k/M, k = -M/2..M/2-1.
ta  = [intval(zeros(M/2-N,1));a;intval(zeros(M/2-N-1,1))];
u_k = real(ifftshift(verifyfft(fftshift(ta),-1))*M);
t   = (2*intval('pi')/M)*(-M/2:M/2-1)';

% Galerkin matrix of L1 (see F_SDDE_rig_Cs); its zeroth row is replaced below.
L1 = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*intval(j));
    L1(N+1,j+N+1)   = (-1)^(j+1)/(1i*intval(j));
end

% ---------- L2 : multiplication by w = ubar'(t-1-eps2 ubar(t)), Prop. 1(c)
% w(t_k) = sum_n (i n) Ubar_n exp(i n theta(t_k)).
V = intval(zeros(M,1));
for k = 0:M-1
    V(k+1) = real(sum((1i*intval(n)).*a.*exp(1i*intval(n)*(t(k+1)-1-eps(2)*u_k(k+1)))));
end
% Wbar_n: M-point DFT of w (the variable V holds w here).
Vbar = fftshift(verifyfft(ifftshift(V))/intval(M));
n2 = (-2*N:2*N)';
% Aliasing error calK_s(w) E_s(n,M) for |n| <= 2N, minimized over s >= 2.
err = inf(4*N+1,1);
for s = 2:smax
    err = min(err, sup(Kw_c(s+1)*alias_bound(n2,M,s)));
end
% The positive Nyquist frequency shares a DFT entry with the negative one
% when M=4N. Wrap indices; each true coefficient keeps its own alias disc.
% (For the parameters of the paper, M > 4N and no wrapping occurs.)
indices=mod(n2+M/2,M)+1;
V = Vbar(indices,1) + complex_error_disc(err);       % W_n, |n| <= 2N
L2 = fourier_multiplication_matrix(V,N);

% ---------- L3 : columns e_k = e^{ik theta(t)}, Prop. 1(d)
% Kcol(s+1,k+N+1) = calK_s(e_k) = sum_{j=1}^s |k|^j B_{s,j}(alpha) (Lemma 3),
% used only for s >= 2.
Kcol = intval(zeros(smax+1,2*N+1));
for s = 2:smax
    for j = 1:s
        Kcol(s+1,:) = Kcol(s+1,:) + abs(intval(n')).^j*Bs(s+1,j+1);
    end
end
% Es(:,s+1) >= E_s(n,M), |n| <= N (Lemma 2).
Es = intval(zeros(2*N+1,smax+1));
for s = 2:smax, Es(:,s+1) = alias_bound(n,M,s); end
L3 = intval(zeros(2*N+1));
for k = -N:N
    if k == 0
        % e_0 = 1: the column k = 0 is exact (calK_s(e_0) = 0 for s >= 1).
        L3(:,N+1) = intval(zeros(2*N+1,1));
        L3(N+1,N+1) = 1;
        continue
    end
    % (L3bar)_{n,k}: M-point DFT of e_k.
    col = exp(1i*intval(k)*(t - 1 - eps(2)*u_k));
    col_bar = fftshift(verifyfft(ifftshift(col))/intval(M));
    % (L3)_{n,k} in (L3bar)_{n,k} + calK_s(e_k) E_s(n,M) B, best s >= 2.
    err = inf(2*N+1,1);
    for s = 2:smax
        err = min(err, sup(Kcol(s+1,k+N+1)*Es(:,s+1)));
    end
    L3(:,k+N+1) = col_bar(M/2-N+1:M/2+N+1,1) + complex_error_disc(err);
end

% dg(ubar) = -I + eps1*eps2*L2 - eps1*L3 and the function block I - L1*dg(ubar).
DGu = -intval(eye(2*N+1)) + eps(1)*eps(2)*L2 - eps(1)*L3;
DFu = intval(eye(2*N+1)) - L1*DGu;

% ---------- zeroth row of L1*DGu, Proposition 1(e):
% (L1 h)_0 = sum_{j~=0} (-1)^{j+1} H_j/(ij),  H_j = coefficients of
% h = dg(ubar) e^{ik.} = -e^{ik.} + eps1*eps2*w*e^{ik.} - eps1*e_k.
% It is approximated by the sum over 0 < |j| < N0/2 of the N0-point DFT
% coefficients Hbar_j; c0 holds the weights (-1)^{j+1}/(ij).
j0 = (-N0/2:N0/2-1)';
c0 = (-1).^(j0+1)./(1i*intval(j0+(j0==0))); c0([1 N0/2+1]) = 0;
% Both Nyquist endpoints are omitted: the remainder below includes |j|>=N0/2.
% S(s+1) = S_s(N0) of Prop. 1(e):
%   S_s(N0) = 2*sum_{j=1}^{N0/2-1} E_s(j,N0)/j + 2*((N0/2)^{-s-1} + (N0/2)^{-s}/s),
% where the last term bounds 2*sum_{j>=N0/2} j^{-s-1} (omitted modes).
S = intval(zeros(smax+1,1)); jj = (1:N0/2-1)'; J = intval(N0/2);
for s = 2:smax
    S(s+1) = 2*sum(alias_bound(jj,N0,s)./intval(jj)) + 2*(1/J^(s+1) + 1/(s*J^s));
end
% Values of ubar and w on the N0-point grid.
t0 = (2*intval('pi')/N0)*j0;
u0 = real(ifftshift(verifyfft(fftshift([intval(zeros(N0/2-N,1));a;intval(zeros(N0/2-N-1,1))]),-1))*N0);
V0 = intval(zeros(N0,1));
for m = 1:N0
    V0(m) = real(sum((1i*intval(n)).*a.*exp(1i*intval(n)*(t0(m)-1-eps(2)*u0(m)))));
end
for k = -N:N
    % Error bound of Prop. 1(e), minimized over s >= 2:
    %   (|eps1 eps2| sum_q binom(s,q)|k|^{s-q} calK_q(w) + |eps1| calK_s(e_k)) S_s(N0),
    % where the Leibniz rule bounds K_s(w e^{ik.}). For k = 0 the sum
    % reduces to calK_s(w), and calK_s(e_0) = 0.
    % K_s of the non-polynomial part eps1*eps2*w*e^{ikt} - eps1*e^{ik theta}
    err = inf;
    for s = 2:smax
        if k == 0
            Kh = Kw_c(s+1);
        else
            Kh = intval(0);
            for i = 0:s
                Kh = Kh + nchoosek(s,i)*abs(intval(k))^(s-i)*Kw_c(i+1);
            end
        end
        Kh = ae1*ae2*Kh + ae1*Kcol(s+1,k+N+1);
        err = min(err, sup(Kh*S(s+1)));
    end
    % The polynomial -exp(i*k*t) has no DFT aliasing or omitted modes
    % because |k|<N0/2. Kh therefore bounds only the other two terms.
    % h(t0) = e^{ikt}(eps1*eps2*w - 1) - eps1*e^{ik theta}, and its DFT G = Hbar.
    g = exp(1i*intval(k)*t0).*(eps(1)*eps(2)*V0 - 1) - eps(1)*exp(1i*intval(k)*(t0 - 1 - eps(2)*u0));
    G = fftshift(verifyfft(ifftshift(g))/intval(N0));
    % Zeroth row of I - L1*dg(ubar): delta_{k0} - (L1 h)_0, with its error disc.
    DFu(N+1,k+N+1) = (k==0) - sum(c0.*G) + complex_error_disc(err);
end

% Assemble pi^{<=N} DF(xbar) pi^{<=N}: first row (0, M dg(ubar)) = (0, zeroth
% row of DGu), first column (0, -e_0) from the derivative in beta.
e0 = intval(zeros(2*N+1,1)); e0(N+1) = 1;
DF = [0 DGu(N+1,:); -e0 DFu];
end
