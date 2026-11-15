function [Y, Z1, Z2, info] = computation_bounds_Cs(x,eps,N_FFT,M,N0,Kv_c,Kw_c,Bs,Kv_g,N_Z1,resopt)
%COMPUTATION_BOUNDS_CS Bounds Y, Z1, Z2 of Theorem 1 for the forced model (3).
%   [Y,Z1,Z2,info] = COMPUTATION_BOUNDS_CS(x,eps,N_FFT,M,N0,Kv_c,Kw_c,Bs,Kv_g,N_Z1,resopt)
%
%   Inputs
%     x      = (betabar, Ubar_{-N}, ..., Ubar_N), N >= N_u + 1 (zero padded)
%     eps    = (eps1, eps2)
%     N_FFT  : grid of the C^0 operator norms (comp_mat_norm, [4, Sect. 2.6.6])
%     M      : DFT grid for V_n, W_n, (L3)_{n,k}  (Prop. 1(a),(c),(d))
%     N0     : DFT grid for the zeroth row of L1 dg(ubar)  (Prop. 1(e))
%     Kv_c, Kw_c, Bs : a priori constants of Lemma 3 (Cs_crude_SDDE)
%     Kv_g(s) >= tildeK_s(v) >= K_s(v) : sharp constants of Lemma 4 (Cs_grid_SDDE)
%     N_Z1   : size of the block on which A is a numerical inverse (default N)
%     resopt : struct with fields a (unpadded coefficients of ubar), Mp (= M')
%              and K, for the residual-based bounds (Lemma 5 and (9)).
%
%   Outputs: Y >= ||A F(xbar)||_H, Z1 >= ||I - A DF(xbar)||, Z2 (global
%   Lipschitz bound on A DF), and a struct info with all intermediate bounds.
%
%   Y is the smallest of three bounds (see the proof of Theorem 2):
%     (i)   Y_Fourier, equation (8), from Proposition 1(a),(b);
%     (ii)  Y_res, the first estimate of Lemma 5;
%     (iii) Y_pw, the bound (9), which involves no sum over the Fourier
%           coefficients of F(xbar). In both runs of the paper (iii) is used.
%   Z1 = Z_{1,0} + Z_{1,1} + Z_{1,2} and Z2 follow [4, Sect. 2.6.2-2.6.4],
%   with the Galerkin matrix enclosed by Proposition 1(c),(d),(e).
%   Only Y and Z_{1,0} use the algebraic Fourier estimates of this paper;
%   Z_{1,1}, Z_{1,2} and Z2 only involve ||ubar'||, ||ubar''|| and theta_N.

N = (length(x)-2)/2;
if nargin < 10 || isempty(N_Z1), N_Z1 = N; end
if nargin < 11, resopt = []; end
assert(N_Z1>=1 && N_Z1<=N && N_Z1==fix(N_Z1));
ind = N-N_Z1+1 : N+N_Z1+1;   io = [1, 1+ind];

%% Approximate inverse A = A^{<=N} + pi^{>N} (Sect. 2).
% A^{<=N} is the floating-point inverse of the numerical Galerkin Jacobian.
% It is a fixed matrix; its quality is measured rigorously by Z_{1,0}.
DF = DF_SDDE([x(1); x(1+ind)],eps,max(2^nextpow2(4*N_Z1),min(M,2^14))); % heuristic grid for A only
A = intval(eye(2*N+2));
Anum = inv(DF);  % A fixed point matrix, not an interval inverse.
% Enforce exact conjugate symmetry: A maps real Fourier series to real ones
% (Sect. 2: A must preserve real-valued functions). perm maps mode n to -n.
perm = [1, 2*N_Z1+2:-1:2];
Anum = (Anum + conj(Anum(perm,perm)))/2;
for ii = 1:size(Anum,1)
    for jj = 1:size(Anum,2)
        ip = perm(ii); jp = perm(jj);
        if ii == ip && jj == jp
            Anum(ii,jj) = real(Anum(ii,jj));
        elseif sub2ind(size(Anum),ii,jj) < sub2ind(size(Anum),ip,jp)
            Anum(ip,jp) = conj(Anum(ii,jj));
        end
    end
end
assert(isequal(Anum,conj(Anum(perm,perm))));
A(io,io) = intval(Anum);

eps = intval(eps);
a = intval(x(2:end));
n = intval((-N:N)');
ae1 = abs(eps(1)); ae2 = abs(eps(2));

%% (i) Y_Fourier: equation (8), Proposition 1(a,b).
% For |n| > N the coefficients of the function component of F(xbar) are
% eps1*V_n/(in), and A acts as the identity there. By Prop. 1(b),
% tail_l1 bounds |eps1| sum_{|n|>N}|V_n|/|n|  >= C^0 norm of this tail,
% tail_der bounds |eps1| sum_{|n|>N}|V_n|     >= C^0 norm of its derivative.
% Taking only s>=2 and s>=3 respectively (the paper allows s>=1 and s>=2)
% is conservative; it omits otherwise valid lower orders without changing
% the mathematical bounds.
disp('F_SDDE_rig_Cs... ')
F = F_SDDE_rig_Cs(x,eps,Kv_c,M);
smax_g = length(Kv_g);
tail_l1 = inf; tail_der = inf; Nt = intval(N);
for s = 2:smax_g
    tail_l1 = min(tail_l1, sup(2*ae1*Kv_g(s)/(s*Nt^s)));
    if s >= 3
        tail_der = min(tail_der, sup(2*ae1*Kv_g(s)/((s-1)*Nt^(s-1))));
    end
end
Fhat = F;                                         % Galerkin part, zeroth row truncated at N
% The zeroth coefficient of L1 g(ubar) also contains
% sum_{|j|>N} (-1)^{j+1} G_j/(ij) with G_j = -eps1*V_j, of modulus <= tail_l1.
F(N+2) = F(N+2) + infsup(-1,1)*tail_l1;          % tail of the zeroth row of L1 G
AF = A*F;
% Y0 = ||pi^{<=N} A F(xbar)||_H <= max(|beta|, sum|u_n|, sum|n u_n|).
Y0 = max([abs(AF(1)) norm(AF(2:end),1) norm(n.*AF(2:end),1)]);
% H-norm of the tail: max of the C^0 norms of the tail and of its derivative.
Y1 = max(tail_l1,tail_der);
Y  = Y0 + Y1;
disp(['Y <= ',num2str(sup(Y)),'   (Galerkin part ',num2str(sup(Y0)),', tail ',num2str(Y1),')'])

%% (ii) residual-based Y: first estimate of Lemma 5,
%   ||A F(xbar)||_H <= ||A^{<=N} Fhat||_H + rho*(1 + theta_N ||A^{<=N} 1||_H),
%   rho >= ||r||_{C^0} + sum_{0<|n|<=N} |n fhat_n|,  r = ubar' - g(ubar) + M g(ubar).
% theta_N is the Kolmogorov constant (Sect. 2, [4, Sect. 2.6.3]).
thetaN = kolmogorov_theta_N(N);
Yres = inf; Ypw = inf; rinfo = [];
if ~isempty(resopt)
    disp('Pointwise residual... ')
    % rnorm >= ||r||_{C^0} (second estimate of Lemma 5), Mg_res ~ M g(ubar).
    [rnorm,Mg_res,rinfo] = residual_norm_Cs(resopt.a,eps,Kv_g,resopt.Mp,resopt.K);
    % e = 1 = (0,1) in H: the constant function 1 (mode 0).
    e  = intval(zeros(2*N+2,1)); e(N+2) = 1;
    AFh = A*Fhat; Ae = A*e;
    Y0h = max([abs(AFh(1)) norm(AFh(2:end),1) norm(n.*AFh(2:end),1)]);
    normAe = max([abs(Ae(1)) norm(Ae(2:end),1) norm(n.*Ae(2:end),1)]);
    rho  = rnorm + norm(n.*Fhat(2:end),1);
    Yres = sup(Y0h + rho*(1 + thetaN*normAe));
    disp(['   ||r||_C0 <= ',num2str(rnorm),' (max|r(t_k)| = ',num2str(rinfo.rmax), ...
          ', max rad = ',num2str(rinfo.rad_r),'),  Y_res <= ',num2str(Yres)])
end
Y = min(sup(Y),Yres);
disp(['Y <= ',num2str(Y)])

%% Operator norms and Z2: [4, Sect. 2.6.2 and 2.6.4]
% ("from infinite to finite" and "Z2 forcing" lemmas of the preprint).
% Ddag = D^dagger(beta,u) = (beta, L1 u) and D(beta,u) = (beta,u') (Sect. 2).
% normAD/normDAD bound A^<=N Ddag and D A^<=N Ddag on H0.
% Their theta_N*||A e|| terms restore the infinite zeroth row of Ddag.
% norm_ADdag bounds max(||A Ddag||_H0, ||D A Ddag||_H0),
% equivalently ||A Ddag|| from H0 to H. The +theta_N and +1 terms
% account for the identity extension of A on the Fourier tail.
L1 = intval(zeros(2*N+1));
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/intval(1i*j);
    L1(N+1,j+N+1)   = (-1)^(j+1)/intval(1i*j);
end
Ddag = [intval(1) intval(zeros(1,2*N+1)); intval(zeros(2*N+1,1)) L1];
D = diag([1 ; 1i*n]);
if sup(thetaN) > 1, error('theta_N > 1'), end
e = intval(zeros(2*N+2,1)); e(N+2) = 1;
DADdag  = D*A*Ddag;
normAD  = comp_mat_norm(A*Ddag,N_FFT) + thetaN*norm(A*e,1);
normDAD = comp_mat_norm(DADdag,N_FFT) + thetaN*norm(D*(A*e),1);
PiJ = diag([intval(0) ; intval(ones(2*N+1,1))]);
norm_ADdag = max(normAD + thetaN, comp_mat_norm(DADdag-PiJ,N_FFT) + thetaN*norm(D*(A*e),1) + 1);
norm_der_u  = norm(n.*a,1);                       % >= ||ubar'||_{C^0}
norm_der2_u = norm(n.^2.*a,1);                    % >= ||ubar''||_{C^0}
% Global Z2 bound of [4, Sect. 2.6.2]; it holds for every r_*, so r_* = r0.
Z2 = norm_ADdag*ae1*ae2*(ae2*norm_der2_u + 2);
disp(['Z2 <= ',num2str(sup(Z2))])

%% (iii) fully pointwise Y: equation (9) of the paper,
%   A F = (A D^dagger)(M g(ubar), r) + f(-pi) A 1,  f(-pi) = ubar(-pi) - betabar,
%   ||A F(xbar)||_H <= ||A D^dagger||_{H0->H} max(|M g(ubar)|, ||r||_{C^0})
%                      + |ubar(-pi) - betabar| ||A^{<=N} 1||_H.
% No sum over Fourier coefficients of F(xbar) is involved.
if ~isempty(resopt)
    Ypw = sup(norm_ADdag*max(sup(abs(Mg_res)),rnorm) + abs(rinfo.upi - x(1))*normAe);
    disp(['   fully pointwise:  Y_pw <= ',num2str(Ypw)])
    Y = min(Y,Ypw);
    disp(['Y <= ',num2str(Y)])
end

%% Z1 = Z_{1,0} + Z_{1,1} + Z_{1,2}: [4, Sect. 2.6.4]
% ("Z1 forcing" and "Z1_12 forcing" lemmas of the preprint).
% Z10 is the finite defect after factoring Ddag. z10 restores the
% evaluation-at-minus-pi term: ||B||_(H->H0) <= ||B Ddag||_H0
% +(1+theta_N)||B e||_H0, for B=defect and D*defect.
% Z11 controls low-output/high-input coupling; Z12 controls high output.
% Thus source Z_{1,0} is bounded here by Z10+z10, not Z10 alone
% (this is the "finite-projection remainder" mentioned in Sect. 4).
norm_dug   = 1 + ae1*(1 + ae2*norm_der_u);
norm_dtdug = 1 + ae1*(1 + ae2*norm_der_u)*(1 + ae2*norm_der2_u) + ae1*ae2*norm_der_u;
Z11 = max(normAD,normDAD)*norm_dug*thetaN;
Z12 = norm_dtdug*thetaN;
disp('DF_SDDE_rig_Cs... ')
% Interval enclosure of pi^{<=N} DF(xbar) pi^{<=N} (Prop. 1(c),(d),(e)).
DF = DF_SDDE_rig_Cs(x,eps,Kw_c,Bs,M,N0);
defect = intval(eye(2*N+2)) - A*DF;
Z10 = max(comp_mat_norm(defect*Ddag,N_FFT), comp_mat_norm(D*defect*Ddag,N_FFT));
z10 = (1+thetaN)*max(norm(defect*e,1), norm(D*(defect*e),1));
Z1  = Z10 + z10 + Z11 + Z12;
disp(['   Z_{1,0} <= ',num2str(sup(Z10+z10)),' (matrix ',num2str(sup(Z10)), ...
      ', evaluation ',num2str(sup(z10)),') , Z_{1,1} <= ',num2str(sup(Z11)), ...
      ' , Z_{1,2} <= ',num2str(sup(Z12))])
disp(['Z1 <= ',num2str(sup(Z1))])
info = struct('Yfourier',sup(Y0+Y1),'Ypointwise',Ypw,'Y',Y,'Yres',Yres,'rinfo',rinfo,'Y0',sup(Y0),'Y1',Y1,'Z10',sup(Z10),'z10',sup(z10),'Z11',sup(Z11), ...
              'Z12',sup(Z12),'normAD',sup(normAD),'normDAD',sup(normDAD),'thetaN',sup(thetaN),'normADdag',sup(norm_ADdag));
end
