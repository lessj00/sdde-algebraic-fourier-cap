function [f,details] = comp_mat_norm(B,M)
%COMP_MAT_NORM Upper bound on the operator norm of a Galerkin matrix on H0.
%   [f,details] = COMP_MAT_NORM(B,M) returns an upper bound f for the norm
%   of the finite-dimensional operator with matrix B acting on
%   H0 = R x C^0_per with norm max(|beta|, ||u||_{C^0}) (Sect. 2 of the
%   paper). B is (2N+2)x(2N+2): coordinate 1 is beta, coordinates 2..end are
%   the Fourier modes -N..N. This is the "C^0 operator-norm estimate" of the
%   companion preprint [4, Sect. 2.6.6], used here for Z2, Z_{1,0}, Z_{1,1}
%   and for the norm of A*Ddag in (9). The grid size M is N_FFT of Sect. 5.
%   Method ([4, Sect. 2.6.6]): the function block acts as an integral
%   operator with kernel Q(x,y) = sum_{n,k} B(n,k) exp(i*n*x - i*k*y). Its
%   C^0 -> C^0 norm is bounded by the maximal row sum of the kernel sampled
%   on an M-point grid, multiplied by the constants C_infinity*C_1 of the
%   Marcinkiewicz-Zygmund-type inequalities of [4, Sect. 2.6.6].
%   The two FFTs enclose Q(x_j,y_l)/M without forming its M-by-M matrix.
%   Coefficients are ordered -N:N; samples are ordered 0:M-1, a permutation
%   of the centered grid in [4]. beta is coordinate 1, not Fourier n=0.
N = (size(B,1)-2)/2; d=2*N+1;
assert(size(B,2)==d+1 && N>=1 && N==fix(N));
assert(M==2^nextpow2(M) && 4*N<=M);
B=intval(B); ind=mod((-N:N),M)+1;
H=intval(zeros(M,d)); H(ind,:)=B(2:end,2:end);
H=verifyfft(H,-1)*M; % H(x,k)=sum_n B(n,k) exp(i n x).
% Maximal row sum of the sampled kernel, computed by blocks of rows.
mx=0; block=32;
for first=1:block:M
    rows=first:min(first+block-1,M);
    C=intval(zeros(M,length(rows)));
    C(ind,:)=H(rows,:).';
    C=verifyfft(C)/M; % sampled integral kernel; no conjugate transpose.
    mx=max(mx,max(sup(sum(abs(C),1))));
end
clear H C
% The constant C_infinity*C_1 depends only on (N,M) and is cached.
persistent cacheN cacheM cacheC
if isempty(cacheN) || cacheN~=N || cacheM~=M
    % De la Vallee Poussin kernel V has degree L and coefficients 1
    % for |k|<=N, tapering to 0 at |k|=M-N ([4, Sect. 2.6.6]).
    % To bound its normalized L1 norm, reconstruct V on the much finer
    % P-grid. The Fejer-kernel triangle bound gives C_1(P,L)<=P/(P-2L).
    % Thus normalized ||V||_1 <= P/(P-2L) * mean_j |V(2*pi*j/P)|.
    % This factor comes from Fourier reconstruction, not a Riemann error.
    % The factor 1/cos(pi*N/M) is the sampling constant C_infinity.
    L=M-N-1; P=256*M; k=(-L:L)';
    Vh=min(intval(1),(intval(L+1)-abs(intval(k)))/intval(L+1-N));
    vv=intval(zeros(P,1)); vv(mod(k,P)+1)=Vh;
    V=verifyfft(vv,-1)*P;
    cacheC=sup(intval(P)/intval(P-2*L)*sum(abs(V))/P/cos(intval('pi')*N/M));
    cacheN=N; cacheM=M;
end
% |u_n|<=||u||_C0 bounds the scalar row; the beta-to-function column
% uses a coefficient l1 bound. max combines the two H0 output components.
f=max(norm(B(1,:),1),intval(cacheC)*intval(mx)+norm(B(2:end,1),1));
if nargout>1
    details=struct('sampledKernel',mx,'samplingFactor',cacheC);
end
end
