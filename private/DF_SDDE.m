function DF = DF_SDDE(x,eps,N_FFT)
%DF_SDDE Floating-point Galerkin Jacobian, used ONLY to choose the matrix A.
%   DF = DF_SDDE(x,eps,N_FFT) is a non-rigorous (double precision) version
%   of DF_SDDE_rig_Cs. Its numerical inverse gives the finite block
%   A^{<=N} of the approximate inverse A = A^{<=N} + pi^{>N} (Sect. 2).
%   No bound or inverse property of this numerical matrix is assumed:
%   DF_SDDE_rig_Cs independently encloses the true Galerkin derivative, and
%   computation_bounds_Cs validates the defect I - A*DF with interval
%   arithmetic (bound Z_{1,0}).
%   Fourier modes are ordered -N:N; the scalar beta is the first coordinate.

a = x(2:end) ;
N = (length(a)-1)/2; n = (-N:N)';

if abs(log2(N_FFT)-ceil(log2(N_FFT)))>0 || N>N_FFT/4
    error('DF_SDDE:grid','Need a power-of-two grid with N_FFT >= 4N.');
end

% Zero padding to size N_FFT
ta = [zeros(N_FFT/2-N,1);a;zeros(N_FFT/2-N-1,1)];
% Values of ubar at the grid points (inverse FFT)
u_k = ifftshift(ifft(fftshift(ta)))*N_FFT;
% Grid t_k = 2*pi*k/N_FFT, k = -N_FFT/2..N_FFT/2-1
t = (2*pi/N_FFT)*(-N_FFT/2:N_FFT/2-1)';

%% Galerkin matrix of the periodic primitive L1 (Sect. 2)
L1 = zeros(2*N+1);
for j = [-N:-1 1:N]
    L1(j+N+1,j+N+1) = 1/(1i*j) ;
end
% zeroth row: (L1 h)_0 = sum_{j~=0} (-1)^{j+1} H_j/(ij), truncated to |j| <= N
for j = [-N:-1 1:N]
    L1(N+1,j+N+1) = (-1)^(j+1)/(1i*j) ;
end

% operator L2: multiplication by w = ubar'(t-1-eps_2 ubar(t)) (DFT coefficients)
V = zeros(N_FFT,1);
for k = 0:N_FFT-1
    V(k+1) = real(sum((1i*n).*a.*exp(1i*n*(t(k+1)-1-eps(2)*u_k(k+1)))));
end
V = fftshift(fft(ifftshift(V))/N_FFT);
% Wrap the positive Nyquist mode if N_FFT=4N; DFT indices are periodic.
indices=mod((-2*N:2*N)'+N_FFT/2,N_FFT)+1;
V = V(indices,1);

L2 = fourier_multiplication_matrix(V,N);

% operator L3: column k holds the DFT coefficients of e_k = exp(i*k*theta)
L3 = zeros(2*N+1);

for k =-N:N
    col = exp(1i*k*(t - 1 - eps(2)*u_k));
    col = fftshift(fft(ifftshift(col))/N_FFT);
    L3(:,k+N+1) = col(N_FFT/2-N+1:N_FFT/2+N+1,1);
end

% dg(ubar) = -I + eps1*eps2*L2 - eps1*L3, function block I - L1*dg(ubar)
DGu = -eye(2*N+1) + eps(1)*eps(2)*L2 - eps(1)*L3;
DFu = eye(2*N+1) - L1*DGu;

% the zeroth row of L1*DGu involves all the Fourier coefficients of the
% functions DGu e_k, which we compute with N0 points (cf. Prop. 1(e))
N0 = 4*N_FFT; j0 = (-N0/2:N0/2-1)';
c0 = (-1).^(j0+1)./(1i*(j0+(j0==0))); c0([1 N0/2+1]) = 0;
t0 = (2*pi/N0)*j0;
u0 = ifftshift(ifft(fftshift([zeros(N0/2-N,1);a;zeros(N0/2-N-1,1)])))*N0;
V0 = zeros(N0,1);
for m = 1:N0
    V0(m) = real(sum((1i*n).*a.*exp(1i*n*(t0(m)-1-eps(2)*u0(m)))));
end
for k =-N:N
    % h = -e^{ikt} + eps1*eps2*w*e^{ikt} - eps1*e^{ik theta} on the N0-grid
    g = exp(1i*k*t0).*(eps(1)*eps(2)*V0 - 1) - eps(1)*exp(1i*k*(t0 - 1 - eps(2)*u0));
    DFu(N+1,k+N+1) = (k==0) - sum(c0.*fftshift(fft(ifftshift(g))/N0));
end


e0 = zeros(2*N+1,1) ; e0(N+1) = 1;

% Assemble [0 M dg(ubar); -e_0 I - L1 dg(ubar)]
DF = [0  DGu(N+1,:);
      -e0  DFu   ];

end
