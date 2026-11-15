function L2 = fourier_multiplication_matrix(W,N)
%FOURIER_MULTIPLICATION_MATRIX Galerkin matrix of multiplication by w.
%   L2 = FOURIER_MULTIPLICATION_MATRIX(W,N), where W stores the Fourier
%   coefficients W_j, j = -2N..2N, of w = ubar' o theta (equation (4)).
%   In Fourier coordinates (w*h)_n = sum_k W_{n-k} h_k, hence the
%   (2N+1)x(2N+1) Toeplitz matrix L2(n,k) = W_{n-k}, |n|,|k| <= N.
%   This is why Proposition 1(c) needs W_n for |n| <= 2N.
% Integer indexing preserves interval entries without any floating-point work.
assert(iscolumn(W) && length(W)==4*N+1);
modes=(-N:N)';
indices=modes-modes.'+2*N+1;               % position of W_{n-k} in W
L2=reshape(W(indices(:)),2*N+1,2*N+1);
end
