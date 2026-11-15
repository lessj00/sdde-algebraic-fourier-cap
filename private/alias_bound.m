function E = alias_bound(n,M,s)
%ALIAS_BOUND Rigorous upper bound for the aliasing sum of Lemma 2.
%   E = ALIAS_BOUND(n,M,s) returns an interval upper bound for
%       E_s(n,M) = sum_{l ~= 0} |n + l*M|^{-s},   |n| < M,  s >= 2,
%   namely (Lemma 2 of the paper)
%       E_s(n,M) <= phi_s(M-|n|) + phi_s(M+|n|),
%       phi_s(c)  = c^{-s} + 1/((s-1)*M*c^{s-1}).
%   By Lemma 2, |F_n - Fbar_n| <= K_s(f)*E_s(n,M), where Fbar_n is the
%   M-point DFT coefficient of f. n may be a vector.
n = abs(intval(n)); M = intval(M);
% phi_s(c) = c^{-s} + c^{1-s}/((s-1)M): the first term of sum_{j>=0}(c+jM)^{-s}
% plus the integral bound for the remaining terms.
phi = @(c) 1./c.^s + 1./(c.^(s-1)*((s-1)*M));
E = phi(M-n) + phi(M+n);
end
