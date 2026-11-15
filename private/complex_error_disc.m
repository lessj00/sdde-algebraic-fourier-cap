function z = complex_error_disc(r)
%COMPLEX_ERROR_DISC Complex interval containing the closed disc r*B.
%   z = COMPLEX_ERROR_DISC(r) returns a complex INTLAB interval that
%   contains every complex number of modulus at most r, that is r*B with B
%   the closed unit disc of Proposition 1. Adding z to a DFT coefficient
%   gives the enclosures of Proposition 1(a),(c),(d),(e).
% Explicit enclosing complex disc. Subtracting i keeps this valid even if
% a zero midpoint would otherwise cause conversion to a real interval.
r = sup(intval(r));
assert(all(isfinite(r(:))) && all(r(:)>=0));
z = midrad(1i*ones(size(r)),r) - intval(1i);
end
