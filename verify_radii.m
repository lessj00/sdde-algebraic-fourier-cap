function certificate = verify_radii(Y,Z1,Z2)
%VERIFY_RADII Verify the Newton-Kantorovich inequalities of Theorem 1.
%   certificate = VERIFY_RADII(Y,Z1,Z2) looks for r0 > 0 such that
%       p(r0) = Y + Z1*r0 + Z2*r0^2/2 - r0 <= 0   and   q(r0) = Z1 + Z2*r0 < 1,
%   and checks both inequalities in interval arithmetic. Theorem 1 of the
%   paper then gives a unique zero of F within distance r0 of xbar in H.
%   With a global Z2 bound one may take r_* = r0 (proof of Theorem 2).
%   certificate.success is true only if both interval tests pass.
% Verified Newton--Kantorovich inequalities. No proof claim on failed checks.
values=[sup(Y),sup(Z1),sup(Z2)];
assert(all(isfinite(values)) && all(values>=0),'Invalid bound.');
% Work with the upper bounds only (p and q are increasing in Y, Z1, Z2).
Y=intval(values(1)); Z1=intval(values(2)); Z2=intval(values(3));
certificate=struct('success',false,'Y',Y,'Z1',Z1,'Z2',Z2, ...
    'r0',NaN,'radiiPolynomial',intval(inf),'contraction',intval(inf));
if sup(Z1)>=1, return; end
% p(r) = 0 has real roots iff its discriminant D = (1-Z1)^2 - 2*Y*Z2 > 0.
D=(1-Z1)^2-2*Y*Z2;
if inf(D)<=0, return; end
% Smaller root of p: r = (1-Z1-sqrt(D))/Z2 = 2Y/(1-Z1+sqrt(D)).
% Stable smaller root avoids subtracting two nearly equal numbers.
root=sup(2*Y/(1-Z1+sqrt(D)));
if root==0, root=realmin; end
% Try radii slightly above the smaller root; accept the first one for which
% both inequalities are verified in interval arithmetic.
for margin=[1e-6,1e-4,1e-2,0.1]
    r=root*(1+margin);
    if ~isfinite(r) || r<=0, continue; end
    p=Y+Z1*intval(r)+Z2*intval(r)^2/2-intval(r);
    q=Z1+Z2*intval(r);
    if sup(p)<=0 && sup(q)<1
        certificate.success=true; certificate.r0=r;
        certificate.radiiPolynomial=p; certificate.contraction=q; return
    end
end
end
