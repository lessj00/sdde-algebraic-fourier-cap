function results = selftest_intlab()
%SELFTEST_INTLAB Small real-INTLAB regression tests of the proof code.
%   results = SELFTEST_INTLAB() runs quick checks of the building blocks on
%   small examples where the answer is known independently:
%     - complex error discs (Proposition 1) and the aliasing bound (Lemma 2),
%     - Faa di Bruno and midpoint bounds on sin(t) (Lemmas 3 and 4),
%     - the C^0 operator-norm bound ([4, Sect. 2.6.6]) against a direct kernel,
%     - residual, Jacobian and the full Y, Z1, Z2 pipeline on the exact
%       solution u = sin(t) of the linear case eps1 = 0,
%     - nonlinear residual/Jacobian enclosures (Prop. 1) and the midpoint and
%       pointwise residual bounds (Lemmas 4 and 5) against independent
%       Taylor/Fourier enclosures,
%     - Fourier multiplication indices, Kolmogorov constants theta_N, and the
%       Newton-Kantorovich radii test (Theorem 1).
%   This is not a certification of the saved orbits (see run_proof).
environment=check_intlab();
passed={};
% Complex error enclosure must include errors in every direction.
z=complex_error_disc(0.25);
assert(sup(abs(z-intval(0.25i)))>=0.25); % enclosure sanity
assert(abs(mid(z)-0.25i)<=rad(z));
passed{end+1}='complex Fourier error disc';
% Aliasing estimate dominates a finite exact partial sum.
for s=[2,3,6]
    E=alias_bound(3,32,s); partial=intval(0);
    for j=[-100:-1,1:100], partial=partial+1/abs(intval(3+32*j))^s; end
    assert(inf(E)>sup(partial));
end
passed{end+1}='algebraic alias bound';
% Exact Fourier mode sin(t): midpoint bound encloses mean |sin|=2/pi.
a=[1i/2;0;-1i/2];
[Kv,Kw,Bs]=Cs_crude_SDDE(a,0,6);
assert(all(sup(Kv)>=1) && all(sup(Kw)>=1));
kg=Cs_grid_SDDE(a,0,4,64);
assert(all(sup(kg)>=sup(2/intval('pi'))));
passed{end+1}='Faà di Bruno and midpoint bounds';
% Compare blocked FFT kernel with direct verified exponential matrices.
N=2; M=16; n=(-N:N); t=2*intval('pi')*intval((0:M-1)')/M;
B=intval(reshape(sin(1:36),6,6)/10);
E=exp(1i*t*intval(n));
direct=max(sup(sum(abs(E*B(2:end,2:end)*E'/M),2)));
% FFT kernel alone enclosed directly: the final norm includes a factor >=1.
[f,normDetails]=comp_mat_norm(B,M);
assert(abs(normDetails.sampledKernel-direct)<1e-10,'Blocked FFT normalization disagrees with the direct kernel.');
fejerUpper=intval(M)/intval(M-2*N)/cos(intval('pi')*N/M);
assert(normDetails.samplingFactor<sup(fejerUpper));
assert(sup(f)>=direct && sup(f)>=sup(norm(B(1,:),1)));
passed{end+1}='blocked Fourier operator bound';
% Exact linear forced solution sin(t), epsilon_1=0.
x=[0;a]; F=F_SDDE_rig_Cs(x,intval([0;2]),Kv,32);
assert(all(inf(real(F))<=0) && all(sup(real(F))>=0));
assert(all(abs(imag(mid(F)))<=rad(F)));
J=DF_SDDE_rig_Cs(x,intval([0;2]),Kw,Bs,32,64);
assert(all(isfinite(sup(abs(J(:))))));
[r,Mg]=residual_norm_Cs(a,intval([0;2]),kg,64,8);
assert(isfinite(r) && r<1e-10 && inf(Mg)<=0 && sup(Mg)>=0);
passed{end+1}='linear exact solution residual and Jacobian';
% End-to-end interval bounds on an exact linear forced solution.
% This exercises the same bound pipeline as the large saved approximations.
N=32; x=[0;zeros(N-1,1);a;zeros(N-1,1)];
[Kv2,Kw2,Bs2]=Cs_crude_SDDE(a,2,6);
kg2=Cs_grid_SDDE(a,2,4,64);
[Y,Z1,Z2,details]=computation_bounds_Cs(x,[0;2],128,256,512, ...
    Kv2,Kw2,Bs2,kg2,[],struct('a',a,'Mp',64,'K',8));
linearCertificate=verify_radii(Y,Z1,Z2);
assert(linearCertificate.success,'Linear full-pipeline certificate failed.');
assert(all(isfinite([sup(Y),sup(Z1),sup(Z2),details.thetaN])));
passed{end+1}='full linear-example interval certificate';
% Independent nonlinear reference: Taylor polynomials of exp(-ik*eps2*u).
% This uses coefficient convolution and a uniform analytic remainder, not
% a sampled DFT or the Bell/aliasing implementation being tested.
nonlinear_reference_test(32);
nonlinear_reference_test(8); % M=4N: positive Nyquist alias wraps to index 1.
passed{end+1}='nonlinear residual/Jacobian versus independent Taylor enclosure';
nonlinear_midpoint_residual_test();
passed{end+1}='nonlinear midpoint and pointwise residual versus independent Fourier bounds';
% Preserve signed index orientation and interval entries in the Toeplitz map.
W=midrad(complex((1:9)',(-4:4)'),(1:9)'/1000);
T=fourier_multiplication_matrix(W,2);
for n=-2:2
    for k=-2:2
        assert(isequal(mid(T(n+3,k+3)),mid(W(n-k+5))));
        assert(isequal(rad(T(n+3,k+3)),rad(W(n-k+5))));
    end
end
% For both production defaults, wrapped selection is exactly the old
% contiguous selection; the boundary fix changes no production arithmetic.
for proofN=[820,1320]
    for grid=[2^14,2^15] % numerical inverse grid and rigorous coefficient grid
        assert(isequal(mod((-2*proofN:2*proofN)'+grid/2,grid)+1, ...
            (grid/2-2*proofN+1:grid/2+2*proofN+1)'));
    end
end
passed{end+1}='signed Fourier multiplication indices and interval radii';
% theta_1 is known explicitly from its root enclosure formula; also check
% the large N used for the two proof cases against a separate safe bound.
for Ntheta=[1,2,820,1320]
    th=kolmogorov_theta_N(Ntheta);
    assert(inf(th)>0 && sup(th)<1);
    if Ntheta>=2, assert(sup(th)<sup(sqrt(2/intval(Ntheta)))); end
end
passed{end+1}='Kolmogorov enclosures at small and production truncations';
% Successful and failed radii tests must be distinguished.
c=verify_radii(1e-10,0.2,100);
assert(c.success && sup(c.radiiPolynomial)<=0 && sup(c.contraction)<1);
cFail=verify_radii(1,1.1,10); assert(~cFail.success);
passed{end+1}='verified radii success/failure';
results=struct('environment',environment,'passed',{passed},'linearCertificate',linearCertificate,'success',true);
fprintf('Passed %d real-INTLAB checks. This is not an orbit certificate.\n',length(passed));
end

function [referenceF,referenceJ] = nonlinear_reference_test(M)
% ubar=sin(t), eps2=1/4 genuinely changes the argument with the state.
N=2; Q=18; L=Q+2*N+1; modes=(-L:L)';
a=intval([1i/2;0;-1i/2]); parameters=intval([1/2;1/4]);
x=[0;0;a;0];
[Kv,Kw,Bs]=Cs_crude_SDDE(a,parameters(2),8);
actualF=F_SDDE_rig_Cs(x,parameters,Kv,M);
actualJ=DF_SDDE_rig_Cs(x,parameters,Kw,Bs,M,64);
assert(all(all(isfinite(DF_SDDE(mid(x),[1/2;1/4],M)))));
V=intval(zeros(2*L+1,1)); W=V; RV=intval(0); RW=RV;
for m=[-1,1]
    [e,R]=taylor_column(m,parameters(2),Q,L);
    V=V+a(m+2)*e; W=W+1i*m*a(m+2)*e;
    RV=RV+abs(a(m+2))*R; RW=RW+abs(m*a(m+2))*R;
end
forcing=intval(zeros(2*L+1,1));
forcing(L)=intval(1)/2+1i*intval(1)/2;
forcing(L+2)=intval(1)/2-1i*intval(1)/2;
u=intval(zeros(2*L+1,1)); u(L:L+2)=a;
G=forcing-u-parameters(1)*V; RG=abs(parameters(1))*RV;
weights=intval(zeros(2*L+1,1)); nz=modes~=0;
weights(nz)=(-1).^(modes(nz)+1)./(1i*intval(modes(nz)));
referenceF=intval(zeros(2*N+2,1));
referenceF(1)=G(L+1)+complex_error_disc(RG);
for n=-N:N
    if n==0
        referenceF(n+N+2)=-sum(weights.*G.*intval(abs(modes)<=N))+complex_error_disc(4*intval('pi')*RG);
    else
        referenceF(n+N+2)=u(n+L+1)-G(n+L+1)/(1i*intval(n))+complex_error_disc(RG/abs(n));
    end
end
referenceJ=intval(zeros(2*N+2)); referenceJ(N+2,1)=-1;
for k=-N:N
    [e,Re]=taylor_column(k,parameters(2),Q,L);
    shiftedW=intval(zeros(2*L+1,1));
    valid=find(abs(modes-k)<=L);
    shiftedW(valid)=W(modes(valid)-k+L+1);
    h=parameters(1)*parameters(2)*shiftedW-parameters(1)*e;
    h(k+L+1)=h(k+L+1)-1;
    Rh=abs(parameters(1)*parameters(2))*RW+abs(parameters(1))*Re;
    referenceJ(1,k+N+2)=h(L+1)+complex_error_disc(Rh);
    for n=-N:N
        if n==0
            referenceJ(n+N+2,k+N+2)=(k==0)-sum(weights.*h)+complex_error_disc(4*intval('pi')*Rh);
        else
            referenceJ(n+N+2,k+N+2)=(n==k)-h(n+L+1)/(1i*intval(n))+complex_error_disc(Rh/abs(n));
        end
    end
end
% The independent reference discs are much narrower than the deliberately
% coarse FFT alias bounds. Require containment, not merely overlap.
neededF=sup(abs(intval(mid(actualF))-intval(mid(referenceF))))+rad(referenceF);
assert(all(neededF<=rad(actualF)), 'Independent nonlinear F enclosure not contained: excess %.3e', max(neededF-rad(actualF)));
neededJ=sup(abs(intval(mid(actualJ))-intval(mid(referenceJ))))+rad(referenceJ);
assert(all(all(neededJ<=rad(actualJ))), 'Independent nonlinear Jacobian enclosure not contained: excess %.3e',max(max(neededJ-rad(actualJ))));
end

function [coefficients,remainder]=taylor_column(k,eps2,Q,L)
% e^{ik(t-1-eps2*sin(t))}; truncate exp(-ik*eps2*sin(t)) after Q.
a=intval([1i/2;0;-1i/2]);
series=intval(zeros(2*Q+1,1)); series(Q+1)=1; term=intval(1);
for j=1:Q
    next=intval(zeros(length(term)+2,1));
    for r=1:3
        next(r:r+length(term)-1)=next(r:r+length(term)-1)+a(r)*term;
    end
    term=(-1i*intval(k)*eps2/j)*next;
    series(Q+1-j:Q+1+j)=series(Q+1-j:Q+1+j)+term;
end
coefficients=intval(zeros(2*L+1,1));
coefficients(L+1-Q+k:L+1+Q+k)=exp(-1i*intval(k))*series;
z=abs(intval(k)*eps2); remainder=exp(z)*z^(Q+1);
for j=1:Q+1, remainder=remainder/j; end
end


function nonlinear_midpoint_residual_test()
% Independently bound v=sin(t-1-eps2*sin(t)) by coefficient convolutions.
% A Cauchy estimate on a unit complex circle bounds derivatives of the
% exponential-series remainder. This reference does not use Bell polynomials.
a=intval([1i/2;0;-1i/2]); eps2=intval(1)/4;
Q=18; L=Q+3; n=intval((-L:L)');
P=intval(zeros(2*L+1,1)); R=intval(0);
for m=[-1,1]
    [col,~]=taylor_column(m,eps2,Q,L);
    P=P+a(m+2)*col;
    z=abs(m*eps2)*cosh(intval(1));
    term=exp(abs(intval(m)))*exp(z)*z^(Q+1);
    for j=1:Q+1, term=term/j; end
    R=R+abs(a(m+2))*term;
end
Kg=Cs_grid_SDDE(a,eps2,4,64);
M=4096; h=2*intval('pi')/M;
ind=mod((-L:L)',M)+1;
for s=1:4
    coeff=(1i*n).^s.*P;
    padded=intval(zeros(M,1)); padded(ind)=coeff.*exp(1i*n*h/2);
    values=verifyfft(padded,-1)*M;
    % Uniform Fourier derivative bound controls the quadrature remainder.
    reference=sum(abs(values))/M + intval('pi')/(2*M)*sum(abs((1i*n).*coeff));
    reference=reference+factorial(s)*R;
    assert(sup(reference)<sup(Kg(s)), ...
        'Nonlinear midpoint bound is below independent Fourier upper bound.');
end
% For ubar=sin(t), r=eps1*(v-mean(v)). The Taylor remainder in this
% mean-adjusted expression is at most 2*R (use the unit-strip bound on R).
coeff=P; coeff(L+1)=0;
padded=intval(zeros(M,1)); padded(ind)=coeff;
values=verifyfft(padded,-1)*M;
referenceR=(max(abs(values))/cos(intval('pi')*L/M)+2*R)/2;
[rnorm,~,~]=residual_norm_Cs(a,intval([1/2;1/4]),Kg,256,12);
assert(sup(referenceR)<rnorm, ...
    'Pointwise residual bound is below independent Fourier upper bound.');
end
