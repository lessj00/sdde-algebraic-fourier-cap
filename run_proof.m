function report = run_proof(Nu, options)
%RUN_PROOF Computer-assisted proof of Theorem 2 of the paper.
%
%   J.-P. Lessard, "Algebraic Fourier Estimates for Computer-Assisted Proofs
%   in State-Dependent Delay Equations" (AMMCS 2026 proceedings).
%
%   Usage (after initializing INTLAB with startintlab):
%       report = run_proof(800);       % or run_proof(1300)
%       report = run_proof(800,struct('OutputDir','my_results'));
%
%   The equation is (3) of the paper,
%       u'(t) = -u(t) - eps1*u(t-1-eps2*u(t)) + cos(t) + sin(t),
%   with (eps1,eps2) = (0.5,2). The zero-finding map is F(beta,u) of Sect. 2,
%       F(beta,u) = ( M g(u) , u - beta - L1 g(u) ),
%   and the proof checks the hypotheses of the Newton-Kantorovich theorem
%   (Theorem 1) at xbar = (betabar, ubar):
%       Y + Z1*r0 + Z2*r0^2/2 <= r0   and   Z1 + Z2*r0 < 1.
%
%   Steps (paper references in brackets):
%     1. Read betabar and the Fourier coefficients of ubar from
%        data/ubar_<Nu>.txt and enforce conjugate symmetry exactly.
%     2. A priori bounds calK_s(v), calK_s(w) and the Bell table (Lemma 3,
%        equation (7)), used for the aliasing errors of Proposition 1.
%     3. Sharp bounds tildeK_s(v) by the midpoint enclosure (Lemma 4), used
%        for the Fourier tails in (8) and in Lemma 5.
%     4. Y, Z1, Z2 in computation_bounds_Cs (equations (8), (9), Lemma 5,
%        Proposition 1, and Sect. 2.6 of the companion preprint [4]).
%     5. Verified radius r0 in verify_radii (Theorem 1).
%   A log (.log) and an interval report (.mat) are written to OutputDir.
%
%   Default parameters are those of Sect. 5 of the paper:
%     N    = Nu+20        Galerkin truncation (Table 2)
%     N_FFT= 2^12 / 2^13  grid of the C^0 operator norms ([4, Sect. 2.6.6])
%     M    = 2^15         DFT grid for V_n, W_n, (L3)_{n,k} (Prop. 1(a),(c),(d))
%     N0   = 2^16         DFT grid for the zeroth row (Prop. 1(e))
%     smax_c = 30         largest s in the a priori bounds (Lemma 3)
%     smax_g = 12, M_g = 2^15   largest s and number of cells in Lemma 4
%     Mp   = 2^15 (= M'), K = 4000   sampling grid and cutoff in Lemma 5
%
% Approximate coefficients are exact stored doubles, not interval solutions.
% Any settings change is recorded; success is determined only by interval tests.
if nargin<1, Nu=800; end
if nargin<2, options=struct; end
root=fileparts(mfilename('fullpath'));
% Abort unless real INTLAB with outward rounding is active.
environment=check_intlab();
assert(ismember(Nu,[800,1300]),'Supplied coefficient orders are 800 and 1300.');
% N_FFT is the smallest power of two with N_FFT >= 4N (2^12 for N=820,
% 2^13 for N=1320), as stated in Sect. 5.
N=Nu+20; nf=2^nextpow2(4*N);
cfg=struct('N',N,'N_FFT',nf,'M',2^15,'N0',2^16, ...
    'smax_c',30,'smax_g',12,'M_g',2^15,'Mp',2^15,'K',4000, ...
    'N_Z1',N,'OutputDir',fullfile(root,'results'));
% Optional user overrides (all recorded in report.settings).
keys=fieldnames(options);
for j=1:length(keys)
    assert(isfield(cfg,keys{j}),'Unknown setting: %s',keys{j});
    cfg.(keys{j})=options.(keys{j});
end
% Hypotheses on the grids used in the paper:
%   N >= N_u+1 (Prop. 1), M > 2N and N0 > 2N (Prop. 1), 2N_u < M_g (Lemma 4),
%   N_u <= K < M'/2 (Lemma 5). The code uses the stronger 4N <= M, which
%   allows the coefficients W_n, |n| <= 2N, to be read from one DFT.
numbers=[cfg.N,cfg.N_FFT,cfg.M,cfg.N0,cfg.smax_c,cfg.smax_g,cfg.M_g,cfg.Mp,cfg.K,cfg.N_Z1];
assert(all(isfinite(numbers)) && all(numbers==fix(numbers)) && all(numbers>=1));
assert(cfg.N>Nu && 4*cfg.N<=cfg.N_FFT && 4*cfg.N<=cfg.M ...
    && 2*cfg.N<cfg.N0 && 2*Nu<cfg.M_g && Nu<=cfg.K && 2*cfg.K<cfg.Mp);
assert(cfg.N_Z1<=cfg.N && cfg.smax_c>=2 && cfg.smax_c<=30 && cfg.smax_g>=3 && cfg.smax_g<=30);
% Lemma 5: the phases n*t_k on the grid t_k = 2*pi*k/M' are reduced
% modulo 2*pi in exact integer arithmetic, which needs n*M' < 2^53.
assert(Nu*cfg.Mp+cfg.Mp/2<2^53,'Residual phase reduction must use exact integers.');
for g=[cfg.N_FFT,cfg.M,cfg.N0,cfg.M_g,cfg.Mp]
    assert(g==2^nextpow2(g),'FFT grids must be powers of two.');
end
if ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end
% Unique file name: validation_<Nu>_<date-time>_<random tag>.
[~,tag]=fileparts(tempname(cfg.OutputDir));
base=fullfile(cfg.OutputDir,sprintf('validation_%d_%s_%s',Nu,datestr(now,30),tag));
diary([base '.log']); cleanup=onCleanup(@()diary('off')); 
report=struct('status','running','environment',environment,'settings',cfg,'Nu',Nu);
% A report records the exact implementation used, not just a folder name.
sourceFiles=[dir(fullfile(root,'*.m'));dir(fullfile(root,'private','*.m'))];
report.source=struct('name',{},'text',{});
for j=1:length(sourceFiles)
    sourcePath=fullfile(sourceFiles(j).folder,sourceFiles(j).name);
    report.source(j).name=erase(sourcePath,[root filesep]);
    report.source(j).text=fileread(sourcePath);
end
clockStart=tic;
try
    % ---- Step 1: approximate solution xbar = (betabar, ubar).
    % Data file: 2*Nu+2 rows (real part, imaginary part). Row 1 is betabar,
    % rows 2..end are Ubar_n, n = -Nu..Nu, so ubar(t) = sum Ubar_n e^{int}.
    file=fullfile(root,'data',sprintf('ubar_%d.txt',Nu));
    raw=load(file);
    assert(isequal(size(raw),[2*Nu+2,2]) && all(isfinite(raw(:))));
    x=complex(raw(:,1),raw(:,2));
    assert(imag(x(1))==0);
    a=x(2:end);
    % ubar must be real valued: Ubar_{-n} = conj(Ubar_n) (Sect. 2).
    % Explicit pair assignments preserve the real subspace exactly.
    for k=1:Nu
        z=(a(Nu+1+k)+conj(a(Nu+1-k)))/2;
        a(Nu+1+k)=z; a(Nu+1-k)=conj(z);
    end
    a(Nu+1)=real(a(Nu+1));
    assert(isequal(a,conj(flipud(a))));
    report.inputFile=file; report.inputData=raw;
    report.beta=x(1); report.coefficients=a; report.parameters=[0.5;2];
    % Auxiliary check, not used in the proof of Theorem 2: a rigorous lower
    % bound for 1 + 2*pi + eps2*ubar(t) >= 1 + 2*pi - 2*sum|Ubar_n|, valid
    % for all t. verify_paper_bounds checks that it stays positive after
    % the perturbation |eps2*(u-ubar)| <= 2*r0.
    report.shiftedLagLower=1+2*intval('pi')-2*sum(abs(intval(a)));
    save([base '.mat'],'report','-v7.3');
    fprintf('MATLAB %s; INTLAB %s; N_u=%d; N=%d\n',environment.matlab,environment.intlab,Nu,cfg.N);
    % ---- Step 2: a priori constants of Lemma 3, equation (7):
    % Kvc(s+1) >= calK_s(v), Kwc(s+1) >= calK_s(w), Bs(s+1,j+1) >= B_{s,j}(alpha).
    fprintf('Crude derivative bounds...\n');
    [Kvc,Kwc,Bs]=Cs_crude_SDDE(a,2,cfg.smax_c);
    % ---- Step 3: sharp constants of Lemma 4: Kvg(s) >= tildeK_s(v) >= K_s(v).
    fprintf('Midpoint derivative bounds...\n');
    [Kvg,KvgMid,KvgRemainder]=Cs_grid_SDDE(a,2,cfg.smax_g,cfg.M_g);
    report.KvgMid=KvgMid; report.KvgRemainder=KvgRemainder;
    report.Kvc=Kvc; report.Kwc=Kwc; report.Bs=Bs; report.Kvg=Kvg;
    save([base '.mat'],'report','-v7.3');
    % ---- Step 4: pad ubar with zeros to the Galerkin size N (modes -N..N).
    % The vector x is (betabar, Ubar_{-N}, ..., Ubar_N).
    x=[x(1);zeros(cfg.N-Nu,1);a;zeros(cfg.N-Nu,1)];
    % Options for the residual-based bounds (Lemma 5 and equation (9)).
    resopt=struct('a',a,'Mp',cfg.Mp,'K',cfg.K);
    [Y,Z1,Z2,info]=computation_bounds_Cs(x,[0.5;2],cfg.N_FFT,cfg.M,cfg.N0,Kvc,Kwc,Bs,Kvg,cfg.N_Z1,resopt);
    % ---- Step 5: Newton-Kantorovich inequalities of Theorem 1.
    report.bounds=info; report.certificate=verify_radii(Y,Z1,Z2);
    report.status='not_validated';
    if report.certificate.success
        report.status='validated';
        fprintf('VERIFIED: both Newton--Kantorovich inequalities hold.\n');
        % Decimal upper bounds, rounded upward (these are the values of Table 2).
        for key={'Y','Z1','Z2','r0'}
            fprintf('%s <= %s\n',key{1},upper_decimal(report.certificate.(key{1})));
        end
    else
        fprintf('NOT VALIDATED: interval inequalities did not close.\n');
    end
    report.elapsedSeconds=toc(clockStart);
    save([base '.mat'],'report','-v7.3');
    fprintf('Saved exact interval report: %s.mat\n',base);
catch problem
    % A failed run is saved as such, then the error is rethrown.
    report.status='error'; report.error=getReport(problem,'extended','hyperlinks','off');
    report.elapsedSeconds=toc(clockStart); save([base '.mat'],'report','-v7.3');
    fprintf('%s\n',report.error); rethrow(problem)
end
end
