function environment = check_intlab()
%CHECK_INTLAB Make sure that genuine INTLAB interval arithmetic is active.
%   environment = CHECK_INTLAB() aborts unless MATLAB (not Octave) runs with
%   a real INTLAB installation whose outward rounding works, and returns the
%   MATLAB/INTLAB versions and paths, which are stored in every report.
%   All bounds of the paper (Sect. 5) are computed in interval arithmetic, so
%   a double-precision stand-in for intval would invalidate the proof.
% Require actual INTLAB, not an Octave/double compatibility shim.
assert(exist('OCTAVE_VERSION','builtin')==0,'Use MATLAB with real INTLAB.');
assert(exist('intval','file')~=0,'Initialize INTLAB with startintlab first.');
x=intval(1);
assert(isa(x,'intval'),'intval is shadowed by a non-interval replacement.');
assert(~intvalinit('RoundingFailure'),'INTLAB reports a rounding failure.');
% 1 + 2^-53 is not a double: a correct enclosure must have sup > 1.
x=intval(1)+intval(2^-53);
assert(inf(x)<=1 && sup(x)>1,'Outward rounding probe failed.');
% The verified FFT is used for all DFT enclosures (Proposition 1).
y=verifyfft(intval([1;0;0;0]));
assert(isa(y,'intval') && all(inf(real(y))<=1) && all(sup(real(y))>=1));
environment=struct('matlab',version,'computer',computer, ...
    'intlab',intvalinit('version',0),'intvalPath',which('intval'), ...
    'verifyfftPath',which('verifyfft'),'timestamp',datestr(now,30));
end
