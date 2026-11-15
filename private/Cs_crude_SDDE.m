function [Kv,Kw,Bs] = Cs_crude_SDDE(a,eps2,smax)
%CS_CRUDE_SDDE A priori derivative bounds of Lemma 3 (equation (7)).
%   [Kv,Kw,Bs] = CS_CRUDE_SDDE(a,eps2,smax)
%   With theta(t) = t - 1 - eps2*ubar(t) and the compositions (4) of the paper,
%       v = ubar o theta,   w = ubar' o theta,   e_k = exp(i*k*theta),
%   this returns, for s = 0..smax,
%       Kv(s+1) = calK_s(v) >= ||v^{(s)}||_{C^0} >= K_s(v),
%       Kw(s+1) = calK_s(w) >= ||w^{(s)}||_{C^0} >= K_s(w),
%       Bs(s+1,j+1) = B_{s,j}(alpha_1,...,alpha_{s-j+1}),
%   where (Lemma 3)
%       mu_q    = sum_n |n|^q |Ubar_n|,
%       alpha_1 = 1 + |eps2|*mu_1,   alpha_q = |eps2|*mu_q  (q >= 2),
%   so that |theta^{(q)}(t)| <= alpha_q by (5), and
%       calK_s(Phi o theta) = sum_{j=0}^{s} (sum_m |m|^j |c_m|) B_{s,j}(alpha).
%   For e_k this gives calK_s(e_k) = sum_{j=1}^s |k|^j Bs(s+1,j+1) (s >= 1),
%   which is assembled in DF_SDDE_rig_Cs.
%   These constants control the aliasing errors in Proposition 1(a),(c),(d),(e).
%   a = (Ubar_n)_{n=-N_u}^{N_u}, Fourier coefficients of ubar. Requires INTLAB.
N  = (length(a)-1)/2;
n  = abs(intval((-N:N)'));
A  = abs(intval(a));
% nj(j+1) = mu_j = sum |n|^j |Ubar_n|, j = 0..smax+1 (with |0|^0 = 1).
nj = intval(zeros(smax+2,1));
for j = 0:smax+1
    nj(j+1) = sum(n.^j.*A);
end
ae2 = abs(intval(eps2));
% Arguments of the Bell polynomials: x{q} = alpha_q (bounds (5) on theta^{(q)}).
x = cell(max(smax,1),1);
x{1} = 1 + ae2*nj(2);
for i = 2:smax
    x{i} = ae2*nj(i+1);
end
B  = bell_table(x,smax);
Kv = intval(zeros(smax+1,1)); Kw = intval(zeros(smax+1,1));
Bs = intval(zeros(smax+1,smax+1));
% s = 0: ||v||_{C^0} <= sum|Ubar_n| and ||w||_{C^0} <= sum|n||Ubar_n|.
Kv(1) = nj(1); Kw(1) = nj(2); Bs(1,1) = 1;
for s = 1:smax
    for j = 1:s
        % Phi = ubar : sum_m |m|^j |c_m| = mu_j.
        Kv(s+1) = Kv(s+1) + nj(j+1)*B{s+1,j+1};
        % Phi = ubar': c_m = i*m*Ubar_m, so sum_m |m|^j |c_m| = mu_{j+1}.
        Kw(s+1) = Kw(s+1) + nj(j+2)*B{s+1,j+1};
        Bs(s+1,j+1) = B{s+1,j+1};
    end
end
end
