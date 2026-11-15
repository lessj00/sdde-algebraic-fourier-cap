function B = bell_table(x,smax)
%BELL_TABLE Partial (exponential) Bell polynomials B_{s,j}.
%   B = BELL_TABLE(x,smax)
%   x : cell array, x{q} = q-th argument x_q (scalar or column vector;
%       double or intval), q = 1,...,smax.
%   B{s+1,j+1} = B_{s,j}(x_1,...,x_{s-j+1}),  0 <= j <= s <= smax,
%   computed with the recurrence of Sect. 3 of the paper,
%       B_{s,j} = sum_{q=1}^{s-j+1} binom(s-1,q-1) x_q B_{s-q,j-1},
%       B_{0,0} = 1,  B_{s,0} = 0 (s >= 1).
%   These polynomials appear in Faa di Bruno's formula (6) (Lemma 3).
%   All coefficients are nonnegative integers, so B_{s,j} is nondecreasing
%   in each argument on [0,inf); this gives the bound (7) of Lemma 3.
%   (In the code below the summation index q of the paper is called i.)
z = 0*x{1};
B = cell(smax+1,smax+1);
B{1,1} = z + 1;                      % B_{0,0} = 1
for s = 1:smax
    B{s+1,1} = z;                    % B_{s,0} = 0 for s >= 1
    for j = 1:s
        acc = z;
        for i = 1:s-j+1
            % binom(s-1,i-1) * x_i * B_{s-i,j-1}
            acc = acc + nchoosek(s-1,i-1)*(x{i}.*B{s-i+1,j});
        end
        B{s+1,j+1} = acc;
    end
end
end
