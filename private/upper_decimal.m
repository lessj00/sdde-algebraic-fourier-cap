function s = upper_decimal(x)
%UPPER_DECIMAL Decimal string that is a verified upper bound of x.
%   s = UPPER_DECIMAL(x) returns a 17-digit decimal string s with s >= sup(x),
%   checked by converting s back to an interval. It is used to print the
%   values of Y, Z1, Z2 and r0 in the log files (Table 2 of the paper).
% A decimal string verified to lie above the exact binary upper bound x.
x=sup(intval(x)); assert(isscalar(x) && isfinite(x) && x>=0);
y=x;
for j=1:10
    s=sprintf('%.17g',y);
    if inf(intval(s))>=x, return; end
    % Increase by one unit in the last place and try again.
    y=y+eps(y);
end
error('Could not format a certified decimal upper bound.');
end
