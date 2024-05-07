function X = fixq( s, bit, rmode, lmode )
%Investigate quantizing with Q# of bits. and see signal loss 
%fixq      simulated fixed-point arithmetic
%-------
%   Usage:  X = fixq( s, bit, rmode, lmode )
%
%       returns the input signal S reduced to a word-length
%       of BIT bits and limited to the range [-1,1). The type of
%       word-length reduction and limitation may be chosen with
%       RMODE:  'round'         rounding to nearest level
%               'trunc'         2's complement truncation
%               'magn'          magnitude truncation
%       LMODE:  'sat'           saturation limiter
%               'overfl'        2's complement overflow
%               'none'          no limiter

if nargin ~= 4;
        error('usage: fixq( s, bit, rmode, lmode ).');
end;
if (bit <= 0 | abs(rem(bit,1)) > eps)
        error('wordlength must be positive integer.');
end;

Plus1 = 2^(bit-1);

X = s * Plus1;

if     strcmp(rmode, 'round');          X = round(X);
elseif strcmp(rmode, 'trunc');          X = floor(X);
elseif strcmp(rmode, 'ceil');           X = ceil(X);
elseif strcmp(rmode, 'magn');           X = fix(X);
else   error('unknown wordlength reduction spec.');
end;

if     strcmp(lmode, 'sat');
        X = min(Plus1 - 1,X);
        X = max(-Plus1,X);
elseif strcmp(lmode, 'overfl');
        X = X + Plus1 * ( 1 - 2*floor((min(min(X),0))/2/Plus1) );
        X = rem(X,2*Plus1) - Plus1;
elseif strcmp(lmode, 'none');           %  limiter switched off
else    error('unknown limiter specifier');
end;

X = X / Plus1;