function r = normrnd(mu,sigma,m,n)
%NORMRND Random matrices from normal distribution.
%   R = NORMRND(MU,SIGMA) returns a matrix of random numbers chosen   
%   from the normal distribution with parameters MU and SIGMA.
%
%   The size of R is the common size of MU and SIGMA if both are matrices.
%   If either parameter is a scalar, the size of R is the size of the other
%   parameter. Alternatively, R = NORMRND(MU,SIGMA,M,N) returns an M by N  
%   matrix.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.11 $  $Date: 2002/03/31 22:26:56 $

r = randn(m,n) .* sigma + mu;

% Fill in elements corresponding to illegal parameter values
if prod(size(sigma)) > 1
    r(sigma < 0) = NaN;
elseif sigma < 0
    r(:) = NaN;
end
