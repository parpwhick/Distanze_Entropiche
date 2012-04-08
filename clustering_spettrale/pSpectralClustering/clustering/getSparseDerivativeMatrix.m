function W3u=getSparseDerivativeMatrix(u,W)
% Returns the upper triangular matrix only.
%
% Usage: W3u=getSparseDerivativeMatrix(u,W)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    Wu=triu(W);
    [i,j]=find(Wu);
    s=size(W,1);
    W3u=sparse(i,j,u(j)-u(i),s,s);

end