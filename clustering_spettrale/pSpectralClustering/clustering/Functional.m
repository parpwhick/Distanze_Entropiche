function result = Functional(u,W,p,normalized,deg)
% Computes the value of the functional F introduced in 
% "Spectral Clustering based on the graph p-Laplacian".
%
% Usage: result = Functional(u,W,p,normalized,deg)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    W3=getSparseDerivativeMatrix(u,W);
    W3=computeAbsPower(W3,p);
    W4=sparse(W.*W3);
    
    enum = full(sum(sum(W4)));

    if (normalized)
        denom=pNormPowDegParallel(u,p,deg);
    else
        denom=pNormPowParallel(u,p);
    end
    
    result = enum/denom;
    
end

