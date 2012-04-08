function gradient = Gradient(u,W,p,normalized,deg,functional)
% Computes the gradient of the Functional F
%
% Usage: gradient = Gradient(u,W,p,normalized,deg,functional)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    W3u = getSparseDerivativeMatrix(u,W);
    W4u=computeAbsPower(W3u,p-1);
    W4ub=W4u.*sign(W3u);
    W5u= sparse(W.*W4ub);
    W5= W5u-W5u';
    
    left=sum(W5)';

    u2=abs(u).^(p-1).*sign(u);
    
    if (normalized)
        denom=pNormPowDegParallel(u,p,deg);
        u2=u2.*deg';
    else
        denom = pNormPowParallel(u,p);
    end
    
    gradient= p/denom * (left-functional*u2);
    
end