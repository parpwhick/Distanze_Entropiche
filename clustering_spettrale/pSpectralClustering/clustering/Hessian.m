function hessian = Hessian(u,W,p,normalized,deg)
% Computes the hessian of the Functional F
%
% hessian = Hessian(u,W,p,normalized,deg)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    if (normalized)
        denom=pNormPowDegParallel(u,p,deg);
    else
        denom=pNormPowParallel(u,p);
    end

	W3u=getSparseDerivativeMatrix(u,sparse(W));
    W4u=computeAbsPower(W3u,p-2);
    W5u= sparse(W.*W4u);
    W5= W5u+W5u';
   
    D=diag(sum(W5));
   
    hessian = p*(p-1)/denom * (D-W5);

end