function [vmin,fmin]= SecondEigenvector(W,normalized)
% Computes the second eigenvector of the standard graph Laplacian.
%
% Usage: [vmin,fmin]= SecondEigenvector(W,normalized)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    options.disp=0;
  
    D=DegreeMatrix(W);
	Lp=D-W;
 
    if(normalized)
        [eigVecs,eigVals] =eigs(Lp,D,2,'SA',options);
    else
        [eigVecs,eigVals] =eigs(Lp,2,'SA',options);
    end
    
	vmin=eigVecs(:,2);
    vmin=vmin/norm(vmin);
    fmin=eigVals(2,2);
    
end