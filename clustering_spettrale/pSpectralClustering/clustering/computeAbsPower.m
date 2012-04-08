function W4=computeAbsPower(W3,p)
% Applies the power function to the absolute value of every component of a sparse matrix.
%
% Usage: W1 = computeAbsPower(W,p)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de
	
    [width,height]=size(W3);
    [i,j,v]=find(W3);
    v1=abs(v).^p;
    W4=sparse(i,j,v1,width,height);
    
end