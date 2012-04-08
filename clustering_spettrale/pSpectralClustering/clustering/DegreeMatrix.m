function D = DegreeMatrix(W)
% Compute the Degree Matrix d_ij = sum(w_ij), i=j , 0 , else
%
% Usage: D = DegreeMatrix(W)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

	D = sparse(diag(sum(W)));
    
end