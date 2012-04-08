function pnorm = pNormPowParallel(u,p)
% Computes the p-norm of u to the power p.
%
% Usage: pnorm = pNormPowParallel(u,p)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

	pnorm= sum(abs(u).^p);
	
end