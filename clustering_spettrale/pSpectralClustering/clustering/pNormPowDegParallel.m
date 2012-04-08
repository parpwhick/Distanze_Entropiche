function pnormdeg = pNormPowDegParallel(u,p,deg)
% Computes the weighted p-norm of u to the power p.
%
% Usage: pnormdeg = pNormPowDegParallel(u,p,deg)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

	pnormdeg = deg*abs(u).^p;
 	
end