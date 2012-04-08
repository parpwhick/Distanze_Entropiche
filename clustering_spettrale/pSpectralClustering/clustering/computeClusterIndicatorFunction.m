function u= computeClusterIndicatorFunction(umin,threshold)
% Computes a cluster indicator function from a vector u
%
% Usage: u= computeClusterIndicatorFunction(umin,threshold)
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de
    
    u= (umin>threshold);
	
end
