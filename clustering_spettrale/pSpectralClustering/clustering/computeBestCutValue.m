function [cut,cheeger,clusters,threshold] = computeBestCutValue(u,W,normalized)
% Creates a cluster indicator function with minimal Cheeger cut value.
%
% Usage: [cut,cheeger,clusters,threshold] = computeBestCutValue(u,W,normalized)
% 
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

	cut=inf;
	cheeger=inf;
	for k=1:size(u,1)
		threshold1=u(k);
		if threshold1~=min(u) && threshold1~= max(u)
			clusters1= computeClusterIndicatorFunction(u,threshold1);
			[cutpart1,cutpart2]=computeCutValue(clusters1,W,normalized);
			
			cut1=cutpart1+cutpart2;
			cheeger1=max(cutpart1,cutpart2);
			
			if cheeger1<cheeger 
				cut=cut1;
				cheeger=cheeger1;
				clusters=clusters1;
				threshold=threshold1;
			end
		end
	end
	
end