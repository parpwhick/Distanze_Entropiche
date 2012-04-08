function [cut, clustersize_max,imbalance] = computeUnbalancedCutWeight(clusters,W)
% Computes the value of the cutweight, the maximal cluster size and the imbalance. 
%
% Usage: [cut, clustersize_max,imbalance] = computeUnbalancedCutWeight(clusters,W)
%
% The imbalance i is defined as i = clustersize_max-clustersize_optimal -1.
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    labels = unique(clusters);
    numClusters=size(labels,1);
    
    cut=0;
    clustersize_max=0;
    
    for k=1:numClusters
        label=labels(k);
        W2 = W(clusters==label,clusters~=label);
        cut=cut + full(sum(sum(W2)));
         
        clustersize=sum(clusters==label);
        if clustersize>clustersize_max
            clustersize_max=clustersize;
        end
    end
    cut=cut/2;
    
    clustersize_opt= ceil(size(W,1)/numClusters);
    
    imbalance=clustersize_max/ clustersize_opt -1;

end