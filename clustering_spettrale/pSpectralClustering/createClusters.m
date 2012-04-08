function [allClusters, cut,cheeger,cutPart1,cutPart2,threshold] =  createClusters(vmin,W,normalized,threshold_type,criterion)
% Transforms an eigenvector into a cluster indicator function by thresholding.
% 
% Usage:	[allClusters, cut,cheeger,cutPart1,cutPart2,threshold] 
%			= createClusters(vmin,W,normalized,threshold_type,criterion)
%
% vmin: The eigenvector.
% W: Weight matrix.
% normalized: True for Ncut/NCC, false for Rcut/RCC.
% threshold_type: 0: zero, 1: median, 2: mean, -1: best
% criterion: thresholding criterion. 1: Ratio/Normalized Cut 2: Ratio/Normalized Cheeger Cut
%
% allClusters: Obtained clustering after thresholding.
% cut: Value of the Normalized/Ratio Cut.
% cheeger: Value of the Normalized/Ratio Cheeger Cut.
% cutpart1,cutpart2: The two components of Ratio/Normalized Cut and Ratio/Normalized Cheeger Cut.
% threshold: The threshold used to obtain the partitioning.
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    if ~isempty(find(diag(W)~=0,1))
        error('Graph contains self loops. W has to have zero diagonal.');
	end
	
	if threshold_type>=0
            threshold= determineThreshold(threshold_type,vmin);
            allClusters= computeClusterIndicatorFunction(vmin,threshold);
            [cutPart1,cutPart2] = computeCutValue(allClusters,W,normalized); %cutPart1: vmin<threshold, cutPart2: vmin>threshold
            cut=cutPart1+cutPart2;
            cheeger=max(cutPart1,cutPart2);
    else

            deg=sum(W);
            [vmin_sorted, index]=sort(vmin);
            W_sorted=W(index,index);

            % sum of all degrees in the cluster minus weights within cluster
            volumes_threshold=cumsum(deg(index));
            tempcuts_threshold=volumes_threshold - 2*cumsum(full(sum(triu(W_sorted))));

            % divide by volume/size
            if(normalized)
                cutparts1_threshold=tempcuts_threshold(1:end-1)./volumes_threshold(1:end-1);
                cutparts2_threshold=tempcuts_threshold(1:end-1)./(volumes_threshold(end)-volumes_threshold(1:end-1));
            else
                sizes_threshold=cumsum(ones(1,size(vmin,1)-1));
                cutparts1_threshold=tempcuts_threshold(1:end-1)./sizes_threshold;
                cutparts2_threshold=tempcuts_threshold(1:end-1)./(size(vmin,1)-sizes_threshold);
            end

            % calculate cuts/cheegers
            cuts_threshold=cutparts1_threshold+cutparts2_threshold;
            cheegers_threshold=max(cutparts1_threshold,cutparts2_threshold);

            % find best cut/cheeger
            if(criterion==1)
                [cut,threshold_index]=min(cuts_threshold);
                cheeger=cheegers_threshold(threshold_index);
     		elseif(criterion==3)
                imbalance=0.05;
                totalsize=size(vmin,1);
                maxsize=floor(ceil(totalsize/2)*(1+imbalance));
                [unbalancedcut,threshold_index]=min(tempcuts_threshold(totalsize-maxsize:maxsize));
                threshold_index=totalsize-maxsize - 1 + threshold_index;
                cut=cuts_threshold(threshold_index);
				cheeger=cheegers_threshold(threshold_index);
            else
                [cheeger,threshold_index]=min(cheegers_threshold);
                cut=cuts_threshold(threshold_index);
            end

            % update
            cutPart1=cutparts1_threshold(threshold_index);
            cutPart2=cutparts2_threshold(threshold_index);

            allClusters= computeClusterIndicatorFunction(vmin,vmin_sorted(threshold_index));
            
            threshold=vmin_sorted(threshold_index);

    end
    
end