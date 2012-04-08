function threshold = determineThreshold(threshold_type,u)
% Select the treshold for the cluster indicator function
%
% Usage: threshold = determineThreshold(threshold_type,u)
%
% Different threshold types: 0: zero 1: median 2: mean
%
% Written by Thomas Bühler and Matthias Hein
% Machine Learning Group, Saarland University
% http://www.ml.uni-saarland.de

    assert(threshold_type==0 || threshold_type==1 || threshold_type==2);
    
    switch threshold_type
        case 0
            threshold=0;
        case 1
            threshold = median(u);
        case 2
            threshold = mean(u);
    end
       
end