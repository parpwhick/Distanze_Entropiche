function [V_out,clusters_out]=...
    spectral_clustering(K,number_of_clusters,kmeans_or_not,plots_or_not)

%function [V_out,clusters_out]=...
%    spectral_clustering(K,number_of_clusters,...
%    kmeans_or_not,plots_or_not)
%
%Performs standard spectral clustering. The rows of V_out are normalized
%after removing the constant column. The user can choose between additive
%or divisive normalization.
%
%INPUTS:
% K                     = the affinity matrix
% number_of_clusters    = how many clusters do you want to find
% kmeans_or_not         = boolean indicating if kmeans needs to be
%                         performed or not
% plots_or_not          = boolean indicating if plots need to be outputted or not
%
%OUTPUTS:
% V_out        = the eigenvectors output, AFTER removing the constant eigenvector,
%                and WITH normalizing their rows
% clusters_out = the clusters output, if requested by putting kmeans_or_not=1
%
%

%Author: Tijl De Bie, december 2003.

%%%%%%%%%%%%%
% Parameters
restarts_for_Kmeans=100;
%%%%%%%%%%%%%

size_K=size(K,1);

if issparse(K)
    D=spdiags(sum(K,2),0,size_K,size_K);
else
    D=diag(sum(K,2));
end
% Compute matrix for which eigenvectors need to be computed

        LH=D-K;
        RH=D;

% Compute eigenvectors, and normalize
opts.disp=0;
[V,eigenvalues]=eigs(LH,RH,number_of_clusters+1,'SM',opts);
[mdismiss,dismiss]=min(abs(diag(eigenvalues)));
V(:,dismiss)=[];

for i=1:size(V,2)
    V(:,i)=V(:,i)/norm(V(:,i));
end
%dv=sqrt(diag(V*V'))+1e-8;
%V_out=diag(1./dv)*V;
V_out=V;


% Perform K-means if requested
if kmeans_or_not
    % Generate initial conditions
    ic=zeros(number_of_clusters,number_of_clusters,restarts_for_Kmeans);
    for ct1=1:restarts_for_Kmeans
        new_ind=ceil(rand*size_K);
        ic_ct1=V_out(new_ind,:);
        for ct2=2:number_of_clusters
            inner_products=V_out*ic_ct1';
            worst_inner_products=max(abs(inner_products),[],2);
            new_ind=find(worst_inner_products==min(worst_inner_products));
            new_ind=new_ind(1);
            ic_ct1=[ic_ct1 ; V_out(new_ind,:)];
        end
        ic(:,:,ct1)=ic_ct1;
    end
    % Try to do kmeans with these initial conditions 
    try
        clusters_out=kmeans(V_out,number_of_clusters,'Start',ic,'EmptyAction','singleton');
    catch
        clusters_out=0;
    end
else
    clusters_out=[];
end

% Give plot output if requested
if plots_or_not
    figure
    plot(V_out)
    grid on
    if kmeans_or_not
        hold on
        plot(clusters_out-number_of_clusters/2-0.5)
    end
end
