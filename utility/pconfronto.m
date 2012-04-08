cvp=zeros(1,100);
Z=linkage(squareform(dist_scelta),'complete');
for i=1:100
cvp(i)=populationEntropy(cluster(Z,'maxclust',i)');
end
fattore=1/mean(cvp(2:100)./log(2:100));
indici=1:50;
plot(indici,cvp(indici),'o-',indici, log(indici)/fattore,'--')
legend('Clustering gerarchico','logaritmo naturale riscalato')
xlabel('p clusters')
ylabel('entropia della popolazione dei cluster')
