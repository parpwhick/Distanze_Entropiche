function [c,numerosita] = draw_ready_clustering(c,n)

N=length(c);

sequenze_date=importdata('sequenze_date.txt','/');
label_anni=sequenze_date(:,1)+(sequenze_date(:,2)-1)/12+(sequenze_date(:,3)-1)/365;
% attenzione al giusto numero di label
if(size(label_anni,1) ~= N)
    label_anni=label_anni(1:N,:);
end

%vettore di prime apparizioni di un cluster,
%in modo da ordinarli per data
first_appearance=zeros(1,n);
for i=1:n
    first_appearance(i)=min(label_anni(c==i));
end
[~,order1]=sort(first_appearance);

%inversione della permutazione(del sorting)
%necessaria a riordinare i cluster
order2=order1;
for i=1:n
    order2(order1(i))=i;
end

%rinominiamo i cluster con il loro label riordinato
c=order2(c);
%calcolo del numero di elementi per cluster
numerosita=1:n;
for i=1:n
    numerosita(i)=sum(c==i);
end

%colori forniti dall'utente, per mostrare come si sono spostati
%i vecchi cluster
if(nargin<3)
    colori=c;
else
    colori=old_clusters;
    provenienza=zeros(n);
     for i=1:n
        for j=1:n
            provenienza(i,j)=sum(c==i & old_clusters==j);
        end
     end
    figure;
    normalizz=ones(n,1);
    normalizz=sum(provenienza,2);
    %normalizz=normalizz ./ max(1,log10(normalizz));
    subplot(2,1,1);
    bar(provenienza ./ (normalizz * ones(1,n)),'stacked')
    xlabel('Nuovi cluster');
    ylabel('Logaritmo del numero di elementi nel cluster');
    title('Composizione percentuale dei nuovi cluster risp. ai vecchi');
    subplot(2,1,2);
    bar(diag(sum(provenienza,2)),'stacked');
    xlabel('Nuovi cluster');
    ylabel('Numero di elementi nel cluster');
    
end
figure;
%scatter plot, usando x=data, y=cluster label, colori=cluster label,
%dimensione dei punti 50
subplot(1,2,1);
scatter(label_anni,c,30,colori,'filled');
ylim([0.5,n+0.5]);
xlim([floor(min(label_anni)),ceil(max(label_anni))]);
subplot(1,2,2);
barh(numerosita);
ylim([0.5,n+0.5]);

%scatterhist(label_anni,c);
%colormap('lines');


% %cerchiamo i cluster piu numerosi!
% l=find(numerosita>100);
% for cur_cluster=l
%     siti=find(c==cur_cluster);
%     da_labellare=siti(ceil(rand(1,10)*length(siti)));
% 
%     for i=1:10
%     text(label_anni(da_labellare(i)),c(da_labellare(i)),num2str(cur_cluster),'EdgeColor','black','BackgroundColor','white');
%     end
% end

%cerchiamo i cluster piu numerosi!
if(nargin<3)
    return;
end

subplot(1,2,1);
da_labellare=ceil(rand(1,100)*N);

for i=da_labellare;
    text(label_anni(i),c(i),num2str(colori(i)),'EdgeColor','black','BackgroundColor','white');
end
    
end

