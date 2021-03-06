function cluster_stagionali(distanza)
[N_tot,m]=size(distanza);
assert(N_tot==m, 'La matrice non e` quadrata');

figura_conta=figure;
figura_scatter=figure;
pause

fid=fopen('vaccine_data.txt');
vaccine = textscan(fid, '%d %s %s');

effettivi=ones(N_tot,1);
effettivi(vaccine{1})=0;
effettivi=logical(effettivi);
n_good=sum(effettivi);
disp(['Lette ' num2str(N_tot) ' sequenze, ' num2str(n_good) 'utili']);
sequenze_date=importdata('sequenze_date.txt','/');
label_anni=sequenze_date(:,1)+(sequenze_date(:,2)-1)/12+(sequenze_date(:,3)-1)/365;
label_anni(isnan(label_anni))=sequenze_date(isnan(label_anni),1);
% attenzione al giusto numero di label
if(size(label_anni,1) ~= N_tot)
    disp('Errato numero di etichette anni');
    return
    %label_anni=label_anni(1:N,:);
end

for anno=1996:2013
    draw_to_season(anno+0.585)
end


    function draw_to_season(fino_a)
        

        scelti = effettivi & (label_anni < fino_a);
        N=sum(scelti);
        if(N<50)
            disp(['Warning, only ' num2str(N) ' sequences until ' num2str(fino_a)])
            return
        end
        
        Z=linkage(squareform(distanza(scelti,scelti)),'complete');
        
        %calcolo del vettore entropie h
        n=choose_cluster_p(Z,N,figura_conta);
        %n=14;
        
        c=cluster(Z,'maxclust',n);
        
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
        colori=c;
        
        
        %scatter plot, usando x=data, y=cluster label, colori=cluster label,
        %dimensione dei punti 80
        figure(figura_scatter);
        title(['Dal 1992 al ' num2str(ceil(fino_a))]);
        ggg=subplot(1,2,1);
        scatter(label_anni(scelti),c,120,colori, ... 
                    'filled','MarkerEdgeColor','black');
        set(ggg,'YTick',1:1:n,'Box','on','YGrid','on');
        ylim([0.5,n+0.5]);
        xlim([floor(min(label_anni(scelti))),ceil(max(label_anni(scelti))+0.5)]);
        da=floor(min(label_anni(scelti)));
        for end_season=(da:floor(fino_a))+0.585
            line([end_season,end_season],[0.5,n+0.5],'Color','Blue','LineStyle','--')
        end
        ggg=subplot(1,2,2);
        barh(numerosita);
        ylim([0.5,n+0.5]);
        set(ggg,'YTick',1:1:n,'Box','on','YGrid','on');
        set(get(ggg,'XLabel'),'Rotation',30)
        
        pause(1)
        if(~isempty(vaccine))
            draw_vaccine(c,scelti);
        end
    end

    function draw_vaccine(c,scelti)
        % now draw vaccine strain info
        subplot(1,2,1);
        for i=1:length(vaccine{1})
            sigla=vaccine{2}{i};
            posizione=vaccine{1}(i);
            colore_testo=vaccine{3}{i};
            if(isempty(colore_testo))
                colore_testo='b';
            end
            %                 for j=1:n
            %                     distanze_da_cluster(j)=mean(distanza(c==j,posizione));
            %                 end
            %                 [~,appartenenza]=min(distanze_da_cluster);
            [~,appartenenza]=min(distanza(scelti,posizione));
            appartenenza=c(appartenenza);
            
            text(label_anni(posizione),appartenenza+0.5,sigla, 'EdgeColor', ...
                'black','BackgroundColor','white','Color',colore_testo);
        end        
    end
end

function n=choose_cluster_p(Z,N,figura_conta)
global nome;
M=35;
indici=1:M;
h=zeros(1,M);
for j=indici
    c=cluster(Z,'maxclust',j);
    for i=1:max(c)
        mu=sum(c==i);
        h(j)=h(j)+mu*log(mu);
    end
    h(j)=-h(j)/N+log(N);
end

fattore=1/mean(h(2:M)./log(2:M));

figure(figura_conta)
plot(indici,h(indici),'o-',indici, log(indici)/fattore,'--')

if(~isempty(whos('global','nome')))
    aggiungere=nome;
else
    aggiungere='';
end

legend(['Clustering gerarchico',aggiungere],'logaritmo naturale riscalato','Location','SouthEast')
xlabel('p clusters')
ylabel('entropia della popolazione dei cluster')

n=input('Numero di cluster? ');
end


