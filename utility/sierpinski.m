function adiacenza=sierpinski(potenza)

if(nargin<1)
    potenza=12;
end

N=2^potenza;
siti=3^potenza;

disp(['Triangolo sierpinski, generazione ' num2str(potenza) ', siti ' num2str(siti)])
disp(['Elementi nonnulli matrice di adiacenza(simm) ' num2str((siti-1)*3)])
% l1: riga 'precedente' nell'iterazione
% l2: riga 'corrente'
l1=zeros(N,1);
l2=zeros(N,1);
l1(1)=1;
l2(1)=1;
% etichette degli elementi nonnulli interazione precedente
etichette1=l1;

%matrice di adiacenza tridiagonale
adiacenza=zeros(siti,3);

%disegno=zeros(N);
%disegno(1,:)=l1;
nonnulli=1;
for i=2:N
    l2(2:i)=mod(l1(2:i)+l1(1:(i-1)) , 2);
    da_numerare=find(l2);
    quanti=length(da_numerare);
    
    % etichette iterazione corrente
    etichette2=zeros(N,1);
    %si assegnano etichette crescenti ai nonnulli
    new_labels=(nonnulli+1):(nonnulli+quanti);
    etichette2(da_numerare)=new_labels;
    nonnulli=nonnulli+quanti;
    
    %controllare i 3 vicini precedenti, per avere la matrice di
    %adiacenza triang. inferiore
    
    % escludiamo il primo elemento, che e' 1, non avendo vicini a sinistra
    % lo trattiamo separatamente 
    da_numerare2=da_numerare(2:end);
    adiacenza(new_labels(1),:)=[etichette1(1),0,0];
    
    % i possibili vicini sono quelli a sinistra, in alto, in alto a
    % sinistra
    adiacenza(new_labels(2:end),:)= [etichette1(da_numerare2-1), ...
                                etichette1(da_numerare2), ...
                                etichette2(da_numerare2-1)];
    
    %disegno(i,:)=l2;
    
    %iterazione corrente diventa precedente
    l1=l2;
    etichette1=etichette2;
end
%disp(nonnulli)
%imagesc(disegno);
disp(' ')
disp('Calcolo completato, compilazione matrice sparse');
% costruzione della matrice di adiacenza propriamente detta
[indici_riga,~,indici_colonna]=find(adiacenza);
clear adiacenza;
adiacenza=sparse(indici_riga,indici_colonna,1,siti,siti,(siti-1)*3);

% commentare la riga se si vuole solo una matrice lower triangular
adiacenza=adiacenza+adiacenza';
end