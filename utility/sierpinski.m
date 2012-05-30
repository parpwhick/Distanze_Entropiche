function [adiacenza,nn]=sierpinski(gen)

if(nargin<1)
    gen=10;
end

total_size = 6;
size = 6;

for g=2:gen
    total_size = 3 * total_size - 3;
end

disp(['Generazione ' num2str(gen) ', size ' num2str(total_size)])

z = zeros(total_size,1);
nn = zeros(total_size,4);

% vettore di numeri di coordinazione
% e struttura dei vicini
% per il triangolo di generazione 1
z(1:6)=[2,4,4,2,4,2];
nn(1,:)=[2,3,0,0];
nn(2,:)=[3,1,4,5];
nn(3,:)=[1,2,5,6];
nn(4,:)=[2,5,0,0];
nn(5,:)=[4,2,3,6];
nn(6,:)=[3,5,0,0];

%angoli
a0 = 1;
a1 = 4;
a2 = 6;

% costruzione iterativa gen successive
for g = 2:gen
    
    %triangolo 1
    %nulla da copiare
    
    %triangolo 2
    i=2:(size-1);
    nn(i+size-1,:) = nn(i,:) + size-1;
    z(i+size-1)=z(i);
    %correzioni
    [corr_i,corr_j]=find(nn(1:size,:)==a0);
    for i=1:length(corr_i)
        nn(corr_i(i)+size-1,corr_j(i))=a1;
    end
    [corr_i,corr_j]=find(nn(1:size,:)==a2);
    for i=1:length(corr_i)
        nn(corr_i(i)+size-1,corr_j(i))=a1+2*size-3;
    end
    
    %triangolo 3
    i=2:size;
    z(i+2*size-3)=z(i);
    nn(i+2*size-3,:) = nn(i,:) + 2*size - 3;
    %correzioni
    [corr_i,corr_j]=find(nn(1:size,:)==a0);
    for i=1:length(corr_i)
        nn(corr_i(i)+2*size-3,corr_j(i))=a2;
    end
    
    %correzioni dei bordi incollati
    z(a1)=4;
    nn(a1,3:4)= [size+1, size+2];
    z(a1+2*size-3)=4;
    nn(a1+2*size-3,3:4) = nn(a2,1:2) + size -1;
    z(a2)=4;
    nn(a2,3:4)=[2*size-1,2*size];
    
    %nuovi siti di angolo
    a1 = a1 + size - 1;
    a2 = a2 + 2* size -3;
    
    %aggiornamento dimensione prossima iterazione
    size= 3*size-3;
    
end

%rimozione link inesistenti, prima di costruire la matrice
nn(z==2,3:4)=0;

%estrazioni vettori di indici dei siti nonnulli
[indici_riga,~,indici_colonna]=find(nn);

%costruzione matrice vera e propria
adiacenza=sparse(indici_riga,indici_colonna,1,total_size,total_size);
end