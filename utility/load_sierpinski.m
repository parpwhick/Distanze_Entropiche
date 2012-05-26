indici_riga=fread(fopen('vector1.bin','r'),inf,'int');
indici_colonna=fread(fopen('vector2.bin','r'),inf,'int');
N=max(max(indici_riga),max(indici_colonna));
adiacenza=sparse(indici_riga,indici_colonna,1,N,N,(N-1)*3);