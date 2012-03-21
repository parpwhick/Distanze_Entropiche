function disegna_laplaciano(distanza,discreto)

if(nargin<2)
    dist_max=max(distanza(:));
    distanza=dist_max-distanza;
end

zI=(diag(sum(distanza,2)));
L=zI-distanza;
v=eig(L);
plot(v(1:100),'-o');
xlabel('Indice autovalori');
ylabel('Autovalori dello pseudolaplaciano');
end
