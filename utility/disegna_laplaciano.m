function [v,indici]=disegna_laplaciano(distanza)


zI=(diag(sum(distanza,2)));
L=zI-distanza;
v=real(eig(L));
indici=1:100;
plot(v(indici),indici,'-o');
% xlabel('Indice autovalori');
% ylabel('Autovalori dello pseudolaplaciano');
end
