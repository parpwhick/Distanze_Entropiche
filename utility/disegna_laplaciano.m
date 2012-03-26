function disegna_laplaciano(distanza)


zI=(diag(sum(distanza,2)));
L=zI-distanza;
v=real(eig(L));
plot(v(1:100),'-o');
xlabel('Indice autovalori');
ylabel('Autovalori dello pseudolaplaciano');
end
