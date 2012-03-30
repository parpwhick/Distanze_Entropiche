function L=laplacian(distanza)

zI=(diag(sum(distanza,2)));
L=zI-distanza;
end
