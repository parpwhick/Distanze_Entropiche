function show_flippable(conf,energies)

[L,~] = size(conf);
[L2,~] = size(energies);
if(L ~= L2)
    disp('Wrong size')
    return
end
N = L^2;


%nearest neighbors to fill up the matrix
nnr = reshape(1:N,L,L);
nnd = circshift(nnr,-1);
nnu = circshift(nnr,1);
nnl = (circshift(nnr',1))';
nnr = (circshift(nnr',-1))';

dH = - 2 * conf .* (conf(nnu)+conf(nnd)+conf(nnl)+conf(nnr));

imagesc(dH <= energies)
end