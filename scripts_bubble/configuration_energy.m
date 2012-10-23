function e_bubble=configuration_energy(t,J,s)

[L,~] = size(s);
N = L^2;

if(max(s)==1)
    s = s / 2;
end

lambda = 600;
%nearest neighbors to fill up the matrix
nnr = reshape(1:N,L,L);
nnd = circshift(nnr,-1);
nnu = circshift(nnr,1);
nnl = (circshift(nnr',1))';
nnr = (circshift(nnr',-1))';
colonne = repmat((1:N)',5,1);
righe = [nnl(:);nnr(:);nnd(:);nnu(:);(1:N)'];
%sum of the neighboring spins for every site
neigh_sum = sum(s([nnl(:),nnr(:),nnd(:),nnu(:)]),2);

% dopon spin up
values=([t*ones(4*N,1);lambda * (+1/2*s(:)+1/4) + J*0.5*neigh_sum]);
H_up= sparse(righe,colonne,values);

% dopon spin down
values=([t*ones(4*N,1);lambda * (-1/2*s(:)+1/4) - J*0.5*neigh_sum]);
H_down= sparse(righe,colonne,values);

opts.issym = true;
E_down = eigs(H_down,1,'SA',opts);
E_up = eigs(H_up,1,'SA',opts);

e_mag = sum(sum(s(nnd) .* s + s(nnl) .* s ));
e_mag = (e_mag + 1/4 * 2*N) * J;

e_bubble = min(E_down,E_up) + e_mag;

end