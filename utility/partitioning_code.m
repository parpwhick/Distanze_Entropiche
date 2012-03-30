H = @(a) - sum(a .* log(a)) / sum(a) + log(sum(a));

makebinpart=@(a) [1, a(2:end) ~= a(1:end-1)]

primi = primes(2000);

binary_part = @(e) e(2:end) ~= e(1:end-1);
xorand = @(a,b) bitxor(a,bitand(a,b));


Hbin  = @(e) H(diff(find([1,e,1])));
Hpart = @(e) H(diff(find([1,binary_part(e),1])));

Hprod = @(a,b) 2*Hbin(bitor(binary_part(a), binary_part(b))) - Hbin(binary_part(a)) - Hbin(binary_part(b));
Hprodbin = @(a,b) 2*Hbin(bitor(a,b)) - Hbin(a) - Hbin(b);
Hprodbinrid =  @(a,b) 2*Hbin(bitxor(a, b)) - Hbin(xorand(a,b)) - Hbin(xorand(b,a));
Hprodrid = @(a,b) Hprodbinrid(binary_part(a),binary_part(b));

prod_gen = @(a,b) primi(a+1) .* primi(b+length(a) +2);
Hprod_gen = @(a,b) 2*Hpart(sort(prod_gen(a,b))) - Hpart(sort(a)) - Hpart(sort(b));

Htopbin = @(e) log(sum(e)+1);

% tic
% for i=1:822
% for j=i+1:822
% d(i,j) = Hprodbinrid(binpart(i,:), binpart(j,:));
% end
% end
% toc