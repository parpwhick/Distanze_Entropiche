fclose('all');
clear fid*;

r1=dir('output-distr.bin');
have_simple= ~isempty(r1);
if(~have_simple)
	disp 'Nessun file presente'
	return
end

n= round(sqrt(r1(1).bytes/8));

fid2=fopen('output-distr.bin','r');

dist_r=fread(fid2,[n,n],'double');

dist_r=dist_r+dist_r';
