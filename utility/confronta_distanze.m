files=dir('output-fuzzy-*');
n=822;
L=length(files);


for i=1:L
    fid1=fopen(files(i).name,'r');
    dist_f=fread(fid1,[n,n],'double');
    [hist_f,labels]=hist(dist_f(dist_f>0),200);
    istogrammi(i).labels=labels;
    istogrammi(i).freq=hist_f;
    istogrammi(i).fuzzy=str2num(files(i).name([14:16]));
end
fclose(fid1);

files=dir('output-fuzzy_t-*');
n=822;
L=length(files);

for i=1:L
    fid1=fopen(files(i).name,'r');
    dist_f=fread(fid1,[n,n],'double');
    [hist_f,labels]=hist(dist_f(dist_f>0),200);
    istogrammi_t(i).labels=labels;
    istogrammi_t(i).freq=hist_f;
    istogrammi_t(i).fuzzy=str2num(files(i).name([14:16]));
end
fclose(fid1);

clear labels;
clear hist_f;
clear fid1;
clear dist_f;
clear i;

h1=plot(istogrammi(1).labels,istogrammi(1).freq,'EraseMode','normal');
legend('Fuzzyness: 1');
for i=1:L
set(h1,'Ydata',istogrammi(i).freq);
legend(['Fuzziness: ',num2str(istogrammi(i).fuzzy)]);
pause(1);
end

