%clear dist*;
fclose('all');
clear fid*;

bins=200;
r1=dir('output*');
have_simple= ~isempty(r1);


if(have_simple)
    n= round(sqrt(r1(2).bytes/8));
    fid1=fopen('output-shan.bin','r');
    fid2=fopen('output-shan_r.bin','r');
    fid3=fopen('output-top.bin','r');
    fid4=fopen('output-top_r.bin','r');
    
    %distanze shannon, normali
    dist_n=fread(fid1,[n,n],'double');
    %distanze shannon, ridotte
    dist_r=fread(fid2,[n,n],'double');
    %distanze topologiche, normali
    dist_t=fread(fid3,[n,n],'double');
    %distanze topologiche, ridotte
    dist_t_r=fread(fid4,[n,n],'double');
    
    [hist_n,labels1]=hist(dist_n(dist_n>0),bins);
    [hist_r,labels2]=hist(dist_r(dist_r>0),bins);
    [hist_t,labels3]=hist(dist_t(dist_t>0),bins);
    [hist_t_r,labels4]=hist(dist_t_r(dist_t_r>0),bins);
    
    dist_n=dist_n+dist_n';
    dist_r=dist_r+dist_r';
    dist_t=dist_t+dist_t';
    dist_t_r=dist_t_r+dist_t_r';
    
    subplot(2,2,1), plot(labels1,hist_n,'-*r',labels3,hist_t,'-*b');
    legend('Shannon','Topologica');
    xlabel('Distanza');
    ylabel('Frequenza');
    title('Distanze normali');
    subplot(2,2,2), plot(labels2,hist_r,'-*r',labels4,hist_t_r,'-*b');
    legend('Shannon rid.','Topologica rid.');
    xlabel('Distanza');
    ylabel('Frequeza');
    title('Distanze ridotte');
end

try 
    fid5=fopen('output-hamm.bin','r');
    dist_ham=fread(fid5,[n,n],'double');
    [hist_h,labels5]=hist(dist_ham(dist_ham>0),100);
    dist_ham=dist_ham+dist_ham';
catch ME
    % do nothing
end
    

if(have_fuzzy)
    n= round(sqrt(r2(1).bytes/8));
    fid6=fopen('output-fuzzy.bin','r');
    fid7=fopen('output-fuzzyt.bin','r');
    fid8=fopen('output-fuzzyr.bin','r');
    fid9=fopen('output-fuzzyrt.bin','r');
    
    %distanze fuzzy
    dist_f=fread(fid6,[n,n],'double');
    %distanze fuzzy topologiche
    dist_f_t=fread(fid7,[n,n],'double');
    %distanze fuzzy topologiche
    dist_f_r=fread(fid8,[n,n],'double');
    %distanze fuzzy topologiche
    dist_f_r_t=fread(fid9,[n,n],'double');
    
    [hist_f,labels6]=hist(dist_f(dist_f>0),100);
    [hist_f_t,labels7]=hist(dist_f_t(dist_f_t>0),100);
    [hist_f_r,labels8]=hist(dist_f_r(dist_f_r>0),100);
    [hist_f_r_t,labels9]=hist(dist_f_r_t(dist_f_r_t>0),100);
    
    %[hist_h,labels5]=hist(dist_ham(dist_ham>0),100);
    %dist_ham=dist_ham+dist_ham';
    
    dist_f=dist_f+dist_f';
    dist_f_t=dist_f_t+dist_f_t';
    dist_f_r=dist_f_r+dist_f_r';
    dist_f_r_t=dist_f_r_t+dist_f_r_t';
    
    subplot(2,2,3);
    plot(labels6,hist_f,'-*r',labels7,hist_f_t,'-*b');
    legend('Salto','Salto top.');
    xlabel('Distanza');
    ylabel('Frequenza');
    title('Partizioni con salto');
    
    subplot(2,2,4);
    plot(labels8,hist_f_r,'-*r',labels9,hist_f_r_t,'-*b');
    legend('Salto rid.','Salto rid. top.');
    xlabel('Distanza');
    ylabel('Frequenza');
    title('Partizioni con salto ridotte');    
end
