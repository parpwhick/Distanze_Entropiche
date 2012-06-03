function save_adiacenza(A)

[j,i]=find(A);

fwrite(fopen('vector1.bin','wb'),i,'int32');
fwrite(fopen('vector2.bin','wb'),j,'int32');
fclose('all');

end