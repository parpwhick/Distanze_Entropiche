stringa1=[1,0,0,0,1,0,0,1,0,0,1,1,1]
stringa2=[1,0,1,0,0,1,0,1,0,0,1,0,1];
comune=bitand(stringa1,stringa2)
ridotto2=bitxor(stringa2,comune)
ridotto1=bitxor(stringa1,comune)
ridotto1(end)=1;
ridotto1(1)=1;
ridotto2(1)=1;
ridotto2(end)=1;

dove=5;
multiplier=40;
sep=1;

a=zeros(60,1000);

for i=1:length(stringa1)
width=2*stringa1(i)+sep;
a(dove:dove+5,(i*multiplier-width):(i*multiplier+width))=1;
end
dove=dove+10;

for i=1:length(stringa1)
width=2*stringa2(i)+sep;
a(dove:dove+5,(i*multiplier-width):(i*multiplier+width))=1;
end
dove=dove+10;

for i=1:length(stringa1)
width=2*comune(i)+sep;
a(dove:dove+5,(i*multiplier-width):(i*multiplier+width))=1;
end
dove=dove+10;

for i=1:length(stringa1)
width=2*ridotto1(i)+sep;
a(dove:dove+5,(i*multiplier-width):(i*multiplier+width))=1;
end
dove=dove+10;

for i=1:length(stringa1)
width=2*ridotto2(i)+sep;
a(dove:dove+5,(i*multiplier-width):(i*multiplier+width))=1;
end
done=dove+10;


imagesc(a)
colormap([1,1,1;0,0,0])
