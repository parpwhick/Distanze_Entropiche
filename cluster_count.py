from pylab import *
import numpy
import scipy.cluster.hierarchy as cluster
from scipy.spatial.distance import squareform
from os import chdir

#base_dir=
chdir('E:\Dawid-software\distanze\H3N2')

distr=numpy.fromfile('output-distr.bin',dtype=double)
N=sqrt(len(distr))

if N != round(N):
        print "Dimensioni sbagliate"
        raise

N=int(N)
n=14
distr=distr.reshape((N,N))
distr=distr+distr.transpose()

Z=cluster.linkage(squareform(distr),method='complete')


mask=set((int(lines.split()[0])  for lines in  open('vaccine_data.txt')))
efficaci=array([line for line in arange(0,N) if line not in mask])
date=[lines.strip().split('/')  for lines in  open('sequenze_date.txt')]
date=[map(int,tripletta) for tripletta in date]

n_max=50;
h=zeros(n_max+1)
x=range(1,n_max+1)
for n in x:
        c=cluster.fcluster(Z,n,criterion='maxclust')
        for i in range(1,n+1):
                mu=sum(c==i)
                h[n]+=mu*log(mu)
        h[n]=-h[n]/N+log(N)

h=h[1:]
fattore=1/mean(h[arange(2,n_max)] / log(arange(2,n_max)))
plot(x,h,'o-',x,log(x)/fattore,'--',linewidth=2)
xlabel('p clusters')
ylabel('entropia della popolazione dei cluster')

label_anni=[]
for i,data in enumerate(date):
        if i not in efficaci:
                continue
        label_anni.append(data[0]+(data[1]-1)/12.0+data[2]/360.0)
label_anni=array(label_anni)

n=int(raw_input("Inserisci numero cluster: "))
close()
c=cluster.fcluster(Z,n,criterion='maxclust')-1


first_appearance=[]
for i in range(0,n):
    first_appearance.append(min(label_anni[c==i]))

order1=[index for key,index in sorted(zip(first_appearance,range(0,n)))]
order2=[index for key,index in sorted(zip(order1,range(0,n)))]
order=array(order2)

c=order[c]
scatter(date_efficaci,c,s=50,c=c)

grid(b=True,axis='y')
yticks(range(0,14))
xlim((min(label_anni)-0.5,max(label_anni)+0.5))
