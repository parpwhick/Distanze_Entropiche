from pylab import *
import numpy
from scipy.cluster.hierarchy import fcluster
from fastcluster import linkage
from scipy.spatial.distance import squareform
from os import chdir
#from hcluster import *

def choose_p(Z):
        n_max = 50
        h = zeros(n_max + 1)
        x = range(1, n_max + 1)
        N = len(Z) + 1
        for n in x:
                c = fcluster(Z, n, criterion='maxclust')
                for i in range(1, n + 1):
                        mu = sum(c == i)
                        if mu:
                            h[n] += mu * log(mu)
                h[n] = -h[n] / N + log(N)

        h = h[1:]
        fattore=1/mean(h[arange(2,n_max)] / log(arange(2,n_max)))
        plot(x,h,'o-',x,log(x)/fattore,'--',linewidth=2)
        xlabel('p clusters')
        ylabel('entropia della popolazione dei cluster')
        print h
        show()

        n = int(raw_input("Inserisci numero cluster: "))
        return(n)

#chdir('E:\Dawid-software\distanze\H3N2')
distr = numpy.fromfile('output-distr.bin', dtype=double)
N = sqrt(len(distr))

if N != round(N):
        print "Dimensioni sbagliate"
        raise

N = int(N)
distr = distr.reshape((N, N))
distr = distr + distr.transpose()



#read dates and vaccines
mask=set((int(lines.split()[0])  for lines in  open('vaccine_data.txt')))
efficaci=array([line for line in arange(0,N) if line not in mask])
date=[lines.strip().split('/')  for lines in  open('sequenze_date.txt')]
date=[map(int,tripletta) for tripletta in date]

#2008,feb,7 -> 2008.10278
all_labels=[]
for i,data in enumerate(date):
        if i not in efficaci:
                continue
        all_labels.append(data[0]+(data[1]-1)/12.0+data[2]/365.0)
all_labels=array(all_labels)

#choose how many clusters and make final vector c
dist=distr[ix_(efficaci,efficaci)]


def hist_per_stagione(start=1992, end=2012):
	stagione=(all_labels > start) & (all_labels < end)
	dist_selected=dist[ix_(stagione,stagione)]
	Z=linkage(squareform(dist_selected),method='complete')
	n=choose_p(Z)
	c=fcluster(Z,n,criterion='maxclust')-1

	label_anni=all_labels[stagione]
	#order by first appearance!
	first_appearance=[]
	for i in range(0,n):
	    first_appearance.append(min(label_anni[c==i]))

	order1=[index for key,index in sorted(zip(first_appearance,range(0,n)))]
	order2=[index for key,index in sorted(zip(order1,range(0,n)))]
	order=array(order2)
	c=order[c]

	#draw scatter plot
	scatter(label_anni,c,s=100,c=c)
	#grid(b=True,axis='y')
	yticks(range(0,n+1))
	xlim((min(label_anni)-0.5,max(label_anni)+0.5))
	ax=gca()
	for i in range(1993,2011+1):
		ax.add_line(Line2D([i+7./12,i+7./12],[0,n+1],linestyle='--'))
	show()

for stagione in [2013]:
	hist_per_stagione(end=stagione)
