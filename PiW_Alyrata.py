from random import *
from numpy import *
from itertools import *
import os
import tabix

#import pysam
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help="VCF.gz, tabixed", default="")
parser.add_option("-d", "--dir", dest="dir", help="input directory", default="")
parser.add_option("-D", "--outdir", dest="outdir", help="output directory", default="")
parser.add_option("-n", "--num", dest="num", help="interval number (line number from intervals file)", default="")
parser.add_option("-C", "--chr", dest="chr", help="scaffold name", default="")
parser.add_option("-s", "--st", dest="st", help="start coordinate", default="")
parser.add_option("-e", "--end", dest="end", help="end coordinate", default="")
parser.add_option("-c", "--colnames", dest="colnames", help="file with column names", default="")
(options, args) = parser.parse_args()

for line in open('%s'%options.colnames):
	names=array(line.strip().split())
pop1,pop2=[],[]
for line in open ('%s.selfing'%options.colnames):
    pop1=array(line.strip().split())
for line in open('%s.outcrossing'%options.colnames):
    pop2=array(line.strip().split())
ind_pop1=in1d(names,pop1)
ind_pop2=in1d(names,pop2)

#gene=options.gene.strip().split('.')
#transcript='%s.t1.%s.%s'%(gene[0], gene[1], gene[2])
##########
L,Lal,Pb,Pw1,Pw2=0,0,0,0,0
#for line in open('%s'%options.ann):
#    if '#' in line:
#        pass
#    else:
#        line=line.strip().split('\t')
#        if transcript in line[8] and line[2]=='CDS':
#            chr,st,end=line[0],min([int(line[3]),int(line[4])]), max([int(line[3]),int(line[4])])
chr,st,end=options.chr, int(options.st), int(options.end)
L+=(end-st)
tb = tabix.open('%s/%s'%(options.dir,options.vcf))
records = tb.query(chr, st, end)
for record in records:
    gen=array([x.split(':')[0] for x in record])
    ref,alt=record[3],record[4]
    if len(ref)==1 and len(alt)==1 and alt!='*':
        gen1, gen2=gen[ind_pop1], gen[ind_pop2]
        bad=['./.', './././.']
        gen1f, gen2f=[x for x in gen1 if x not in bad], [x for x in gen2 if x not in bad]
        if len(gen1f)>=0.5*len(gen1) and len(gen2f)>=0.5*len(gen2): #mis <50%
            gen1f=[item for sublist in [list(x[0::2]) for x in gen1f] for item in sublist]
            gen2f=[item for sublist in [list(x[0::2]) for x in gen2f] for item in sublist]
            gen1f,gen2f=[int(x) for x in gen1f],[int(x) for x in gen2f]
            Lal+=1
            if len(set(gen1f+gen2f))==2:
                between=[(i, j) for i in gen1f for j in gen2f]
                P=0
                for i in between:
                    if len(set(list(i)))>1:
                        P+=1
                Pb+=float(P)/len(between)
		if len(set(gen1f))==2:
		    Ncomp,P=0,0
		    for s1,s2 in combinations(gen1f,2):
		        Ncomp+=1
		        if s1!=s2:
		            P+=1
		    Pw1+=(float(P)/Ncomp)
		if len(set(gen2f))==2:
		    Ncomp,P=0,0
		    for s1,s2 in combinations(gen2f,2):
		        Ncomp+=1
		        if s1!=s2:
		            P+=1
		    Pw2+=(float(P)/Ncomp)
print (L, Lal)
out=open('%s/%s.Pi'%(options.outdir, options.num), 'w')
out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(chr, st, end, L,Lal,Pb,Pw1,Pw2))
out.close()
