#Coge to circular plot 

cd /netscratch/dep_mercier/grp_novikova/Neobatrachus/Jozefien_scratch/CoGe_analysis
#Step 1: change Coge into synteny blocks
#######################################


file='/netscratch/dep_mercier/grp_novikova/Neobatrachus/Jozefien_scratch/CoGe_analysis/52485_59218.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.ks.txt' # what we downloaded from coge
#file is from "	Results with synonymous/non-synonymous rate values" in download section of synmap
grep -v "##" $file > t
grep -v "#This" t > t2
mv t2 t


grep -A 1 "Ks" t > t.starts
grep -B 1 "Ks" t > t.ends
grep -v "#Ks" t.starts > x.starts
perl -i -p -e 's/-//g' x.starts

grep -v "#Ks" t.ends > x.ends
perl -i -p -e 's/-//g' x.ends
wc -l x.starts;wc -l x.ends
#Important x.ends is always one less than start because we miss the last line for t.ends by grep based on '#Ks'
#So we need to simply add this line back in
tail -1 $file


#copy and paste the last line into x.ends --> or do this in R manually like I did
#now x.starts and x.ends are the same length
#lets cbind them in R



R
library(dplyr)
