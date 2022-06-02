#Coge to circular plot 

#######################################


file='/path/to/coge_outputfile' # what we downloaded from coge
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


#copy and paste the last line into x.ends 

#now x.starts and x.ends are the same length
#lets cbind them in R


