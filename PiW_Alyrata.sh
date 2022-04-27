#! /bin/bash
#BSUB -J pi[1-25000:100]. #go over jobindices in increments of 100, here 25000 in total, since we had 24028 intervals
#BSUB -q normal
#BSUB -M 2000 -R "rusage[mem=2000] "
#BSUB -n 1
#BSUB -o pi_%I_%J_output.txt
#BSUB -e pi_%I_%J_error.txt

inputdir='/path/to/folder/with/vcf'
vcf='example.vcf.gz'
intervals='/path/to/bedfile/containing/10kb/intervals/' 
#example of bedfile
#scaffold_1_RagTag	0	10000
#scaffold_1_RagTag	10000	20000
#scaffold_1_RagTag	20000	30000
#scaffold_1_RagTag	30000	40000
#...

out=${inputdir}'/Pi'

for i in {0..99}
do
 num=$(($i+${LSB_JOBINDEX})) #use jobindex as numerator
 export l=$num\p
 chr=`sed -n $l $intervals | cut -f1`
 st=`sed -n $l $intervals | cut -f2`
 end=`sed -n $l $intervals | cut -f3`
 python2 /path/to/pythonscript/PiW_Alyrata.py  -v $vcf -d $inputdir -D $out -C $chr -s $st -e $end -n $num -c ${inputdir}/${vcf}.colnames
done

