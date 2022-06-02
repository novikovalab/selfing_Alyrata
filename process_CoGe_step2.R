#merge x.starts and x.ends file together
####################################


mystarts=read.delim(file='x.starts', header=F)
myends=read.delim(file='x.ends', header=F) #manually pasted extra line omn command line

mystarts=mystarts[,c(3,5,7,9)]
myends=myends[,c(3,6,7,10)]
mystart_ends=cbind(mystarts, myends)
mystart_ends_blocks=mystart_ends[,c(1,2,6,3,4,8)]
colnames(mystart_ends_blocks)=c( "newMN47_chr", "newMN47_start", "newMN47_end","newNT1_chr", "newNT1_start", "newNT1_end") #name columns according to assembly
mystart_ends_blocks$newNT1_chr=gsub("b63282", "newNT1", mystart_ends_blocks$newNT1_chr)
mystart_ends_blocks$newMN47_chr=gsub("a62989_", "newMN47_", mystart_ends_blocks$newMN47_chr)
mystart_ends_blocks$newMN47_chr=gsub("scafold", "scaffold", mystart_ends_blocks$newMN47_chr)
write.table(mystart_ends_blocks, file='synteny.newMN47_newNT1_nonmerge.txt', col.names=T, row.names=F, quote=F, sep='\t') #save as a new table, containg all colinear and inversed regions between the two assemblies

#########################
#now do this for each combination of assemblies
