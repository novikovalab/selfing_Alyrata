library(circlize)
library(RColorBrewer)
library(stringr)
library(GISTools)

##########Fig 1A: three way- oldMN47 centered

#synteny blocks

oldMN47_Asue<-read.table('NT1_030222/synteny.oldMN47_Asue_nonmerge.txt', header=T)
oldMN47_NT1<-read.table('NT1_030222/synteny.oldMN47_newNT1_nonmerge.txt', header=T)
#cytobands
cytoband_oldMN47 = read.table("MN47/MN47_cytobands_all.txt", colClasses = c("character", "numeric", "numeric"), sep = " ")
cytoband_oldMN47$V1=gsub("scaffold_", "oldMN47_scaffold_", cytoband_oldMN47$V1)
cytoband_oldMN47_uniq<-cytoband_oldMN47[which(cytoband_oldMN47$V1 %in% unique(oldMN47_NT1$oldMN47_chr)),]
cytoband_oldMN47_uniq_rev<-cytoband_oldMN47_uniq[dim(cytoband_oldMN47_uniq)[1]:1,]

cytoband_Asue = read.table("suecica/cytobands_suecica.txt", colClasses = c("character", "numeric", "numeric"), sep = "\t")
cytoband_Asue_uniq<-cytoband_Asue[which(cytoband_Asue$V1 %in% unique(oldMN47_Asue$Asue_chr)),]
cytoband_Asue_uniq_rev<-cytoband_Asue_uniq[dim(cytoband_Asue_uniq)[1]:1,] #reverse cytobands of suecica as well, so that they are shown counter clock wise

cytoband_NT1 = read.table("newMN47/NT1_cytobands.txt", colClasses = c("character", "numeric", "numeric"), sep = " ")
cytoband_NT1$V1=gsub("scaffold_", "newNT1_scaffold_", cytoband_NT1$V1)
cytoband_NT1_uniq<-cytoband_NT1[which(cytoband_NT1$V1 %in% unique(oldMN47_NT1$newNT1_chr)),]
cytoband_NT1_uniq_rev<-cytoband_NT1_uniq[dim(cytoband_NT1_uniq)[1]:1,]


#step 3: merge all cytobands for one big sugergenome
library(data.table)
cytoband <- rbindlist(list(  cytoband_NT1_uniq, cytoband_oldMN47_uniq_rev, cytoband_Asue_uniq), use.names=T)
cytoband[[1]] = factor(cytoband[[1]], levels = cytoband$V1)


#step 4: inverse syntenyblocks in chr6,7,9 of Asue
data=oldMN47_Asue
for (i in 1:dim(data)[1]){
  if (data$Asue_chr[i]=='Asue_scaffold6'){
    
    data[i,5]=cytoband_Asue_uniq[1,3]-oldMN47_Asue[i,5]
    data[i,6]=cytoband_Asue_uniq[1,3]-oldMN47_Asue[i,6]
  }
  if (data$Asue_chr[i]=='Asue_scaffold7'){
    print(i)
    print(cytoband_Asue_uniq[2,3]-oldMN47_Asue[i,5])
    data[i,5]=cytoband_Asue_uniq[2,3]-oldMN47_Asue[i,5]
    data[i,6]=cytoband_Asue_uniq[2,3]-oldMN47_Asue[i,6]
  }
  if (data$Asue_chr[i]=='Asue_scaffold9'){
    data[i,5]=cytoband_Asue_uniq[4,3]-oldMN47_Asue[i,5]
    data[i,6]=cytoband_Asue_uniq[4,3]-oldMN47_Asue[i,6]
  }
}
oldMN47_Asue=data

#step 4: make the figure
facing <- c(rep("bending.inside", 8), rep("clockwise", 1), rep("bending.outside", 8), rep("bending.inside", 8))
padding_y <-c(rep(0.15, 8),rep(0.35, 1),  rep(0.15, 8),rep(0.15, 8))

sc_labels <- c(cytoband_NT1_uniq$V1, cytoband_oldMN47_uniq_rev$V1,  cytoband_Asue_uniq$V1)
sc_labels <- c(paste("Chr", 1:8 ,sep=''), "Scaffold9", paste("Chr", 8:1, sep=''), paste("Chr", 6:13, sep='') )

#iniate circos cirle and cytobands
colors <-brewer.pal(5, "BrBG")
#colors <-brewer.pal(5, "RdYlBu")

circos.clear()

circos.par("start.degree" = 90,canvas.ylim=c(-1.1,1.1),canvas.xlim=c(-1.1,1.1),clock.wise=T,
           track.margin=c(0,0),track.height=0.1, "gap.degree"= c(rep(2,7),15,rep(2,8),15,rep(2,7),15))
circos.genomicInitialize(cytoband, plotType = NULL,tickLabelsStartFromZero = TRUE,
                         axis.labels.cex = 0.3*par("cex"),track.height=0.1)
#circos.track(ylim = c(0, 0.05))
circos.track(ylim = c(0, 0.05),
             bg.col = c(rep(colors[2],8), rep(colors[5],9), rep(colors[1],8)),
             bg.border = NA, track.height = 0.05, panel.fun = function(x, y) {
               chr =get.cell.meta.data("sector.index")
               xlim =get.cell.meta.data("xlim")
               ylim =get.cell.meta.data("ylim")
               idx=match(get.cell.meta.data("sector.index"), cytoband$V1)
               circos.text(mean(xlim), mean(ylim)+padding_y[idx], sc_labels[idx], cex = 0.9, facing = facing[idx], niceFacing = F)
               #circos.text(mean(xlim), mean(ylim)+padding_y[idx], sc_labels[idx], cex = 0.6, facing = 'clockwise', niceFacing = T)
               #circos.text( sc_labels[idx], cex = 0.6, niceFacing = TRUE)
               #circos.axis(major.tick=TRUE, minor.ticks = 0,labels = c("30MB", "20MB","10MB", ""),labels.cex = 0.5)
               
             })

rev_x = function(x, xrange = CELL_META$xlim) {
  as.numeric(xrange[2]) - x + as.numeric(xrange[1])
}
color_synteny_first<-function(relations, cytobands_firstgenome, inv_color){
  colnames(relations)=c("V1", "V2", "V3", "V4", "V5", "V6")
  relations.rev = relations
  for ( i in 1:dim(cytobands_firstgenome)[1]){ #reverse synteny blocks
    inv <- cytobands_firstgenome[i,1]
    l = relations.rev$V1 == inv; relations.rev$V2=as.numeric(relations.rev$V2); relations.rev$V3=as.numeric(relations.rev$V3);
    relations.rev$V3[l] = rev_x(relations$V2[l], get.cell.meta.data("xlim",sector.index = inv))
    relations.rev$V2[l] = rev_x(relations$V3[l], get.cell.meta.data("xlim",sector.index = inv))
    relations.rev = relations.rev
  }
  data=relations.rev
  for (i in 1:dim(cytobands_firstgenome)[1]){ #plot synteny blocks per scaffold, can modify color per scaffold this way
    data_Scaffolds<-data[which(data$V1==cytobands_firstgenome[i,1]),]
    for (j in 1:length(data_Scaffolds[,1])){
      #print(i)
      #print(j)
      #MN47_idx=match( data_Scaffolds[j,1], cytoband_MN47_uniq$V1)
      #circos.genomicLink(data_Scaffolds[j,1:3], data_Scaffolds[j,4:6], col = mycolors[3], border = mycolors[3])
      if (as.numeric(data_Scaffolds[j,6])<as.numeric(data_Scaffolds[j,5])){
        print(data_Scaffolds[j,])
        print(abs(as.numeric(data_Scaffolds[j,6])-as.numeric(data_Scaffolds[j,5])))
        print(abs(as.numeric(data_Scaffolds[j,6])-as.numeric(data_Scaffolds[j,5]))>1000000)
        if (abs(as.numeric(data_Scaffolds[j,6])-as.numeric(data_Scaffolds[j,5]))>1000000){
          circos.genomicLink(data_Scaffolds[j,1:3], data_Scaffolds[j,4:6], col = inv_color, border = NA)
        } else {
          circos.genomicLink(data_Scaffolds[j,1:3], data_Scaffolds[j,4:6], col = inv_color, border = NA)
        }
        
      } else {
        circos.genomicLink(data_Scaffolds[j,1:3], data_Scaffolds[j,4:6], col = add.alpha("#D3D3D3",0.5), border = NA)
      }
    }
  }
}


#color_synteny(oldMN47_NT1, cytoband_NT1_uniq, colors[5])
color_synteny_first(oldMN47_NT1, cytoband_oldMN47_uniq, colors[5])
color_synteny_first(oldMN47_Asue, cytoband_oldMN47_uniq, colors[5])

##########Fig 1B: two way - NT1 newMN47

newNT1_newMN47<-read.table('NT1_030222/synteny.newMN47_newNT1_nonmerge.txt', header=T)
cytoband_NT1 = read.table("newMN47/NT1_cytobands.txt", colClasses = c("character", "numeric", "numeric"), sep = " ")
cytoband_NT1$V1=gsub("scaffold_", "newNT1_scaffold_", cytoband_NT1$V1)
cytoband_NT1_uniq<-cytoband_NT1[which(cytoband_NT1$V1 %in% unique(newNT1_newMN47$newNT1_chr)),]
cytoband_NT1_uniq_rev<-cytoband_NT1_uniq[dim(cytoband_NT1_uniq)[1]:1,]
cytoband_newMN47 = read.table("newMN47/newMN47_cytobands.txt", colClasses = c("character", "numeric", "numeric"), sep = " ")
cytoband_newMN47$V1=gsub("scaffold_", "newMN47_scaffold_", cytoband_newMN47$V1)
cytoband_newMN47$V1=gsub("scafold_", "newMN47_scaffold_", cytoband_newMN47$V1)
cytoband_newMN47_uniq<-cytoband_newMN47[which(cytoband_newMN47$V1 %in% unique(newNT1_newMN47$newMN47_chr)),]
cytoband_newMN47_uniq_rev<-cytoband_newMN47_uniq[dim(cytoband_newMN47_uniq)[1]:1,] #reverse cytobands of newMN47 as well, so that they are shown counter clock wise

cytoband <- rbindlist(list( cytoband_newMN47_uniq_rev, cytoband_NT1_uniq), use.names=T)
cytoband[[1]] = factor(cytoband[[1]], levels = cytoband$V1)
#step 4: make the figure

facing <- c( rep("bending.outside", 8), rep("bending.inside", 8))
padding_y <-c(rep(0.15, 8), rep(0.15, 8))
sc_labels <- c(paste("Chr", 8:1, sep=''), paste("Chr", 1:8, sep='') )

#iniate circos cirle and cytobands
colors <-brewer.pal(5, "RdYlBu")
colors <-brewer.pal(5, "BrBG")
#suecica, NT1, fill, newMN47, old MN47
circos.clear()

circos.par("start.degree" = 0,canvas.ylim=c(-1.1,1.1),canvas.xlim=c(-1.1,1.1),clock.wise=T,
           track.margin=c(0,0),track.height=0.1, "gap.degree"= c(rep(2,7),15,rep(2,7),15))
circos.genomicInitialize(cytoband, plotType = NULL,tickLabelsStartFromZero = TRUE,
                         axis.labels.cex = 0.3*par("cex"),track.height=0.1)
#circos.track(ylim = c(0, 0.05))
circos.track(ylim = c(0, 0.05),
             bg.col = c(rep(colors[4],8), rep(colors[2],8)),
             bg.border = NA, track.height = 0.05, panel.fun = function(x, y) {
               chr =get.cell.meta.data("sector.index")
               xlim =get.cell.meta.data("xlim")
               ylim =get.cell.meta.data("ylim")
               idx=match(get.cell.meta.data("sector.index"), cytoband$V1)
               circos.text(mean(xlim), mean(ylim)+padding_y[idx], sc_labels[idx], cex = 0.9, facing = facing[idx], niceFacing = F)
               
               #circos.text( sc_labels[idx], cex = 0.6, niceFacing = TRUE)
               #circos.axis(major.tick=TRUE, minor.ticks = 0,labels = c("", "10MB","20MB", "30MB"),labels.cex = 0.5)
               
             })


#color_synteny(newNT1_newMN47, cytoband_NT1_uniq_rev, colors[4])
color_synteny_first(newNT1_newMN47, cytoband_newMN47, colors[4])

pdf("testfig_threeway2.pdf")
dev.off()

pdf('legend.pdf')
colors <-brewer.pal(5, "RdYlBu")
plot(1,1)
legend(1, 1 ,legend = c('A. lyrata - NT1','A. lyrata - MN47 updated', 'A. suecica', 'A. lyrata - MN47 v1' ),
       col=c(colors[2], colors[4], colors[1], colors [5]), lty=1, lwd = 7, cex=0.6)
legtext=c('Colinear','A. suecica specific inversion', 'original MN47 specific inversion', 'updated MN47 specific inversion' , 'NT1 specific inversion')
xcoords <- c(0, 0.07, 0.30, 0.50, 0.70)
secondvector <- (1:length(legtext))-1
textwidths <- xcoords/secondvector # this works for all but the first element
textwidths[1] <- 0 # 

legend(1,1 ,legend = legtext, title= 'Syntenic blocks',
       col=c( '#D3D3D3', colors[1], colors[5], colors[4], colors [2]), cex=0.8, horiz=F, border="white", pt.cex=1.5, pch=15)

dev.off()
