library(RColorBrewer)
args= commandArgs(TRUE)
name=args[1]
data<-read.delim(name,header=T,check.names=F,sep="\t")
org<-data
data<-t(org)
haha=data.frame(data[2:nrow(data),],rownames(data)[2:nrow(data)])
colnames(haha)<-c("percentage","term")
haha[,1]<-as.vector(haha[,1])
a<-grep("%",haha[,1])
for (i in a){
          haha[i,1]<-as.numeric(gsub("%","",haha[i,1]))
}
other<-100-as.numeric(haha["CDS",][1])-as.numeric(haha["UTR",][1])-as.numeric(haha["intron",][1])
intergenic<-data.frame(percentage=other,term = "intergenic",row.names="intergenic")
#rd<-rbind(haha["CDS",],haha["UTR",],haha["intron",],intergenic,haha["rRNA",],haha["tRNA",],haha["MTRNA",],haha["ERCC",])
rd<-rbind(haha["CDS",],haha["UTR",],haha["intron",],intergenic)
rd[,1]<-as.numeric(rd[,1])
rd$term<-factor(rd$term,levels=c("CDS","UTR","intron","intergenic"))
total_mapping<-as.numeric(haha["total%_mapping",1][1])
total_mapping<-paste("Total Mapping: ",total_mapping,"%",sep="")
#colors<-c("darkgoldenrod1","gray75","lightskyblue","plum2","mediumorchid1","khaki1","cornflowerblue","mistyrose","palegreen1")
x<-as.numeric(rd[,1])
labels<-paste(rd$term,": ",rd$percentage,"%",sep="")
colors=brewer.pal(8,"Pastel2")

#junk_RNA<-rbind(haha["MTRNA",][1],haha["rRNA",][1],haha["tRNA",][1],haha["ERCC",][1])
#junk_RNA<-data.frame(junk_RNA,term=rownames(junk_RNA))
b<-c()
b<-grep("*RNA",rownames(haha))
b<-c(b,grep("ERCC",rownames(haha)))
junk_RNA<-haha[b,]
junk_RNA[,1]<-as.numeric(junk_RNA[,1])
other_RNA<-100-sum(junk_RNA[,1])
other_RNA<-data.frame(percentage=other_RNA,term="others")
#other_RNA<-c(100-sum(junk_RNA[,1]),"Others")
junk_RNA<-rbind(junk_RNA,others=other_RNA)
vals<-junk_RNA$percentage
names(vals)<-junk_RNA$term
#pie(x,labels = labels,col=colors[as.factor(labels)],main=total_mapping)
pie(x,labels = labels,col=colors[as.factor(labels)],main="Reads Distribution on different features",cex.main=1.5)
legend("topleft",total_mapping,text.col = "red",bty = "n",cex = 1.2)

mp<-barplot(vals, ylim = c(0, 100),width = 0.1, space=0.7,ylab="percentage(%)",main="Reads Distribution on different types of RNA",cex.main=1.5)
text(mp, vals, labels = paste(vals,"%"), pos = 3)
