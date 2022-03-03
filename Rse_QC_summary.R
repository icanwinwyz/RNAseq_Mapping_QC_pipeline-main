#rm(list=ls())
#setwd("/Users/joewang/")
#table1<-read.delim("~/01-C1-48h_S1.final_results_summary.txt", check.names=F,row.names=NULL, quote="", stringsAsFactors=FALSE)
#table2<-read.delim("~/02-C2-48h_S2.final_results_summary.txt", check.names=F,row.names=NULL, quote="", stringsAsFactors=FALSE)

library(reshape2)
library(ggplot2)
library(plyr)
library(scales)

args=commandArgs(TRUE)
name=args[length(args)]
data<-data.frame()
for (i in 1:(length(args)-1)){
	  data<-rbind(data,read.delim(args[i], check.names=F,row.names=NULL, quote="", stringsAsFactors=FALSE))
}

#print(dat)
#print(length(args))
#data<-merge(table1,table2,all=TRUE)

if(length(grep("read1",colnames(data),ignore.case=T))==1){
haha<-data.frame(samples=data[,grep("Sample",colnames(data))],data[,c(grep("CDS",colnames(data)),grep("UTR",colnames(data)),grep("intron",colnames(data)))])

try<-data.frame(samples=haha[,1],sapply(haha[,2:4], FUN = function(x) as.numeric(gsub("%", "", x))))
try$intergenic<-100-try[,2]-try[,3]-try[,4]
try_format<-ddply(melt(try, id.vars = 'samples'), .(samples), mutate, prop = value / sum(value))
colnames(try_format)[2:4]<-c("Features","value","Percentage")
write.table(data,paste(name,"_QC.txt",sep=""),col=NA,sep="\t",quote = F)


million_label<-function(x){
	paste(x/1000000,"M",sep="")
}


cal<-data.frame(samples = data[,1],data[,5:6],data[,8:9],unmapped=(data[,2]+data[,3]-data[,5]-data[,6]-data[,8]-data[,9]),check.names = F)
cal_format<-ddply(melt(cal, id.vars = 'samples'), .(samples), mutate,prop = value / sum(value))

b<-c()
b<-c(b,grep("#_raw_Read1",colnames(data)))
b<-c(b,grep("#_raw_Read2",colnames(data)))
total<-data.frame(samples=data[,1],data[,b],check.names = F)
total<-ddply(melt(total, id.vars = 'samples'), .(samples), mutate)


b<-c()
b<-grep("*RNA",colnames(data))
b<-c(b,grep("ERCC",colnames(data)))
junk_RNA<-data.frame(samples=data[,1],data[,b])
junk_RNA<-data.frame(samples=junk_RNA[,1],sapply(junk_RNA[,2:ncol(junk_RNA)], FUN = function(x) as.numeric(gsub("%", "", x))))
junk_RNA$other<-100-rowSums(junk_RNA[,2:ncol(junk_RNA)])
junk_RNA_format<-ddply(melt(junk_RNA, id.vars = 'samples'), .(samples), mutate, prop = value / sum(value))
colnames(junk_RNA_format)[2:4]<-c("RNA","value","Percentage")

#print(total)
#print(cal_format)
#print(try_format)
#print(junk_RNA_format)

pdf_name=paste(name,"_QC.pdf",sep="")

pdf(file=pdf_name)
print(ggplot(total,aes(x=samples,y=value,fill=variable))+geom_bar(position="dodge",stat="identity",width=0.6)+scale_y_continuous(label=million_label)+theme_bw()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)))
print(ggplot(cal_format, aes(x = samples, y = prop, fill = variable))+ theme_bw() +theme(panel.grid.major = element_blank())+geom_bar(stat = 'identity',width = 0.5)+scale_y_continuous(labels = percent_format())+coord_equal(1/0.2)+scale_fill_manual(values=c("royalblue1","blue","orange1","orange3","indianred1"))+ coord_flip()+xlab("")+ylab("")+guides(fill=guide_legend(title=NULL)))
print(ggplot(try_format, aes(x = samples, y = Percentage, fill = Features))+ theme_bw() +theme(panel.grid.major = element_blank())+geom_bar(stat = 'identity',width = 0.5)+scale_y_continuous(labels = percent_format())+coord_equal(1/0.2)+scale_fill_manual(values=c("darkgoldenrod1","gray45","lightskyblue","orchid"))+ coord_flip()+xlab("")+ylab("")+guides(fill=guide_legend(title=NULL)))
print(ggplot(junk_RNA_format, aes(x = samples, y = Percentage, fill = RNA))+ theme_bw() +theme(panel.grid.major = element_blank())+geom_bar(stat = 'identity',width = 0.6)+scale_y_continuous(labels = percent_format())+coord_equal(1/0.2)+scale_fill_manual(values=c("lightslategray","navajowhite4","lightskyblue","mistyrose3","olivedrab3"))+ coord_flip()) 

dev.off()

}else{
haha<-data.frame(samples=data[,grep("Sample",colnames(data))],data[,c(grep("CDS",colnames(data)),grep("UTR",colnames(data)),grep("intron",colnames(data)))])
try<-data.frame(samples=haha[,1],sapply(haha[,2:4], FUN = function(x) as.numeric(gsub("%", "", x))))
try$intergenic<-100-try[,2]-try[,3]-try[,4]
try_format<-ddply(melt(try, id.vars = 'samples'), .(samples), mutate, prop = value / sum(value))
colnames(try_format)[2:4]<-c("Features","value","Percentage")
#write.table(data,"all_samples_stat.txt",col=NA,sep="\t",quote = F)
write.table(data,paste(name,"_QC.txt",sep=""),col=NA,sep="\t",quote = F)


cal<-data.frame(samples = data[,1],unique=data[,3],multi=data[,5],unmapped=data[,2]-data[,3]-data[,5],check.names = F)
colnames(cal)<-c("samples","unique_mapped","multiple_mapped","unmapped")
cal_total<-data.frame(samples=data[,1],total=data[,2],check.names = F)
cal_format<-ddply(melt(cal, id.vars = 'samples'), .(samples), mutate, prop = value / sum(value))

million_label<-function(x){
	paste(x/1000000,"M",sep="")
}


b<-c()
b<-grep("*RNA",colnames(data))
b<-c(b,grep("ERCC",colnames(data)))
junk_RNA<-data.frame(samples=data[,1],data[,b])
junk_RNA<-data.frame(samples=junk_RNA[,1],sapply(junk_RNA[,2:ncol(junk_RNA)], FUN = function(x) as.numeric(gsub("%", "", x))))
junk_RNA$other<-100-rowSums(junk_RNA[,2:ncol(junk_RNA)])
junk_RNA_format<-ddply(melt(junk_RNA, id.vars = 'samples'), .(samples), mutate, prop = value / sum(value))
colnames(junk_RNA_format)[2:4]<-c("RNA","value","Percentage")


pdf_name=paste(name,"_QC.pdf",sep="")

pdf(file=pdf_name)

#pdf("all_samples_reads_dist_single_end.pdf",16,13)
print(ggplot(cal_total,aes(x=samples,y=total))+geom_bar(stat="identity",fill="lightblue")+scale_y_continuous(label=million_label)+theme_bw()+xlab("samples")+ylab("total#reads")+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)))
print(ggplot(cal_format, aes(x = samples, y = prop, fill = variable))+ theme_bw() +theme(panel.grid.major = element_blank())+geom_bar(stat = 'identity',width = 0.5)+scale_y_continuous(labels = percent_format())+coord_equal(1/0.2)+scale_fill_manual(values=c("royalblue1","orange","indianred1"))+ coord_flip()+xlab("")+ylab("")+guides(fill=guide_legend(title=NULL)))
print(ggplot(try_format, aes(x = samples, y = Percentage, fill = Features))+ theme_bw() +theme(panel.grid.major = element_blank())+geom_bar(stat = 'identity',width = 0.5)+scale_y_continuous(labels = percent_format())+coord_equal(1/0.2)+scale_fill_manual(values=c("darkgoldenrod1","gray45","lightskyblue","orchid"))+ coord_flip()+xlab("")+ylab("")+guides(fill=guide_legend(title=NULL)))
print(ggplot(junk_RNA_format, aes(x = samples, y = Percentage, fill = RNA))+ theme_bw() +theme(panel.grid.major = element_blank())+geom_bar(stat = 'identity',width = 0.6)+scale_y_continuous(labels = percent_format())+coord_equal(1/0.2)+scale_fill_manual(values=c("lightslategray","navajowhite4","lightskyblue","mistyrose3","olivedrab3"))+ coord_flip()) 
dev.off()
}
