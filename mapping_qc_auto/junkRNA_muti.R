args=commandArgs(TRUE)
total=as.numeric(args[1])
name=args[2]
#path=paste("./",name,"/",name,"_temp/junk_RNA_tempfile.txt",sep="")
#path=paste("./",name,"/junk_RNA_tempfile.txt",sep="")

#print(total)
#print(name)
#print(path)

data<-read.table("junk_RNA_tempfile.txt")
data<-data[,-3]
#data
for(i in 1:nrow(data)){
total<-total-data[i,][2]
}
final<-rbind(data,data.frame(V1="Other_polyA_RNA",V2=as.vector(total[1,])))
write.table(final,paste(name,"junk_RNA_muti.txt",sep="."),row.names=F,col.names=F,quote=F)


