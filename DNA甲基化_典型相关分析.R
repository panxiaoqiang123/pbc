dat<-read.csv("????????????.csv")
id<-read.csv("?????????gene_id.csv")
gid<-id$genename
gid<-unique(gid)
gid<-as.data.frame(gid)
result<-array(0,c((nrow(gid)*nrow(gid)),3))
count<-1
for(i in c(1:nrow(gid)))
{ for(j in c(1:nrow(gid)))
{
  id1<-id[id$genename==as.character(gid[i,1]),]
  dat1<-semi_join(dat,id1,by=c("X"="id"))
  dat1<-dat1[,-1]
  id2<-id[id$genename==as.character(gid[j,1]),]
  dat2<-semi_join(dat,id2,by=c("X"="id"))
  dat2<-dat2[,-1]
  dat1<-as.matrix(dat1)
  dat2<-as.matrix(dat2)
  dat1<-t(dat1)
  dat2<-t(dat2)
  ca<-cancor(dat1,dat2)
  result[count,1]<-as.character(gid[i,1])
  result[count,2]<-as.character(gid[j,1])
  result[count,3]<-max(ca$cor)
  
  count<-count+1
}
  print(i)
  print(result[count-1,1])
  
}