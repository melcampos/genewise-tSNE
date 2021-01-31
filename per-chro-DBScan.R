library(RColorBrewer)
library(dbscan)
library(ggplot2)

# Load the main file
file <- read.table("tsne1-2_chro2.txt", sep="\t", header=T)
head(file)

# Subset to genes on chromosome 2L only. 
chro2L <- file[which(file$chro == '2L'),]
cols <- c("#e41a1c","#377eb8")
cols_t1<-cols[as.factor(chro2L$chro3)]

# Plot for Figure showing collinear and inverted regions on chromosome 2L.
plot(file$tsne1,file$tsne2,main="", col=alpha("#d9d9d9", 0.5), cex=0.6, pch=16,axes=F,xlab="",ylab="")
par(new=T)
plot(chro2L$tsne1, chro2L$tsne2,main="", col=alpha(cols_t1, 0.5), cex=0.5, pch=16,xlab="",ylab="",axes=F,
     xlim=range(file$tsne1),ylim=range(file$tsne2))
legend(x=-15,y=0,legend=paste(c( "collinear","a")),
       col=cols,pch=16,bty="n",ncol=1,cex=0.6,pt.cex=1,xpd=TRUE)

# Performing DBScan for collinear regions on chromosome 2L
chro2L <- file[which(file$chro3 == '2L'),]
dim(chro2L)
tsne <- chro2L[,c(3,4)]

dbs = dbscan(tsne, eps=0.6, minPts=15)
ggplot(tsne, aes(tsne1, tsne2)) + geom_point(size=1, aes(color=as.factor(dbs$cluster))) + scale_fill_manual(values=col_vector)

cols<-brewer.pal(n=8,name="Dark2")
cols <- c("darkgray", "#1B9E77", "#D95F02", "#7570B3", "#E7298A" ,"#66A61E" ,"#E6AB02" ,
          "#A6761D", "#666666")
plot(1:8,1:8,col=cols,cex=3,pch=16)
cols_t1<-cols[as.factor(dbs$cluster)]
cols_t1

plot(file$tsne1,file$tsne2,main="", col=alpha("#d9d9d9",0.5), cex=0.6, pch=16,axes=F,xlab="",ylab="")
par(new=T)
plot(tsne$tsne1, tsne$tsne2,main="", col=alpha(cols_t1,0.5), cex=0.6, pch=16,axes=F,ylab="",xla="",ylim=range(file$tsne2),xlim=range(file$tsne1))
legend(x=-15,y=-1,legend=paste(c("NA","2L-i","2L-ii","2L-iii","2L-iv","2L-v")),
       col=cols,pch=16,bty="n",ncol=1,cex=0.6,pt.cex=1,xpd=TRUE)


