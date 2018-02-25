anno$len=anno$End-anno$Start

promoters=read.table("Promoters_final.txt",sep="\t",header=T)
terminator=read.table("4tts.txt",sep="\t",header=T)


#######################
##calculate promoters
#######################
promoters$Start=promoters$TSS-1000
promoters$End=promoters$TSS+1000
prom=promoters
colnames(prom)[1]='Chr'


snp_custom_annotation(snps,prom)

minus=subset(snp_annotated_list,STRAND=="-")
plus=subset(snp_annotated_list,STRAND=="+")

#обрабатываем плюс-цепь
temp=plus
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
temp_graph_df_plus=data.frame()
temp_df_plus=data.frame()
for (i in 0:2000){
temp_df=subset(temp,V1==i)
temp_df[1,2]=sum(temp_df$V2)
temp_df=temp_df[1,]
temp_graph_df_plus=rbind(temp_graph_df_plus,temp_df)}
plot(temp_graph_df_plus)
temp_graph_df_plus$V1=temp_graph_df_plus$V1-1000
colnames(temp_graph_df_plus)=c('snp_position','plus')


#обрабатываем минус-цепь
temp=minus
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
#x=data.frame()
temp_graph_df_minus=data.frame()
temp_df=data.frame()
for (i in 0:2000){
temp_df=subset(temp,V1==i)
temp_df[1,2]=sum(temp_df$V2)
temp_df=temp_df[1,]
temp_graph_df_minus=rbind(temp_graph_df_minus,temp_df)}
plot(temp_graph_df_minus)
temp_graph_df_minus$V2=rev(temp_graph_df_minus$V2) #reverse minus-strand
temp_graph_df_minus$V1=temp_graph_df_minus$V1-1000
colnames(temp_graph_df_minus)=c("V1","V2")


graph_df=temp_graph_df_plus
graph_df[,2]=graph_df[,2]+temp_graph_df_minus[,2]
graph_df[,2]=graph_df[,2]/20367*1000

promoter_plot=graph_df #create plot for transcription promoter



#######################
##calculate terminators
#######################
plus=subset(terminator,STRAND=="+")
minus=subset(terminator,STRAND=="-")
minus$Start=minus$FRAG_START-1000
minus$End=minus$FRAG_START+1000
plus$Start=plus$FRAG_STOP-1000
plus$End=plus$FRAG_STOP+1000
ter=rbind(plus,minus)
colnames(ter)[1]='Chr'


snp_custom_annotation(snps,ter)

minus=subset(snp_annotated_list,STRAND=="-")
plus=subset(snp_annotated_list,STRAND=="+")

#обрабатываем плюс-цепь
temp=plus
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
temp_graph_df_plus=data.frame()
temp_df_plus=data.frame()
for (i in 0:2000){
temp_df=subset(temp,V1==i)
temp_df[1,2]=sum(temp_df$V2)
temp_df=temp_df[1,]
temp_graph_df_plus=rbind(temp_graph_df_plus,temp_df)}
plot(temp_graph_df_plus)
temp_graph_df_plus$V1=temp_graph_df_plus$V1-1000
colnames(temp_graph_df_plus)=c('snp_position','plus')


#обрабатываем минус-цепь
temp=minus
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
#x=data.frame()
temp_graph_df_minus=data.frame()
temp_df=data.frame()
for (i in 0:2000){
temp_df=subset(temp,V1==i)
temp_df[1,2]=sum(temp_df$V2)
temp_df=temp_df[1,]
temp_graph_df_minus=rbind(temp_graph_df_minus,temp_df)}
plot(temp_graph_df_minus)
temp_graph_df_minus$V2=rev(temp_graph_df_minus$V2) #reverse minus-strand
temp_graph_df_minus$V1=temp_graph_df_minus$V1-1000
colnames(temp_graph_df_minus)=c("V1","V2")


graph_df=temp_graph_df_plus
graph_df[,2]=graph_df[,2]+temp_graph_df_minus[,2]
graph_df[,2]=graph_df[,2]/20367*1000

terminator_plot=graph_df #create plot for transcription terminator

#############
##Build a graph
#############


#tiff("prom_trancr_ter_coverage_filtered1.tiff",height = 18, width = 12, units = 'in',res=600)
dev.off()
par(mfrow=c(1,2))
par(mar = c(5,4,4,4) + 0.1)
plot(promoter_plot,ylab="Number of SNPs per 1000 genes",
xlab="Distance from TSS")
y=(density(subset(subset(anno,len<1000),Frag=='five_prime_UTR')$len))$y*2000+10
x=(density(subset(subset(anno,len<1000),Frag=='five_prime_UTR')$len))$x
x=x[19:500]
y=y[19:500]
lines(x,y,col='blue',lwd=3)
axis(4,at=seq(0,0.005,by=0.001)*1300+3,label=seq(0,0.005,by=0.001))
mtext("Distribution of 5'UTR lengths", side=4, line=2)
lines(predict(smooth.spline(tss)),col='red',lwd=3)
legend("topright", c("Average number of SNPs","Distribution of 5'UTR lengths"),fill=c("red","blue"))
mtext("A)", adj=0, line=1, font=2,cex=1.5)

plot(terminator_plot,ylab="Number of SNPs per 1000 genes"
,xlab="Distance from TTS")
y=(density(subset(subset(anno,len<1000),Frag=='three_prime_UTR')$len))$y*2500+10
x=(density(subset(subset(anno,len<1000),Frag=='three_prime_UTR')$len))$x
x=x[28:500]
y=y[28:500]
lines(predict(smooth.spline(terminator_plot,df=10)),col='red',lwd=3)
lines(predict(smooth.spline(x*(-1),y,df=7)),col='blue',lwd=3)
axis(4,at=seq(0,0.005,by=0.001)*1300+3,label=seq(0,0.005,by=0.001))
mtext("Distribution of 3'UTR lengths", side=4, line=2)
legend("topright", c("Average number of SNPs","Distribution of 3'UTR lengths"),fill=c("red","blue"))
mtext("B)", adj=0, line=1, font=2,cex=1.5)
