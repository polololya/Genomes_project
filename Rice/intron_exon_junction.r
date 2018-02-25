colnames(snp)=c("Chr_snp","global","coords","from","to")

plus=subset(anno,Strand=="+")
minus=subset(anno,Strand=="-")

exoin=subset(plus,Frag=='exon')
exoin$Start=exoin$End-100
exoin$End=exoin$End+100
inexo=subset(plus,Frag=='exon')
inexo$End=inexo$Start+100
inexo$Start=inexo$Start-100

snp_custom_annotation(snps,exoin)

temp=snp_annotated_list
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
temp_graph_df_exoin=data.frame()
temp_df_plus=data.frame()
for (i in 0:200){
temp_df=subset(temp,V1==i)
temp_df[1,2]=sum(temp_df$V2)
temp_df=temp_df[1,]
temp_graph_df_exoin=rbind(temp_graph_df_exoin,temp_df)}
plot(temp_graph_df_exoin)
temp_graph_df_exoin$V1=temp_graph_df_exoin$V1-100
colnames(temp_graph_df_exoin)=c('V1','V2')


snp_custom_annotation(snps,inexo)

temp=snp_annotated_list
temp=cbind(temp$Start,temp$snp_position)
temp=as.data.frame(temp)
temp$V1=temp$V2-temp$V1
temp$V2=1
temp=na.omit(temp)
temp=temp[order(temp),]
temp_graph_df_inexo=data.frame()
temp_df_plus=data.frame()
for (i in 0:200){
temp_df=subset(temp,V1==i)
temp_df[1,2]=sum(temp_df$V2)
temp_df=temp_df[1,]
temp_graph_df_inexo=rbind(temp_graph_df_inexo,temp_df)}
plot(temp_graph_df_inexo)
temp_graph_df_inexo$V1=temp_graph_df_inexo$V1-100
colnames(temp_graph_df_inexo)=c('V1','V2')

temp_graph_df_exoin$V3=c(1,2,3)
temp_graph_df_inexo$V3=c(1,2,3)

#temp_graph_df_exoin=x
#temp_graph_df_exoin$V2=temp_graph_df_exoin$V2*100
#temp_graph_df_exoin$V1=temp_graph_df_exoin$V1-100
#temp_graph_df_exoin$V2=temp_graph_df_exoin$V2/10339
#temp_graph_df_exoin$V3=c(1,2,3)
#temp_graph_df_exoin=temp_graph_df_inexo

library(scales)
#tiff("splice_combined_1.tiff",height = 7, width = 14, units = 'in',res=600)
par(mfrow=c(1,2))
plot(type='n',temp_graph_df_exoin$V1,temp_graph_df_exoin$V2,xlab="Distanse from splicing site",
ylab="Number of SNPs per 1000 genes",
main="Number of SNPs near splicing site exon -> intron")
x2=subset(temp_graph_df_exoin,V3=="1")
points(x2$V1,x2$V2,col=alpha('red', 0.5),pch=20)
lines(x2$V1,x2$V2,col=alpha("red", 0.5))
x2=subset(temp_graph_df_exoin,V3=="2")
lines(x2$V1,x2$V2,col=alpha('darkgreen', 0.5))
points(x2$V1,x2$V2,col=alpha('darkgreen', 0.5),pch=20)
x2=subset(temp_graph_df_exoin,V3=="3")
lines(x2$V1,x2$V2,col=alpha('blue', 0.3))
points(x2$V1,x2$V2,col=alpha('blue', 0.3),pch=20)
legend("topright", c("1st base","2st base","3st base"),fill=c("red","blue","darkgreen"))
#dev.off()

#tiff("splice_exo2intro.tiff",height = 7, width = 7, units = 'in',res=600)
plot(type='n',temp_graph_df_inexo$V1,temp_graph_df_inexo$V2,xlab="Distanse from splicing site",
ylab=NA,
main="Number of SNPs near splicing site intron -> exon")
x2=subset(temp_graph_df_inexo,V3=="3")
points(x2$V1,x2$V2,col=alpha('red', 0.5),pch=20)
lines(x2$V1,x2$V2,col=alpha("red", 0.5))
x2=subset(temp_graph_df_inexo,V3=="1")
lines(x2$V1,x2$V2,col=alpha('darkgreen', 0.5))
points(x2$V1,x2$V2,col=alpha('darkgreen', 0.5),pch=20)
x2=subset(temp_graph_df_inexo,V3=="2")
lines(x2$V1,x2$V2,col=alpha('blue', 0.3))
points(x2$V1,x2$V2,col=alpha('blue', 0.3),pch=20)
legend("bottomright", c("1st base","2st base","3st base"),fill=c("red","blue","darkgreen"))
#dev.off()
