snps=read.table('snp_strict_filter.txt',header=T)
colnames(snps)[c(2,3)]=c('snp_position','snp_chrom')

#подгружаем дату, которую сгенерила Таня
atg=read.table("data/ATGs.txt",sep="\t",header=T)
atg$Start=atg$ATG-300
atg$End=atg$ATG+300
atg$X=NULL
colnames(atg)[1]='Chr'

snp_custom_annotation(snps,atg)
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
for (i in 0:600){
temp_df=subset(temp,V1==i)
temp_df[1,2]=sum(temp_df$V2)
temp_df=temp_df[1,]
temp_graph_df_plus=rbind(temp_graph_df_plus,temp_df)}
plot(temp_graph_df_plus)
temp_graph_df_plus$V1=temp_graph_df_plus$V1-300
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
for (i in 0:600){
temp_df=subset(temp,V1==i)
temp_df[1,2]=sum(temp_df$V2)
temp_df=temp_df[1,]
temp_graph_df_minus=rbind(temp_graph_df_minus,temp_df)}
plot(temp_graph_df_minus)
temp_graph_df_minus$V1=temp_graph_df_minus$V1-300
temp_graph_df_minus$V2=rev(temp_graph_df_minus$V2) #reverse minus-strand
colnames(temp_graph_df_minus)=c("V1","minus")


graph_df=temp_graph_df_plus
graph_df[,2]=graph_df[,2]+temp_graph_df_minus[,2]
graph_df[,2]=graph_df[,2]/20367*1000
plot(graph_df)
#pdf("ATG_normalized_fin.pdf")
plot(graph_df,main="Number of SNPs near ATG",xlab="Distance from ATG",
ylab="Number of SNPs per 1000 genes")
#dev.off()
nucl_count=c(1,2,3)
plot_df_final=as.data.frame(cbind(graph_df$snp_position,graph_df$plus,nucl_count))
colnames(plot_df_final)=c('V1','V2','V3')

library(scales)
plot(type='n',plot_df_final$V1,plot_df_final$V2,xlab="Distanse from ATG",ylab="Number of SNPs per 1000 genes",main="Number of SNPs near ATG")
curren_nucl=subset(plot_df_final,V3=="2")
points(curren_nucl$V1,curren_nucl$V2,col=alpha('red', 0.5),pch=20)
lines(curren_nucl$V1,curren_nucl$V2,col=alpha("red", 0.5))
curren_nucl=subset(plot_df_final,V3=="3")
lines(curren_nucl$V1,curren_nucl$V2,col=alpha('darkgreen', 0.5))
points(curren_nucl$V1,curren_nucl$V2,col=alpha('darkgreen', 0.5),pch=20)
curren_nucl=subset(plot_df_final,V3=="1")
lines(curren_nucl$V1,curren_nucl$V2,col=alpha('blue', 0.3))
points(curren_nucl$V1,curren_nucl$V2,col=alpha('blue', 0.3),pch=20)
legend("bottomright", c("1st base","2nd base","3st base"),fill=c("red","blue","darkgreen"))