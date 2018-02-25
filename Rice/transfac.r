data=read.table('./data/transfac.txt',header=T)


#tiff("transfac_ultimate.tiff",height = 7, width = 7, units = 'in',res=600)
plot(data,
main="Transcription Factor Binding Sites (TRANSFAC database)",
ylab='Fraction of promoters with TFBS',
xlab='Distance from TSS, nt')
lines(x,y,col='blue',lwd=3)
lines(predict(smooth.spline(data,df=10)),col='red',lwd=3)
abline(v=0)
#dev.off()
