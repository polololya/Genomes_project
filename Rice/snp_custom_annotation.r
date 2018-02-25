#Обязательно должен быть установлен пакет R data.table
#На вход 2 дата фрейма. data- данные со снипами, обязательны колонки snp_chrom- хромосома снипа, snp_position- позиция снипа на хромосоме
#annotation- любая кастомная аннотация, обязательны колонки chromosome (которая содержит данные по хромосоме в форме Chr#),
#Start- начало локуса, End- конец локуса. 
#в переменной chrom_vector (ex chrom_vector=c(3,6:8), будут обсчитаны 1,3,6,7,8 хромосомы; chrom_vector=0, будет обсчитана только 1-я хромосома) можно выбрать, какие хромосомы, помимо первой, будут обсчитаны. Первая хромосома будет обсчитана в любом случае.
#переменная enable_XY (обозначения в аннотации дожны быть ChrX, ChrY, в снипах X, Y) и enable_MT (в аннотации обозначения должны быть ChrMT, в снипах- MT)могут принимать значение TRUE, в этом случае
#в отчет будут добавлены снипы, локализованные на половых и митохондриальной хромосомах
#
#
#
snp_custom_annotation=function(data,annotation,chrom_vector=NULL,enable_XY=NULL,enable_MT=NULL){
require(data.table)
temp_vec=vector()
if(is.null(chrom_vector)){ #Так делать, конечно, плохо, но ничего лучше не придумал. Тащим из аннотации только циферные значения хромосом
options(warn=-1)
temp_vec = as.numeric(unique(data$snp_chrom))
temp_vec  = temp_vec [!is.na(temp_vec)]
chrom_vector  = temp_vec[-1]
options(warn=0)
}
chr_annotation=subset(annotation,Chr=="Chr1") #Аннотируем первую хромосому, как тестовую, т.к. она есть у всех
data_chr_current=subset(data,snp_chrom==1)
print(paste('working on 1 chromosome of ',length(chrom_vector)+1,sep=''))
setDT(data_chr_current)
setDT(chr_annotation)
setkeyv(chr_annotation,c("Start","End"))
data_chr_current[,V11:=snp_position]
temp <- foverlaps(
  x=data_chr_current,y=chr_annotation,
  by.x=c("V11","snp_position"),
  type="within")
temp[,V11:=NULL]
chr_out=na.omit(temp)
if (length(chrom_vector)>0){
for (i in chrom_vector){ #Аннотируем остальные хромосомы
	chr_annotation=subset(annotation,Chr==paste('Chr',i,sep=''))
		data_chr_current=subset(data,snp_chrom==i)
		setDT(data_chr_current)
	print(paste('working on ',i,' chromosome of ',length(chrom_vector)+1,sep=''))
		setDT(chr_annotation)
		setkeyv(chr_annotation,c("Start","End"))
		data_chr_current[,V11:=snp_position]
		temp <- foverlaps(
			x=data_chr_current,y=chr_annotation,
			by.x=c("V11","snp_position"),
			type="within")
		temp[,V11:=NULL]
		temp=na.omit(temp)
	chr_out=rbind(chr_out,temp)}
	}
	
if(is.null(enable_XY)==FALSE){ #Если поставили enable XY=TRUE- теперь аннотируем половые хромосомы. Сначала X
	chr_annotation=subset(annotation,Chr=="ChrX")
	print('working on X chromosome')
		data_chr_current=subset(data,snp_chrom=="X")
		setDT(data_chr_current)
		setDT(chr_annotation)
		setkeyv(chr_annotation,c("Start","End"))
		data_chr_current[,V11:=snp_position]
		temp <- foverlaps(
			x=data_chr_current,y=chr_annotation,
			by.x=c("V11","snp_position"),
			type="within")
		temp[,V11:=NULL]
		temp=na.omit(temp)
	chr_out=rbind(chr_out,temp)

	chr_annotation=subset(annotation,Chr=="ChrY") #Теперь Y
	print('working on Y chromosome')
		data_chr_current=subset(data,snp_chrom=="Y")
		setDT(data_chr_current)
		setDT(chr_annotation)
		setkeyv(chr_annotation,c("Start","End"))
		data_chr_current[,V11:=snp_position]
		temp <- foverlaps(
			x=data_chr_current,y=chr_annotation,
			by.x=c("V11","snp_position"),
			type="within")
		temp[,V11:=NULL]
		temp=na.omit(temp)
	chr_out=rbind(chr_out,temp)}
if(is.null(enable_MT)==FALSE){ #если enable_MT=TRUE, аннотируем митохондриальную хромосому.
	print('working on MT chromosome, almost done')
	chr_annotation=subset(annotation,Chr=="ChrMT")
		data_chr_current=subset(data,snp_chrom=="MT")
		setDT(data_chr_current)
		setDT(chr_annotation)
		setkeyv(chr_annotation,c("Start","End"))
		data_chr_current[,V11:=snp_position]
		temp <- foverlaps(
			x=data_chr_current,y=chr_annotation,
			by.x=c("V11","snp_position"),
			type="within")
		temp[,V11:=NULL]
		temp=na.omit(temp)
	chr_out=rbind(chr_out,temp)}
snp_annotated_list<<-unique(na.omit(rbind(chr_out,temp))) #Вывод и уборка за функцией
rm(temp,chr_out,chrom_vector,data_chr_current,chr_annotation,temp_vec)
print(paste('finished annotating ',dim(snp_annotated_list)[1],' SNPs. Data contained in data.frame called snp_annotated_list',sep=''))
}

