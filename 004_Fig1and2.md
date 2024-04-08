### FIG1
##### \[5.\]lncRNA的数量统计
```r
awk '{print $4}' /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.ref_lncRNA_exon.bed|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.3Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.2Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.ref_lncRNA_exon.bed|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.3Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.2Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.ref_lncRNA_exon.bed|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.3Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.2Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.ref_lncRNA_exon.bed|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.3Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.2Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.ref_lncRNA_exon.bed|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.2Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.ref_lncRNA_exon.bed|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.2Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.ref_lncRNA_exon.bed|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.2Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l
awk '{print $4}' /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.ref_lncRNA_exon.bed|sort|uniq|wc -l
awk '{print $4}' "/home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.2Gnovel_lncRNA_exon.bed"|sort|uniq|wc -l

#简化结果
sp,source,num
Human,ENSEMBL,18864
Human,Novel,20597
Mouse,ENSEMBL,11527
Mouse,Novel,18911
Chicken,ENSEMBL,11707
Chicken,Novel,9726
Zebrafish,ENSEMBL,2202
Zebrafish,Novel,6800
Macaque,ENSEMBL,4596
Macaque,Novel,13076
Rat,ENSEMBL,2465
Rat,Novel,21844
Rabbit,ENSEMBL,3429
Rabbit,Novel,12878
Opossum,ENSEMBL,10371
Opossum,Novel,17435





#结果
sp,source,num
Human,ENSEMBL,18864
Human,Novel (from long-read RNA-seq),510
Human,Novel (from short-read RNA-seq),20087
Mouse,ENSEMBL,11527
Mouse,Novel (from long-read RNA-seq),420
Mouse,Novel (from short-read RNA-seq),18491





lncRNA_source <- read.csv("/home/zhzhang/PG/1Identification/SP.lncRNA_source.csv")
lncRNA_source$sp <- factor(lncRNA_source$sp,levels = rev(c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish")))
lncRNA_source$source <- factor(lncRNA_source$source,levels = c("Novel","ENSEMBL"))
#plot
p2 <- ggplot(data = lncRNA_source,aes(y=sp,x=num/1000))+
  geom_col(width = 0.8,alpha=0.7,aes(fill=source,color=source))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  theme_half_open()+
  labs(y =NULL , x ="No. of lncRNAs (×10E+3)",fill=NULL,color=NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,41))+
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = c(0.78, 0.1))+
  theme(legend.text = element_text(size = 15))
  
ggsave("/home/zhzhang/PG/1Identification/sp_lncRNAsource_num.png",
       p2,width =9, height =4.5,dpi=1200, units = "in", device='png',bg = "transparent")




#pglnc数量展示
#
pglncnum <- read.csv("~/PG/1Identification/SP.pglncnum.tj.txt", sep="")
pglncnum$sp <- factor(pglncnum$sp,levels = rev(c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish")))

p2 <- ggplot(data = pglncnum,aes(y=sp,x=num))+
  geom_point(aes(fill=sp,color=sp,size=num))+
  geom_col(width = 0.01,alpha=0.7,aes(fill=sp,color=sp))+
  geom_text(aes(label=num),color="white")+
  scale_size(range = c(9,20))+
  scale_x_continuous(expand = c(0,0),limits = c(0,1300))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  theme_half_open()+
  labs(y =NULL , x ="No. of pseudogene-derived lncRNAs",fill=NULL,color=NULL) +
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = "none")+
  theme(legend.text = element_text(size = 15))
ggsave("/home/zhzhang/PG/1Identification/sp_pglncnum.pdf",
       p2,width =6, height =5,dpi=1200, units = "in", device='pdf',bg = "transparent")


#pglnc来源展示
#pglnc source
splist1 <- c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")
splist2 <- c("Danio_rerio","Gallus_gallus","Monodelphis_domestica","Oryctolagus_cuniculus","Rattus_norvegicus","Mus_musculus","Macaca_mulatta","Homo_sapiens")
all <- data.frame()
for (i in 1:8) {
assign("pglnc",read.delim(paste("/home/zhzhang/PG/lncRNA_class/",splist2[i],".lncRNA_class.txt",sep=""))%>%
  filter(type=="Pseudogene-derived lncRNA")%>%
  mutate(type="ENSEMBL")%>%
  mutate(sp=splist1[i]))
rownovel <- c(grep(".",pglnc$geneid,ignore.case=F,fixed=T,value=F),
              grep("chr",pglnc$geneid,ignore.case=F,fixed=T,value=F),
              grep("XLOC",pglnc$geneid,ignore.case=F,fixed=T,value=F),
              grep("PB",pglnc$geneid,ignore.case=F,fixed=T,value=F))
pglnc$type[rownovel] <- "Novel"
gc()
all <- rbind(all,pglnc)
}
#
tj <- group_by(all,sp,type)%>%
  summarise(num=n())
tjsp <- group_by(all,sp)%>%
  summarise(all=n())
tjmerge <- left_join(tj,tjsp,by="sp")%>%
  mutate(percentage=num*100/all)
data.table::fwrite(tjmerge,
                   file ="/home/zhzhang/PG/1Identification/SP.pglncsource.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
tjmerge$sp <- factor(tjmerge$sp,levels = rev(c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish")))
tjmerge$type <- factor(tjmerge$type,levels = c("Novel","ENSEMBL"))
p2 <- ggplot(data = tjmerge,aes(y=sp,x=percentage))+
  geom_col(width = 0.6,alpha=0.7,aes(fill=type,color=type))+
  geom_text(data = filter(tjmerge,type=="ENSEMBL"),x=13,color="white",size=5,
            aes(label=paste(round(percentage,1),"%",sep = "")))+
  geom_text(data = filter(tjmerge,type=="Novel"),x=87,color="white",size=5,
            aes(label=paste(round(percentage,1),"%",sep = "")))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  scale_x_continuous(labels = c("0","25","50","75","100"),expand = c(0,0),limits = c(0,105))+
  theme_half_open()+
  labs(y =NULL , x ="Percentage (%)",fill=NULL,color=NULL) +
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 0),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = "none")+
  theme(legend.text = element_text(size = 15)) +
  theme(axis.ticks.y = element_line(size = 0)) + 
  theme(axis.ticks.y = element_line(colour = "white"))
ggsave("/home/zhzhang/PG/1Identification/sp_pglncsource.pdf",
       p2,width =3.3, height =5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```








##### \[6.\]假基因数量,密度，比例展示
```r
#假基因组成
hshpg <- read.delim("~/PG/pg_message/Homo_sapiens_hpg.txt")%>%
  select(geneid,type)
mmhpg <- read.delim("/home/zhzhang/PG/pg_message/Mus_musculus_hpg.txt")%>%
  select(geneid,type)
gghpg <- read.delim("/home/zhzhang/PG/pg_message/Gallus_gallus_hpg.txt")%>%
  select(geneid,type)
drhpg <- read.delim("/home/zhzhang/PG/pg_message/Danio_rerio_hpg.txt")%>%
  select(geneid,type)
machpg <- read.delim("/home/zhzhang/PG/pg_message/Macaca_mulatta_hpg.txt")%>%
  select(geneid,type)
rathpg <- read.delim("/home/zhzhang/PG/pg_message/Rattus_norvegicus_hpg.txt")%>%
  select(geneid,type)
ochpg <- read.delim("/home/zhzhang/PG/pg_message/Oryctolagus_cuniculus_hpg.txt")%>%
  select(geneid,type)
ophpg <- read.delim("/home/zhzhang/PG/pg_message/Monodelphis_domestica_hpg.txt")%>%
  select(geneid,type)
#合并
hpg <- rbind(mutate(hshpg,sp="Human"),
             mutate(machpg,sp="Macaque"),
             mutate(mmhpg,sp="Mouse"),
             mutate(rathpg,sp="Rat"),
             mutate(ochpg,sp="Rabbit"),
             mutate(ophpg,sp="Opossum"),
             mutate(gghpg,sp="Chicken"),
             mutate(drhpg,sp="Zebrafish")
             )
#tj
tj <- group_by(hpg,sp,type)%>%
  summarise(num=n())
tj$sp <- factor(tj$sp,levels = rev(c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish")))
tj$type <- factor(tj$type,levels = c("Duplicated","Processed","Fragment"))
data.table::fwrite(tj,
                   file ="/home/zhzhang/PG/1Identification/SP.pgnum.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
p1 <- ggplot(data = tj,aes(y=sp,x=num/1000))+
  geom_col(width = 0.8,alpha=1,aes(fill=type,color=type))+
  ggsci::scale_fill_jco()+
  ggsci::scale_color_jco()+
  theme_half_open()+
  labs(y =NULL , x ="No. of Pseudogenes (×10E+3)",fill=NULL,color=NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,21))+
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = c(0.75, 0.15))+
  theme(legend.text = element_text(size = 15))
ggsave("/home/zhzhang/PG/1Identification/sp_pgnum.png",
       p1,width = 9, height = 4.5,dpi=1200, units = "in", device='png',bg = "transparent")




#密度
#导入假基因数量和基因组大小
pgnum <- read.delim("~/PG/1Identification/SP.pgnum.tj.txt")
genome <- read.csv("~/PG/1Identification/sp_allgenome_size.txt", sep="")%>%
  mutate(gs=gs/1000000)
merge <- left_join(pgnum,genome,by="sp")%>%
  mutate(des=num/gs)%>%
  select(-num,-gs)
#
merge$sp <- factor(merge$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
merge$type <- factor(merge$type,levels = c("Duplicated","Processed","Fragment"))
data.table::fwrite(merge,
                   file ="/home/zhzhang/PG/1Identification/SP.pgdensity.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
p1 <- ggplot(data = merge,aes(x=sp,y=des))+
  geom_col(position = "dodge",width = 0.6,alpha=1,aes(fill=type),color="black")+
  ggsci::scale_fill_jco()+
  theme_half_open()+
  labs(y ="Pseudogene density\n(Pseudogene/Mb)" , x =NULL,fill=NULL,color=NULL) +
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = "top")+
  theme(legend.text = element_text(size = 15))
ggsave("/home/zhzhang/PG/1Identification/sp_pgdensity.pdf",
       p1,width = 9, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")




#三类假基因比例
allpgnum <- group_by(pgnum,sp)%>%
  summarise(all=sum(num))
pgpro <- left_join(pgnum,allpgnum,by="sp")%>%
  mutate(pro=num*100/all)
#
pgpro$sp <- factor(pgpro$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
pgpro$type <- factor(pgpro$type,levels = c("Duplicated","Processed","Fragment"))
#plot
p2 <- ggplot(data = pgpro,aes(x=sp,y=pro))+
  geom_col(position = "dodge",width = 0.6,alpha=1,aes(fill=type),color="black")+
  ggsci::scale_fill_jco()+
  theme_half_open()+
  labs(y ="Percentage (%)" , x =NULL,fill=NULL,color=NULL) +
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = "top")+
  theme(legend.text = element_text(size = 15))
ggsave("/home/zhzhang/PG/1Identification/sp_pgpercentage.pdf",
       p2,width = 9, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")











#假基因长度对比
hshpg <- read.delim("~/PG/pg_message/Homo_sapiens_hpg.txt")%>%
  mutate(len=len*3)%>%
  select(geneid,type,len)
mmhpg <- read.delim("/home/zhzhang/PG/pg_message/Mus_musculus_hpg.txt")%>%
  mutate(len=len*3)%>%
  select(geneid,type,len)
gghpg <- read.delim("/home/zhzhang/PG/pg_message/Gallus_gallus_hpg.txt")%>%
  mutate(len=len*3)%>%
  select(geneid,type,len)
drhpg <- read.delim("/home/zhzhang/PG/pg_message/Danio_rerio_hpg.txt")%>%
  mutate(len=len*3)%>%
  select(geneid,type,len)
machpg <- read.delim("/home/zhzhang/PG/pg_message/Macaca_mulatta_hpg.txt")%>%
  mutate(len=len*3)%>%
  select(geneid,type,len)
rathpg <- read.delim("/home/zhzhang/PG/pg_message/Rattus_norvegicus_hpg.txt")%>%
  mutate(len=len*3)%>%
  select(geneid,type,len)
ochpg <- read.delim("/home/zhzhang/PG/pg_message/Oryctolagus_cuniculus_hpg.txt")%>%
  mutate(len=len*3)%>%
  select(geneid,type,len)
ophpg <- read.delim("/home/zhzhang/PG/pg_message/Monodelphis_domestica_hpg.txt")%>%
  mutate(len=len*3)%>%
  select(geneid,type,len)
#合并
hpg <- rbind(mutate(hshpg,sp="Human"),
             mutate(machpg,sp="Macaque"),
             mutate(mmhpg,sp="Mouse"),
             mutate(rathpg,sp="Rat"),
             mutate(ochpg,sp="Rabbit"),
             mutate(ophpg,sp="Opossum"),
             mutate(gghpg,sp="Chicken"),
             mutate(drhpg,sp="Zebrafish")
)
#tj
tj <- group_by(hpg,sp,type)%>%
  summarise(median=median(len))
data.table::fwrite(tj,
                   file ="/home/zhzhang/PG/1Identification/SP.pglen.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
hpg$sp <- factor(hpg$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
hpg$type <- factor(hpg$type,levels = c("Duplicated","Processed","Fragment"))
p1 <- ggplot(data = filter(hpg,type!="Fragment"),aes(y=log10(len),x=sp))+
  geom_boxplot(width = 0.6,outlier.alpha = 0,notch = T,aes(fill=type))+
  geom_signif(annotations=c(rep("***",times=8)),
              y_position=c(5,5,5.15,5,4.4,5.1,4.7,5.3),tip_length = 0,
              xmin = c(0.85,1.85,2.85,3.85,4.85,5.85,6.85,7.85),
              xmax = c(1.15,2.15,3.15,4.15,5.15,6.15,7.15,8.15))+
  ggsci::scale_fill_jco()+
  ggsci::scale_color_jco()+
  theme_half_open()+
  labs(y =expression("P"*"s"*"e"*"u"*"d"*"o"*"g"*"e"*"n"*"e"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[10]*"("*"b"*"p"*")"*")") ,
       x =NULL,fill=NULL) +
  coord_cartesian(ylim = c(2, 5.4))+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 18)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(legend.position = "top", legend.direction = "horizontal")
ggsave("/home/zhzhang/PG/1Identification/sp_pglen.png",
       p1,width = 8, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")


```


```r
library(ggsci)
library(cowplot, lib.loc = "/usr/local/lib/R/site-library")
#各物种母基因数量PLOT
#导入假基因信息
Danio_rerio_hpg <- read.delim("~/PG/pg_message/Danio_rerio_hpg.txt")%>%
  mutate(sp="Zebrafish")%>%
  mutate(num=1)%>%
  distinct(parentgene_id,.keep_all = T)
Gallus_gallus_hpg <- read.delim("~/PG/pg_message/Gallus_gallus_hpg.txt")%>%
  mutate(sp="Chicken")%>%
  mutate(num=1)%>%
  distinct(parentgene_id,.keep_all = T)
Homo_sapiens_hpg <- read.delim("~/PG/pg_message/Homo_sapiens_hpg.txt")%>%
  mutate(sp="Human")%>%
  mutate(num=1)%>%
  distinct(parentgene_id,.keep_all = T)
Mus_musculus_hpg <- read.delim("~/PG/pg_message/Mus_musculus_hpg.txt")%>%
  mutate(sp="Mouse")%>%
  mutate(num=1)%>%
  distinct(parentgene_id,.keep_all = T)
Oryzias_latipes_hpg <- read.delim("~/PG/pg_message/Oryzias_latipes_hpg.txt")%>%
  mutate(sp="Medaka")%>%
  mutate(num=1)%>%
  distinct(parentgene_id,.keep_all = T)
#统计各物种母基因数量
ALL_PARG <- rbind(Danio_rerio_hpg,Gallus_gallus_hpg,Homo_sapiens_hpg,Mus_musculus_hpg,Oryzias_latipes_hpg)%>%
  group_by(sp)%>%
  summarise(num=sum(num))
ALL_PARG$sp <- factor(ALL_PARG$sp,levels = c("Zebrafish","Medaka","Chicken","Mouse","Human"))
#各物种母基因数量PLOT
p1 <- ggplot(data = ALL_PARG,aes(x=num/1000,y=sp))+
  geom_col(width = 0.7,alpha=0.7,aes(fill=sp,color=sp))+
  geom_text(size=7,aes(x=(num/1000)+0.5,label=num))+
  scale_fill_npg()+
  scale_color_npg()+
  theme_half_open()+
  labs(x = "No. of Parent Genes (×10E+3)", y = NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,4.5),breaks = c(0,1,2,3,4))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 15),axis.ticks.y=element_blank())
ggsave("/home/zhzhang/PG/1Identification/sp_parentgenenum.png", p1,width = 5, height = 8,dpi=1200, units = "in", device='png',bg = "transparent")










#各物种假基因数量PLOT
#导入假基因信息
Danio_rerio_hpg <- read.delim("~/PG/pg_message/Danio_rerio_hpg.txt")%>%
  mutate(sp="Zebrafish")%>%
  mutate(num=1)
Gallus_gallus_hpg <- read.delim("~/PG/pg_message/Gallus_gallus_hpg.txt")%>%
  mutate(sp="Chicken")%>%
  mutate(num=1)
Homo_sapiens_hpg <- read.delim("~/PG/pg_message/Homo_sapiens_hpg.txt")%>%
  mutate(sp="Human")%>%
  mutate(num=1)
Mus_musculus_hpg <- read.delim("~/PG/pg_message/Mus_musculus_hpg.txt")%>%
  mutate(sp="Mouse")%>%
  mutate(num=1)
Oryzias_latipes_hpg <- read.delim("~/PG/pg_message/Oryzias_latipes_hpg.txt")%>%
  mutate(sp="Medaka")%>%
  mutate(num=1)
#统计各物种两类假基因数量
ALL_PG <- rbind(Danio_rerio_hpg,Gallus_gallus_hpg,Homo_sapiens_hpg,Mus_musculus_hpg,Oryzias_latipes_hpg)%>%
  group_by(sp)%>%
  summarise(num=sum(num))
ALL_PG$sp <- factor(ALL_PG$sp,levels = c("Zebrafish","Medaka","Chicken","Mouse","Human"))
#各物种假基因数量PLOT
p1 <- ggplot(data = ALL_PG,aes(x=num/1000,y=sp))+
  geom_col(width = 0.7,alpha=0.7,aes(fill=sp,color=sp))+
  geom_text(size=7,aes(x=(num/1000)+1.4,label=num))+
  scale_fill_npg()+
  scale_color_npg()+
  theme_half_open()+
  labs(x = "No. of Pseudogenes (×10E+3)", y = NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,15),breaks = c(0,2,4,6,8,10,12))+
  theme(legend.position = "none")+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 15),axis.ticks.y=element_blank())
ggsave("/home/zhzhang/PG/1Identification/sp_pgnum.png", p1,width = 5, height = 8,dpi=1200, units = "in", device='png',bg = "transparent")





#各物种假基因数量PLOT(分为两类)
#导入假基因信息
Danio_rerio_hpg <- read.delim("~/PG/pg_message/Danio_rerio_hpg.txt")%>%
  mutate(sp="Zebrafish")%>%
  mutate(num=1)
Gallus_gallus_hpg <- read.delim("~/PG/pg_message/Gallus_gallus_hpg.txt")%>%
  mutate(sp="Chicken")%>%
  mutate(num=1)
Homo_sapiens_hpg <- read.delim("~/PG/pg_message/Homo_sapiens_hpg.txt")%>%
  mutate(sp="Human")%>%
  mutate(num=1)
Mus_musculus_hpg <- read.delim("~/PG/pg_message/Mus_musculus_hpg.txt")%>%
  mutate(sp="Mouse")%>%
  mutate(num=1)
Oryzias_latipes_hpg <- read.delim("~/PG/pg_message/Oryzias_latipes_hpg.txt")%>%
  mutate(sp="Medaka")%>%
  mutate(num=1)
#统计各物种两类假基因数量
ALL_PG <- rbind(Danio_rerio_hpg,Gallus_gallus_hpg,Homo_sapiens_hpg,Mus_musculus_hpg,Oryzias_latipes_hpg)%>%
  group_by(sp,type)%>%
  summarise(num=sum(num))
ALL_PG$sp <- factor(ALL_PG$sp,levels = c("Zebrafish","Medaka","Chicken","Mouse","Human"))
#数量标签
ALL_PG_label <- mutate(ALL_PG,weizhi=num)
ALL_PG_label[1,4] <- 1980/2+819
ALL_PG_label[3,4] <- 2575/2+8144
ALL_PG_label[5,4] <- 4039/2+5630
ALL_PG_label[7,4] <- 2053/2+9904
ALL_PG_label[9,4] <- 2091/2+916
ALL_PG_label[2,4] <- 819/2
ALL_PG_label[4,4] <- 8144/2
ALL_PG_label[6,4] <- 5630/2
ALL_PG_label[8,4] <- 9904/2
ALL_PG_label[10,4] <- 916/2
#各物种假基因数量PLOT(分为两类)
p1 <- ggplot(data = ALL_PG,aes(x=num/1000,y=sp))+
  geom_col(width = 0.7,alpha=0.7,aes(fill=type,color=type))+
  geom_text(data=ALL_PG_label,aes(x=weizhi/1000,label=num))+
  scale_fill_manual(values = c("#5EB47D","#F9E654"),limits=c("Duplicated", "Processed"))+
  scale_color_manual(values = c("#5EB47D","#F9E654"),limits=c("Duplicated", "Processed"))+
  theme_half_open()+
  guides(color="none")+
  labs(x = "No. of Pseudogenes (×10E+3)", y = NULL, fill = NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,12.5),breaks = c(0,2,4,6,8,10,12))+
  theme(legend.position = c(0.6, 0.1))+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 15),axis.ticks.y=element_blank())
ggsave("/home/zhzhang/PG/1Identification/sp_pgnum.png", p1,width = 5, height = 8,dpi=1200, units = "in", device='png',bg = "transparent")


#各物种假基因比对覆盖率图
ALL_PG_alc <- select(rbind(Danio_rerio_hpg,Gallus_gallus_hpg,Homo_sapiens_hpg,Mus_musculus_hpg,Oryzias_latipes_hpg),18,6)
ALL_PG_alc$sp <- factor(ALL_PG_alc$sp,levels = c("Zebrafish","Medaka","Chicken","Mouse","Human"))
p2 <- ggplot(data = ALL_PG_alc,aes(x=frac,y=sp))+
  geom_boxplot(notch=T,width = 0.7,alpha=0.7,outlier.alpha = 0,aes(fill=sp))+
  scale_fill_npg()+
  theme_half_open()+
  guides(fill="none")+
  labs(x = "Alignment Coverage", y = NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,1.05),labels = c("0","0.25","0.5","0.75","1")) + 
  theme(axis.title = element_text(size = 20),axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 15),axis.ticks.y=element_blank())
ggsave("/home/zhzhang/PG/1Identification/sp_pg_alignmentcover.png", p2,width = 5, height = 8,dpi=1200, units = "in", device='png',bg = "transparent")


#各物种基因组总大小图
sp_allgenome_size <- read.csv("~/PG/1Identification/sp_allgenome_size.txt", sep="")
sp_allgenome_size$sp <- factor(sp_allgenome_size$sp,levels = c("Zebrafish","Medaka","Chicken","Mouse","Human"))
p3 <- ggplot(data = sp_allgenome_size,aes(x=genome_size/1000000000,y=sp))+
  geom_col(width = 0.7,alpha=0.7,aes(fill=sp,color=sp))+
  scale_fill_npg()+
  scale_color_npg()+
  theme_half_open()+
  guides(fill="none",color="none")+
  labs(x = "Genome Size (Gb)", y = NULL) +
  scale_x_continuous(expand = c(0,0)) +
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 23),
        axis.text.x = element_text(size = 15))
ggsave("/home/zhzhang/PG/1Identification/sp_genome_size.png", p3,width = 6, height = 8,dpi=1200, units = "in", device='png',bg = "transparent")


#各物种coding gene数量
sp_pcg_num <- read.csv("~/PG/1Identification/sp_pcg_num.txt", sep="")
sp_pcg_num$sp <- factor(sp_pcg_num$sp,levels = c("Zebrafish","Medaka","Chicken","Mouse","Human"))
p4 <- ggplot(data = sp_pcg_num,aes(x=pcg_num/10000,y=sp))+
  geom_col(width = 0.7,alpha=0.7,aes(fill=sp,color=sp))+
  scale_fill_npg()+
  scale_color_npg()+
  theme_half_open()+
  guides(fill="none",color="none")+
  labs(x = "No. of Coding Genes (×10E+4)", y = NULL) +
  scale_x_continuous(expand = c(0,0),breaks = c(0,1,2,3),limits = c(0,3)) +
  theme(axis.title = element_text(size = 20),axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 15),axis.ticks.y=element_blank())
ggsave("/home/zhzhang/PG/1Identification/sp_pcg_num.png", p4,width = 5, height = 8,dpi=1200, units = "in", device='png',bg = "transparent")


#各物种假基因候选数量
sp_pg_candidates_num <- read.csv("~/PG/1Identification/sp_pg_candidates_num.txt", sep="")
sp_pg_candidates_num$sp <- factor(sp_pg_candidates_num$sp,levels = c("Zebrafish","Medaka","Chicken","Mouse","Human"))
p5 <- ggplot(data = sp_pg_candidates_num,aes(x=pg_candidates_num/10000,y=sp))+
  geom_col(width = 0.7,alpha=0.7,aes(fill=sp,color=sp))+
  scale_fill_npg()+
  scale_color_npg()+
  theme_half_open()+
  guides(fill="none",color="none")+
  labs(x = "No. of Pseudogene candidates (×10E+4)", y = NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,8)) +
  theme(axis.title = element_text(size = 20),axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 15),axis.ticks.y=element_blank())
ggsave("/home/zhzhang/PG/1Identification/sp_pgcandidates_num.png", p5,width = 5, height = 8,dpi=1200, units = "in", device='png',bg = "transparent")



#各物种WGD来源的假基因
sp_WGDpg_num <- read.csv("~/PG/1Identification/sp_WGDpg_num.txt", sep="")
sp_WGDpg_num$sp <- factor(sp_WGDpg_num$sp,levels = c("Zebrafish","Medaka","Chicken","Mouse","Human"))
p6 <- ggplot(data = sp_WGDpg_num,aes(x=WGDpg_num,y=sp))+
  geom_col(width = 0.7,alpha=0.7,aes(fill=sp,color=sp))+
  geom_text(size=7,aes(x=WGDpg_num+8,label=WGDpg_num))+
  scale_fill_npg()+
  scale_color_npg()+
  theme_half_open()+
  guides(fill="none",color="none")+
  labs(x = "No. of WGD-derived Pseudogenes", y = NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,140),breaks = c(0,25,50,75,100,125)) +
  theme(axis.title = element_text(size = 20),axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 15),axis.ticks.y=element_blank())
ggsave("/home/zhzhang/PG/1Identification/sp_WGDpg_num.png", p6,width = 5, height = 8,dpi=1200, units = "in", device='png',bg = "transparent")





#各物种lncrna假基因数量
sp_lncrnapg_num <- read.csv("~/PG/1Identification/sp_lncrnapg_num.txt", sep="")
sp_lncrnapg_num$sp <- factor(sp_lncrnapg_num$sp,levels = c("Zebrafish","Medaka","Chicken","Mouse","Human"))
p7 <- ggplot(data = sp_lncrnapg_num,aes(x=lncpg_num,y=sp))+
  geom_col(width = 0.7,alpha=0.7,aes(fill=sp,color=sp))+
  geom_text(size=7,aes(x=lncpg_num-90,label=lncpg_num))+
  scale_fill_npg()+
  scale_color_npg()+
  theme_half_open()+
  guides(fill="none",color="none")+
  labs(x = "Pseudogenes Encoding LncRNA", y = NULL) +
  theme(axis.title = element_text(size = 20),axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 15),axis.ticks.y=element_blank())+
  scale_x_continuous(expand = c(0,0))
ggsave("/home/zhzhang/PG/1Identification/sp_lncrnapg_num.png", 
       p7,width = 5, height = 8,dpi=1200, units = "in", device='png',bg = "transparent")


```


##### \[9.\]参与lncRNA起源的假基因的类型分布
```r
#探究各类假基因参与lncRNA起源的比例
#pgdlncRNA_pg输入pgdlnc及其对应的假基因，hpg输入高质量假基因信息，函数返回每种类型假基因中衍生出lncrna的假基因数
tjdlnctyped <- function(pgdlncRNA_pg,hpg) {
  #导入pgdlnc及其对应的假基因
  pgdlncRNA_pg <- read.delim(pgdlncRNA_pg)
  colnames(pgdlncRNA_pg)[2] <- "geneid"
  #导入假基因信息
  pg_parent <- read.delim(hpg)
  #添加是否衍生lncRNA的属性
  pg_parent <- left_join(pg_parent,mutate(distinct(pgdlncRNA_pg,geneid),lnc=1),by="geneid")
  pg_parent$lnc[is.na(pg_parent$lnc)==T] <- 0
  #统计不同类型的假基因衍生lncrna的情况
  type_lnc <- group_by(pg_parent,type)%>%
    summarise(allnum=n(),dlncnum=sum(lnc))%>%
    mutate(nodlncnum=allnum-dlncnum,dlncratio=dlncnum*100/allnum)
  return(type_lnc)
}
#人
hs_type_lnc <- tjdlnctyped("~/PG/lncRNA_class/Homo_sapiens.pgdlncRNA_pg.txt",
                           "~/PG/pg_message/Homo_sapiens_hpg.txt")%>%
  mutate(sp="Human")
#小鼠
mm_type_lnc <- tjdlnctyped("~/PG/lncRNA_class/Mus_musculus.pgdlncRNA_pg.txt",
                           "~/PG/pg_message/Mus_musculus_hpg.txt")%>%
  mutate(sp="Mouse")
#鸡
gg_type_lnc <- tjdlnctyped("~/PG/lncRNA_class/Gallus_gallus.pgdlncRNA_pg.txt",
                           "~/PG/pg_message/Gallus_gallus_hpg.txt")%>%
  mutate(sp="Chicken")
#斑马鱼
dr_type_lnc <- tjdlnctyped("~/PG/lncRNA_class/Danio_rerio.pgdlncRNA_pg.txt",
                           "~/PG/pg_message/Danio_rerio_hpg.txt")%>%
  mutate(sp="Zebrafish")
#猕猴
mac_type_lnc <- tjdlnctyped("~/PG/lncRNA_class/Macaca_mulatta.pgdlncRNA_pg.txt",
                           "~/PG/pg_message/Macaca_mulatta_hpg.txt")%>%
  mutate(sp="Macaque")
#大鼠
rat_type_lnc <- tjdlnctyped("~/PG/lncRNA_class/Rattus_norvegicus.pgdlncRNA_pg.txt",
                            "~/PG/pg_message/Rattus_norvegicus_hpg.txt")%>%
  mutate(sp="Rat")
#兔子
oc_type_lnc <- tjdlnctyped("~/PG/lncRNA_class/Oryctolagus_cuniculus.pgdlncRNA_pg.txt",
                            "~/PG/pg_message/Oryctolagus_cuniculus_hpg.txt")%>%
  mutate(sp="Rabbit")
#负鼠
op_type_lnc <- tjdlnctyped("~/PG/lncRNA_class/Monodelphis_domestica.pgdlncRNA_pg.txt",
                           "~/PG/pg_message/Monodelphis_domestica_hpg.txt")%>%
  mutate(sp="Opossum")
#合并
type_lnc <- rbind(hs_type_lnc,mac_type_lnc,mm_type_lnc,rat_type_lnc,oc_type_lnc,op_type_lnc,gg_type_lnc,dr_type_lnc)
#参与lncRNA起源的假基因类型分布
typedis <- select(type_lnc,1,3,6)%>%
  select(3,1,2)
colnames(typedis)[3] <- "num"
#导入假基因类型分布
SPpgnum <- read.delim("~/PG/1Identification/SP.pgnum.tj.txt")
#合并
he <- rbind(mutate(SPpgnum,class="Pseudogene"),
            mutate(typedis,class="Pseudogenelnc"))
#转换为百分比
hetj <- group_by(he,sp,class)%>%
  summarise(allnum=sum(num))
he <- left_join(mutate(he,forhe=paste(sp,class)),
                select(mutate(hetj,forhe=paste(sp,class)),forhe,allnum),
                by="forhe")%>%
  select(1,2,3,4,7)
colnames(he)[1] <- "sp"
he <- mutate(he,percent=num*100/allnum)%>%
  mutate(nonum=allnum-num)
#检验参与lncRNA衍生的假基因中DUP比例显著高于全部假基因中DUP比例
#p-value < 2.2e-16
fisher.test(filter(he,sp=="Human" & type=="Duplicated")[,c(3,7)])
#p-value < 2.2e-16
fisher.test(filter(he,sp=="Macaque" & type=="Duplicated")[,c(3,7)])
#p-value < 2.2e-16
fisher.test(filter(he,sp=="Mouse" & type=="Duplicated")[,c(3,7)])
#p-value < 2.2e-16
fisher.test(filter(he,sp=="Rat" & type=="Duplicated")[,c(3,7)])
#p-value < 2.2e-16
fisher.test(filter(he,sp=="Rabbit" & type=="Duplicated")[,c(3,7)])
#p-value < 2.2e-16
fisher.test(filter(he,sp=="Opossum" & type=="Duplicated")[,c(3,7)])
#p-value = 0.00000000006354
fisher.test(filter(he,sp=="Chicken" & type=="Duplicated")[,c(3,7)])
#p-value < 2.2e-16
fisher.test(filter(he,sp=="Zebrafish" & type=="Duplicated")[,c(3,7)])
#plot
he$class <- factor(he$class,levels = c("Pseudogene","Pseudogenelnc"))
he$type <- factor(he$type,levels = c("Duplicated","Processed","Fragment"))
he$sp <- factor(he$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
sp <- ggplot(data = he,aes(x=class,y=percent))+
  geom_signif(annotations=c("***"),
              y_position=c(100),tip_length = 0,size=0,textsize=5,
              xmin = c(1),
              xmax = c(2))+
  geom_col(aes(fill=type,color=type),width = 0.5)+
  ggalluvial::geom_alluvium(aes(alluvium=type,fill=type,color=type),
                            width = 0.5,curve_type="linear",alpha=0.3)+
  ggsci::scale_fill_jco()+
  ggsci::scale_color_jco()+
  theme_half_open()+
  labs(x = NULL, y ="Percentage (%)",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 0)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "none", legend.direction = "horizontal")+
  facet_wrap(~sp,nrow=1,strip.position="top",as.table=F)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))
ggsave("/home/zhzhang/PG/1Identification/SP.dlncpg_typedistribution.pdf", 
       sp,width = 8, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


##### \[10.\]不同期望概率区间内，对比DUP/PSSD假基因参与lncRNA衍生的比例
```r
#排除长度因素，证明dup比pssd更有利于lnc起源
#pgdlncRNA_pg输入pglnc及其对应的假基因，hpg输入假基因信息，
#shufflepgdlncRNA_pg输入shufflepglnc及其对应的shuffle假基因，shufflehpg输入shuffle假基因
#函数输出假基因以及其是否参与lncRNA起源，预期参与概率（根据shuffle计算）
pglendelnc <- function(pgdlncRNA_pg,hpg,shufflepgdlncRNA_pg,shufflehpg,pgexp){
  #导入pgdlnc及其对应的假基因
  pgdlncRNA_pg <- read.delim(pgdlncRNA_pg)
  colnames(pgdlncRNA_pg)[2] <- "geneid"
  #导入假基因
  pg <- read.delim(hpg)%>%
    select(geneid,type)
  #添加是否衍生lncRNA的属性
  pg <- left_join(pg,mutate(distinct(pgdlncRNA_pg,geneid),lnc=1),by="geneid")
  pg$lnc[is.na(pg$lnc)==T] <- 0
  #导入shufflepgdlnc及其对应的shuffle假基因
  shufflepgdlncRNA_pg <- read.delim(shufflepgdlncRNA_pg)
  colnames(shufflepgdlncRNA_pg)[2] <- "geneid"
  #导入shuffle假基因信息
  shufflepg <- read.delim(shufflehpg,header = F)%>%
    select(V7)
  colnames(shufflepg)[1] <- "geneid"
  #添加是否衍生lncRNA的属性
  shufflepg <- left_join(shufflepg,mutate(distinct(shufflepgdlncRNA_pg,geneid),lnc=1),by="geneid")
  shufflepg$lnc[is.na(shufflepg$lnc)==T] <- 0
  #计算每个假基因在1000次shuffle中，预期参与lncRNA衍生的概率
  shufflepg_exp <- separate(shufflepg,geneid,into = c("one","two","three"),sep = "_")%>%
    unite("geneid",one,two,three,sep = "_")%>%
    group_by(geneid)%>%
    summarise(exp=sum(lnc)/1000)
  #合并
  he <- left_join(pg,shufflepg_exp,by="geneid")
  #输出
  data.table::fwrite(he,
                     file =pgexp,
                     sep = '\t',row.names = F,quote = F,col.names = T)
}
#人
pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Homo_sapiens.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Homo_sapiens_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Homo_sapiens.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Homo_sapiens_hpg_shuffle.sort.bed",
                    pgexp="/home/zhzhang/PG/1Identification/Homo_sapiens.hpg_dlncexp.txt")
#小鼠
pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Mus_musculus.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Mus_musculus_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Mus_musculus.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Mus_musculus_hpg_shuffle.sort.bed",
                    pgexp="/home/zhzhang/PG/1Identification/Mus_musculus.hpg_dlncexp.txt")
#鸡
pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Gallus_gallus.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Gallus_gallus_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Gallus_gallus.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Gallus_gallus_hpg_shuffle.sort.bed",
                    pgexp="/home/zhzhang/PG/1Identification/Gallus_gallus.hpg_dlncexp.txt")
#斑马鱼
pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Danio_rerio.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Danio_rerio_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Danio_rerio.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Danio_rerio_hpg_shuffle.sort.bed",
                    pgexp="/home/zhzhang/PG/1Identification/Danio_rerio.hpg_dlncexp.txt")
#猕猴
pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Macaca_mulatta.pgdlncRNA_pg.txt",
           hpg="~/PG/pg_message/Macaca_mulatta_hpg.txt",
           shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Macaca_mulatta.pgdlncRNA_pg.txt",
           shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Macaca_mulatta_hpg_shuffle.sort.bed",
           pgexp="/home/zhzhang/PG/1Identification/Macaca_mulatta.hpg_dlncexp.txt")
#大鼠
pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Rattus_norvegicus.pgdlncRNA_pg.txt",
           hpg="~/PG/pg_message/Rattus_norvegicus_hpg.txt",
           shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Rattus_norvegicus.pgdlncRNA_pg.txt",
           shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Rattus_norvegicus_hpg_shuffle.sort.bed",
           pgexp="/home/zhzhang/PG/1Identification/Rattus_norvegicus.hpg_dlncexp.txt")
#兔子
pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Oryctolagus_cuniculus.pgdlncRNA_pg.txt",
           hpg="~/PG/pg_message/Oryctolagus_cuniculus_hpg.txt",
           shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Oryctolagus_cuniculus.pgdlncRNA_pg.txt",
           shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Oryctolagus_cuniculus_hpg_shuffle.sort.bed",
           pgexp="/home/zhzhang/PG/1Identification/Oryctolagus_cuniculus.hpg_dlncexp.txt")
#负鼠
pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Monodelphis_domestica.pgdlncRNA_pg.txt",
           hpg="~/PG/pg_message/Monodelphis_domestica_hpg.txt",
           shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Monodelphis_domestica.pgdlncRNA_pg.txt",
           shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Monodelphis_domestica_hpg_shuffle.sort.bed",
           pgexp="/home/zhzhang/PG/1Identification/Monodelphis_domestica.hpg_dlncexp.txt")


#人
hs_hpg_dlncexp <- read.delim("~/PG/1Identification/Homo_sapiens.hpg_dlncexp.txt")%>%
  mutate(sp="Human")
#小鼠
mm_hpg_dlncexp <- read.delim("~/PG/1Identification/Mus_musculus.hpg_dlncexp.txt")%>%
  mutate(sp="Mouse")
#鸡
gg_hpg_dlncexp <- read.delim("~/PG/1Identification/Gallus_gallus.hpg_dlncexp.txt")%>%
  mutate(sp="Chicken")
#斑马鱼
dr_hpg_dlncexp <- read.delim("~/PG/1Identification/Danio_rerio.hpg_dlncexp.txt")%>%
  mutate(sp="Zebrafish")
#猕猴
mac_hpg_dlncexp <- read.delim("~/PG/1Identification/Macaca_mulatta.hpg_dlncexp.txt")%>%
  mutate(sp="Macaque")
#大鼠
rat_hpg_dlncexp <- read.delim("~/PG/1Identification/Rattus_norvegicus.hpg_dlncexp.txt")%>%
  mutate(sp="Rat")
#兔子
oc_hpg_dlncexp <- read.delim("~/PG/1Identification/Oryctolagus_cuniculus.hpg_dlncexp.txt")%>%
  mutate(sp="Rabbit")
#负鼠
op_hpg_dlncexp <- read.delim("~/PG/1Identification/Monodelphis_domestica.hpg_dlncexp.txt")%>%
  mutate(sp="Opossum")
#根据预期概率划分类型
hpg_dlncexp <- rbind(hs_hpg_dlncexp,mac_hpg_dlncexp,mm_hpg_dlncexp,rat_hpg_dlncexp,oc_hpg_dlncexp,op_hpg_dlncexp,gg_hpg_dlncexp,dr_hpg_dlncexp)%>%
  mutate(exptype="un")
hpg_dlncexp <- hpg_dlncexp[is.na(hpg_dlncexp$exp)==F,]

hpg_dlncexp$exptype[hpg_dlncexp$exp < 0.005] <- "[Min, 5‰)"
hpg_dlncexp$exptype[hpg_dlncexp$exp < 0.01 & hpg_dlncexp$exp >= 0.005] <- "[5‰, 10‰)"
hpg_dlncexp$exptype[hpg_dlncexp$exp < 0.015 & hpg_dlncexp$exp >= 0.01] <- "[10‰, 15‰)"
hpg_dlncexp$exptype[hpg_dlncexp$exp < 0.02 & hpg_dlncexp$exp >= 0.015] <- "[15‰, 20‰)"
hpg_dlncexp$exptype[hpg_dlncexp$exp >= 0.02 ] <- "[20‰, Max]"

#统计不同预期概率的三类基因，参与lncRNA起源的比例
tj <- group_by(hpg_dlncexp,sp,exptype,type)%>%
  summarise(dlncnum=sum(lnc),allnum=n(),expmedian=median(exp))%>%
  mutate(nonum=allnum-dlncnum,delncratio=dlncnum*100/allnum)
tj$type <- factor(tj$type,levels = c("Duplicated","Processed","Fragment"))
tj$exptype <- factor(tj$exptype,levels = c("[Min, 5‰)","[5‰, 10‰)","[10‰, 15‰)","[15‰, 20‰)","[20‰, Max]"))
data.table::fwrite(tj,
                   file ="/home/zhzhang/PG/1Identification/SP.exptype.type.PG.dlncratio.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
tj <- filter(tj,type=="Duplicated" |type=="Processed" )
#plot
spplot <- "Macaque" #Human Macaque Mouse Rat Rabbit Opossum Zebrafish
tjforsp <- filter(tj,sp==spplot)
tjforsp <- arrange(tjforsp,exptype)
anno <- c()
for (i in c("[Min, 5‰)","[5‰, 10‰)","[10‰, 15‰)","[15‰, 20‰)","[20‰, Max]")) {
  pvalue <- fisher.test(filter(tjforsp,exptype==i)[,c(4,7)])[["p.value"]]
  if (pvalue > 0.05) {
    pvalue <- "N.S."
  }
  if (pvalue < 0.05) {
    pvalue <- "*"
  }
  if (pvalue < 0.01) {
    pvalue <- "**"
  }
  if (pvalue < 0.001) {
    pvalue <- "***"
  }
  anno <- c(anno,pvalue)
}
pp <- ggplot(data = filter(tj,sp==spplot),aes(x=exptype,y=delncratio))+
  geom_col(aes(fill=type),width = 0.6,position = "dodge")+
  geom_signif(annotations=anno,
              y_position=filter(tjforsp,type=="Duplicated")$delncratio+1.5,tip_length = 0,
              xmin = c(0.85,1.85,2.85,3.85,4.85),xmax = c(1.15,2.15,3.15,4.15,5.15))+
  ggsci::scale_fill_jco()+
  theme_half_open()+
  labs(x = "Expected probability", y ="Pseudogenes\ninvolved in lncRNA origins (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "top", legend.direction = "horizontal")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(paste("/home/zhzhang/PG/1Identification/",spplot,".exptype.2typepg_dlncratio.pdf",sep = ""), 
       pp,width = 6.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```



##### \[7.\]证明蛋白编码基因的遗物序列(假基因)有利于lncRNA起源
```r
#将假基因在基因组上重排到非假基因以及非蛋白编码基因外显子区域(排除蛋白编码基因遗物影响，同时遵循假基因的位置分布准则，不改变每个假基因长度)
#人
#合并蛋白编码基因外显子区域+假基因
awk '{print $1"\t"$2"\t"$3}' /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg_pcgexon_forshuffle.bed
awk '{print $2"\t"$3"\t"$4}' "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/Homo_sapiens.GRCh38.108.exLocs.txt" >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg_pcgexon_forshuffle.bed
#重排
for i in $(seq 1 1000)
do
bedtools shuffle -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -excl /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg_pcgexon_forshuffle.bed -noOverlapping|awk -v time="${i}" '{print $0"_"time}' >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg_shuffle.bed
done
#排序添加新ID
bedtools sort -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg_shuffle.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg_shuffle.sort.bed
#小鼠
#合并蛋白编码基因外显子区域+假基因
awk '{print $1"\t"$2"\t"$3}' /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg_pcgexon_forshuffle.bed
awk '{print $2"\t"$3"\t"$4}' "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.exLocs.txt" >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg_pcgexon_forshuffle.bed
#重排
for i in $(seq 1 1000)
do
bedtools shuffle -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg.bed -g /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chrsize.txt -excl /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg_pcgexon_forshuffle.bed -noOverlapping |awk -v time="${i}" '{print $0"_"time}' >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg_shuffle.bed
done
#排序添加新ID
bedtools sort -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg_shuffle.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg_shuffle.sort.bed
#鸡
#合并蛋白编码基因外显子区域+假基因
awk '{print $1"\t"$2"\t"$3}' /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg_pcgexon_forshuffle.bed
awk '{print $2"\t"$3"\t"$4}' "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.exLocs.txt" >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg_pcgexon_forshuffle.bed
#重排
for i in $(seq 1 1000)
do
bedtools shuffle -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg.bed -g "/home/zhzhang/PG/refgenome/Gallus_gallus.GRCg7b.dna.chrsize.txt" -excl /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg_pcgexon_forshuffle.bed -noOverlapping |awk -v time="${i}" '{print $0"_"time}' >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg_shuffle.bed
done
#排序添加新ID
bedtools sort -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg_shuffle.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg_shuffle.sort.bed
#斑马鱼
#合并蛋白编码基因外显子区域+假基因
awk '{print $1"\t"$2"\t"$3}' /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg_pcgexon_forshuffle.bed
awk '{print $2"\t"$3"\t"$4}' "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/Danio_rerio.GRCz11.108.exLocs.txt" >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg_pcgexon_forshuffle.bed
#重排
for i in $(seq 1 1000)
do
bedtools shuffle -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg.bed -g "/home/zhzhang/PG/refgenome/Danio_rerio.GRCz11.dna.chrsize.txt" -excl /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg_pcgexon_forshuffle.bed -noOverlapping |awk -v time="${i}" '{print $0"_"time}' >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg_shuffle.bed
done
#排序添加新ID
bedtools sort -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg_shuffle.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg_shuffle.sort.bed
#猕猴
#合并蛋白编码基因外显子区域+假基因
awk '{print $1"\t"$2"\t"$3}' /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg_pcgexon_forshuffle.bed
awk '{print $2"\t"$3"\t"$4}' "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/Macaca_mulatta.Mmul_10.108.exLocs.txt" >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg_pcgexon_forshuffle.bed
#重排
for i in $(seq 1 1000)
do
bedtools shuffle -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg.bed -g "/home/zhzhang/PG/refgenome/Macaca_mulatta.Mmul_10.dna.chrsize.txt" -excl /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg_pcgexon_forshuffle.bed -noOverlapping |awk -v time="${i}" '{print $0"_"time}' >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg_shuffle.bed
echo "${i}"
done
#排序添加新ID
bedtools sort -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg_shuffle.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg_shuffle.sort.bed
#大鼠
#合并蛋白编码基因外显子区域+假基因
awk '{print $1"\t"$2"\t"$3}' /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg_pcgexon_forshuffle.bed
awk '{print $2"\t"$3"\t"$4}' "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/Rattus_norvegicus.mRatBN7.2.108.exLocs.txt" >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg_pcgexon_forshuffle.bed
#重排
for i in $(seq 1 1000)
do
bedtools shuffle -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg.bed -g "/home/zhzhang/PG/refgenome/Rattus_norvegicus.mRatBN7.2.dna.chrsize.txt" -excl /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg_pcgexon_forshuffle.bed -noOverlapping |awk -v time="${i}" '{print $0"_"time}' >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg_shuffle.bed
echo "${i}"
done
#排序添加新ID
bedtools sort -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg_shuffle.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg_shuffle.sort.bed
#兔子
#合并蛋白编码基因外显子区域+假基因
awk '{print $1"\t"$2"\t"$3}' /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg_pcgexon_forshuffle.bed
awk '{print $2"\t"$3"\t"$4}' "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/Oryctolagus_cuniculus.OryCun2.0.108.exLocs.txt" >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg_pcgexon_forshuffle.bed
#重排
for i in $(seq 1 1000)
do
bedtools shuffle -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg.bed -g "/home/zhzhang/PG/refgenome/Oryctolagus_cuniculus.OryCun2.0.dna.chrsize.txt" -excl /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg_pcgexon_forshuffle.bed -noOverlapping |awk -v time="${i}" '{print $0"_"time}' >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg_shuffle.bed
echo "${i}"
done
#排序添加新ID
bedtools sort -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg_shuffle.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg_shuffle.sort.bed
#负鼠
#合并蛋白编码基因外显子区域+假基因
awk '{print $1"\t"$2"\t"$3}' /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg_pcgexon_forshuffle.bed
awk '{print $2"\t"$3"\t"$4}' "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/Monodelphis_domestica.ASM229v1.108.exLocs.txt" >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg_pcgexon_forshuffle.bed
#重排
for i in $(seq 1 1000)
do
bedtools shuffle -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg.bed -g "/home/zhzhang/PG/refgenome/Monodelphis_domestica.ASM229v1.dna.chrsize.txt" -excl /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg_pcgexon_forshuffle.bed -noOverlapping |awk -v time="${i}" '{print $0"_"time}' >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg_shuffle.bed
echo "${i}"
done
#排序添加新ID
bedtools sort -i /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg_shuffle.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg_shuffle.sort.bed





#使用重排的假基因重新对lncRNA分类，找到衍生了lncRNA的重排假基因
#获取全部lncRNA外显子与重拍假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/）
cat "/home/zhzhang/PG/sp.txt" |while read i
do
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/${i}/${i}.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg_shuffle.sort.bed -wo > /home/zhzhang/PG/lncRNA_class/${i}/${i}.alllncRNAexon_intersect_shufflepg.bed
done

bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg_shuffle.sort.bed -wo > /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.alllncRNAexon_intersect_shufflepg.bed
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg_shuffle.sort.bed -wo > /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.alllncRNAexon_intersect_shufflepg.bed
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg_shuffle.sort.bed -wo > /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.alllncRNAexon_intersect_shufflepg.bed
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg_shuffle.sort.bed -wo > /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.alllncRNAexon_intersect_shufflepg.bed
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg_shuffle.sort.bed -wo > /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.alllncRNAexon_intersect_shufflepg.bed
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg_shuffle.sort.bed -wo > /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.alllncRNAexon_intersect_shufflepg.bed
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg_shuffle.sort.bed -wo > /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.alllncRNAexon_intersect_shufflepg.bed





#lncRNA基因分类，确定参与lncrna产生的shuffle假基因
#a输入假基因与lncRNA外显子交集文件，b输入lncRNA geneid列表，c输入假基因信息
#LncRNA_class输出LncRNA_pg输出pglnc对应的假基因
lncrnaclass <- function(a,b,c,LncRNA_pg){
  #导入lncRNA外显子与假基因交集文件
  alllncRNAexon_intersect_pg <- read.delim(a, header=FALSE)
  colnames(alllncRNAexon_intersect_pg)[c(4,5,13,14)] <- c("lncgid","lnctid","pgid","overlap_len")
  #导入lncRNA geneid列表
  all_lncRNA_geneid <- read.table(b, quote="\"", comment.char="")
  colnames(all_lncRNA_geneid) <- c("lncgid")
  #导入假基因及其母基因信息
  pg_parent <- read.delim(c)%>%
    select(16,17)%>%
    separate(parentgene_id,c("trans_acting_target"))
  colnames(pg_parent)[1] <- "pgid"
  #计算每个lncRNA基因每个转录本的全部外显子与假基因在相同链上的总重叠长度,并筛选出总重叠长度>=200的转录本（即与假基因显著交集的转录本）
  lncRNAPG <- filter(alllncRNAexon_intersect_pg,V6==V12)%>%
    group_by(lncgid,lnctid,pgid)%>%
    summarise("total_overlap_len"=sum(overlap_len))%>%
    filter(total_overlap_len>=200)
  #获得假基因来源的lncRNA基因的id
  pgdlncRNAid <- data.frame(lncRNAPG$lncgid)%>%
    distinct(lncRNAPG.lncgid)
  colnames(pgdlncRNAid) <- "lncgid"
  pgdlncRNAid <- mutate(pgdlncRNAid,type="Pseudogene-derived lncRNA")
  #获得与假基因有交集但不显著的lncRNA基因的id（Interference）【全部交集的lncid去除显著交集的lncid即为不显著的lncRNA基因的id（Interference）】
  InterferencelncRNAid <- distinct(alllncRNAexon_intersect_pg,lncgid)%>%
    dplyr::setdiff(select(pgdlncRNAid,lncgid))%>%
    mutate(type="Interference lncRNA")
  #lncRNA基因分类文件
  lncRNA_class <- left_join(all_lncRNA_geneid,pgdlncRNAid,by="lncgid")%>%
    left_join(InterferencelncRNAid,by="lncgid")%>%
    mutate(type="un")
  lncRNA_class$type[lncRNA_class$type.x=="Pseudogene-derived lncRNA"] <- "Pseudogene-derived lncRNA"
  lncRNA_class$type[lncRNA_class$type.y=="Interference lncRNA"] <- "Interference lncRNA"
  lncRNA_class$type[lncRNA_class$type=="un"] <- "Non-pseudogene-derived lncRNA"
  lncRNA_class <- select(lncRNA_class,lncgid,type)
  colnames(lncRNA_class)[1] <- "geneid"
  #确定假基因来源的lncRNA基因 对应的假基因
  #并根据假基因-母基因，确定假基因来源的lncRNA基因潜在的反式调控靶标基因
  pgdlncRNA_PG <- select(data.frame(lncRNAPG),-2) %>%
    arrange(lncgid,pgid,desc(total_overlap_len))%>%
    distinct(lncgid,pgid,.keep_all = T)%>%
    select(-3)%>%
    left_join(pg_parent,by="pgid")
  pgdlncRNA_Pse <- select(pgdlncRNA_PG,-3)%>%
    distinct(lncgid,pgid)
  #输出参与的假基因
  data.table::fwrite(pgdlncRNA_Pse,file =LncRNA_pg,sep = '\t',row.names = F,quote = F,col.names = T)
}

#小鼠
lncrnaclass("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Mus_musculus.alllncRNAexon_intersect_shufflepg.bed",
            "~/PG/lncRNA_class/Mus_musculus.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Mus_musculus_hpg.txt",
            "/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Mus_musculus.pgdlncRNA_pg.txt")
#人
lncrnaclass("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Homo_sapiens.alllncRNAexon_intersect_shufflepg.bed",
            "~/PG/lncRNA_class/Homo_sapiens.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Homo_sapiens_hpg.txt",
            "/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Homo_sapiens.pgdlncRNA_pg.txt")
#鸡
lncrnaclass("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Gallus_gallus.alllncRNAexon_intersect_shufflepg.bed",
            "~/PG/lncRNA_class/Gallus_gallus.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Gallus_gallus_hpg.txt",
            "/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Gallus_gallus.pgdlncRNA_pg.txt")
#斑马鱼
lncrnaclass("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Danio_rerio.alllncRNAexon_intersect_shufflepg.bed",
            "~/PG/lncRNA_class/Danio_rerio.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Danio_rerio_hpg.txt",
            "/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Danio_rerio.pgdlncRNA_pg.txt")
#猕猴
lncrnaclass("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Macaca_mulatta.alllncRNAexon_intersect_shufflepg.bed",
            "~/PG/lncRNA_class/Macaca_mulatta.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Macaca_mulatta_hpg.txt",
            "/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Macaca_mulatta.pgdlncRNA_pg.txt")
#大鼠
lncrnaclass("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Rattus_norvegicus.alllncRNAexon_intersect_shufflepg.bed",
            "~/PG/lncRNA_class/Rattus_norvegicus.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Rattus_norvegicus_hpg.txt",
            "/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Rattus_norvegicus.pgdlncRNA_pg.txt")
#兔子
lncrnaclass("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Oryctolagus_cuniculus.alllncRNAexon_intersect_shufflepg.bed",
            "~/PG/lncRNA_class/Oryctolagus_cuniculus.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Oryctolagus_cuniculus_hpg.txt",
            "/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Oryctolagus_cuniculus.pgdlncRNA_pg.txt")
#负鼠
lncrnaclass("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Monodelphis_domestica.alllncRNAexon_intersect_shufflepg.bed",
            "~/PG/lncRNA_class/Monodelphis_domestica.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Monodelphis_domestica_hpg.txt",
            "/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Monodelphis_domestica.pgdlncRNA_pg.txt")




###shuffle和正常假基因 产生lncRNA的假基因比例对比
#证明蛋白编码基因的遗物序列，有利于lncRNA起源
#pgdlncRNA_pg输入pglnc及其对应的假基因，hpg输入假基因信息，
#shufflepgdlncRNA_pg输入shufflepglnc及其对应的shuffle假基因，shufflehpg输入shuffle假基因
#函数返回假基因和每次shuffle假基因，参与衍生lncrna的数量比例
pglendelnc <- function(pgdlncRNA_pg,hpg,shufflepgdlncRNA_pg,shufflehpg){
  #导入pgdlnc及其对应的假基因
  pgdlncRNA_pg <- read.delim(pgdlncRNA_pg)
  colnames(pgdlncRNA_pg)[2] <- "geneid"
  #导入假基因
  pg <- read.delim(hpg)%>%
    mutate(type="Pseudogene")%>%
    select(geneid,type)
  #添加是否衍生lncRNA的属性
  pg <- left_join(pg,mutate(distinct(pgdlncRNA_pg,geneid),lnc=1),by="geneid")
  pg$lnc[is.na(pg$lnc)==T] <- 0
  pg <- mutate(pg,times=1)%>%
    select(times,type,lnc)
  #导入shufflepgdlnc及其对应的shuffle假基因
  shufflepgdlncRNA_pg <- read.delim(shufflepgdlncRNA_pg)
  colnames(shufflepgdlncRNA_pg)[2] <- "geneid"
  #导入shuffle假基因信息
  shufflepg <- read.delim(shufflehpg,header = F)%>%
    select(V7)%>%
    mutate(type="Pseudogene (shuffled)")
  colnames(shufflepg)[1] <- "geneid"
  #添加是否衍生lncRNA的属性
  shufflepg <- left_join(shufflepg,mutate(distinct(shufflepgdlncRNA_pg,geneid),lnc=1),by="geneid")
  shufflepg$lnc[is.na(shufflepg$lnc)==T] <- 0
  #统计每次shuffle结果
  shufflepg_time <- separate(shufflepg,geneid,c("one","two","three","times"),sep="_",remove = T)%>%
    select(times,type,lnc)
  #合并
  he <- rbind(pg,shufflepg_time)
  #统计参与衍生lncRNA的假基因比例
  tj <- group_by(he,type,times)%>%
    summarise(delncnum=sum(lnc),allnum=n())%>%
    mutate(nonum=allnum-delncnum,delncratio=delncnum*100/allnum)
}
#人
hs_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Homo_sapiens.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Homo_sapiens_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Homo_sapiens.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Homo_sapiens_hpg_shuffle.sort.bed")%>%
  mutate(sp="Human")
#小鼠
mm_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Mus_musculus.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Mus_musculus_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Mus_musculus.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Mus_musculus_hpg_shuffle.sort.bed")%>%
  mutate(sp="Mouse")
#鸡
gg_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Gallus_gallus.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Gallus_gallus_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Gallus_gallus.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Gallus_gallus_hpg_shuffle.sort.bed")%>%
  mutate(sp="Chicken")
#斑马鱼
dr_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Danio_rerio.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Danio_rerio_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Danio_rerio.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Danio_rerio_hpg_shuffle.sort.bed")%>%
  mutate(sp="Zebrafish")
#猕猴
mac_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Macaca_mulatta.pgdlncRNA_pg.txt",
                     hpg="~/PG/pg_message/Macaca_mulatta_hpg.txt",
                     shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Macaca_mulatta.pgdlncRNA_pg.txt",
                     shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Macaca_mulatta_hpg_shuffle.sort.bed")%>%
  mutate(sp="Macaque")
#大鼠
rat_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Rattus_norvegicus.pgdlncRNA_pg.txt",
                     hpg="~/PG/pg_message/Rattus_norvegicus_hpg.txt",
                     shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Rattus_norvegicus.pgdlncRNA_pg.txt",
                     shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Rattus_norvegicus_hpg_shuffle.sort.bed")%>%
  mutate(sp="Rat")
#兔子
oc_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Oryctolagus_cuniculus.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Oryctolagus_cuniculus_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Oryctolagus_cuniculus.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Oryctolagus_cuniculus_hpg_shuffle.sort.bed")%>%
  mutate(sp="Rabbit")
#负鼠
op_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Monodelphis_domestica.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Monodelphis_domestica_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Monodelphis_domestica.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Monodelphis_domestica_hpg_shuffle.sort.bed")%>%
  mutate(sp="Opossum")

#合并
pgratio <- rbind(hs_pg,mm_pg,gg_pg,dr_pg,mac_pg,rat_pg,oc_pg,op_pg)
#统计shuffle 1000次衍生lnc比例的均值和置信度区间
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(1024)
tj1 <- group_by(filter(pgratio,type=="Pseudogene (shuffled)"),sp,type)%>%
  summarise(mean=mean(delncratio),median=median(delncratio),
            confmin=boot::boot.ci(boot::boot(delncratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(delncratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][5])
tj2 <- group_by(filter(pgratio,type=="Pseudogene"),sp,type)%>%
  summarise(mean=mean(delncratio),median=median(delncratio))%>%
  mutate(confmin=mean,confmax=mean)
tj <- rbind(tj1,tj2)
data.table::fwrite(tj,
                   file ="/home/zhzhang/PG/1Identification/SP.pg_shufflepg_dlncratio.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#1000次shuffle结果合计与正常进行fishertest检验
fortest <- group_by(pgratio,sp,type)%>%
  summarise(delncnum=sum(delncnum),nonum=sum(nonum))
#检验人类假基因比shuffle衍生lnc的比例更高（p-value < 2.2e-16）
fisher.test(fortest[c(3,4),c(3,4)])
#检验小鼠假基因比shuffle衍生lnc的比例更高（p-value < 2.2e-16）
fisher.test(fortest[c(7,8),c(3,4)])
#检验鸡假基因比shuffle衍生lnc的比例更高（p-value = 0.0000000002598）
fisher.test(fortest[c(1,2),c(3,4)])
#检验斑马鱼假基因比shuffle衍生lnc的比例更高（p-value < 2.2e-16）
fisher.test(fortest[c(15,16),c(3,4)])
#检验猕猴假基因比shuffle衍生lnc的比例更高（p-value = 0.0000001584）
fisher.test(fortest[c(5,6),c(3,4)])
#检验大鼠假基因比shuffle衍生lnc的比例更高（p-value < 2.2e-16）
fisher.test(fortest[c(13,14),c(3,4)])
#检验兔子假基因比shuffle衍生lnc的比例更高（p-value < 2.2e-16）
fisher.test(fortest[c(11,12),c(3,4)])
#检验负鼠假基因比shuffle衍生lnc的比例更高（p-value < 2.2e-16）
fisher.test(fortest[c(9,10),c(3,4)])
#plot
tj$type <- factor(tj$type,levels = c("Pseudogene","Pseudogene (shuffled)"))
tj$sp <- factor(tj$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
pp <- ggplot(data = tj,aes(x=sp,y=mean))+
  geom_col(aes(fill=type),width = 0.6,position = "dodge")+
  geom_errorbar(data = filter(tj,type=="Pseudogene (shuffled)"),
                aes(ymin=confmin,ymax=confmax),width=0.15,position = position_nudge(x=0.15))+
  geom_signif(annotations=c(rep("***",times=6),"***","***"),
              y_position=c(8.3,1.5,4,2.5,2,2.6,6.3,5.3),tip_length = 0,
              xmin = c(0.85,1.85,2.85,3.85,4.85,5.85,6.85,7.85),
              xmax = c(1.15,2.15,3.15,4.15,5.15,6.15,7.15,8.15))+
  scale_fill_manual(values = c("#B75458","#FAA467"))+
  scale_y_continuous(breaks = c(0,2,4,6,8))+
  theme_half_open()+
  labs(x = NULL, y ="Pseudogenes\ninvolved in lncRNA origins (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 18)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(legend.position = "top", legend.direction = "horizontal")
ggsave("/home/zhzhang/PG/1Identification/SP.pg_shufflepg_dlncratio.pdf", 
       pp,width = 6.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 阴性对照：证明假基因不利于/无关于基因间区（无功能区域）起源
```r
#基因间区（人，小鼠）
cp /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed /home/zhzhang/PG/1Identification/intergenic/
cp /home/zhzhang/PG/Evolution/conserve/Mus_musculus.20000random_3kb_intergenic.bed /home/zhzhang/PG/1Identification/intergenic/

##获取基因间区bed文件
#获取全部基因bed文件
#其他
for i in $(seq 1 6)
do
sp=`sed -n ${i}p "/home/zhzhang/PG/sp2.txt"`
gtf=`find /home/zhzhang/PG/RNAseqdata/newGTF/${sp}*.rmpg.novellncRNA.gtf`
tail -n +6 ${gtf} |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/1Identification/intergenic/${sp}.pre_allgene.bed
done
#其他
sp <- c("Monodelphis_domestica",
"Rattus_norvegicus",
"Oryctolagus_cuniculus",
"Macaca_mulatta",
"Gallus_gallus",
"Danio_rerio")
for (i in 1:6) {
  forsp=sp[i]
a <- paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".pre_allgene.bed",sep="")
b <- paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".allgene.bed",sep="")
#导入gene外显子bed
pre_allgene <- read.delim(a, header=FALSE)
#生成gene bed
allgene <- group_by(pre_allgene,V1,V4,V6)%>%
  summarise(start=min(V2),end=max(V3))%>%
  mutate(an="gene")%>%
  select(1,4,5,2,6,3)
#导出genebed
data.table::fwrite(allgene,
                   file =b,
                   sep = '\t',row.names = F,quote = F,col.names = F)
}

#获取全部基因启动子bed文件
for i in $(seq 1 6)
do
sp=`sed -n ${i}p "/home/zhzhang/PG/sp2.txt"`
chrsize=`find /home/zhzhang/PG/refgenome/${sp}.*.dna.chrsize.txt`
cat "/home/zhzhang/PG/1Identification/intergenic/${sp}.allgene.bed"|awk '$6=="+" {print $1"\t"$2-2000"\t"$2"\t"$4"\tpromoter\t"$6} $6=="-" {print $1"\t"$3"\t"$3+2000"\t"$4"\tpromoter\t"$6}'|awk '$2<0 {print $1"\t0\t"$3"\t"$4"\t"$5"\t"$6} $2>=0 {print $0}'|bedtools sort > /home/zhzhang/PG/1Identification/intergenic/${sp}.allgenepromoter.bed
#合并基因区和启动子区域bed文件
cp /home/zhzhang/PG/1Identification/intergenic/${sp}.allgene.bed /home/zhzhang/PG/1Identification/intergenic/${sp}.allgeneApromoter.bed
cat /home/zhzhang/PG/1Identification/intergenic/${sp}.allgenepromoter.bed >> /home/zhzhang/PG/1Identification/intergenic/${sp}.allgeneApromoter.bed
bedtools sort -i /home/zhzhang/PG/1Identification/intergenic/${sp}.allgeneApromoter.bed|bedtools merge > /home/zhzhang/PG/1Identification/intergenic/${sp}.allgeneApromoter.sort.bed
rm /home/zhzhang/PG/1Identification/intergenic/${sp}.allgeneApromoter.bed
#补集获得基因间区bed文件
bedtools complement -i /home/zhzhang/PG/1Identification/intergenic/${sp}.allgeneApromoter.sort.bed -g ${chrsize} > /home/zhzhang/PG/1Identification/intergenic/${sp}.intergenic.bed
#生成随机的20000个3kb区域bed文件
bedtools random -n 20000 -l 3000 -seed 1024 -g ${chrsize} > /home/zhzhang/PG/1Identification/intergenic/${sp}.20000random_3kb_region.bed
#生成随机的20000个3kb基因间区bed文件
bedtools shuffle -i /home/zhzhang/PG/1Identification/intergenic/${sp}.20000random_3kb_region.bed -g ${chrsize} -incl /home/zhzhang/PG/1Identification/intergenic/${sp}.intergenic.bed -noOverlapping|bedtools sort|awk '{print $1"\t"$2"\t"$3"\tintergenic_"$4"\tintergenic\t"$6}' > /home/zhzhang/PG/1Identification/intergenic/${sp}.20000random_3kb_intergenic.bed
done

#使用假基因/重排的假基因对基因间区分类，找到衍生了基因间区基因间区的假基因/重排假基因
#获取全部基因间区与假基因重排假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/1Identification/intergenic/）
cat "/home/zhzhang/PG/sp.txt" |while read i
do
bedtools intersect -a /home/zhzhang/PG/1Identification/intergenic/${i}.20000random_3kb_intergenic.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg.bed -wo > /home/zhzhang/PG/1Identification/intergenic/${i}.intergenic_intersect_pg.bed
bedtools intersect -a /home/zhzhang/PG/1Identification/intergenic/${i}.20000random_3kb_intergenic.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg_shuffle.sort.bed -wo > /home/zhzhang/PG/1Identification/intergenic/${i}.intergenic_intersect_shufflepg.bed
done

#基因间区id列表
cat "/home/zhzhang/PG/sp.txt" |while read i
do
awk '{print $4}' "/home/zhzhang/PG/1Identification/intergenic/${i}.20000random_3kb_intergenic.bed" > /home/zhzhang/PG/1Identification/intergenic/${i}.intergenic.id.txt
done






#基因间区分类，确定参与基因间区产生的假基因/shuffle假基因
#a输入假基因与基因间区交集文件，b输入基因间区id列表，c输入假基因信息
#LncRNA_pg输出参与基因间区衍生的假基因
lncrnaclass <- function(a,b,c,LncRNA_pg){
  #导入基因间区与假基因交集文件
  alllncRNAexon_intersect_pg <- read.delim(a, header=FALSE)
  colnames(alllncRNAexon_intersect_pg)[c(4,13,14)] <- c("lncgid","pgid","overlap_len")
  #导入基因间区id列表
  all_lncRNA_geneid <- read.table(b, quote="\"", comment.char="")
  colnames(all_lncRNA_geneid) <- c("lncgid")
  #导入假基因及其母基因信息
  pg_parent <- read.delim(c)%>%
    select(16,17)%>%
    separate(parentgene_id,c("trans_acting_target"))
  colnames(pg_parent)[1] <- "pgid"
  #计算每个lncRNA基因每个转录本的全部外显子与假基因在相同链上的总重叠长度,并筛选出总重叠长度>=200的转录本（即与假基因显著交集的转录本）
  lncRNAPG <- filter(alllncRNAexon_intersect_pg,V6==V12)%>%
    group_by(lncgid,pgid)%>%
    summarise("total_overlap_len"=sum(overlap_len))%>%
    filter(total_overlap_len>=200)
  #获得假基因来源的lncRNA基因的id
  pgdlncRNAid <- data.frame(lncRNAPG$lncgid)%>%
    distinct(lncRNAPG.lncgid)
  colnames(pgdlncRNAid) <- "lncgid"
  pgdlncRNAid <- mutate(pgdlncRNAid,type="Pseudogene-derived lncRNA")
  #获得与假基因有交集但不显著的lncRNA基因的id（Interference）【全部交集的lncid去除显著交集的lncid即为不显著的lncRNA基因的id（Interference）】
  InterferencelncRNAid <- distinct(alllncRNAexon_intersect_pg,lncgid)%>%
    dplyr::setdiff(select(pgdlncRNAid,lncgid))%>%
    mutate(type="Interference lncRNA")
  #lncRNA基因分类文件
  lncRNA_class <- left_join(all_lncRNA_geneid,pgdlncRNAid,by="lncgid")%>%
    left_join(InterferencelncRNAid,by="lncgid")%>%
    mutate(type="un")
  lncRNA_class$type[lncRNA_class$type.x=="Pseudogene-derived lncRNA"] <- "Pseudogene-derived lncRNA"
  lncRNA_class$type[lncRNA_class$type.y=="Interference lncRNA"] <- "Interference lncRNA"
  lncRNA_class$type[lncRNA_class$type=="un"] <- "Non-pseudogene-derived lncRNA"
  lncRNA_class <- select(lncRNA_class,lncgid,type)
  colnames(lncRNA_class)[1] <- "geneid"
  #确定假基因来源的lncRNA基因 对应的假基因
  #并根据假基因-母基因，确定假基因来源的lncRNA基因潜在的反式调控靶标基因
  pgdlncRNA_PG <- data.frame(lncRNAPG)%>%
    arrange(lncgid,pgid,desc(total_overlap_len))%>%
    distinct(lncgid,pgid,.keep_all = T)%>%
    select(-3)%>%
    left_join(pg_parent,by="pgid")
  pgdlncRNA_Pse <- select(pgdlncRNA_PG,-3)%>%
    distinct(lncgid,pgid)
  #输出参与的假基因
  data.table::fwrite(pgdlncRNA_Pse,file =LncRNA_pg,sep = '\t',row.names = F,quote = F,col.names = T)
}
sp <- c("Monodelphis_domestica",
        "Rattus_norvegicus",
        "Oryctolagus_cuniculus",
        "Macaca_mulatta",
        "Gallus_gallus",
        "Danio_rerio",
        "Homo_sapiens",
        "Mus_musculus")
#假基因循环
for (i in 1:8) {
  forsp <- sp[i]
  lncrnaclass(paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".intergenic_intersect_pg.bed",sep=""),
              paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".intergenic.id.txt",sep=""),
              paste("/home/zhzhang/PG/pg_message/",forsp,"_hpg.txt",sep=""),
              paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".pgdintergenic_pg.txt",sep=""))
}
#shuffle假基因循环
for (i in 1:8) {
  forsp <- sp[i]
  lncrnaclass(paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".intergenic_intersect_shufflepg.bed",sep=""),
              paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".intergenic.id.txt",sep=""),
              paste("/home/zhzhang/PG/pg_message/",forsp,"_hpg.txt",sep=""),
              paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".pgdintergenic_shufflepg.txt",sep=""))
}






###shuffle和正常假基因 产生基因间区的假基因比例对比
#证明蛋白编码基因的遗物序列，无关于于基因间区起源
#pgdlncRNA_pg输入pgd基因间区及其对应的假基因，hpg输入假基因信息，
#shufflepgdlncRNA_pg输入shufflepg基因间区及其对应的shuffle假基因，shufflehpg输入shuffle假基因
#函数返回假基因和每次shuffle假基因，参与衍生基因间区的数量比例
pglendelnc <- function(pgdlncRNA_pg,hpg,shufflepgdlncRNA_pg,shufflehpg){
  #导入pgdlnc及其对应的假基因
  pgdlncRNA_pg <- read.delim(pgdlncRNA_pg)
  colnames(pgdlncRNA_pg)[2] <- "geneid"
  #导入假基因
  pg <- read.delim(hpg)%>%
    mutate(type="Pseudogene")%>%
    select(geneid,type)
  #添加是否衍生lncRNA的属性
  pg <- left_join(pg,mutate(distinct(pgdlncRNA_pg,geneid),lnc=1),by="geneid")
  pg$lnc[is.na(pg$lnc)==T] <- 0
  pg <- mutate(pg,times=1)%>%
    select(times,type,lnc)
  #导入shufflepgdlnc及其对应的shuffle假基因
  shufflepgdlncRNA_pg <- read.delim(shufflepgdlncRNA_pg)
  colnames(shufflepgdlncRNA_pg)[2] <- "geneid"
  #导入shuffle假基因信息
  shufflepg <- read.delim(shufflehpg,header = F)%>%
    select(V7)%>%
    mutate(type="Pseudogene (shuffled)")
  colnames(shufflepg)[1] <- "geneid"
  #添加是否衍生lncRNA的属性
  shufflepg <- left_join(shufflepg,mutate(distinct(shufflepgdlncRNA_pg,geneid),lnc=1),by="geneid")
  shufflepg$lnc[is.na(shufflepg$lnc)==T] <- 0
  #统计每次shuffle结果
  shufflepg_time <- separate(shufflepg,geneid,c("one","two","three","times"),sep="_",remove = T)%>%
    select(times,type,lnc)
  #合并
  he <- rbind(pg,shufflepg_time)
  #统计参与衍生lncRNA的假基因比例
  tj <- group_by(he,type,times)%>%
    summarise(delncnum=sum(lnc),allnum=n())%>%
    mutate(nonum=allnum-delncnum,delncratio=delncnum*100/allnum)
}
sp <- c("Monodelphis_domestica",
        "Rattus_norvegicus",
        "Oryctolagus_cuniculus",
        "Macaca_mulatta",
        "Gallus_gallus",
        "Danio_rerio",
        "Homo_sapiens",
        "Mus_musculus")
#循环
pgratio <- data.frame()
for (i in 1:8) {
  forsp <- sp[i]
  hs_pg <- pglendelnc(pgdlncRNA_pg=paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".pgdintergenic_pg.txt",sep=""),
                    hpg=paste("/home/zhzhang/PG/pg_message/",forsp,"_hpg.txt",sep=""),
                    shufflepgdlncRNA_pg=paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".pgdintergenic_shufflepg.txt",sep=""),
                    shufflehpg=paste("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/",forsp,"_hpg_shuffle.sort.bed",sep=""))%>%
  mutate(sp=forsp)
  #合并
  pgratio <- rbind(pgratio,hs_pg)
}
#统计shuffle 1000次衍生基因间区比例的均值和置信度区间
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(1024)
tj1 <- group_by(filter(pgratio,type=="Pseudogene (shuffled)"),sp,type)%>%
  summarise(mean=mean(delncratio),median=median(delncratio),
            confmin=boot::boot.ci(boot::boot(delncratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(delncratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][5])
tj2 <- group_by(filter(pgratio,type=="Pseudogene"),sp,type)%>%
  summarise(mean=mean(delncratio),median=median(delncratio))%>%
  mutate(confmin=mean,confmax=mean)
tj <- rbind(tj1,tj2)
data.table::fwrite(tj,
                   file ="/home/zhzhang/PG/1Identification/intergenic/SP.pg_shufflepg_dintergenicratio.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)

#1000次shuffle结果合计与正常进行fishertest检验
fortest <- group_by(pgratio,sp,type)%>%
  summarise(delncnum=sum(delncnum),nonum=sum(nonum))
#检验人类假基因比shuffle衍生lnc的比例更高（p-value = 0.764）
fisher.test(fortest[c(5,6),c(3,4)])
#检验小鼠假基因比shuffle衍生lnc的比例更高（p-value = 0.1945）
fisher.test(fortest[c(11,12),c(3,4)])
#检验鸡假基因比shuffle衍生lnc的比例更高（p-value = 0.000000000002816）
fisher.test(fortest[c(3,4),c(3,4)])
#检验斑马鱼假基因比shuffle衍生lnc的比例更高（p-value = 0.0261）
fisher.test(fortest[c(1,2),c(3,4)])
#检验猕猴假基因比shuffle衍生lnc的比例更高（p-value = 0.8197）
fisher.test(fortest[c(7,8),c(3,4)])
#检验大鼠假基因比shuffle衍生lnc的比例更高（p-value = 0.6652）
fisher.test(fortest[c(15,16),c(3,4)])
#检验兔子假基因比shuffle衍生lnc的比例更高（p-value = 0.1324）
fisher.test(fortest[c(13,14),c(3,4)])
#检验负鼠假基因比shuffle衍生lnc的比例更高（p-value = 0.07009）
fisher.test(fortest[c(9,10),c(3,4)])
#plot
tj$sp[tj$sp=="Homo_sapiens"] <- "Human"
tj$sp[tj$sp=="Danio_rerio"] <- "Zebrafish"
tj$sp[tj$sp=="Gallus_gallus"] <- "Chicken"
tj$sp[tj$sp=="Macaca_mulatta"] <- "Macaque"
tj$sp[tj$sp=="Monodelphis_domestica"] <- "Opossum"
tj$sp[tj$sp=="Mus_musculus"] <- "Mouse"
tj$sp[tj$sp=="Oryctolagus_cuniculus"] <- "Rabbit"
tj$sp[tj$sp=="Rattus_norvegicus"] <- "Rat"
tj$type <- factor(tj$type,levels = c("Pseudogene","Pseudogene (shuffled)"))
tj$sp <- factor(tj$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
pp <- ggplot(data = tj,aes(x=sp,y=mean))+
  geom_col(aes(fill=type),width = 0.6,position = "dodge")+
  geom_errorbar(data = filter(tj,type=="Pseudogene (shuffled)"),
                aes(ymin=confmin,ymax=confmax),width=0.15,position = position_nudge(x=0.15))+
  geom_signif(annotations=c(rep("N.S.",times=6),"***","*"),
              y_position=c(2.5,2.5,2.5,2.5,2.5,2.5,8,5.2),tip_length = 0,
              xmin = c(0.85,1.85,2.85,3.85,4.85,5.85,6.85,7.85),
              xmax = c(1.15,2.15,3.15,4.15,5.15,6.15,7.15,8.15))+
  scale_fill_manual(values = c("#4F5F87","#ADC0F0"))+
  scale_y_continuous(breaks = c(0,2,4,6,8))+
  theme_half_open()+
  labs(x = NULL, y ="Pseudogenes involved in\nintergenic region origins (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 18)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(legend.position = "top", legend.direction = "horizontal")
ggsave("/home/zhzhang/PG/1Identification/intergenic/SP.pg_shufflepg_dintergenicratio.pdf", 
       pp,width = 6.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```


##### \[8.\]证明每种类型的假基因都有利于lncRNA起源
```r
###shuffle和正常假基因 产生lncRNA的假基因比例对比
#证明蛋白编码基因的遗物序列，有利于lncRNA起源
#pgdlncRNA_pg输入pglnc及其对应的假基因，hpg输入假基因信息，
#shufflepgdlncRNA_pg输入shufflepglnc及其对应的shuffle假基因，shufflehpg输入shuffle假基因
#函数返回假基因和每次shuffle假基因中,三类假基因参与衍生lncrna的数量比例
pglendelnc <- function(pgdlncRNA_pg,hpg,shufflepgdlncRNA_pg,shufflehpg){
  #导入pgdlnc及其对应的假基因
  pgdlncRNA_pg <- read.delim(pgdlncRNA_pg)
  colnames(pgdlncRNA_pg)[2] <- "geneid"
  #导入假基因
  pg <- read.delim(hpg)%>%
    select(geneid,type)%>%
    mutate(times=1)
  #添加是否衍生lncRNA的属性
  pg <- left_join(pg,mutate(distinct(pgdlncRNA_pg,geneid),lnc=1),by="geneid")
  pg$lnc[is.na(pg$lnc)==T] <- 0
  #导入shufflepgdlnc及其对应的shuffle假基因
  shufflepgdlncRNA_pg <- read.delim(shufflepgdlncRNA_pg)
  colnames(shufflepgdlncRNA_pg)[2] <- "geneid"
  #导入shuffle假基因信息
  shufflepg <- read.delim(shufflehpg,header = F)%>%
    select(V7)%>%
    separate(V7,into = c("o","tw","th","times"),sep = "_",remove = F)%>%
    unite("geneid",o,tw,th)%>%
    left_join(select(pg,geneid,type),by="geneid")%>%
    select(V7,type,times)%>%
    mutate(type=paste(type,"(shuffled)",sep = " "))
  colnames(shufflepg)[1] <- "geneid"
  #添加是否衍生lncRNA的属性
  shufflepg <- left_join(shufflepg,mutate(distinct(shufflepgdlncRNA_pg,geneid),lnc=1),by="geneid")
  shufflepg$lnc[is.na(shufflepg$lnc)==T] <- 0
  #合并
  he <- rbind(pg,shufflepg)
  #统计参与衍生lncRNA的假基因比例
  tj <- group_by(he,type,times)%>%
    summarise(delncnum=sum(lnc),allnum=n())%>%
    mutate(nonum=allnum-delncnum,delncratio=delncnum*100/allnum)
}
#人
hs_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Homo_sapiens.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Homo_sapiens_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Homo_sapiens.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Homo_sapiens_hpg_shuffle.sort.bed")%>%
  mutate(sp="Human")
#小鼠
mm_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Mus_musculus.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Mus_musculus_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Mus_musculus.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Mus_musculus_hpg_shuffle.sort.bed")%>%
  mutate(sp="Mouse")
#鸡
gg_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Gallus_gallus.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Gallus_gallus_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Gallus_gallus.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Gallus_gallus_hpg_shuffle.sort.bed")%>%
  mutate(sp="Chicken")
#斑马鱼
dr_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Danio_rerio.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Danio_rerio_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Danio_rerio.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Danio_rerio_hpg_shuffle.sort.bed")%>%
  mutate(sp="Zebrafish")
#猕猴
mac_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Macaca_mulatta.pgdlncRNA_pg.txt",
                     hpg="~/PG/pg_message/Macaca_mulatta_hpg.txt",
                     shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Macaca_mulatta.pgdlncRNA_pg.txt",
                     shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Macaca_mulatta_hpg_shuffle.sort.bed")%>%
  mutate(sp="Macaque")
#大鼠
rat_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Rattus_norvegicus.pgdlncRNA_pg.txt",
                     hpg="~/PG/pg_message/Rattus_norvegicus_hpg.txt",
                     shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Rattus_norvegicus.pgdlncRNA_pg.txt",
                     shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Rattus_norvegicus_hpg_shuffle.sort.bed")%>%
  mutate(sp="Rat")
#兔子
oc_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Oryctolagus_cuniculus.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Oryctolagus_cuniculus_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Oryctolagus_cuniculus.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Oryctolagus_cuniculus_hpg_shuffle.sort.bed")%>%
  mutate(sp="Rabbit")
#负鼠
op_pg <- pglendelnc(pgdlncRNA_pg="~/PG/lncRNA_class/Monodelphis_domestica.pgdlncRNA_pg.txt",
                    hpg="~/PG/pg_message/Monodelphis_domestica_hpg.txt",
                    shufflepgdlncRNA_pg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Monodelphis_domestica.pgdlncRNA_pg.txt",
                    shufflehpg="/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/Monodelphis_domestica_hpg_shuffle.sort.bed")%>%
  mutate(sp="Opossum")
#合并
pgratio <- rbind(hs_pg,mac_pg,mm_pg,rat_pg,oc_pg,op_pg,gg_pg,dr_pg)
#统计shuffle 1000次衍生lnc比例的均值和置信度区间
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(1024)
tj1 <- group_by(filter(pgratio,type=="Duplicated (shuffled)"|type=="Processed (shuffled)"|type=="Fragment (shuffled)"),sp,type)%>%
  summarise(mean=mean(delncratio),median=median(delncratio),
            confmin=boot::boot.ci(boot::boot(delncratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(delncratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][5])
tj2 <- group_by(filter(pgratio,type=="Duplicated"|type=="Processed"|type=="Fragment"),sp,type)%>%
  summarise(mean=mean(delncratio),median=median(delncratio))%>%
  mutate(confmin=mean,confmax=mean)
tj <- rbind(tj1,tj2)
data.table::fwrite(tj,
                   file ="/home/zhzhang/PG/1Identification/SP.3typepg_3typeshufflepg_dlncratio.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#检验假基因比shuffle衍生lnc的比例更高
#1000次shuffle结果合计与正常进行fishertest检验
fortest <- group_by(pgratio,sp,type)%>%
  summarise(delncnum=sum(delncnum),nonum=sum(nonum))
fortest$sp <- factor(fortest$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
fortest <- arrange(fortest,sp)
fisherdata <- data.frame("row"=seq(1,48,by=2),"type"=rep(c("Duplicated","Fragment",
                                                           "Processed"),4))
for (i in seq(1,48,by=2)) {
  fisherfor <- fisher.test(fortest[c(i,i+1),c(3,4)])[["p.value"]]
  fisherdata[fisherdata$row==i,3] <- fisherfor
}
fisherdata <- mutate(fisherdata,anno="un")
fisherdata$anno[fisherdata$V3>0.05] <- "N.S."
fisherdata$anno[fisherdata$V3>0.01 & fisherdata$V3<0.05] <- "*"
fisherdata$anno[fisherdata$V3>0.001 & fisherdata$V3<0.01] <- "**"
fisherdata$anno[fisherdata$V3<0.001] <- "***"
#plot
tj$type <- factor(tj$type,levels = c("Duplicated","Duplicated (shuffled)",
                                               "Processed","Processed (shuffled)",
                                               "Fragment","Fragment (shuffled)"))
tj$sp <- factor(tj$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
pp <- ggplot(data = tj,aes(x=sp,y=mean))+
  geom_col(aes(fill=type),width = 0.6,position = "dodge")+
  geom_errorbar(data = filter(tj,type=="Duplicated (shuffled)"),
                aes(ymin=confmin,ymax=confmax),
                width=0.05,position = position_nudge(x=-0.15))+
  geom_errorbar(data = filter(tj,type=="Processed (shuffled)"),
                aes(ymin=confmin,ymax=confmax),
                width=0.05,position = position_nudge(x=0.05))+
  geom_errorbar(data = filter(tj,type=="Fragment (shuffled)"),
                aes(ymin=confmin,ymax=confmax),
                width=0.05,position = position_nudge(x=0.25))+
  geom_signif(annotations=filter(fisherdata,type=="Duplicated")$anno,
              y_position=c(28,4.5,16,7,6,7.2,12,11),tip_length = 0,
              xmin = c(0.75,1.75,2.75,3.75,4.75,5.75,6.75,7.75),
              xmax = c(0.85,1.85,2.85,3.85,4.85,5.85,6.85,7.85))+
  geom_signif(annotations=filter(fisherdata,type=="Processed")$anno,
              y_position=c(5,1.8,3.5,2.5,2,2.6,4.5,4.5),tip_length = 0,
              xmin = c(0.95,1.95,2.95,3.95,4.95,5.95,6.95,7.95),
              xmax = c(1.05,2.05,3.05,4.05,5.05,6.05,7.05,8.05))+
  geom_signif(annotations=filter(fisherdata,type=="Fragment")$anno,
              y_position=c(5.5,2,2.5,2.6,2.3,2.9,5.5,3.5),tip_length = 0,
              xmin = c(1.15,2.15,3.15,4.15,5.15,6.15,7.15,8.15),
              xmax = c(1.25,2.25,3.25,4.25,5.25,6.25,7.25,8.25))+
  scale_fill_manual(values = c("#0073C2","#7FB8E0","#EFC000",
                               "#F7DF7F","#868686","#C2C2C2"))+
  scale_y_continuous(breaks = c(0,5,10,15,20,25))+
  theme_half_open()+
  labs(x = NULL, y ="Pseudogenes\ninvolved in lncRNA origins (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 18)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(legend.position = "top", legend.direction = "horizontal")
ggsave("/home/zhzhang/PG/1Identification/SP.3typepg_3typeshufflepg_dlncratio.pdf", 
       pp,width = 9, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
##### 阴性对照
```r
###shuffle和正常假基因 产生基因间区的假基因比例对比
#pgdlncRNA_pg输入pgd基因间区及其对应的假基因，hpg输入假基因信息，
#shufflepgdlncRNA_pg输入shufflepg基因间区及其对应的shuffle假基因，shufflehpg输入shuffle假基因
#函数返回假基因和每次shuffle假基因，三类假基因参与衍生基因间区的数量比例
pglendelnc <- function(pgdlncRNA_pg,hpg,shufflepgdlncRNA_pg,shufflehpg){
  #导入pgdlnc及其对应的假基因
  pgdlncRNA_pg <- read.delim(pgdlncRNA_pg)
  colnames(pgdlncRNA_pg)[2] <- "geneid"
  #导入假基因
  pg <- read.delim(hpg)%>%
    select(geneid,type)%>%
    mutate(times=1)
  #添加是否衍生lncRNA的属性
  pg <- left_join(pg,mutate(distinct(pgdlncRNA_pg,geneid),lnc=1),by="geneid")
  pg$lnc[is.na(pg$lnc)==T] <- 0
  #导入shufflepgdlnc及其对应的shuffle假基因
  shufflepgdlncRNA_pg <- read.delim(shufflepgdlncRNA_pg)
  colnames(shufflepgdlncRNA_pg)[2] <- "geneid"
  #导入shuffle假基因信息
  shufflepg <- read.delim(shufflehpg,header = F)%>%
    select(V7)%>%
    separate(V7,into = c("o","tw","th","times"),sep = "_",remove = F)%>%
    unite("geneid",o,tw,th)%>%
    left_join(select(pg,geneid,type),by="geneid")%>%
    select(V7,type,times)%>%
    mutate(type=paste(type,"(shuffled)",sep = " "))
  colnames(shufflepg)[1] <- "geneid"
  #添加是否衍生lncRNA的属性
  shufflepg <- left_join(shufflepg,mutate(distinct(shufflepgdlncRNA_pg,geneid),lnc=1),by="geneid")
  shufflepg$lnc[is.na(shufflepg$lnc)==T] <- 0
  #合并
  he <- rbind(pg,shufflepg)
  #统计参与衍生lncRNA的假基因比例
  tj <- group_by(he,type,times)%>%
    summarise(delncnum=sum(lnc),allnum=n())%>%
    mutate(nonum=allnum-delncnum,delncratio=delncnum*100/allnum)
}
#
sp <- c("Monodelphis_domestica",
        "Rattus_norvegicus",
        "Oryctolagus_cuniculus",
        "Macaca_mulatta",
        "Gallus_gallus",
        "Danio_rerio",
        "Homo_sapiens",
        "Mus_musculus")
pgratio <- data.frame()
for (i in 1:8) {
  forsp <- sp[i]
  hs_pg <- pglendelnc(pgdlncRNA_pg=paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".pgdintergenic_pg.txt",sep=""),
                      hpg=paste("/home/zhzhang/PG/pg_message/",forsp,"_hpg.txt",sep=""),
                      shufflepgdlncRNA_pg=paste("/home/zhzhang/PG/1Identification/intergenic/",forsp,".pgdintergenic_shufflepg.txt",sep=""),
                      shufflehpg=paste("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/",forsp,"_hpg_shuffle.sort.bed",sep=""))%>%
    mutate(sp=forsp)
  #合并
  pgratio <- rbind(pgratio,hs_pg)
}

#统计shuffle 1000次衍生基因间区比例的均值和置信度区间
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(1024)
tj1 <- group_by(filter(pgratio,type=="Duplicated (shuffled)"|type=="Processed (shuffled)"|type=="Fragment (shuffled)"),sp,type)%>%
  summarise(mean=mean(delncratio),median=median(delncratio),
            confmin=boot::boot.ci(boot::boot(delncratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(delncratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][5])
tj2 <- group_by(filter(pgratio,type=="Duplicated"|type=="Processed"|type=="Fragment"),sp,type)%>%
  summarise(mean=mean(delncratio),median=median(delncratio))%>%
  mutate(confmin=mean,confmax=mean)
tj <- rbind(tj1,tj2)
data.table::fwrite(tj,
                   file ="/home/zhzhang/PG/1Identification/intergenic/SP.3typepg_3typeshufflepg_dintergenicratio.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)

#检验假基因比shuffle衍生基因间区的比例差异
#1000次shuffle结果合计与正常进行fishertest检验
fortest <- group_by(pgratio,sp,type)%>%
  summarise(delncnum=sum(delncnum),nonum=sum(nonum))
fortest$sp[fortest$sp=="Homo_sapiens"] <- "Human"
fortest$sp[fortest$sp=="Danio_rerio"] <- "Zebrafish"
fortest$sp[fortest$sp=="Gallus_gallus"] <- "Chicken"
fortest$sp[fortest$sp=="Macaca_mulatta"] <- "Macaque"
fortest$sp[fortest$sp=="Monodelphis_domestica"] <- "Opossum"
fortest$sp[fortest$sp=="Mus_musculus"] <- "Mouse"
fortest$sp[fortest$sp=="Oryctolagus_cuniculus"] <- "Rabbit"
fortest$sp[fortest$sp=="Rattus_norvegicus"] <- "Rat"
fortest$sp <- factor(fortest$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
fortest <- arrange(fortest,sp)
fisherdata <- data.frame("row"=seq(1,48,by=2),"type"=rep(c("Duplicated","Fragment",
                                                           "Processed"),4))
for (i in seq(1,48,by=2)) {
  fisherfor <- fisher.test(fortest[c(i,i+1),c(3,4)])[["p.value"]]
  fisherdata[fisherdata$row==i,3] <- fisherfor
}
fisherdata <- mutate(fisherdata,anno="un")
fisherdata$anno[fisherdata$V3>0.05] <- "N.S."
fisherdata$anno[fisherdata$V3>0.01 & fisherdata$V3<0.05] <- "*"
fisherdata$anno[fisherdata$V3>0.001 & fisherdata$V3<0.01] <- "**"
fisherdata$anno[fisherdata$V3<0.001] <- "***"
#plot
tj$sp[tj$sp=="Homo_sapiens"] <- "Human"
tj$sp[tj$sp=="Danio_rerio"] <- "Zebrafish"
tj$sp[tj$sp=="Gallus_gallus"] <- "Chicken"
tj$sp[tj$sp=="Macaca_mulatta"] <- "Macaque"
tj$sp[tj$sp=="Monodelphis_domestica"] <- "Opossum"
tj$sp[tj$sp=="Mus_musculus"] <- "Mouse"
tj$sp[tj$sp=="Oryctolagus_cuniculus"] <- "Rabbit"
tj$sp[tj$sp=="Rattus_norvegicus"] <- "Rat"
tj$type <- factor(tj$type,levels = c("Duplicated","Duplicated (shuffled)",
                                     "Processed","Processed (shuffled)",
                                     "Fragment","Fragment (shuffled)"))
tj$sp <- factor(tj$sp,levels = c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish"))
pp <- ggplot(data = tj,aes(x=sp,y=mean))+
  geom_col(aes(fill=type),width = 0.6,position = "dodge")+
  geom_errorbar(data = filter(tj,type=="Duplicated (shuffled)"),
                aes(ymin=confmin,ymax=confmax),
                width=0.05,position = position_nudge(x=-0.15))+
  geom_errorbar(data = filter(tj,type=="Processed (shuffled)"),
                aes(ymin=confmin,ymax=confmax),
                width=0.05,position = position_nudge(x=0.05))+
  geom_errorbar(data = filter(tj,type=="Fragment (shuffled)"),
                aes(ymin=confmin,ymax=confmax),
                width=0.05,position = position_nudge(x=0.25))+
  geom_signif(annotations=filter(fisherdata,type=="Duplicated")$anno,
              y_position=c(5,5,6,5.7,5.5,5.6,14.5,11),tip_length = 0,
              xmin = c(0.75,1.75,2.75,3.75,4.75,5.75,6.75,7.75),
              xmax = c(0.85,1.85,2.85,3.85,4.85,5.85,6.85,7.85))+
  geom_signif(annotations=filter(fisherdata,type=="Processed")$anno,
              y_position=c(2.2,2.2,2.2,2.2,2.2,2.2,5.5,4.5),tip_length = 0,
              xmin = c(0.95,1.95,2.95,3.95,4.95,5.95,6.95,7.95),
              xmax = c(1.05,2.05,3.05,4.05,5.05,6.05,7.05,8.05))+
  geom_signif(annotations=filter(fisherdata,type=="Fragment")$anno,
              y_position=c(1.8,1.8,1.8,1.8,1.8,1.8,5.5,3.3),tip_length = 0,
              xmin = c(1.15,2.15,3.15,4.15,5.15,6.15,7.15,8.15),
              xmax = c(1.25,2.25,3.25,4.25,5.25,6.25,7.25,8.25))+
  scale_fill_manual(values = c("#0073C2","#7FB8E0","#EFC000",
                               "#F7DF7F","#868686","#C2C2C2"))+
  scale_y_continuous(breaks = c(0,5,10,15,20,25))+
  theme_half_open()+
  labs(x = NULL, y ="Pseudogenes involved in\nintergenic region origins (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 18)) +
  theme(legend.text = element_text(size = 14)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(legend.position = "top", legend.direction = "horizontal")
ggsave("/home/zhzhang/PG/1Identification/intergenic/SP.3typepg_3typeshufflepg_dintergenicratio.pdf", 
       pp,width = 13, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
