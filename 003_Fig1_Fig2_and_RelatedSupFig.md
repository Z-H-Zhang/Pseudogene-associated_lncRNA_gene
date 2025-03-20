# 【UP】1Pg-lnc


### 三.Figure1和Figure2相关结果
##### 1.各种lncRNA基因的数量统计
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






lncRNA_source <- read.csv("/home/zhzhang/PG/1Identification/SP.lncRNA_source.csv")
lncRNA_source$sp <- factor(lncRNA_source$sp,levels = rev(c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish")))
lncRNA_source$source <- factor(lncRNA_source$source,levels = c("Novel","ENSEMBL"))
#plot
p2 <- ggplot(data = lncRNA_source,aes(y=sp,x=num/1000))+
  geom_col(width = 0.8,alpha=0.7,aes(fill=source,color=source))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  theme_half_open()+
  labs(y =NULL , x ="No. of lncRNA genes (×1E+3)",fill=NULL,color=NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,41))+
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = c(0.78, 0.1))+
  theme(legend.text = element_text(size = 15))

ggsave("/home/zhzhang/PG/1Identification/sp_lncRNAsource_num.pdf",
       p2,width =9, height =4.5,dpi=1200, units = "in", device='pdf',bg = "transparent")




#palnc数量展示
#导入lncRNA类型
splist1 <- c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")
splist2 <- c("Danio_rerio","Gallus_gallus","Monodelphis_domestica","Oryctolagus_cuniculus","Rattus_norvegicus","Mus_musculus","Macaca_mulatta","Homo_sapiens")
all <- data.frame()
for (i in 1:8) {
  lncRNA_class <- read.delim(paste("/home/zhzhang/PG/lncRNA_class/new_r1/",splist2[i],".lncRNA_class.txt",sep=""))%>%
    filter(type!="Non-pseudogene-associated lncRNA")%>%
    mutate(sp=splist1[i])
  all=rbind(all,lncRNA_class)
}
#统计PAlnc总数
pglncnum <- group_by(all,sp)%>%
  summarise(num=n())%>%data.frame()
pglncnum$sp <- factor(pglncnum$sp,levels = rev(c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish")))
#展示palnc总数
p2 <- ggplot(data = pglncnum,aes(y=sp,x=num))+
  geom_point(aes(fill=sp,color=sp,size=num))+
  geom_col(width = 0.01,alpha=0.7,aes(fill=sp,color=sp))+
  geom_text(aes(label=num),color="white")+
  scale_size(range = c(9,20))+
  scale_x_continuous(expand = c(0,0),limits = c(0,2700),breaks = c(0,500,1000,1500,2000,2500))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  theme_half_open()+
  labs(y =NULL , x ="No. of pseudogene-associated lncRNA genes",fill=NULL,color=NULL) +
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = "none")+
  theme(legend.text = element_text(size = 15))
ggsave("/home/zhzhang/PG/1Identification/sp_palncnum.pdf",
       p2,width =6.5, height =5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#统计PAlnc两种类型的数量
patypenum <- group_by(all,sp,type)%>%
  summarise(typenum=n())%>%data.frame()
tjmerge <- left_join(patypenum,pglncnum,by="sp")%>%
  mutate(percentage=typenum*100/num)
data.table::fwrite(tjmerge,
                   file ="/home/zhzhang/PG/1Identification/SP.palnctype.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
tjmerge$sp <- factor(tjmerge$sp,levels = rev(c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish")))
tjmerge$type <- factor(tjmerge$type,levels = c("Pseudogene-associated antisense lncRNA","Pseudogene-associated sense lncRNA"))
p2 <- ggplot(data = tjmerge,aes(y=sp,x=percentage))+
  geom_col(width = 0.6,alpha=0.7,aes(fill=type,color=type))+
  geom_text(data = filter(tjmerge,type=="Pseudogene-associated sense lncRNA"),x=13,color="white",size=5,
            aes(label=paste(round(percentage,1),"%",sep = "")))+
  geom_text(data = filter(tjmerge,type=="Pseudogene-associated antisense lncRNA"),x=87,color="white",size=5,
            aes(label=paste(round(percentage,1),"%",sep = "")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-derived lncRNA","Random intergenic"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-derived lncRNA","Random intergenic"))+
  scale_x_continuous(labels = c("0","25","50","75","100"),expand = c(0,0),limits = c(0,105))+
  theme_half_open()+
  labs(y =NULL , x ="Percentage (%)",fill=NULL,color=NULL) +
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 0),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = "none")+
  theme(legend.text = element_text(size = 15)) +
  theme(axis.ticks.y = element_line(size = 0)) + 
  theme(axis.ticks.y = element_line(colour = "white"))
ggsave("/home/zhzhang/PG/1Identification/sp_palnctypenum.pdf",
       p2,width =3, height =5,dpi=1200, units = "in", device='pdf',bg = "transparent")




#palnc来源展示
splist1 <- c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")
splist2 <- c("Danio_rerio","Gallus_gallus","Monodelphis_domestica","Oryctolagus_cuniculus","Rattus_norvegicus","Mus_musculus","Macaca_mulatta","Homo_sapiens")
all <- data.frame()
for (i in 1:8) {
  assign("pglnc",read.delim(paste("/home/zhzhang/PG/lncRNA_class/new_r1/",splist2[i],".lncRNA_class.txt",sep=""))%>%
           filter(type!="Non-pseudogene-associated lncRNA")%>%
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
                   file ="/home/zhzhang/PG/1Identification/SP.palncsource.tj.txt",
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
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = "none")+
  theme(legend.text = element_text(size = 15)) 
ggsave("/home/zhzhang/PG/1Identification/sp_palncsource.pdf",
       p2,width =5, height =5,dpi=1200, units = "in", device='pdf',bg = "transparent")



#palnc比例
ls /home/zhzhang/PG/lncRNA_class/new_r1/|awk -F '.' '{print $1}'|sort -u|while read i; do palnc=`grep -w "Pseudogene-associated" "/home/zhzhang/PG/lncRNA_class/new_r1/${i}.lncRNA_class.txt"|wc -l`; alllnc=`tail -n +2 "/home/zhzhang/PG/lncRNA_class/new_r1/${i}.lncRNA_class.txt"|wc -l`; awk -v sp=$i -v fenzi=$palnc -v fenmu=$alllnc '{print sp"\t"fenzi/fenmu}' "/home/zhzhang/PG/lncRNA_class/new_r1/${i}.lncRNA_class.txt"; done|uniq
#
Danio_rerio     0.0827594
Gallus_gallus   0.0239817
Homo_sapiens    0.0622893
Macaca_mulatta  0.029708
Monodelphis_domestica   0.0359275
Mus_musculus    0.0523687
Oryctolagus_cuniculus   0.0436009
Rattus_norvegicus       0.0407668

```








##### 2.假基因数量,密度和比例展示
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
  labs(y =NULL , x ="No. of Pseudogenes (×1E+3)",fill=NULL,color=NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,21))+
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = c(0.75, 0.15))+
  theme(legend.text = element_text(size = 15))
ggsave("/home/zhzhang/PG/1Identification/sp_pgnum.pdf",
       p1,width = 9, height = 4.5,dpi=1200, units = "in", device='pdf',bg = "transparent")



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


```


##### 3.ISO-seq质量表征，对比Iso-seq和ref annotation的转录本长度以及外显子数量
```r
#提取外显子bed，为求REF转录本注释长度做准备【chr start end geneid transcriptid strand】 
#斑马鱼
tail -n +6 "/home/zhzhang/PG/HPG/GTF_PG/Danio_rerio.GRCz11.108.chr.rmpg.gtf" |grep -w exon|awk '$13=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/ISOseqdata/A_translen/Danio_rerio.ensemble.allgeneexon.bed
#小鼠
tail -n +6 "/home/zhzhang/PG/HPG/GTF_PG/Mus_musculus.GRCm39.108.chr.rmpg.gtf" |grep -w exon|awk '$13=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/ISOseqdata/A_translen/Mus_musculus.ensemble.allgeneexon.bed
#鸡
tail -n +6 "/home/zhzhang/PG/HPG/GTF_PG/Gallus_gallus.GRCg7b.108.chr.rmpg.gtf" |grep -w exon|awk '$13=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/ISOseqdata/A_translen/Gallus_gallus.ensemble.allgeneexon.bed


#提取ISo-seq转录本注释的外显子bed
#斑马鱼
cat "/home/zhzhang/PG/ISOseqdata/Danio_rerio/gtf/drALL.collapsed.gff" |grep -w exon|awk '{print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/ISOseqdata/A_translen/Danio_rerio.isoseq.allgeneexon.bed
#小鼠
cat "/home/zhzhang/PG/ISOseqdata/Mus_musculus/gtf/mmALL.collapsed.filtered.gff" |grep -w exon|awk '{print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/ISOseqdata/A_translen/Mus_musculus.isoseq.allgeneexon.bed
#鸡
cat "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/gtf/ggALL.collapsed.filtered.gff" |grep -w exon|awk '{print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/ISOseqdata/A_translen/Gallus_gallus.isoseq.allgeneexon.bed



```
```r
#对比两类annotation转录本长度
#计算三类基因每个转录本的长度/外显子数量
  #导入isoseq外显子bed文件
  splist1 <- c("Zebrafish","Chicken","Mouse")
  splist2 <- c("Danio_rerio","Gallus_gallus","Mus_musculus")
  for (i in 1:3) {
    assign(paste(splist1[i],"iso",sep = "_"),read.delim(paste("/home/zhzhang/PG/isoseq/A_translen/",splist2[i],".isoseq.allgeneexon.bed",sep=""), header=FALSE)%>%
             select(-1,-6)%>%
             mutate(sp=splist1[i])%>%
             mutate(type="PacBio"))
        gc()
  }
  isomerge <- rbind(Zebrafish_iso,Chicken_iso,Mouse_iso)
  colnames(isomerge) <- c("start","end","geneid","transcriptid","sp","type")
  #导入ref外显子bed文件
  for (i in 1:3) {
    assign(paste(splist1[i],"ref",sep = "_"),read.delim(paste("/home/zhzhang/PG/isoseq/A_translen/",splist2[i],".ensemble.allgeneexon.bed",sep=""), header=FALSE)%>%
             select(-1,-6)%>%
             mutate(sp=splist1[i])%>%
             mutate(type="ENSEMBL"))
    gc()
  }
  refmerge <- rbind(Zebrafish_ref,Chicken_ref,Mouse_ref)
  colnames(refmerge) <- c("start","end","geneid","transcriptid","sp","type")
  #
  alll <- rbind(isomerge,refmerge)
  #计算每个外显子长度，再计算每个基因每个转录本长度
  alltran_len <- mutate(alll,len=end-start)%>%
    group_by(sp,type,geneid,transcriptid)%>%
    summarise(len=sum(len),exonnum=n())
#
  alltran_len$type <- factor(alltran_len$type,
                          levels = c("PacBio","ENSEMBL"))
  alltran_len$sp <- factor(alltran_len$sp,
                             levels = c("Zebrafish","Chicken","Mouse"))
#统计储存
tj <- group_by(alltran_len,sp,type)%>%
  summarise(mediannum=median(exonnum),medianlen=median(len))
data.table::fwrite(tj,file ="/home/zhzhang/PG/isoseq/A_translen/3SP_iso_ref_transcript_len_exonnum.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot_translen
isotranslen <- ggplot(data = alltran_len,aes(x=type,y=log10(len)))+
  geom_violin(width=0.9,aes(fill=type,color=type))+
  geom_boxplot(width=0.1,outlier.alpha = 0)+
  geom_signif(map_signif_level=T,y_position = c(5.5,5),
              comparisons=list(c("PacBio","ENSEMBL")))+
  scale_fill_manual(values=c("#D8838C","#F8CA7E"),
                    limits=c("PacBio","ENSEMBL"))+
  scale_color_manual(values=c("#D8838C","#F8CA7E"),
                     limits=c("PacBio","ENSEMBL"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0.5, 6))+
  scale_y_continuous(breaks = c(1,2,3,4,5),labels = c("1","2","3","4","5"))+
  scale_x_discrete(labels = c("PacBio","ENSEMBL"))+
  labs(x = NULL, y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[10]*"("*"b"*"p"*")"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  facet_wrap(~sp,nrow=1,strip.position="top",as.table=F)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))

ggsave("/home/zhzhang/PG/isoseq/A_translen/3SP_iso_ref_transcript_len.pdf", 
       isotranslen,width = 4, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")

#plot_transexonnum
isotransnum <- ggplot(data = alltran_len,aes(x=type,y=log10(exonnum)))+
  geom_violin(width=0.9,aes(fill=type,color=type))+
  geom_boxplot(width=0.1,outlier.alpha = 0)+
  geom_signif(map_signif_level=T,y_position = c(2.7),
              comparisons=list(c("PacBio","ENSEMBL")))+
  scale_fill_manual(values=c("#D8838C","#F8CA7E"),
                    limits=c("PacBio","ENSEMBL"))+
  scale_color_manual(values=c("#D8838C","#F8CA7E"),
                     limits=c("PacBio","ENSEMBL"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 3))+
  scale_x_discrete(labels = c("PacBio","ENSEMBL"))+
  labs(x = NULL, y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"e"*"x"*"o"*"n"~"n"*"u"*"m"*"b"*"e"*"r"~"("*"l"*"o"*"g"[10]*"("*"c"*"o"*"u"*"n"*"t"*")"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  facet_wrap(~sp,nrow=1,strip.position="top",as.table=F)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))

ggsave("/home/zhzhang/PG/isoseq/A_translen/3SP_iso_ref_transcript_exonnum.pdf", 
       isotransnum,width = 4, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 4.LNCBOOK compare
```r
#LNCBOOK数据库gtf-bed(去除chrm的五条)
cat "/share/home/zhzhang24/PG/Z_DF/lncbook/lncRNA_LncBookv2.0_GRCh38.gtf"|awk '$1!="chrM" && $3=="gene" {print $1"\t"$4"\t"$5"\t"$10"\tlncbook\t"$7}'|sed "s/^chr//g;s/\"//g;s/\;//g" > /share/home/zhzhang24/PG/Z_DF/lncbook/lncRNA_LncBook.bed
#our human lnc bed
tail -n +2 "/share/home/zhzhang24/PG/lncRNA_class/new_r1/Homo_sapiens.lncRNA_class.txt"|awk '{print $1}' > /share/home/zhzhang24/PG/Z_DF/lncbook/human.lnrna.geneid.txt
grep -w -Ff "/share/home/zhzhang24/PG/Z_DF/lncbook/human.lnrna.geneid.txt" /share/home/zhzhang24/PG/Evolution/conserve/Homo_sapiens.allgene.bed > /share/home/zhzhang24/PG/Z_DF/lncbook/human.lnrna.bed
#相较于lncbook的新lnc基因
bedtools intersect -a "/share/home/zhzhang24/PG/Z_DF/lncbook/human.lnrna.bed" -b "/share/home/zhzhang24/PG/Z_DF/lncbook/lncRNA_LncBook.bed" -v > /share/home/zhzhang24/PG/Z_DF/lncbook/NovellncRNA_compareLncBook.bed

#统计相较于lncbook新的人类基因的来源
Novel <- read.delim("~/PG/Z_DF/lncbook/NovellncRNA_compareLncBook.bed", header=FALSE)
#三类来源基因数量统计
num1=length(grep("chr",Novel$V4))+length(grep("[.]",Novel$V4))
num2=length(grep("Hum_XLOC",Novel$V4))
num3=length(Novel$V4)-(num1+num2)
tj=data.frame(source=c("ENSEMBL",
                       "Long-read RNA-seq",
                       "Short-read RNA-seq"),
              num=c(num3,num1,num2))%>%
  mutate(label=paste(source," (",num,")",sep=""))
tj$source=factor(tj$source,levels = c("ENSEMBL",
                                      "Long-read RNA-seq",
                                      "Short-read RNA-seq"))
#plot
pp=ggplot(data=tj,aes(x=1,y=num))+
  geom_col(aes(fill=label),position="fill")+
  coord_polar(theta="y")+
  scale_fill_manual(values = c("#F8CA7E","#D8838C","#7DC6EB"))+
  theme_void()+
  labs(fill=NULL,color=NULL) +
  theme(legend.text = element_text(size = 11))
ggsave("/home/zhzhang/PG/Z_DF/lncbook/novel_lncRNA_source_num.pdf",
       pp,width =5, height =4,dpi=1200, units = "in", device='pdf',bg = "transparent")

  


```
##### 5.各种注释来源的lncRNA对比，展示长读长转录组优势
```r
#对比三类不同注释来源的lnc转录本长度
#函数计算三类基因每个转录本的长度
#a输入全部基因外显子bed，b输入基因分类文件,c输入物种名
gettlen <- function(a,b,c){
  #导入外显子bed文件
  allgeneexon <- read.delim(a, header=FALSE)%>%
    select(-1,-6)
  colnames(allgeneexon) <- c("start","end","geneid","transcriptid")
  #导入基因分类文件
  geneid_class <- read.delim(b)
  #计算每个外显子长度，再计算每个基因每个转录本长度
  alltran_len <- mutate(allgeneexon,len=end-start)%>%
    group_by(geneid,transcriptid)%>%
    summarise(len=sum(len))%>%
    data.frame()%>%
    group_by(geneid)%>%
    summarise(len=max(len))%>%
    data.frame()
  #合并分类和转录本长度
  geneclass_tranlen <- left_join(geneid_class,alltran_len,by="geneid")%>%
    mutate(sp=c)
}
#人
hstlen <- gettlen(a="~/PG/Evolution/tran_len/Homo_sapiens.allgeneexon.bed",
                  b="/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.lncRNA_class.txt",
                  c="Human")
#小鼠
mmtlen <- gettlen("~/PG/Evolution/tran_len/Mus_musculus.allgeneexon.bed",
                  "/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.lncRNA_class.txt",
                  "Mouse")
#鸡
ggtlen <- gettlen("/home/zhzhang/PG/Evolution/tran_len/Gallus_gallus.allgeneexon.bed",
                  "/home/zhzhang/PG/lncRNA_class/new_r1/Gallus_gallus.lncRNA_class.txt",
                  "Chicken")
#斑马鱼
drtlen <- gettlen("/home/zhzhang/PG/Evolution/tran_len/Danio_rerio.allgeneexon.bed",
                  "/home/zhzhang/PG/lncRNA_class/new_r1/Danio_rerio.lncRNA_class.txt",
                  "Zebrafish")
#合并不同物种信息
allsp_tlen <- rbind(hstlen,mmtlen,ggtlen,drtlen)
#添加注释来源
allsp_tlen$type="ENSEMBL"
allsp_tlen$type[grep("[.]",allsp_tlen$geneid)]="Long-read RNA-seq"
allsp_tlen$type[grep("PB",allsp_tlen$geneid)]="Long-read RNA-seq"
allsp_tlen$type[grep("chr",allsp_tlen$geneid)]="Long-read RNA-seq"
allsp_tlen$type[grep("XLOC",allsp_tlen$geneid)]="Short-read RNA-seq"
tj=group_by(allsp_tlen,sp,type)%>%
  summarise(median=median(len))
#
allsp_tlen$type <- factor(allsp_tlen$type,
                           levels = c("ENSEMBL","Long-read RNA-seq","Short-read RNA-seq"))
allsp_tlen$sp <- factor(allsp_tlen$sp,
                         levels = rev(c("Zebrafish","Chicken","Mouse","Human")))
#plot_translen
isotranslen <- ggplot(data = allsp_tlen,aes(x=type,y=log10(len)))+
  geom_violin(width=0.9,aes(fill=type,color=type))+
  geom_boxplot(width=0.1,outlier.alpha = 0)+
  ggsignif::geom_signif(map_signif_level=T,y_position = c(5.7,5.4),tip_length = 0.01,
              comparisons=list(c("Long-read RNA-seq","ENSEMBL"),
                               c("Long-read RNA-seq","Short-read RNA-seq")))+
  scale_fill_manual(values = c("#F8CA7E","#D8838C","#7DC6EB"))+
  scale_color_manual(values = c("#F8CA7E","#D8838C","#7DC6EB"))+
  theme_half_open()+
  coord_cartesian(ylim = c(1.5, 6))+
  scale_y_continuous(breaks = c(1,2,3,4,5,6),labels = c("1","2","3","4","5","6"))+
  labs(x = NULL, y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[10]*"("*"b"*"p"*")"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 13),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  facet_wrap(~sp,nrow=1,strip.position="top",as.table=F)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))
ggsave("/home/zhzhang/PG/isoseq/A_translen/4SP_3lncsource_transcript_len.pdf", 
       isotranslen,width = 7.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")








#对比三类lncRNA外显子的保守性
#人类
a <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paalnc_lncRNA_exon.merge.phastCons.txt"
ab <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.npalnc_lncRNA_exon.merge.phastCons.txt"
#小鼠
a1 <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paslnc_lncRNA_exon.merge.phastCons.txt"
b1 <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paalnc_lncRNA_exon.merge.phastCons.txt"
ab1 <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.npalnc_lncRNA_exon.merge.phastCons.txt"
#导入paslnc外显子的保守性打分
paslnc_lncRNA_exon <- rbind(read.delim(a, header=FALSE)%>%mutate(sp="Human"),
                            read.delim(a1, header=FALSE)%>%mutate(sp="Mouse"))%>%
  filter(V3!=0)%>%
  select(1,3,4,7)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(sp,geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,2,5)
#导入paalnc外显子的保守性打分
paalnc_lncRNA_exon <- rbind(read.delim(b, header=FALSE)%>%mutate(sp="Human"),
                            read.delim(b1, header=FALSE)%>%mutate(sp="Mouse"))%>%
  filter(V3!=0)%>%
  select(1,3,4,7)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(sp,geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,2,5)
#导入npa lnc外显子的保守性打分
npalnc_lncRNA_exon <- rbind(read.delim(ab, header=FALSE)%>%mutate(sp="Human"),
                            read.delim(ab1, header=FALSE)%>%mutate(sp="Mouse"))%>%
  filter(V3!=0)%>%
  select(1,3,4,7)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(sp,geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,2,5)
#合并
alldata <- rbind(paslnc_lncRNA_exon,paalnc_lncRNA_exon,npalnc_lncRNA_exon)
#添加注释来源
alldata$type="ENSEMBL"
alldata$type[grep("[.]",alldata$geneid)]="Long-read RNA-seq"
alldata$type[grep("PB",alldata$geneid)]="Long-read RNA-seq"
alldata$type[grep("chr",alldata$geneid)]="Long-read RNA-seq"
alldata$type[grep("XLOC",alldata$geneid)]="Short-read RNA-seq"
tj=group_by(alldata,sp,type)%>%
  summarise(median=median(score))
#PLOT
p1=ggplot(data = alldata,aes(x=type,y=score))+
  geom_boxplot(width=0.5,outlier.alpha = 0,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position = c(0.38,0.45),tip_length = 0.005,
                        comparisons=list(c("Long-read RNA-seq","ENSEMBL"),
                                         c("Long-read RNA-seq","Short-read RNA-seq")))+
  scale_fill_manual(values = c("#F8CA7E","#D8838C","#7DC6EB"))+
  scale_color_manual(values = c("#F8CA7E","#D8838C","#7DC6EB"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 0.5))+
  #scale_y_continuous(breaks = c(1,2,3,4,5,6),labels = c("1","2","3","4","5","6"))+
  labs(x = NULL, y ="PhastCons score",
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 13),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  facet_wrap(~sp,nrow=1,strip.position="top",as.table=F)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))
ggsave("/home/zhzhang/PG/isoseq/A_translen/2SP_3lncsource_phastcons.pdf", 
       p1,width = 4, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")




#对比基因的表达比例
#（a输入基因表达TPM矩阵，b输入基因分类文件，c输入物种名，d输出在物种中表达的基因id和类型文件）
three_gene_exppercent <- function(a,b,c,d){
  #导入基因表达TPM矩阵
  allsample_TPM <- read.delim(a, row.names=1)
  #导入基因ID分类
  geneid_class <- read.delim(b)
  geneid_class$type="ENSEMBL"
  geneid_class$type[grep("[.]",geneid_class$geneid)]="Long-read RNA-seq"
  geneid_class$type[grep("PB",geneid_class$geneid)]="Long-read RNA-seq"
  geneid_class$type[grep("chr",geneid_class$geneid)]="Long-read RNA-seq"
  geneid_class$type[grep("XLOC",geneid_class$geneid)]="Short-read RNA-seq"
  #循环计算不同阈值下，三类基因中表达的基因
  gene_data <- data.frame()
  for (i in c(0.5,1,1.5,2)) {
    #计算基因最大表达量，基因最大表达量TPM>=i视为表达，增加基因表达的定性属性，表达为1，不表达为0
    gene_exp <- data.frame(max=apply(allsample_TPM,1,max))%>%
      rownames_to_column("geneid")%>%
      mutate(exp=0)
    gene_exp$exp[gene_exp$max>=i] <- 1
    #合并基因类型和最大表达量
    gene <- left_join(geneid_class,select(gene_exp,1,3),by="geneid")%>%
      mutate(threshold=i)
    gene_data <- rbind(gene_data,gene)
  }
  #储存在物种中表达的基因id和类型(基因最大表达量(在至少一个样本中)TPM>=1视为表达)
  geneforsave <- filter(gene_data,threshold==1 & exp==1)%>%
    select(geneid,type)
  data.table::fwrite(geneforsave,file =d,sep = '\t',row.names = F,quote = F,col.names = T)
  #统计三类基因的表达比例
  genetj <- group_by(gene_data,threshold,type)%>%
    summarise(expnum=sum(exp),allnum=n())%>%
    mutate(spexpratio=expnum/allnum)
  #添加物种信息
  genetj <- mutate(genetj,sp=c)
  return(genetj)
}
#计算各物种三类基因表达比例
#斑马鱼Zebrafish
dr_spexp <- three_gene_exppercent(a ="/home/zhzhang/PG/RNAseq/Danio_rerio/allsample_TPM.txt",
                                  b = "/home/zhzhang/PG/lncRNA_class/new_r1/Danio_rerio.lncRNA_class.txt",
                                  c = "Zebrafish",
                                  d="/home/zhzhang/PG/isoseq/Danio_rerio.spexp_geneid_class.txt")
#human
hs_spexp <- three_gene_exppercent(a ="/home/zhzhang/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                                  b = "/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.lncRNA_class.txt",
                                  c = "Human",
                                  d="/home/zhzhang/PG/isoseq/Homo_sapiens.spexp_geneid_class.txt")
#mm
mm_spexp <- three_gene_exppercent(a ="/home/zhzhang/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                                  b ="/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.lncRNA_class.txt",
                                  c = "Mouse",
                                  d="/home/zhzhang/PG/isoseq/Mus_musculus.spexp_geneid_class.txt")
#gg
gg_spexp <- three_gene_exppercent(a ="/home/zhzhang/PG/RNAseq/Gallus_gallus/allsample_TPM.txt",
                                  b ="/home/zhzhang/PG/lncRNA_class/new_r1/Gallus_gallus.lncRNA_class.txt",
                                  c = "Chicken",
                                  d="/home/zhzhang/PG/isoseq/Gallus_gallus.spexp_geneid_class.txt")
#合并不同物种信息
allsp_tlen <- rbind(hs_spexp,mm_spexp,gg_spexp,dr_spexp)
allsp_tlen$type <- factor(allsp_tlen$type,
                          levels = c("ENSEMBL","Long-read RNA-seq","Short-read RNA-seq"))
allsp_tlen$sp <- factor(allsp_tlen$sp,
                        levels = rev(c("Zebrafish","Chicken","Mouse","Human")))
#画图
ph <- ggplot(data =allsp_tlen,aes(x=threshold,y=spexpratio*100))+
  geom_point(aes(color=type),size=2)+
  geom_line(aes(group=type,color=type))+
  scale_color_manual(values = c("#F8CA7E","#D8838C","#7DC6EB"))+
  theme_half_open()+
  labs(x = "Expression cutoff (TPM)", y ="Expressed genes (%)",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "top", legend.direction = "horizontal")+
  facet_wrap(~sp,nrow=1,strip.position="top",as.table=F,scales="free_y")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))
ggsave("/home/zhzhang/PG/isoseq/A_translen/4SP_3lncsource_expratio.pdf", 
       ph,width = 8, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


##画图对比各物种三类lnc基因的最大表达量
#统计三类基因（蛋白编码，两类lncrna基因）在物种内的最大表达量
#函数（a输入基因表达TPM矩阵，b输入在物种中表达的基因的分类文件,c输入物种名）
tgenemaxTPM <- function(a,b,c){
  #导入基因表达TPM矩阵，计算出每个基因在所有样本中的最大表达量
  allsample_TPM <- read.delim(a, row.names=1)
  gene_max_TPM <- apply(allsample_TPM,1,max)%>%
    data.frame()
  colnames(gene_max_TPM) <- "maxTPM"
  gene_max_TPM <- rownames_to_column(gene_max_TPM,"geneid")
  #基因分类和maxtpm合并
  geneid_class <- read.delim(b)%>%
    left_join(gene_max_TPM,by="geneid")%>%
    mutate(sp=c)
  return(geneid_class)
}
#小鼠
mm_TPM <- tgenemaxTPM(a = "/home/zhzhang/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/isoseq/Mus_musculus.spexp_geneid_class.txt",
                      c = "Mouse")
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/isoseq/Homo_sapiens.spexp_geneid_class.txt",
                      c = "Human")
#鸡
gg_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Gallus_gallus/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/isoseq/Gallus_gallus.spexp_geneid_class.txt",
                      c = "Chicken")
#斑马鱼
dr_TPM <- tgenemaxTPM(a = "/home/zhzhang/PG/RNAseq/Danio_rerio/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/isoseq/Danio_rerio.spexp_geneid_class.txt",
                      c = "Zebrafish")
#不同物种三类基因表达量数据合并
allsp_tlen <- rbind(hs_TPM,mm_TPM,gg_TPM,dr_TPM)
tj=group_by(allsp_tlen,sp,type)%>%
  summarise(median=median(maxTPM))
#
allsp_tlen$type <- factor(allsp_tlen$type,
                          levels = c("ENSEMBL","Long-read RNA-seq","Short-read RNA-seq"))
allsp_tlen$sp <- factor(allsp_tlen$sp,
                        levels = rev(c("Zebrafish","Chicken","Mouse","Human")))
phs <- ggplot(data = allsp_tlen,aes(x=type,y=log2(maxTPM)))+
  geom_boxplot(fatten = 3,notch = T,width=0.5,outlier.alpha = 0,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position = c(9.2,8.6),tip_length = 0.005,
                        comparisons=list(c("Long-read RNA-seq","ENSEMBL"),
                                         c("Long-read RNA-seq","Short-read RNA-seq")))+
  scale_fill_manual(values = c("#F8CA7E","#D8838C","#7DC6EB"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 10.5))+
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10,12.5),labels = c("0","2.5","5","7.5","10","12.5"))+
  labs(x = NULL, y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        legend.position = "none", legend.direction = "horizontal")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  facet_wrap(~sp,nrow=1,strip.position="top",as.table=F)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))
ggsave("/home/zhzhang/PG/isoseq/A_translen/4SP_3lncsource_maxexp.pdf", 
       phs,width = 8, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")




#三代注释的lnc数量
geneid_class <- rbind(read.delim("/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.lncRNA_class.txt")%>%mutate(sp="Human"),
                      read.delim("/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.lncRNA_class.txt")%>%mutate(sp="Mouse"),
                      read.delim("/home/zhzhang/PG/lncRNA_class/new_r1/Gallus_gallus.lncRNA_class.txt")%>%mutate(sp="Chicken"),
                      read.delim("/home/zhzhang/PG/lncRNA_class/new_r1/Danio_rerio.lncRNA_class.txt")%>%mutate(sp="Zebrafish"))
geneid_class$type="ENSEMBL"
geneid_class$type[grep("[.]",geneid_class$geneid)]="Long-read RNA-seq"
geneid_class$type[grep("PB",geneid_class$geneid)]="Long-read RNA-seq"
geneid_class$type[grep("chr",geneid_class$geneid)]="Long-read RNA-seq"
geneid_class$type[grep("XLOC",geneid_class$geneid)]="Short-read RNA-seq"
tj=filter(geneid_class,type=="Long-read RNA-seq")%>%
  group_by(sp)%>%
  summarise(num=n())%>%
  data.frame()
tj$sp=factor(tj$sp,levels = c("Human","Mouse","Chicken","Zebrafish"))
#plot
p2 <- ggplot(data = tj,aes(y=num,x=sp))+
  geom_point(aes(fill=sp,color=sp,size=num))+
  geom_col(width = 0.01,alpha=0.7,aes(fill=sp,color=sp))+
  geom_text(aes(label=num),color="white")+
  scale_size(range = c(9,20))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1100))+
  ggsci::scale_fill_aaas()+
  ggsci::scale_color_aaas()+
  theme_half_open()+
  labs(x =NULL , y ="No. of lncRNA genes",fill=NULL,color=NULL) +
  theme(axis.title = element_text(size = 20),axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 15)) +
  theme(legend.position = "none")+
  theme(legend.text = element_text(size = 15))+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/isoseq/A_translen/4SP_longreadLNCRNA_num.pdf",
       p2,width =4, height =5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
##### 6.物种lncRNA库的大小是影响假基因相关lncRNA基因数量的主要因素
```r
#lncRNA总数与palnc数量相关性
#lnc总数
lncRNA_source <- read.csv("/home/zhzhang/PG/1Identification/SP.lncRNA_source.csv")%>%
  group_by(sp)%>%
  summarise(allnum=sum(num))%>%
  data.frame()
#导入lncRNA类型
splist1 <- c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")
splist2 <- c("Danio_rerio","Gallus_gallus","Monodelphis_domestica","Oryctolagus_cuniculus","Rattus_norvegicus","Mus_musculus","Macaca_mulatta","Homo_sapiens")
all <- data.frame()
for (i in 1:8) {
  lncRNA_class <- read.delim(paste("/home/zhzhang/PG/lncRNA_class/new_r1/",splist2[i],".lncRNA_class.txt",sep=""))%>%
    filter(type!="Non-pseudogene-associated lncRNA")%>%
    mutate(sp=splist1[i])
  all=rbind(all,lncRNA_class)
}
#统计PAlnc总数
pglncnum <- group_by(all,sp)%>%
  summarise(num=n())%>%data.frame()
pglncnum$sp <- factor(pglncnum$sp,levels = rev(c("Human","Macaque","Mouse","Rat","Rabbit","Opossum","Chicken","Zebrafish")))
patypenum <- group_by(all,sp,type)%>%
  summarise(typenum=n())%>%data.frame()%>%
  mutate(type=case_when(type=="Pseudogene-associated antisense lncRNA" ~ "PAA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS"))%>%
  spread(type,typenum)
#合并
he=left_join(lncRNA_source,pglncnum,by="sp")%>%
  left_join(patypenum,by="sp")
#plot
p1 <- ggplot(data=he,aes(x=allnum,y=num))+
  geom_point(color="#99CCFF",size=3)+
  ggrepel::geom_text_repel(aes(label=sp),size=5)+
  geom_smooth(method="lm",alpha=0.1,color="#3A528A")+
  theme_half_open()+
  ggpubr::stat_cor(size=5,label.x=0.25*max(he$allnum),label.y =0.9*max(he$num))+
  ggpmisc::stat_poly_eq(size=5,label.x=0.4*max(he$allnum),label.y =0.96*max(he$num),geom="text",aes(label = ..eq.label..))+
  labs(x = "No. of lncRNA genes", y = "No. of pseudogene-associated lncRNA genes")+
  theme(axis.title = element_text(size = 18),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),plot.margin = margin(l=10,r = 15))
p2 <- ggplot(data=he,aes(x=allnum,y=PAS))+
  geom_point(color="#7197AD",size=3)+
  ggrepel::geom_text_repel(aes(label=sp),size=5)+
  geom_smooth(method="lm",alpha=0.1,color="#3A528A")+
  theme_half_open()+
  ggpubr::stat_cor(size=5,label.x=0.25*max(he$allnum),label.y =0.9*max(he$PAS))+
  ggpmisc::stat_poly_eq(size=5,label.x=0.4*max(he$allnum),label.y =0.96*max(he$PAS),geom="text",aes(label = ..eq.label..))+
  labs(x = "No. of lncRNA genes", y = "No. of PAS lncRNA genes")+
  theme(axis.title = element_text(size = 18),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),plot.margin = margin(l=10,r = 15))
p3 <- ggplot(data=he,aes(x=allnum,y=PAA))+
  geom_point(color="#628255",size=3)+
  ggrepel::geom_text_repel(aes(label=sp),size=5)+
  geom_smooth(method="lm",alpha=0.1,color="#3A528A")+
  theme_half_open()+
  ggpubr::stat_cor(size=5,label.x=0.25*max(he$allnum),label.y =0.9*max(he$PAA))+
  ggpmisc::stat_poly_eq(size=5,label.x=0.4*max(he$allnum),label.y =0.96*max(he$PAA),geom="text",aes(label = ..eq.label..))+
  labs(x = "No. of lncRNA genes", y = "No. of PAA lncRNA genes")+
  theme(axis.title = element_text(size = 18),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),plot.margin = margin(l=10,r = 15))
pp=plot_grid(p1,p2,p3,nrow=1)
ggsave("/home/zhzhang/PG/1Identification/lncRNAnum_cor_palncnum.pdf", pp,width = 18, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 7.证明lncRNA基因富集在假基因上
```r
#正常重排，将假基因在基因组上重排到非假基因以及非蛋白编码基因外显子区域(排除蛋白编码基因遗物影响，同时遵循假基因的位置分布准则，不改变每个假基因长度)
#人
#合并蛋白编码基因外显子区域+假基因
awk '{print $1"\t"$2"\t"$3}' /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg.bed > /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg_pcgexon_forshuffle.bed
awk '{print $2"\t"$3"\t"$4}' "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/Homo_sapiens.GRCh38.108.exLocs.txt" >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg_pcgexon_forshuffle.bed
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
#传至实验室服务/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg_shuffle.sort.bed /home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg_shuffle.sort.bed /home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg_shuffle.sort.bed /home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg_shuffle.sort.bed /home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg_shuffle.sort.bed /home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg_shuffle.sort.bed /home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg_shuffle.sort.bed /home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg_shuffle.sort.bed /home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/




#固定chr重排
cat "/share/home/zhzhang24/PG/sp.txt"|while read iii
do
chrsize=`ls /share/home/zhzhang24/PG/refgenome/${iii}*dna.chrsize.txt`
for i in $(seq 1 1000)
do
micromamba run -n SEQ bedtools shuffle -i "/share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg.bed" -g ${chrsize} -excl /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_pcgexon_forshuffle.bed -noOverlapping -chrom |awk -v time="${i}" '{print $0"_"time}' >> /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_shuffle.keepchr.bed
done
micromamba run -n SEQ bedtools sort -i /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_shuffle.keepchr.bed > /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_shuffle.keepchr.sort.bed
done
#传至实验室服务/home/zhzhang/PG/1Identification/shufflepg_keepchr_lncRNAclass/
cat "/share/home/zhzhang24/PG/sp.txt"|while read iii
do
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_shuffle.keepchr.sort.bed zhzhang@122.205.95.67:/home/zhzhang/PG/1Identification/shufflepg_keepchr_lncRNAclass/
done



#下载重复序列区域ucsc
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz -O /share/home/zhzhang24/PG/refgenome/Homo_sapiens.rmsk.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/database/rmsk.txt.gz -O /share/home/zhzhang24/PG/refgenome/Gallus_gallus.rmsk.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/database/rmsk.txt.gz -O /share/home/zhzhang24/PG/refgenome/Macaca_mulatta.rmsk.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/oryCun2/database/rmsk.txt.gz -O /share/home/zhzhang24/PG/refgenome/Oryctolagus_cuniculus.rmsk.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/rn7/database/rmsk.txt.gz -O /share/home/zhzhang24/PG/refgenome/Rattus_norvegicus.rmsk.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/monDom5/database/rmsk.txt.gz -O /share/home/zhzhang24/PG/refgenome/Monodelphis_domestica.rmsk.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/rmsk.txt.gz -O /share/home/zhzhang24/PG/refgenome/Mus_musculus.rmsk.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/database/rmsk.txt.gz -O /share/home/zhzhang24/PG/refgenome/Danio_rerio.rmsk.txt.gz
#转bed
gzip -d /share/home/zhzhang24/PG/refgenome/*rmsk.txt.gz
cat "/share/home/zhzhang24/PG/sp.txt"|while read i
do
awk '{print $6"\t"$7"\t"$8"\t"$12"\t"$13"\t"$10}' "/share/home/zhzhang24/PG/refgenome/${i}.rmsk.txt" > "/share/home/zhzhang24/PG/refgenome/${i}.rmsk.bed"
done
liftOver "/share/home/zhzhang24/PG/refgenome/Gallus_gallus.rmsk.bed" "/share/home/zhzhang24/software/liftover/galGal6ToGRCg7b.over.chain.gz" /share/home/zhzhang24/PG/refgenome/Gallus_gallus.rmsk.1.bed /dev/null
#NCBI染色体号转正常数字（ENSEMBLE格式）
cat /share/home/zhzhang24/PG/refgenome/Gallus_gallus.rmsk.1.bed|sed "s/NC_052532.1/1/g;s/NC_052533.1/2/g;s/NC_052534.1/3/g;s/NC_052535.1/4/g;s/NC_052536.1/5/g;s/NC_052537.1/6/g;s/NC_052538.1/7/g;s/NC_052539.1/8/g;s/NC_052540.1/9/g;s/NC_052541.1/10/g;s/NC_052542.1/11/g;s/NC_052543.1/12/g;s/NC_052544.1/13/g;s/NC_052545.1/14/g;s/NC_052546.1/15/g;s/NC_052547.1/16/g;s/NC_052548.1/17/g;s/NC_052549.1/18/g;s/NC_052550.1/19/g;s/NC_052551.1/20/g;s/NC_052552.1/21/g;s/NC_052553.1/22/g;s/NC_052554.1/23/g;s/NC_052555.1/24/g;s/NC_052556.1/25/g;s/NC_052557.1/26/g;s/NC_052558.1/27/g;s/NC_052559.1/28/g;s/NC_052560.1/29/g;s/NC_052561.1/30/g;s/NC_052562.1/31/g;s/NC_052563.1/32/g;s/NC_052564.1/33/g;s/NC_052565.1/34/g;s/NC_052566.1/35/g;s/NC_052567.1/36/g;s/NC_052568.1/37/g;s/NC_052569.1/38/g;s/NC_052570.1/39/g;s/NC_052571.1/W/g;s/NC_052572.1/Z/g" > /share/home/zhzhang24/PG/refgenome/Gallus_gallus.rmsk.2.bed
rm /share/home/zhzhang24/PG/refgenome/Gallus_gallus.rmsk.1.bed /share/home/zhzhang24/PG/refgenome/Gallus_gallus.rmsk.bed
mv /share/home/zhzhang24/PG/refgenome/Gallus_gallus.rmsk.2.bed /share/home/zhzhang24/PG/refgenome/Gallus_gallus.rmsk.bed
cat "/share/home/zhzhang24/PG/sp.txt"|while read i
do
cat "/share/home/zhzhang24/PG/refgenome/${i}.rmsk.bed"|sed 's/^chr//g'|awk '$1~/^[0-9]/ {print $0} $1~/^[XYZW]/ {print $0}' > /share/home/zhzhang24/PG/refgenome/${i}.rmsk.ok.bed
rm /share/home/zhzhang24/PG/refgenome/${i}.rmsk.bed
done
#产生掩蔽重复序列shuffle所需的bed
cat "/share/home/zhzhang24/PG/sp.txt"|while read i
do
cat /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg_pcgexon_forshuffle.bed > /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg_pcgexon_forshuffle.rm.bed
awk '{print $1"\t"$2"\t"$3}' /share/home/zhzhang24/PG/refgenome/${i}.rmsk.ok.bed >> /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg_pcgexon_forshuffle.rm.bed
done
#获取与重复序列重叠的假基因/lnc的ID
cat "/share/home/zhzhang24/PG/sp.txt"|while read iii
do
#获取与重复序列重叠的假基因
micromamba run -n SEQ bedtools intersect -a /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg.bed -b /share/home/zhzhang24/PG/refgenome/${iii}.rmsk.ok.bed -wa |awk '{print $7}'|sort -u > /share/home/zhzhang24/PG/lncRNA_class/RE_overlap/${iii}_hpg.overlap_re.id.txt
#获取与重复序列重叠的lnc基因
micromamba run -n SEQ bedtools intersect -a "/share/home/zhzhang24/PG/lncRNA_class/${iii}/${iii}.all_lncRNA_exon.merge.bed" -b /share/home/zhzhang24/PG/refgenome/${iii}.rmsk.ok.bed -wa |awk '{print $4}'|sort -u > /share/home/zhzhang24/PG/lncRNA_class/RE_overlap/${iii}_lncRNA.overlap_re.id.txt
done
#掩蔽重复序列(repetitive element)重排假基因
cat "/share/home/zhzhang24/PG/sp.txt"|while read iii
do
grep -w -v -Ff /share/home/zhzhang24/PG/lncRNA_class/RE_overlap/${iii}_hpg.overlap_re.id.txt /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg.bed > /share/home/zhzhang24/PG/lncRNA_class/RE_overlap/${iii}_hpg.nonoverlap_re.bed
chrsize=`ls /share/home/zhzhang24/PG/refgenome/${iii}*dna.chrsize.txt`
for i in $(seq 1 1000)
do
micromamba run -n SEQ bedtools shuffle -i /share/home/zhzhang24/PG/lncRNA_class/RE_overlap/${iii}_hpg.nonoverlap_re.bed -g ${chrsize} -excl /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_pcgexon_forshuffle.rm.bed -noOverlapping |awk -v time="${i}" '{print $0"_"time}' >> /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_shuffle.rm.bed
done
micromamba run -n SEQ bedtools sort -i /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_shuffle.rm.bed > /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_shuffle.rm.sort.bed
done
#传至实验室服务/home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/
cat "/share/home/zhzhang24/PG/sp.txt"|while read iii
do
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${iii}/pgenes/${iii}_hpg_shuffle.rm.sort.bed zhzhang@122.205.95.67:/home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/
done



#使用重排的假基因重新对lncRNA分类，找到与lncRNA重叠的重排假基因，传至实验室服务器
#获取全部lncRNA merge外显子与重排假基因/固定chr重排假基因重叠情况
cat "/share/home/zhzhang24/PG/sp.txt"|while read i
do
micromamba run -n SEQ bedtools intersect -a "/share/home/zhzhang24/PG/lncRNA_class/${i}/${i}.all_lncRNA_exon.merge.bed" -b "/share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg_shuffle.sort.bed" -wo > /share/home/zhzhang24/PG/1Identification/shufflepg_lncRNAclass/${i}.alllncRNAexonmerge_intersect_shufflepg.bed
micromamba run -n SEQ bedtools intersect -a "/share/home/zhzhang24/PG/lncRNA_class/${i}/${i}.all_lncRNA_exon.merge.bed" -b "/share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg_shuffle.keepchr.sort.bed" -wo > /share/home/zhzhang24/PG/1Identification/shufflepg_keepchr_lncRNAclass/${i}.alllncRNAexonmerge_intersect_shufflepg.bed
done
#lncRNA外显子文件去除与RE重叠的lnc基因
cat "/share/home/zhzhang24/PG/sp.txt"|while read iii
do
grep -w -v -Ff /share/home/zhzhang24/PG/lncRNA_class/RE_overlap/${iii}_lncRNA.overlap_re.id.txt "/share/home/zhzhang24/PG/lncRNA_class/${iii}/${iii}.all_lncRNA_exon.merge.bed" > /share/home/zhzhang24/PG/lncRNA_class/RE_overlap/${iii}.all_lncRNA_exon.merge.nonoverlap_re.bed
done
#获取全部lncRNA merge外显子与rm重排假基因重叠情况
cat "/share/home/zhzhang24/PG/sp.txt"|while read i
do
micromamba run -n SEQ bedtools intersect -a /share/home/zhzhang24/PG/lncRNA_class/RE_overlap/${i}.all_lncRNA_exon.merge.nonoverlap_re.bed -b "/share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg_shuffle.rm.sort.bed" -wo > /share/home/zhzhang24/PG/1Identification/shufflepg_rm_lncRNAclass/${i}.alllncRNAexonmerge_intersect_shufflepg.bed
done
#传至实验室服务/home/zhzhang/PG/1Identification/对应/
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/1Identification/shufflepg_lncRNAclass/ zhzhang@122.205.95.67:/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/1Identification/shufflepg_keepchr_lncRNAclass/ zhzhang@122.205.95.67:/home/zhzhang/PG/1Identification/shufflepg_keepchr_lncRNAclass/
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/1Identification/shufflepg_rm_lncRNAclass/ zhzhang@122.205.95.67:/home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/



#确定与lncrna重叠的shuffle假基因
#a输入假基因与lncRNA外显子交集文件
#pgid输出与lncRNA外显子重叠的所有假基因id
lncrnaclass <- function(a,pgid){
  #导入lncRNA外显子与假基因交集文件
  alllncRNAexon_intersect_pg <- read.delim(a, header=FALSE)
  colnames(alllncRNAexon_intersect_pg)[c(4,5,13,14)] <- c("lncgid","exonlen","pgid","overlap_len")
  #输出与lncRNA外显子重叠的所有假基因id
  allpgid <- rbind(mutate(distinct(filter(alllncRNAexon_intersect_pg,V6==V12),pgid,lncgid),type="sense"),
                   mutate(distinct(filter(alllncRNAexon_intersect_pg,V6!=V12),pgid,lncgid),type="antisense"))%>%
    dplyr::rename(geneid=pgid)
  data.table::fwrite(allpgid,file =pgid,sep = '\t',row.names = F,quote = F,col.names = T)
}
splist <- c("Danio_rerio","Gallus_gallus","Monodelphis_domestica","Oryctolagus_cuniculus","Rattus_norvegicus","Mus_musculus","Macaca_mulatta","Homo_sapiens")
#正常shuffle
for (i in 1:8) {
lncrnaclass(paste("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/",splist[i],".alllncRNAexonmerge_intersect_shufflepg.bed",sep=""),
            paste("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/",splist[i],".shufflepg_overlap_lnc.id.txt",sep=""))
}
#keepchr shuffle
for (i in 1:8) {
  lncrnaclass(paste("/home/zhzhang/PG/1Identification/shufflepg_keepchr_lncRNAclass/",splist[i],".alllncRNAexonmerge_intersect_shufflepg.bed",sep=""),
              paste("/home/zhzhang/PG/1Identification/shufflepg_keepchr_lncRNAclass/",splist[i],".shufflepg_overlap_lnc.id.txt",sep=""))
}
#rm shuffle
for (i in 1:8) {
  lncrnaclass(paste("/home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/",splist[i],".alllncRNAexonmerge_intersect_shufflepg.bed",sep=""),
              paste("/home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/",splist[i],".shufflepg_overlap_lnc.id.txt",sep=""))
}



#产生去除与RE重叠的lncid和pgid的二者重叠文件以及hpg文件，用于计算rm shuffle时的FE
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/lncRNA_class/RE_overlap/ zhzhang@122.205.95.67:/home/zhzhang/PG/lncRNA_class/RE_overlap/
cat "/home/zhzhang/PG/sp.txt"|while read iii
do
grep -w -v -Ff /home/zhzhang/PG/lncRNA_class/RE_overlap/${iii}_lncRNA.overlap_re.id.txt /home/zhzhang/PG/lncRNA_class/new_r1/${iii}.pseudogene_overlap_lnc.id.txt|grep -w -v -Ff /home/zhzhang/PG/lncRNA_class/RE_overlap/${iii}_hpg.overlap_re.id.txt > /home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/${iii}.pseudogene_overlap_lnc.id.txt
grep -w -v -Ff /home/zhzhang/PG/lncRNA_class/RE_overlap/${iii}_hpg.overlap_re.id.txt /home/zhzhang/PG/pg_message/${iii}_hpg.txt > /home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/${iii}_hpg.txt
done




###假基因与shuffle区域相比 lncRNA的富集倍数
#palncRNA_pg输入与lnc重叠的假基因，hpg输入假基因信息，
#shufflepalncRNA_pg输入与lnc重叠的shuffle假基因，shufflehpg输入shuffle假基因
#shuffle_type输入shuffle类型,sptype输入英文物种名
#函数返回假基因和shuffle假基因，与lncrna overlap的比例,FE,经验p值
pglendelnc <- function(palncRNA_pg,hpg,shufflepalncRNA_pg,shufflehpg,shuffle_type,sptype){
  #导入假基因
  pg <- read.delim(hpg)%>%
    mutate(type2="All")%>%
    select(geneid,type,type2)
  #统计假基因总数量
  tjall=rbind(group_by(pg,type2)%>%summarise(total_number=n())%>%data.frame()%>%dplyr::rename(type=type2),
              group_by(pg,type)%>%summarise(total_number=n())%>%data.frame())
  #导入与lnc重叠的假基因
  pgdlncRNA_pg <- read.delim(palncRNA_pg)
  colnames(pgdlncRNA_pg)[2] <- "geneid"
  colnames(pgdlncRNA_pg)[3] <- "Btype"
  pgdlncRNA_pg <- left_join(pgdlncRNA_pg,pg,by="geneid")
  #统计overlap的假基因数量
  tjover=rbind(group_by(pgdlncRNA_pg,type2,Btype)%>%summarise(overlap_number=n())%>%data.frame()%>%dplyr::rename(type=type2),
              group_by(pgdlncRNA_pg,type,Btype)%>%summarise(overlap_number=n())%>%data.frame())
  #生成假基因的重叠情况
  pgme <- left_join(tjall,tjover,by="type")%>%
    mutate(nonoverlap_number=total_number-overlap_number,
           ratio=overlap_number*100/total_number)
  #导入shuffle假基因信息
  shufflepg <- data.table::fread(shufflehpg,header = F)%>%data.frame()%>%
    select(V7)%>%
    tidyr::separate(V7,c("one","two","three","times"),sep="_",remove = T)%>%
    mutate(geneid=paste(one,two,three,sep="_"))%>%
    select(geneid,times)%>%
    left_join(pg,by="geneid")
  colnames(shufflepg)[1] <- "geneid"
  #统计shuffle假基因总数量
  tjshuffleall=rbind(group_by(shufflepg,type2,times)%>%summarise(shuffle_total_number=n())%>%data.frame()%>%dplyr::rename(type=type2),
              group_by(shufflepg,type,times)%>%summarise(shuffle_total_number=n())%>%data.frame())
  #导入与lnc重叠的shuffle假基因
  shufflepgdlncRNA_pg <- data.table::fread(shufflepalncRNA_pg)%>%data.frame()
  colnames(shufflepgdlncRNA_pg)[2] <- "geneid"
  colnames(shufflepgdlncRNA_pg)[3] <- "Btype"
  shufflepgdlncRNA_pg <- tidyr::separate(shufflepgdlncRNA_pg,geneid,c("one","two","three","times"),sep="_",remove = T)%>%
    mutate(geneid=paste(one,two,three,sep="_"))%>%
    select(geneid,times,Btype)%>%
    left_join(pg,by="geneid")
  #统计overlap的shuffle假基因数量
  pre_tjshuffleover=data.frame(type=rep(c("All","Duplicated","Processed","Fragment"),times=1,each=2000),
                               Btype=rep(c("antisense","sense"),times=4,each=1000),
                               times=as.character(rep(c(1:1000),times=8,each=1)))%>%
    mutate(left=paste(type,Btype,times,sep="_"))
  pre2_tjshuffleover=rbind(group_by(shufflepgdlncRNA_pg,type2,Btype,times)%>%summarise(shuffle_overlap_number=n())%>%data.frame()%>%dplyr::rename(type=type2),
               group_by(shufflepgdlncRNA_pg,type,Btype,times)%>%summarise(shuffle_overlap_number=n())%>%data.frame())%>%
    mutate(left=paste(type,Btype,times,sep="_"))%>%
    select(-type,-Btype,-times)
  tjshuffleover=left_join(pre_tjshuffleover,
                          pre2_tjshuffleover,
                          by="left")%>%
    select(-left)
  tjshuffleover$shuffle_overlap_number[is.na(tjshuffleover$shuffle_overlap_number)==T] <- 0
  #pre
  ob1=select(filter(pgme,type=="All" & Btype=="antisense"),ratio)$ratio
  ob2=select(filter(pgme,type=="All" & Btype=="sense"),ratio)$ratio
  ob3=select(filter(pgme,type=="Duplicated" & Btype=="antisense"),ratio)$ratio
  ob4=select(filter(pgme,type=="Duplicated" & Btype=="sense"),ratio)$ratio
  ob5=select(filter(pgme,type=="Processed" & Btype=="antisense"),ratio)$ratio
  ob6=select(filter(pgme,type=="Processed" & Btype=="sense"),ratio)$ratio
  ob7=select(filter(pgme,type=="Fragment" & Btype=="antisense"),ratio)$ratio
  ob8=select(filter(pgme,type=="Fragment" & Btype=="sense"),ratio)$ratio
  #生成shuffle假基因的重叠情况
  shuffle_pgme <- left_join(mutate(tjshuffleall,left=paste(type,times,sep="_"))%>%select(-type,-times),
                    mutate(tjshuffleover,left=paste(type,times,sep="_"))%>%select(-type,-times),
                    by="left")%>%
    mutate(shuffle_nonoverlap_number=shuffle_total_number-shuffle_overlap_number,
           shuffle_ratio=shuffle_overlap_number*100/shuffle_total_number)%>%
    tidyr::separate(left,c("type","times"),sep="_",remove = T)%>%
    mutate(greater=case_when(type=="All" & Btype=="antisense" & shuffle_ratio >= ob1 ~ 1,
                             type=="All" & Btype=="sense" & shuffle_ratio >= ob2 ~ 1,
                             type=="Duplicated" & Btype=="antisense" & shuffle_ratio >= ob3 ~ 1,
                             type=="Duplicated" & Btype=="sense" & shuffle_ratio >= ob4 ~ 1,
                             type=="Processed" & Btype=="antisense" & shuffle_ratio >= ob5 ~ 1,
                             type=="Processed" & Btype=="sense" & shuffle_ratio >= ob6 ~ 1,
                             type=="Fragment" & Btype=="antisense" & shuffle_ratio >= ob7 ~ 1,
                             type=="Fragment" & Btype=="sense" & shuffle_ratio >= ob8 ~ 1,
                             T~0))
  #生成shuffle假基因的重叠情况(平均)
  mean_forboot <- function(data, index) {
    return(mean(data[index]))
  }
  set.seed(1024)
  shuffle_pgme_tj <- group_by(shuffle_pgme,type,Btype)%>%
    summarise(shuffle_total_number=sum(shuffle_total_number),
              shuffle_overlap_number=sum(shuffle_overlap_number),
              shuffle_nonoverlap_number=sum(shuffle_nonoverlap_number),
              shuffle_ratio_mean=mean(shuffle_ratio),
              shuffle_ratio_confmin=boot::boot.ci(boot::boot(shuffle_ratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][4],
              shuffle_ratio_confmax=boot::boot.ci(boot::boot(shuffle_ratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][5],
              forpvalue=sum(greater))%>%
    data.frame()%>%
    mutate(prepvalue=case_when(forpvalue>0 ~ forpvalue,
                               forpvalue==0 ~ 1))%>%
    mutate(pvalue=prepvalue/1000,
           left=paste(type,Btype,sep="_"))%>%
    select(-type,-Btype,-forpvalue,-prepvalue)
  #leftjoin pg和shufflepg的情况
  output=left_join(mutate(pgme,left=paste(type,Btype,sep="_")),
                   shuffle_pgme_tj,by="left")%>%
    select(1,3,2,4:6,8:14)%>%
    mutate(FE=ratio/shuffle_ratio_mean)%>%
    mutate(shuffle_type=shuffle_type,species=sptype)%>%
    mutate(Btype=stringr::str_to_title(Btype))%>%
    select(15,16,1:14)
}
#物种list
splist1 <- c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")
splist2 <- c("Danio_rerio","Gallus_gallus","Monodelphis_domestica","Oryctolagus_cuniculus","Rattus_norvegicus","Mus_musculus","Macaca_mulatta","Homo_sapiens")
all <- data.frame()
for (i in 1:8) {
  forbind1 <- pglendelnc(palncRNA_pg=paste("/home/zhzhang/PG/lncRNA_class/new_r1/",splist2[i],".pseudogene_overlap_lnc.id.txt",sep=""),
                        hpg=paste("/home/zhzhang/PG/pg_message/",splist2[i],"_hpg.txt",sep=""),
                        shufflepalncRNA_pg=paste("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/",splist2[i],".shufflepg_overlap_lnc.id.txt",sep=""),
                        shufflehpg=paste("/home/zhzhang/PG/1Identification/shufflepg_lncRNAclass/",splist2[i],"_hpg_shuffle.sort.bed",sep=""),
                        shuffle_type="Normal",
                        sptype=splist1[i])
  forbind2 <- pglendelnc(palncRNA_pg=paste("/home/zhzhang/PG/lncRNA_class/new_r1/",splist2[i],".pseudogene_overlap_lnc.id.txt",sep=""),
                         hpg=paste("/home/zhzhang/PG/pg_message/",splist2[i],"_hpg.txt",sep=""),
                         shufflepalncRNA_pg=paste("/home/zhzhang/PG/1Identification/shufflepg_keepchr_lncRNAclass/",splist2[i],".shufflepg_overlap_lnc.id.txt",sep=""),
                         shufflehpg=paste("/home/zhzhang/PG/1Identification/shufflepg_keepchr_lncRNAclass/",splist2[i],"_hpg_shuffle.keepchr.sort.bed",sep=""),
                         shuffle_type="Keepchr",
                         sptype=splist1[i])
  forbind3 <- pglendelnc(palncRNA_pg=paste("/home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/",splist2[i],".pseudogene_overlap_lnc.id.txt",sep=""),
                         hpg=paste("/home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/",splist2[i],"_hpg.txt",sep=""),
                         shufflepalncRNA_pg=paste("/home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/",splist2[i],".shufflepg_overlap_lnc.id.txt",sep=""),
                         shufflehpg=paste("/home/zhzhang/PG/1Identification/shufflepg_rm_lncRNAclass/",splist2[i],"_hpg_shuffle.rm.sort.bed",sep=""),
                         shuffle_type="Repeatmasked",
                         sptype=splist1[i])
  all=rbind(all,forbind1,forbind2,forbind3)
}
#储存
data.table::fwrite(all,
                   file ="/home/zhzhang/PG/1Identification/SP.3shuffle.pg_shufflepg_overlaplncratio.tj.tsv",
                   sep = '\t',row.names = F,quote = F,col.names = T)



```
```r
#FE热图
#input
all <- read.delim("~/PG/1Identification/SP.3shuffle.pg_shufflepg_overlaplncratio.tj.tsv")%>%
  mutate(panno=case_when(pvalue<=0.001 ~ "***",
                         pvalue<=0.01 ~ "**",
                         pvalue<=0.05 ~ "*",
                         T~" "))
all$shuffle_type <- factor(all$shuffle_type,levels = c("Normal","Keepchr","Repeatmasked"))
#ALL pg
#提取FE
all_FE <- filter(all,type=="All")%>%
  select(1:4,16)%>%
  mutate(colname=paste(Btype,"-",species,sep=""))%>%
  select(-species,-Btype)%>%
  spread(colname,FE)%>%
  column_to_rownames("shuffle_type")%>%
  select(-type)%>%   ####
  select(paste(rep(c("Sense","Antisense"),each=8),
               rep(rev(c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")),times=2),sep="-"))
#提取pvalue
all_p <- filter(all,type=="All")%>%
  select(1:4,17)%>%
  mutate(colname=paste(Btype,"-",species,sep=""))%>%
  select(-species,-Btype)%>%
  spread(colname,panno)%>%
  column_to_rownames("shuffle_type")%>%
  select(-type)%>%   ####
  select(paste(rep(c("Sense","Antisense"),each=8),
               rep(rev(c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")),times=2),sep="-"))
#plot
bk <- c(seq(0,1,by=0.01),seq(1.1,6,by=0.01))
colnames(all_FE)=rep(rev(c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")),times=2)
ph <- pheatmap::pheatmap(all_FE,
                         scale="none",cluster_cols = F, cluster_rows=F,
                         color = c(colorRampPalette(colors = c("#6082B0","white"))(length(bk)*1/6),colorRampPalette(colors = c("white","#CE6463"))(length(bk)*5/6)),
                         breaks=bk,
                         border_color=NA,angle_col=45,
                         fontsize=15,fontsize_col=12,fontsize_row=12,
                         display_numbers=all_p,
                         gaps_col = c(8),gaps_row = 1:3)
ggsave("/home/zhzhang/PG/1Identification/SP.3shuffle.pg_shufflepg_overlaplncratio.FE.pdf", 
       ph,width = 11, height = 3.5,dpi=1200, units = "in", device='pdf',bg = "transparent")

#3type pg
#提取FE
all_FE <- filter(all,type!="All")%>%
  select(1:4,16)%>%
  mutate(colname=paste(Btype,"-",species,sep=""))%>%
  select(-species,-Btype)%>%
  spread(colname,FE)%>%
  mutate(rowname=paste(shuffle_type,"-",type,sep=""))%>%
  column_to_rownames("rowname")%>%
  select(paste(rep(c("Sense","Antisense"),each=8),
               rep(rev(c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")),times=2),sep="-"))
all_FE=all_FE[paste(rep(c("Normal","Keepchr","Repeatmasked"),each=3),
                    rep(c("Duplicated","Processed","Fragment"),times=3),sep="-"),]
all_FE=log2(all_FE)
#提取pvalue
all_p <- filter(all,type!="All")%>%
  select(1:4,17)%>%
  mutate(colname=paste(Btype,"-",species,sep=""))%>%
  select(-species,-Btype)%>%
  spread(colname,panno)%>%
  mutate(rowname=paste(shuffle_type,"-",type,sep=""))%>%
  column_to_rownames("rowname")%>%
  select(paste(rep(c("Sense","Antisense"),each=8),
               rep(rev(c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")),times=2),sep="-"))
all_p=all_p[paste(rep(c("Normal","Keepchr","Repeatmasked"),each=3),
                    rep(c("Duplicated","Processed","Fragment"),times=3),sep="-"),]
#plot
bk <- c(seq(-1.5,0,by=0.01),seq(0.1,4.5,by=0.01))
colnames(all_FE)=rep(rev(c("Zebrafish","Chicken","Opossum","Rabbit","Rat","Mouse","Macaque","Human")),times=2)
ph <- pheatmap::pheatmap(all_FE,
                         scale="none",cluster_cols = F, cluster_rows=F,
                         color = c(colorRampPalette(colors = c("#6082B0","white"))(length(bk)*1.5/6),colorRampPalette(colors = c("white","#CE6463"))(length(bk)*4.5/6)),
                         breaks=bk,
                         border_color=NA,angle_col=45,
                         fontsize=15,fontsize_col=12,fontsize_row=1,
                         display_numbers=all_p,
                         gaps_col = c(8),gaps_row = c(3,6,9))
ggsave("/home/zhzhang/PG/1Identification/SP.3shuffle.3typepg_shufflepg_overlaplncratio.FE.pdf", 
       ph,width = 11, height = 7,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


##### 8.对比与lncRNA基因重叠的假基因中，与lnc重叠区域和非重叠区域保守性
```r
#产生lncRNA基因的假基因中，与lnc重叠区域和非重叠区域的bed
cat "/share/home/zhzhang24/PG/sp.txt" |while read i
do
#获取与sense lncRNA overlap的假基因bed
awk '$2!="pgid" && $3=="sense"{print $2}' "/share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.pseudogene_overlap_lnc.id.txt"|sort -u > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpgid.sense.txt
grep -w -Ff /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpgid.sense.txt "/share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg.bed"|awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6}' > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.bed
#获取与antisense lncRNA overlap的假基因bed
awk '$2!="pgid" && $3=="antisense"{print $2}' "/share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.pseudogene_overlap_lnc.id.txt"|sort -u > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpgid.antisense.txt
grep -w -Ff /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpgid.antisense.txt "/share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg.bed"|awk '{print $1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6}' > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.bed
done
#获取重叠区域bed
cat "/share/home/zhzhang24/PG/sp.txt" |while read i
do
micromamba run -n SEQ bedtools intersect -a /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.bed -b "/share/home/zhzhang24/PG/lncRNA_class/${i}/${i}.all_lncRNA_exon.merge.bed" -s|awk '{print $1"\t"$2"\t"$3"\t"$4"____"NR}' > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.overlap_region.bed
micromamba run -n SEQ bedtools intersect -a /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.bed -b "/share/home/zhzhang24/PG/lncRNA_class/${i}/${i}.all_lncRNA_exon.merge.bed" -S|awk '{print $1"\t"$2"\t"$3"\t"$4"____"NR}' > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.overlap_region.bed
done
#非重叠区bed
cat "/share/home/zhzhang24/PG/sp.txt" |while read i
do
micromamba run -n SEQ bedtools subtract -a /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.bed -b /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.overlap_region.bed|awk '{print $1"\t"$2"\t"$3"\t"$4"____"NR}' > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.nonoverlap_region.bed
micromamba run -n SEQ bedtools subtract -a /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.bed -b /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.overlap_region.bed|awk '{print $1"\t"$2"\t"$3"\t"$4"____"NR}' > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.nonoverlap_region.bed
done




#phastcons
echo -e "Homo_sapiens\nMus_musculus"|while read i
do
if [ "$i" == "Homo_sapiens" ]
then
bw="/share/home/zhzhang24/PG/Evolution/conserve/hg38.phastCons100way.bw"
fi
if [ "$i" == "Mus_musculus" ]
then
bw="/share/home/zhzhang24/PG/Evolution/conserve/mm39.phastCons35way.bw"
fi
#
awk '{print "chr"$0}' /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.overlap_region.bed > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.overlap_region.ucsc.bed
awk '{print "chr"$0}' /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.nonoverlap_region.bed > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.nonoverlap_region.ucsc.bed
awk '{print "chr"$0}' /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.overlap_region.bed > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.overlap_region.ucsc.bed
awk '{print "chr"$0}' /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.nonoverlap_region.bed > /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.nonoverlap_region.ucsc.bed
#计算phastCons
micromamba run -n SEQ bigWigAverageOverBed ${bw} /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.overlap_region.ucsc.bed /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.overlap_region.phastCons.txt
micromamba run -n SEQ bigWigAverageOverBed ${bw} /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.nonoverlap_region.ucsc.bed /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.sense.nonoverlap_region.phastCons.txt
micromamba run -n SEQ bigWigAverageOverBed ${bw} /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.overlap_region.ucsc.bed /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.overlap_region.phastCons.txt
micromamba run -n SEQ bigWigAverageOverBed ${bw} /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.nonoverlap_region.ucsc.bed /share/home/zhzhang24/PG/Evolution/hpg_conserve/${i}.hpg.antisense.nonoverlap_region.phastCons.txt
done


rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/Evolution/hpg_conserve/ zhzhang@122.205.95.67:/home/zhzhang/PG/Evolution/hpg_conserve/

```
```r
SP <- c("Homo_sapiens","Mus_musculus")
sp <- c("Human","Mouse")
#
zdf <- data.frame()
for (i in 1:length(sp)) {
  #path
  a <- paste("/home/zhzhang/PG/Evolution/hpg_conserve/",SP[i],".hpg.sense.overlap_region.phastCons.txt",sep="")
  b <- paste("/home/zhzhang/PG/Evolution/hpg_conserve/",SP[i],".hpg.sense.nonoverlap_region.phastCons.txt",sep="")
  c <- paste("/home/zhzhang/PG/Evolution/hpg_conserve/",SP[i],".hpg.antisense.overlap_region.phastCons.txt",sep="")
  d <- paste("/home/zhzhang/PG/Evolution/hpg_conserve/",SP[i],".hpg.antisense.nonoverlap_region.phastCons.txt",sep="")
  #导入与sense lnc重叠的假基因区域的保守性打分
  sense_overlapregion <- read.delim(a, header=FALSE)%>%
    filter(V3!=0)%>%
    select(1,3,4)%>%
    separate(V1,c("pgid"),sep="____")%>%
    group_by(pgid)%>%
    summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
    mutate(score=sumscore/sumlen)%>%
    select(1,4)%>%
    mutate(type="Overlapping",lnctype="Pseudogene overlapping with sense lncRNA gene")%>%
    mutate(sp=sp[i])
  #
  sense_nonoverlapregion <- read.delim(b, header=FALSE)%>%
    filter(V3!=0)%>%
    select(1,3,4)%>%
    separate(V1,c("pgid"),sep="____")%>%
    group_by(pgid)%>%
    summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
    mutate(score=sumscore/sumlen)%>%
    select(1,4)%>%
    mutate(type="Non-overlapping",lnctype="Pseudogene overlapping with sense lncRNA gene")%>%
    mutate(sp=sp[i])
  #导入与antisense lnc重叠的假基因区域的保守性打分
  antisense_overlapregion <- read.delim(c, header=FALSE)%>%
    filter(V3!=0)%>%
    select(1,3,4)%>%
    separate(V1,c("pgid"),sep="____")%>%
    group_by(pgid)%>%
    summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
    mutate(score=sumscore/sumlen)%>%
    select(1,4)%>%
    mutate(type="Overlapping",lnctype="Pseudogene overlapping with antisense lncRNA gene")%>%
    mutate(sp=sp[i])
  #
  antisense_nonoverlapregion <- read.delim(d, header=FALSE)%>%
    filter(V3!=0)%>%
    select(1,3,4)%>%
    separate(V1,c("pgid"),sep="____")%>%
    group_by(pgid)%>%
    summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
    mutate(score=sumscore/sumlen)%>%
    select(1,4)%>%
    mutate(type="Non-overlapping",lnctype="Pseudogene overlapping with antisense lncRNA gene")%>%
    mutate(sp=sp[i])
  #
  zdf <- rbind(zdf,sense_overlapregion,sense_nonoverlapregion,antisense_overlapregion,antisense_nonoverlapregion)
}
#
zdf$type <- factor(zdf$type,levels = c("Overlapping","Non-overlapping"))
zdf$lnctype <- factor(zdf$lnctype,levels = c("Pseudogene overlapping with sense lncRNA gene","Pseudogene overlapping with antisense lncRNA gene"))
zdf$sp <- factor(zdf$sp,levels = sp)
#统计
tj <- group_by(zdf,sp,lnctype,type)%>%
  summarise(mean=mean(score),median=median(score))
data.table::fwrite(tj,file ="/home/zhzhang/PG/Evolution/hpg_conserve/2SP.hpg_overlapSaASlnc.overnonover.phast.tj.tsv",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
#Pseudogene overlapping with sense lncRNA gene
p1 <- ggplot(data=filter(zdf,lnctype=="Pseudogene overlapping with sense lncRNA gene"), aes(x=type,y=score))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  ggsignif::geom_signif(test = "t.test",map_signif_level=T,comparisons=list(c("Overlapping","Non-overlapping")))+
  scale_fill_manual(values = c("#99CCFF","#FFCC66"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  labs(y = "PhastCons score", x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 0),legend.position = "right") +
  theme(legend.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  facet_wrap(~sp,nrow=1)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15)) +
  theme(axis.ticks.x = element_line(size =0))
ggsave("/home/zhzhang/PG/Evolution/hpg_conserve/2SP.hpg_overlapSlnc.overnonover.phast.pdf",
       p1,width = 6, height = 4,dpi=1200, units = "in", device='pdf',bg = "transparent")
#Pseudogene overlapping with antisense lncRNA gene
p1 <- ggplot(data=filter(zdf,lnctype=="Pseudogene overlapping with antisense lncRNA gene"), aes(x=type,y=score))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  ggsignif::geom_signif(test = "t.test",map_signif_level=T,comparisons=list(c("Overlapping","Non-overlapping")))+
  scale_fill_manual(values = c("#99CCFF","#FFCC66"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  labs(y = "PhastCons score", x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 0),legend.position = "right") +
  theme(legend.text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  facet_wrap(~sp,nrow=1)+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15)) +
  theme(axis.ticks.x = element_line(size =0))
ggsave("/home/zhzhang/PG/Evolution/hpg_conserve/2SP.hpg_overlapASlnc.overnonover.phast.pdf",
       p1,width = 6, height = 4,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


