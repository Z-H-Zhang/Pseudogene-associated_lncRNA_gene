### 六.RBP结合相关分析以及功能推断
##### 1.m6A修饰相关分析
```r
#REPIC数据库下载 m6A peak位点（人类hg38，小鼠mm10）
wget https://repicmod.uchicago.edu/repic/data/download/m6A=sites=species=human=hg38.txt.gz --no-check-certificate
wget https://repicmod.uchicago.edu/repic/data/download/m6A=sites=species=mouse=mm10.txt.gz --no-check-certificate
gzip -d "/home/zhzhang/PG/m6a/m6A=sites=species=human=hg38.txt.gz"
gzip -d "/home/zhzhang/PG/m6a/m6A=sites=species=mouse=mm10.txt.gz"



#提取简约信息（POS）,转化为非冗余的m6A位点bed（所有潜在的m6A位点）
#人
tail -n +2 "/home/zhzhang/PG/m6a/m6A=sites=species=human=hg38.txt"|awk '{print $1}' > /home/zhzhang/PG/m6a/m6A=sites=human=hg38.txt
grep -w "+" /home/zhzhang/PG/m6a/m6A=sites=human=hg38.txt|sed 's/\[+\]//g'|awk -F ":" '{print $1"\t"$2}'|awk -F "-" '{print $1"\t"$2}'|sed 's/^chr//g'|bedtools sort|bedtools merge|awk '{print $0"\t+"}' > /home/zhzhang/PG/m6a/Homo_sapiens.m6Asites.bed
grep -w "-" /home/zhzhang/PG/m6a/m6A=sites=human=hg38.txt|sed 's/\[-\]//g'|awk -F ":" '{print $1"\t"$2}'|awk -F "-" '{print $1"\t"$2}'|sed 's/^chr//g'|bedtools sort|bedtools merge|awk '{print $0"\t-"}' >> /home/zhzhang/PG/m6a/Homo_sapiens.m6Asites.bed
cat /home/zhzhang/PG/m6a/Homo_sapiens.m6Asites.bed|bedtools sort|awk '{print $1"\t"$2"\t"$3"\tm6Asite_"NR"\tm6Asite\t"$4}' > /home/zhzhang/PG/m6a/Homo_sapiens.m6Asites.sort.bed
rm /home/zhzhang/PG/m6a/Homo_sapiens.m6Asites.bed
#小鼠(外加liftover到grch39)
tail -n +2 "/home/zhzhang/PG/m6a/m6A=sites=species=mouse=mm10.txt"|awk '{print $1}' > /home/zhzhang/PG/m6a/m6A=sites=mouse=mm10.txt
grep -w "+" "/home/zhzhang/PG/m6a/m6A=sites=mouse=mm10.txt"|sed 's/\[+\]//g'|awk -F ":" '{print $1"\t"$2}'|awk -F "-" '{print $1"\t"$2}'|bedtools sort|bedtools merge|awk '{print $0"\t+"}' > /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.bed
grep -w "-" "/home/zhzhang/PG/m6a/m6A=sites=mouse=mm10.txt"|sed 's/\[-\]//g'|awk -F ":" '{print $1"\t"$2}'|awk -F "-" '{print $1"\t"$2}'|bedtools sort|bedtools merge|awk '{print $0"\t-"}' >> /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.bed
cat /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.bed|bedtools sort|awk '{print $1"\t"$2"\t"$3"\tm6Asite_"NR"\tm6Asite\t"$4}' > /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.sort.mm10.bed
rm /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.bed
liftOver /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.sort.mm10.bed /home/zhzhang/software/liftover/mm10ToMm39.over.chain.gz /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.sort.mm39.bed /home/zhzhang/PG/m6a/Mus_musculus.unmap
rm /home/zhzhang/PG/m6a/Mus_musculus.unmap
cat /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.sort.mm39.bed |sed 's/^chr//g'|bedtools sort > /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.sort.bed
rm /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.sort.mm*.bed



```
```r
#非冗余的m6a位点与lncRNA基因外显子交集，获取每个基因外显子具有的m6a位点数量
#人
bedtools intersect -a /home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon.bed -b /home/zhzhang/PG/m6a/Homo_sapiens.m6Asites.sort.bed -c -s > /home/zhzhang/PG/m6a/Homo_sapiens.allgeneexon_intersect_m6A.txt
#小鼠
bedtools intersect -a /home/zhzhang/PG/RBP/Mus_musculus.allgeneexon.bed -b /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.sort.bed -c -s > /home/zhzhang/PG/m6a/Mus_musculus.allgeneexon_intersect_m6A.txt



#m6a位点与基因间区交集
bedtools intersect -a /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed -b /home/zhzhang/PG/m6a/Homo_sapiens.m6Asites.sort.bed -c -s > /home/zhzhang/PG/m6a/Homo_sapiens.intergenic_intersect_m6A.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/conserve/Mus_musculus.20000random_3kb_intergenic.bed -b /home/zhzhang/PG/m6a/Mus_musculus.m6Asites.sort.bed -c -s > /home/zhzhang/PG/m6a/Mus_musculus.intergenic_intersect_m6A.txt


```
##### 对比蛋白编码/lncRNA基因转录本的m6a位点密度
```r
#函数统计每个基因全部外显子&基因间区具有的m6a位点密度
#a输入每个基因全部外显子与m6a位点交集文件，b输入每个基因间区与m6a位点交集文件
#c输入基因分类文件，e输入全部基因外显子bed文件
geneexonM6Ad <- function(a,b,c,e){
  #导入全部外显子与m6a位点交集,获取每个基因全部外显子(全部)中每种m6a位点数量
  allgeneexon_intersect_m6a <- data.table::fread(a, header=FALSE)%>%
    data.frame()%>%
    select(4,7)%>%
    group_by(V4)%>%
    summarise(num=sum(V7))
  colnames(allgeneexon_intersect_m6a) <- c("geneid","num")
  #计算基因全部外显子（全部转录本）长度
  allgeneexon <- read.delim(e, header=FALSE)%>%
    mutate(len=(V3-V2+1)/1000)%>%
    group_by(V4)%>%
    summarise(len=sum(len))
  colnames(allgeneexon) <- c("geneid","len")
  #合并长度和m6a位点数量信息，计算密度（m6a位点数量/基因全部外显子长度）
  allgeneexon_m6a_len <- left_join(allgeneexon_intersect_m6a,allgeneexon,by="geneid")%>%
    mutate(m6adensity=num/len)
  #导入基因分类
  geneid_class <- read.delim(c)
  allgeneexon_m6a_len <- left_join(geneid_class,allgeneexon_m6a_len,by="geneid")%>%
    filter(type!="Interference lncRNA")
  #导入基因间区与RBP交集,获取每个基因间区每种RBP的结合长度比例
  intergenic_intersect_m6a <- data.table::fread(b, header=FALSE)%>%
    data.frame()%>%
    select(4,7)%>%
    group_by(V4)%>%
    summarise(num=sum(V7))%>%
    mutate(len=3,type="Random intergenic")%>%
    mutate(m6adensity=num/len)%>%
    select(1,4,2,3,5)
  colnames(intergenic_intersect_m6a) <- colnames(allgeneexon_m6a_len)
  #合并
  allm6ainfo <- rbind(allgeneexon_m6a_len,intergenic_intersect_m6a)
}

#人
hsm6ad <- geneexonM6Ad(a="/home/zhzhang/PG/m6a/Homo_sapiens.allgeneexon_intersect_m6A.txt",
                       b="/home/zhzhang/PG/m6a/Homo_sapiens.intergenic_intersect_m6A.txt",
                       c="/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                       e="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon.bed")
#合并所有物种结果
ALLSP_m6a <- rbind(mutate(hsm6ad,sp="Human")%>%filter(type!="Random intergenic"))%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
ALLSP_m6a$type <- factor(ALLSP_m6a$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
data.table::fwrite(ALLSP_m6a,file ="/home/zhzhang/PG/m6a/SPhs_3gene_m6Adensity.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#统计
tj <- group_by(ALLSP_m6a,sp,type)%>%
  summarise(mean=mean(m6adensity),median=median(m6adensity))
data.table::fwrite(tj,file ="/home/zhzhang/PG/m6a/SPhs_3gene_m6Adensity.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
phs <- ggplot(data=ALLSP_m6a, aes(x=type,y=m6adensity))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position=c(6.7,3.2,4.2,1),tip_length = 0.005,
              comparisons = list(c("Protein-coding","PAS lncRNA"),
                                 c("PAS lncRNA","PAA lncRNA"),
                                 c("PAS lncRNA","NPA lncRNA"),
                                 c("PAA lncRNA","NPA lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                     limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0, 8))+
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA"))+
  labs(y = "Density of m6A sites",x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave("/home/zhzhang/PG/m6a/Homo_sapiens.3gene_m6Adensity.pdf",
       phs,width = 3, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
##### 严格阈值对m6a位点密度的影响
```r
#函数统计每个基因全部外显子&基因间区具有的m6a位点密度
#a输入每个基因全部外显子与m6a位点交集文件，b输入每个基因间区与m6a位点交集文件
#c输入基因分类文件，e输入全部基因外显子bed文件
geneexonM6Ad <- function(a,c,e){
  #导入全部外显子与m6a位点交集,获取每个基因全部外显子(全部)中每种m6a位点数量
  allgeneexon_intersect_m6a <- data.table::fread(a, header=FALSE)%>%
    data.frame()%>%
    select(4,7)%>%
    group_by(V4)%>%
    summarise(num=sum(V7))
  colnames(allgeneexon_intersect_m6a) <- c("geneid","num")
  #计算基因全部外显子（全部转录本）长度
  allgeneexon <- read.delim(e, header=FALSE)%>%
    mutate(len=(V3-V2+1)/1000)%>%
    group_by(V4)%>%
    summarise(len=sum(len))
  colnames(allgeneexon) <- c("geneid","len")
  #合并长度和m6a位点数量信息，计算密度（m6a位点数量/基因全部外显子长度）
  allgeneexon_m6a_len <- left_join(allgeneexon_intersect_m6a,allgeneexon,by="geneid")%>%
    mutate(m6adensity=num/len)
  #导入基因分类
  geneid_class <- read.delim(c)%>%select(2,4,7)%>%
    dplyr::rename(geneid=lncgid)
  allgeneexon_m6a_len <- left_join(geneid_class,allgeneexon_m6a_len,by="geneid")
}

#人
hsm6ad <- geneexonM6Ad(a="/home/zhzhang/PG/m6a/Homo_sapiens.allgeneexon_intersect_m6A.txt",
                       c="/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt",
                       e="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon.bed")
#
pas <- rbind(mutate(filter(hsm6ad,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(hsm6ad,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(hsm6ad,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(hsm6ad,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(hsm6ad,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>200),cuttype="200"))
paa <- rbind(mutate(filter(hsm6ad,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(hsm6ad,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(hsm6ad,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(hsm6ad,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(hsm6ad,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>200),cuttype="200"))
#合并
alldata <- rbind(pas,paa)%>%
  mutate(type=case_when(type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA"))
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","50","100","150","200"))
#统计储存
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(m6adensity),median=median(m6adensity))
#PAS
p1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=m6adensity))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),test.args = c("greater"),
                        y_position=c(3.9),tip_length = 0.01,size=0.5,textsize=4)+
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0, 4.4))+
  labs(y ="Density of m6A sites",
       x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/m6a/Homo_sapiens.paslnc_m6Adensity.pdf",
       p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=m6adensity))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=F,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),test.args = c("greater"),
                        y_position=c(2),tip_length = 0.005,size=0.5,textsize=4)+
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0, 2.5))+
  labs(y ="Density of m6A sites",
       x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/m6a/Homo_sapiens.paalnc_m6Adensity.pdf",
       p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
##### 对比不同表达水平区间内，三类基因m6a密度
```r
##画图对比不同表达水平区间内，三类基因m6a密度
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
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.spexp_geneid_class.txt",
                      c = "Human")
#不同物种三类基因表达量数据合并
ALLsp_TPM <- hs_TPM
#每个基因全部外显子&基因间区具有的m6a位点密度
ALLSP_m6a <- read.delim("~/PG/m6a/SPhsmm_3gene_m6Adensity.txt")%>%
  filter(type!="Random intergenic" & sp!="Mouse")%>%
  select(1,5)
#合并
he <- left_join(ALLsp_TPM,ALLSP_m6a,by="geneid")%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
he$type <- factor(he$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
#
tj4exp <- group_by(he,sp)%>%
  summarise(fen25=quantile(maxTPM,0.25),fen50=quantile(maxTPM,0.50),
            fen75=quantile(maxTPM,0.75))
#增加表达区间
he <- mutate(he,expinter=case_when(sp=="Human" & maxTPM < tj4exp$fen25[1] ~ "[Min,25%)",
                                   sp=="Human" & maxTPM < tj4exp$fen50[1] ~ "[25%,50%)",
                                   sp=="Human" & maxTPM < tj4exp$fen75[1] ~ "[50%,75%)",
                                   sp=="Human" & maxTPM >= tj4exp$fen75[1] ~ "[75%,Max]"
))
#统计
he$expinter <- factor(he$expinter,levels = c("[Min,25%)","[25%,50%)","[50%,75%)","[75%,Max]"))
tj <- group_by(he,sp,type,expinter)%>%
  summarise(mean=mean(m6adensity),median=median(m6adensity))
data.table::fwrite(tj,file ="/home/zhzhang/PG/m6a/SPhs_3gene_m6Adensity.expinterval.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#每个表达区间内，m6adensity差异性检验
interval <- c("[Min,25%)","[25%,50%)","[50%,75%)","[75%,Max]")
pvaluedata <- data.frame()
for (i in 1:4) {
  interfor <- interval[i]
  pvalue1 <- wilcox.test(filter(he,expinter==interfor & type=="NPA lncRNA")$m6adensity,
                        filter(he,expinter==interfor & type=="PAS lncRNA")$m6adensity)[["p.value"]]
  pvalue2 <- wilcox.test(filter(he,expinter==interfor & type=="NPA lncRNA")$m6adensity,
                         filter(he,expinter==interfor & type=="PAA lncRNA")$m6adensity)[["p.value"]]
  pvalue3 <- wilcox.test(filter(he,expinter==interfor & type=="PAS lncRNA")$m6adensity,
                         filter(he,expinter==interfor & type=="PAA lncRNA")$m6adensity)[["p.value"]]
  pvaluedata[i,1] <- interfor
  pvaluedata[i,2] <- pvalue1
  pvaluedata[i,3] <- pvalue2
  pvaluedata[i,4] <- pvalue3
}
colnames(pvaluedata) <- c("expinter","pvalue1","pvalue2","pvalue3")
pvaluedata <- mutate(pvaluedata,anno1=case_when(pvalue1<0.001 ~ "***",
                                               pvalue1<0.01 ~ "**",
                                               pvalue1<0.05 ~ "*",
                                               pvalue1>=0.05 ~ "N.S."),
                     anno2=case_when(pvalue2<0.001 ~ "***",
                                     pvalue2<0.01 ~ "**",
                                     pvalue2<0.05 ~ "*",
                                     pvalue2>=0.05 ~ "N.S."),
                     anno3=case_when(pvalue3<0.001 ~ "***",
                                     pvalue3<0.01 ~ "**",
                                     pvalue3<0.05 ~ "*",
                                     pvalue3>=0.05 ~ "N.S."))
data.table::fwrite(pvaluedata,file ="/home/zhzhang/PG/m6a/SPhs_3gene_m6Adensity.expinterval.difpvalue.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
phs <- ggplot(data=he, aes(x=expinter,y=m6adensity))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.6,notch=F,aes(fill=type))+
  ggsignif::geom_signif(annotations =pvaluedata$anno1 ,y_position=c(2.5,3.5,5.5,6.5)+1,tip_length = 0.005,
                        xmin = c(0.925:3.925),
                        xmax = c(1.225:4.225))+
  ggsignif::geom_signif(annotations =pvaluedata$anno3 ,y_position=c(2.5,3.5,5.5,6.5),tip_length = 0.005,
                        xmin = c(0.925:3.925),
                        xmax = c(1.075:4.075))+
  ggsignif::geom_signif(annotations =pvaluedata$anno2 ,y_position=c(1,2,4,5),tip_length = 0.005,
                        xmin = c(1.075:4.075),
                        xmax = c(1.225:4.225))+
  ggsignif::geom_signif(annotations =c("***","***","***","***") ,y_position=c(6,6,6.2,8.3),tip_length = 0.005,
                        xmin = c(0.775:3.775),
                        xmax = c(0.925:3.925))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0,8.5))+
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8))+
  scale_x_discrete(labels = c("[Min, 25%)","[25%, 50%)","[50%, 75%)","[75%, Max]"))+
  labs(y = "Density of m6A sites",x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "top",legend.direction = "vertical") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1)) + theme(legend.direction = "horizontal")
ggsave("/home/zhzhang/PG/m6a/Homo_sapiens.3gene_m6Adensity.expinter.pdf",
       phs,width = 5, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


##### 2.RBP相关分析
```r
#POSTAR3数据库 实验检测的RBP结合位点(HG38)
wget https://cloud.tsinghua.edu.cn/seafhttp/files/d05e6a97-0c75-4029-86a8-8b135368aae4/human.txt.gz
wget https://cloud.tsinghua.edu.cn/seafhttp/files/b21ab9c5-3e6a-48f0-a6ff-34a84939dc5a/mouse.txt.gz
gzip -d "/home/zhzhang/PG/RBP/human.txt.gz"
gzip -d "/home/zhzhang/PG/RBP/mouse.txt.gz"
#转bed【chr start end RBPname ID strand】(去掉RBP不明确的条目)
#人
grep -w -v RBP_occupancy /home/zhzhang/PG/RBP/human.txt|awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}'|sed 's/^chr//g'|bedtools sort > /home/zhzhang/PG/RBP/Homo_sapiens.RBP_BS.bed


```
```r
#RBP结合位点信息整合(获得每个RBP非冗余的结合位点，再合并到一起)（所有RBP所有潜在的结合位点）
#人(220种)
cat /home/zhzhang/PG/RBP/Homo_sapiens.RBP_BS.bed|awk '{print $4}'|sort|uniq > /home/zhzhang/PG/RBP/Homo_sapiens.RBPname.txt
################"/home/zhzhang/PG/RBP/hsRBPrmd.sh"内容
#!/bin/bash
cat /home/zhzhang/PG/RBP/Homo_sapiens.RBPname.txt|while read i
do
RBP=${i}
grep -w "${RBP}" /home/zhzhang/PG/RBP/Homo_sapiens.RBP_BS.bed|awk '$6=="+"{print $0}'|bedtools sort|bedtools merge|awk -v RBP="${RBP}" '{print $1"\t"$2"\t"$3"\t"RBP"\tRBPBS\t+"}' > /home/zhzhang/PG/RBP/hstmp/${RBP}.bed
grep -w "${RBP}" /home/zhzhang/PG/RBP/Homo_sapiens.RBP_BS.bed|awk '$6=="-"{print $0}'|bedtools sort|bedtools merge|awk -v RBP="${RBP}" '{print $1"\t"$2"\t"$3"\t"RBP"\tRBPBS\t-"}' >> /home/zhzhang/PG/RBP/hstmp/${RBP}.bed
cat /home/zhzhang/PG/RBP/hstmp/${RBP}.bed >> /home/zhzhang/PG/RBP/hstmp/Homo_sapiens.RBP_BS.NonRedundant.bed
echo "${RBP}" >> /home/zhzhang/PG/RBP/hstmp/Homo_sapiens.allRBP.log
done
bedtools sort -i /home/zhzhang/PG/RBP/hstmp/Homo_sapiens.RBP_BS.NonRedundant.bed > /home/zhzhang/PG/RBP/Homo_sapiens.RBP_BS.NonRedundant.sort.bed
###############################################################
nohup "/home/zhzhang/PG/RBP/hsRBPrmd.sh" &



#非冗余的RBP结合位点与基因外显子交集，获取每个基因外显子具有的RBP情况
#人
#提取所有基因外显子bed【chr start end geneid tid strand】
cat "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf"|awk '$3=="exon" && $13=="transcript_id"{print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $3=="exon" && $11=="transcript_id"{print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g"|bedtools sort > /home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon.bed
#交集
bedtools intersect -a /home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon.bed -b /home/zhzhang/PG/RBP/Homo_sapiens.RBP_BS.NonRedundant.sort.bed -wo -s > /home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon_intersect_RBP.txt
scp -P 5614 /home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon_intersect_RBP.txt zhzhang@47.102.45.227:/home/zhzhang/PG/RBP/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon_intersect_RBP.txt /home/zhzhang/PG/RBP/


#RBP结合位点与基因间区交集
bedtools intersect -a /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed -b /home/zhzhang/PG/RBP/Homo_sapiens.RBP_BS.NonRedundant.sort.bed -wo -s > /home/zhzhang/PG/RBP/Homo_sapiens.intergenic_intersect_RBP.txt



```
##### 统计每个基因全部外显子(全部转录本)&基因间区具有的RBP种类及对应的结合比例(结合长度/全部外显子长度)信息文件
```r
！！！生科院脚本跑
#函数统计每个基因全部外显子&基因间区具有的RBP种类及对应结合长度比例的信息文件
#a输入每个基因全部外显子与RBP交集文件，b输入每个基因间区与RBP交集文件，e输入全部基因外显子bed文件
#c输出每个基因外显子区域&基因间区拥有的RBP种类及其结合长度比例统计结果
geneexonRBP <- function(a,b,c,e){
  #导入启动子与RBP交集,获取每个基因全部外显子(全部)中每种RBP的结合长度
  allgeneexon_intersect_RBP <- data.table::fread(a, header=FALSE)%>%
    data.frame()%>%
    select(4,10,13)%>%
    separate(V10,c("RBP"),remove = T,sep = "_")%>%
    group_by(V4,RBP)%>%
    summarise(len=sum(V13))
  colnames(allgeneexon_intersect_RBP) <- c("geneid","RBP","num")
  #计算基因全部外显子（全部转录本）长度
  allgeneexon <- read.delim(e, header=FALSE)%>%
    mutate(len=(V3-V2+1)/1000)%>%
    group_by(V4)%>%
    summarise(len=sum(len))
  colnames(allgeneexon) <- c("geneid","len")
  #合并长度和RBP信息，计算结合长度占比（RBP总结合长度/基因全部外显子长度）
  allgeneexon_RBP_len <- left_join(allgeneexon_intersect_RBP,allgeneexon,by="geneid")%>%
    mutate(num=num/1000)%>%
    mutate(ratio=num/len)%>%
    select(1,2,5)
  #导入基因间区与RBP交集,获取每个基因间区每种RBP的结合长度比例
  intergenic_intersect_RBP <- data.table::fread(b, header=FALSE)%>%
    data.frame()%>%
    select(4,10,13)%>%
    separate(V10,c("RBP"),remove = T,sep = "_")%>%
    group_by(V4,RBP)%>%
    summarise(num=sum(V13))%>%
    mutate(num=num/1000)%>%
    mutate(ratio=num/3)%>%
    select(1,2,4)
  colnames(intergenic_intersect_RBP) <- c("geneid","RBP","ratio")
  #合并
  allrbpinfo <- rbind(allgeneexon_RBP_len,intergenic_intersect_RBP)
  #储存
  data.table::fwrite(allrbpinfo,file =c,sep = '\t',row.names = F,quote = F,col.names = T)
}

#人
geneexonRBP(a="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon_intersect_RBP.txt",
            b="/home/zhzhang/PG/RBP/Homo_sapiens.intergenic_intersect_RBP.txt",
            c="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.txt",
            e="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon.bed")
#小鼠
geneexonRBP(a="/home/zhzhang/PG/RBP/Mus_musculus.allgeneexon_intersect_RBP.txt",
            b="/home/zhzhang/PG/RBP/Mus_musculus.intergenic_intersect_RBP.txt",
            c="/home/zhzhang/PG/RBP/Mus_musculus.allgeneexonAintergenic.RBP.txt",
            e="/home/zhzhang/PG/RBP/Mus_musculus.allgeneexon.bed")

#传至实验室服务器/home/zhzhang/PG/RBP/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.txt /home/zhzhang/PG/RBP/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/RBP/Mus_musculus.allgeneexonAintergenic.RBP.txt /home/zhzhang/PG/RBP/
```
##### 对比蛋白编码/lncRNA基因转录本以及基因间区的RBP种类数量
```r
#函数根据基因具有的RBP种类及其结合比例信息文件，统计RBP种类数
#a输入基因&基因间区具有的RBP种类及其结合比例信息文件，b输入物种英文名，c输入物种基因分类文件，d输入基因间区文件
getRBP_NUM <- function(a,b,c,d){
  #导入基因具有的RBP种类及其结合比例信息文件,计算RBP种类数
  allgene_RBP <- read.delim(a)%>%
    mutate(ratio=1)%>%
    group_by(geneid)%>%
    summarise(RBPnum=sum(ratio))
  #导入基因分类
  geneid_class <- read.delim(c)
  #导入基因间区分类
  intergenic <- read.delim(d, header=FALSE)%>%
    select(1)%>%
    mutate(type="Random intergenic")
  colnames(intergenic)[1] <- "geneid"
  #合并
  allclass <- rbind(geneid_class,intergenic)
  geneid_class_RBP <- left_join(allclass,allgene_RBP,by="geneid")%>%
    filter(type!="Interference lncRNA")%>%
    mutate(sp=b)
  geneid_class_RBP$RBPnum[is.na(geneid_class_RBP$RBPnum)==T] <- 0
  return(geneid_class_RBP)
}
#人
hs_RBP_NUM <- getRBP_NUM(a="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.txt",
                         b="Human",
                         c="~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                         d="~/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.phastCons.txt")
#合并所有物种结果
ALLSP_RBP_NUM <- filter(hs_RBP_NUM,type!="Random intergenic")%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
ALLSP_RBP_NUM$type <- factor(ALLSP_RBP_NUM$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
data.table::fwrite(ALLSP_RBP_NUM,file ="/home/zhzhang/PG/RBP/SPhs_3gene_RBPrichness.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#统计RBP丰富度(排除未被检测到与RBP结合的基因)
tj <- group_by(filter(ALLSP_RBP_NUM,RBPnum!=0),sp,type)%>%
  summarise(mean=mean(RBPnum),median=median(RBPnum))
data.table::fwrite(tj,file ="/home/zhzhang/PG/RBP/SPhs_3gene_RBPrichness.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot(排除未被检测到与RBP结合的基因)
phs <- ggplot(data=filter(ALLSP_RBP_NUM,RBPnum!=0), aes(x=type,y=log10(RBPnum)))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position=c(2.3,2.1,2.5,2.3),
                        comparisons = list(c("Protein-coding","PAS lncRNA"),
                                           c("PAS lncRNA","PAA lncRNA"),
                                           c("PAS lncRNA","NPA lncRNA"),
                                           c("PAA lncRNA","NPA lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  scale_y_continuous(breaks = c(0,1,2),labels = c("0","1","2"))+
  labs(y = expression("l"*"o"*"g"[10]*"("*"R"*"B"*"P"~"r"*"i"*"c"*"h"*"n"*"e"*"s"*"s"*")"),
       x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.3gene_RBPrichness.pdf",
       phs,width = 3, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")








#
##画图对比不同表达水平区间内，三类基因RBP丰富度
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
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.spexp_geneid_class.txt",
                      c = "Human")
#不同物种三类基因表达量数据合并
ALLsp_TPM <- hs_TPM%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
#每个基因全部外显子&基因间区具有的m6a位点密度
ALLSP_m6a <- read.delim("/home/zhzhang/PG/RBP/SPhs_3gene_RBPrichness.txt")%>%
  select(1,3)
#合并
he <- left_join(ALLsp_TPM,ALLSP_m6a,by="geneid")
he$type <- factor(he$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
#
tj4exp <- group_by(he,sp)%>%
  summarise(fen25=quantile(maxTPM,0.25),fen50=quantile(maxTPM,0.50),
            fen75=quantile(maxTPM,0.75))
#增加表达区间
he <- mutate(he,expinter=case_when(sp=="Human" & maxTPM < tj4exp$fen25[1] ~ "[Min,25%)",
                                   sp=="Human" & maxTPM < tj4exp$fen50[1] ~ "[25%,50%)",
                                   sp=="Human" & maxTPM < tj4exp$fen75[1] ~ "[50%,75%)",
                                   sp=="Human" & maxTPM >= tj4exp$fen75[1] ~ "[75%,Max]"
))
#统计
he$expinter <- factor(he$expinter,levels = c("[Min,25%)","[25%,50%)","[50%,75%)","[75%,Max]"))
tj <- group_by(filter(he,RBPnum!=0),sp,type,expinter)%>%
  summarise(mean=mean(RBPnum),median=median(RBPnum))
data.table::fwrite(tj,file ="/home/zhzhang/PG/RBP/SPhs_3gene_RBPrichness.expinterval.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#每个表达区间内，RBPrichness差异性检验
he=filter(he,RBPnum!=0)
interval <- c("[Min,25%)","[25%,50%)","[50%,75%)","[75%,Max]")
pvaluedata <- data.frame()
for (i in 1:4) {
  interfor <- interval[i]
  pvalue1 <- wilcox.test(filter(he,expinter==interfor & type=="PAS lncRNA")$RBPnum,
                         filter(he,expinter==interfor & type=="NPA lncRNA")$RBPnum,
                         alternative ="greater")[["p.value"]]
  pvalue2 <- wilcox.test(filter(he,expinter==interfor & type=="PAA lncRNA")$RBPnum,
                         filter(he,expinter==interfor & type=="NPA lncRNA")$RBPnum,
                         alternative ="greater")[["p.value"]]
  pvalue3 <- wilcox.test(filter(he,expinter==interfor & type=="PAS lncRNA")$RBPnum,
                         filter(he,expinter==interfor & type=="PAA lncRNA")$RBPnum,
                         alternative ="greater")[["p.value"]]
  pvalue4 <- wilcox.test(filter(he,expinter==interfor & type=="Protein-coding")$RBPnum,
                         filter(he,expinter==interfor & type=="PAS lncRNA")$RBPnum,
                         alternative ="greater")[["p.value"]]
  pvaluedata[i,1] <- interfor
  pvaluedata[i,2] <- pvalue1
  pvaluedata[i,3] <- pvalue2
  pvaluedata[i,4] <- pvalue3
  pvaluedata[i,5] <- pvalue4
}
colnames(pvaluedata) <- c("expinter","pvalue1","pvalue2","pvalue3","pvalue4")
pvaluedata <- mutate(pvaluedata,anno1=case_when(pvalue1<0.001 ~ "***",
                                                pvalue1<0.01 ~ "**",
                                                pvalue1<0.05 ~ "*",
                                                pvalue1>=0.05 ~ "N.S."),
                     anno2=case_when(pvalue2<0.001 ~ "***",
                                     pvalue2<0.01 ~ "**",
                                     pvalue2<0.05 ~ "*",
                                     pvalue2>=0.05 ~ "N.S."),
                     anno3=case_when(pvalue3<0.001 ~ "***",
                                     pvalue3<0.01 ~ "**",
                                     pvalue3<0.05 ~ "*",
                                     pvalue3>=0.05 ~ "N.S."),
                     anno4=case_when(pvalue4<0.001 ~ "***",
                                     pvalue4<0.01 ~ "**",
                                     pvalue4<0.05 ~ "*",
                                     pvalue4>=0.05 ~ "N.S."))
data.table::fwrite(pvaluedata,file ="/home/zhzhang/PG/RBP/SPhs_3gene_RBPrichness.expinterval.difpvalue.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
phs <- ggplot(data=filter(he,RBPnum!=0), aes(x=expinter,y=log10(RBPnum)))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.6,notch=F,aes(fill=type))+
  ggsignif::geom_signif(annotations =pvaluedata$anno1 ,y_position=c(2.5),tip_length = 0.001,
                        xmin =c(0.925:3.925),
                        xmax =c(1.225:4.225))+
  ggsignif::geom_signif(annotations =pvaluedata$anno2 ,y_position=c(2.3),tip_length = 0.001,
                        xmin = c(1.075:4.075),
                        xmax = c(1.225:4.225))+
  ggsignif::geom_signif(annotations =pvaluedata$anno3 ,y_position=c(2.15),tip_length = 0.001,
                        xmin = c(0.925:3.925),
                        xmax = c(1.075:4.075))+
  ggsignif::geom_signif(annotations =pvaluedata$anno4 ,y_position=c(2.7),tip_length = 0.001,
                        xmin = c(0.775:3.775),
                        xmax = c(0.925:3.925))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0,2.7))+
  scale_y_continuous(breaks = c(0,1,2))+
  scale_x_discrete(labels = c("[Min, 25%)","[25%, 50%)","[50%, 75%)","[75%, Max]"))+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "top",legend.direction = "vertical") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  labs(y = expression("l"*"o"*"g"[10]*"("*"R"*"B"*"P"~"r"*"i"*"c"*"h"*"n"*"e"*"s"*"s"*")"),
       x =NULL,fill = NULL,color = NULL)+ theme(legend.direction = "horizontal")
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.3gene_RBPrichness.expinter.pdf",
       phs,width = 5, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")



```


##### 严格阈值对RBP richness影响
```r
#函数根据基因具有的RBP种类及其结合比例信息文件，统计RBP种类数
#a输入基因&基因间区具有的RBP种类及其结合比例信息文件，b输入物种英文名，c输入物种基因分类文件，d输入基因间区文件
getRBP_NUM <- function(a,b,c,d){
  #导入基因具有的RBP种类及其结合比例信息文件,计算RBP种类数
  allgene_RBP <- read.delim(a)%>%
    mutate(ratio=1)%>%
    group_by(geneid)%>%
    summarise(RBPnum=sum(ratio))
  #导入基因分类
  geneid_class <- read.delim(c)%>%select(2,4,7)%>%
    dplyr::rename(geneid=lncgid)
  #合并
  allclass <- geneid_class
  geneid_class_RBP <- left_join(allclass,allgene_RBP,by="geneid")%>%
    mutate(sp=b)
  geneid_class_RBP$RBPnum[is.na(geneid_class_RBP$RBPnum)==T] <- 0
  return(geneid_class_RBP)
}
#人
hs_RBP_NUM <- getRBP_NUM(a="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.txt",
                         b="Human",
                         c="/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt",
                         d="~/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.phastCons.txt")
hs_RBP_NUM=filter(hs_RBP_NUM,RBPnum!=0)
#
pas <- rbind(mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>200),cuttype="200"))
paa <- rbind(mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(hs_RBP_NUM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>200),cuttype="200"))
#合并
alldata <- rbind(pas,paa)%>%
  mutate(type=case_when(type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA"))
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","50","100","150","200"))
#统计储存
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(RBPnum),median=median(RBPnum))
#PAS
p1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=log10(RBPnum)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=log10(mean)),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=log10(mean),group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),test.args = c("greater"),
                        y_position=c(2.1),tip_length = 0.01,size=0.5,textsize=4)+
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0, 2.2))+
  labs(y = expression("l"*"o"*"g"[10]*"("*"R"*"B"*"P"~"r"*"i"*"c"*"h"*"n"*"e"*"s"*"s"*")"),
       x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.paslnc_RBPrichness.pdf",
       p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=log10(RBPnum)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=log10(mean)),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=log10(mean),group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),test.args = c("greater"),
                        y_position=c(2.1),tip_length = 0.01,size=0.5,textsize=4)+
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0, 2.2))+
  labs(y = expression("l"*"o"*"g"[10]*"("*"R"*"B"*"P"~"r"*"i"*"c"*"h"*"n"*"e"*"s"*"s"*")"),
       x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.paalnc_RBPrichness.pdf",
       p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


##### 根据基因具有的RBP种类及其结合长度信息文件，输出基因RBP结合长度矩阵（每行代表基因，每列代表一类RBP）
```r
!!!!!极为耗时，写到脚本里用/usr/bin/Rscript跑！！！！！
#函数根据基因&基因间区具有的RBP种类及其结合长度信息文件，输出基因&基因间区RBP结合长度矩阵（每行代表基因，每列代表一类RBP）
#a输入基因具有的RBP种类及其结合长度信息文件，b输出基因RBP结合长度矩阵文件（每行代表基因，每列代表一类RBP）
getgene_RBP_matrix <- function(a,b){
  #导入基因具有的RBP种类及其结合长度信息文件
  pcgApg_RBP <- read.delim(a)
  #统计全部RBP类型
  RBP_class <- select(pcgApg_RBP,RBP)%>%
    distinct(RBP)
  #统计全部geneid
  geneid <- select(pcgApg_RBP,geneid)%>%
    distinct(geneid)
  #先循环获得每行为一类RBP，每列代表基因的矩阵
  gene_RBP_matrix <- RBP_class
  genesample_num <- nrow(geneid)
  for (i in c(1:nrow(geneid))) {
    #提取出第i个基因的RBP结合长度信息
    genei <- geneid[i,1]#i
    genei_RBP <- filter(pcgApg_RBP,geneid==genei)%>%
      select(-1)
    colnames(genei_RBP)[2] <- genei
    #将基因i的各RBP结合长度信息以列增添的形式循环添加到总矩阵中
    gene_RBP_matrix <- left_join(gene_RBP_matrix,genei_RBP,by="RBP")
    gene_RBP_matrix[is.na(gene_RBP_matrix[,i+1])==T,i+1] <- 0#i+1
    #输出进度
    print(i/genesample_num)
  }
  #转置生成基因RBP结合长度矩阵（每行代表基因，每列代表一类RBP）
  gene_RBP_matrix_t <- column_to_rownames(gene_RBP_matrix,"RBP")%>%
    t()%>%
    data.frame()
  #储存
  data.table::fwrite(gene_RBP_matrix_t,file =b,sep = '\t',row.names = T,quote = F,col.names = T)
}

#小鼠
getgene_RBP_matrix(a="/home/zhzhang/PG/RBP/Mus_musculus.allgeneexonAintergenic.RBP.txt",
                   b="/home/zhzhang/PG/RBP/Mus_musculus.allgeneexonAintergenic.RBP.matrix.txt")
#人
getgene_RBP_matrix(a="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.txt",
                   b="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.matrix.txt")

#传至实验室服务器/home/zhzhang/PG/RBP/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.matrix.txt /home/zhzhang/PG/RBP/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/RBP/Mus_musculus.allgeneexonAintergenic.RBP.matrix.txt /home/zhzhang/PG/RBP/

```


##### 对比蛋白编码/lncRNA基因转录本以及基因间区的RBP多样性diversity
```r
#函数根据RBP结合比例矩阵文件，统计RBP多样性diversity
#a输入基因&基因间区RBP结合比例矩阵文件，b输入物种名，c输入基因分类文件,d输入基因间区文件
getRBP_diversity <- function(a,b,c,d){
  #导入RBP结合比例矩阵
  gene_RBP_matrix <- read.delim(a, row.names=1)
  #计算每个基因的RBP多样性
  gene_RBP_diversity <- data.frame("geneid"=rownames(gene_RBP_matrix),
                                   "simpson"=vegan::diversity(gene_RBP_matrix, index = "simpson"),
                                   "shannon"=vegan::diversity(gene_RBP_matrix, index = "shannon"))
  #导入基因分类文件
  geneid_class <- read.delim(c)
  #导入基因间区分类
  intergenic <- read.delim(d, header=FALSE)%>%
    select(1)%>%
    mutate(type="Random intergenic")
  colnames(intergenic)[1] <- "geneid"
  #合并分类信息
  allclass <- rbind(geneid_class,intergenic)
  #基因分类以及RBP多样性信息合并
  geneid_class_RBP <- left_join(allclass,gene_RBP_diversity,by="geneid")%>%
    filter(type!="Interference lncRNA")%>%
    mutate(sp=b)
  geneid_class_RBP$simpson[is.na(geneid_class_RBP$simpson)==T] <- 0
  geneid_class_RBP$shannon[is.na(geneid_class_RBP$shannon)==T] <- 0
  return(geneid_class_RBP)
}
#人
hs_RBP_diversity <- getRBP_diversity(a="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.matrix.txt",
                                     b="Human",
                                     c="~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                                     d="~/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.phastCons.txt")

#合并所有物种结果
ALLSP_RBP_diversity <- filter(hs_RBP_diversity,type!="Random intergenic")%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
ALLSP_RBP_diversity$type <- factor(ALLSP_RBP_diversity$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
data.table::fwrite(ALLSP_RBP_diversity,file ="/home/zhzhang/PG/RBP/SPhs_3gene_RBPdiversity.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#统计(排除基因间区，以及未被检测到与RBP结合的基因)
tj <- group_by(filter(ALLSP_RBP_diversity,type!="Random intergenic" & shannon!=0),sp,type)%>%
  summarise(shannonmean=mean(shannon),shannonmedian=median(shannon))
data.table::fwrite(tj,file ="/home/zhzhang/PG/RBP/SPhs_3gene_RBPdiversity.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot(排除基因间区，以及未被检测到与RBP结合的基因)
phsshannon <- ggplot(data=filter(filter(ALLSP_RBP_diversity,type!="Random intergenic" & shannon!=0),sp=="Human"), aes(x=type,y=shannon))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position=c(4.7,4.4,5.1,4.7),
                        comparisons = list(c("Protein-coding","PAS lncRNA"),
                                           c("PAS lncRNA","PAA lncRNA"),
                                           c("PAS lncRNA","NPA lncRNA"),
                                           c("PAA lncRNA","NPA lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  scale_y_continuous(breaks = c(0,1,2,3,4,5))+
  labs(y = "Shannon-Wiener index",
       x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =90)) +
  theme(axis.text.x = element_text(hjust=0.5,vjust = 0.5))
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.3gene_RBPdiversity.shannon.pdf",
       phsshannon,width = 3, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")













##画图对比不同表达水平区间内，三类基因RBP多样性
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
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.spexp_geneid_class.txt",
                      c = "Human")
#不同物种三类基因表达量数据合并
ALLsp_TPM <- hs_TPM%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
#每个基因全部外显子&基因间区具有的m6a位点密度
ALLSP_m6a <- read.delim("/home/zhzhang/PG/RBP/SPhs_3gene_RBPdiversity.txt")%>%
  select(1,4)
#合并
he <- left_join(ALLsp_TPM,ALLSP_m6a,by="geneid")
he$type <- factor(he$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
#
tj4exp <- group_by(he,sp)%>%
  summarise(fen25=quantile(maxTPM,0.25),fen50=quantile(maxTPM,0.50),
            fen75=quantile(maxTPM,0.75))
#增加表达区间
he <- mutate(he,expinter=case_when(sp=="Human" & maxTPM < tj4exp$fen25[1] ~ "[Min,25%)",
                                   sp=="Human" & maxTPM < tj4exp$fen50[1] ~ "[25%,50%)",
                                   sp=="Human" & maxTPM < tj4exp$fen75[1] ~ "[50%,75%)",
                                   sp=="Human" & maxTPM >= tj4exp$fen75[1] ~ "[75%,Max]"
))
#统计
he$expinter <- factor(he$expinter,levels = c("[Min,25%)","[25%,50%)","[50%,75%)","[75%,Max]"))
tj <- group_by(filter(he,shannon!=0),sp,type,expinter)%>%
  summarise(shannonmean=mean(shannon),shannonmedian=median(shannon),num=n())
data.table::fwrite(tj,file ="/home/zhzhang/PG/RBP/SPhs_3gene_RBPdiversity.expinterval.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#每个表达区间内，差异性检验
he=filter(he,shannon!=0)
interval <- c("[Min,25%)","[25%,50%)","[50%,75%)","[75%,Max]")
pvaluedata <- data.frame()
for (i in 1:4) {
  interfor <- interval[i]
  pvalue1 <- wilcox.test(filter(he,expinter==interfor & type=="PAS lncRNA")$shannon,
                        filter(he,expinter==interfor & type=="NPA lncRNA")$shannon,
                        alternative = "greater")[["p.value"]]
  pvalue2 <- wilcox.test(filter(he,expinter==interfor & type=="PAA lncRNA")$shannon,
                         filter(he,expinter==interfor & type=="NPA lncRNA")$shannon,
                         alternative ="greater")[["p.value"]]
  pvalue3 <- wilcox.test(filter(he,expinter==interfor & type=="PAS lncRNA")$shannon,
                         filter(he,expinter==interfor & type=="PAA lncRNA")$shannon,
                         alternative ="greater")[["p.value"]]
  pvalue4 <- wilcox.test(filter(he,expinter==interfor & type=="Protein-coding")$shannon,
                         filter(he,expinter==interfor & type=="PAS lncRNA")$shannon,
                         alternative ="greater")[["p.value"]]
  pvaluedata[i,1] <- interfor
  pvaluedata[i,2] <- pvalue1
  pvaluedata[i,3] <- pvalue2
  pvaluedata[i,4] <- pvalue3
  pvaluedata[i,5] <- pvalue4
}
colnames(pvaluedata) <- c("expinter","pvalue1","pvalue2","pvalue3","pvalue4")
pvaluedata <- mutate(pvaluedata,anno1=case_when(pvalue1<0.001 ~ "***",
                                                pvalue1<0.01 ~ "**",
                                                pvalue1<0.05 ~ "*",
                                                pvalue1>=0.05 ~ "N.S."),
                     anno2=case_when(pvalue2<0.001 ~ "***",
                                     pvalue2<0.01 ~ "**",
                                     pvalue2<0.05 ~ "*",
                                     pvalue2>=0.05 ~ "N.S."),
                     anno3=case_when(pvalue3<0.001 ~ "***",
                                     pvalue3<0.01 ~ "**",
                                     pvalue3<0.05 ~ "*",
                                     pvalue3>=0.05 ~ "N.S."),
                     anno4=case_when(pvalue4<0.001 ~ "***",
                                     pvalue4<0.01 ~ "**",
                                     pvalue4<0.05 ~ "*",
                                     pvalue4>=0.05 ~ "N.S."))
data.table::fwrite(pvaluedata,file ="/home/zhzhang/PG/RBP/SPhs_3gene_RBPdiversity.expinterval.difpvalue.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
phs <- ggplot(data=filter(he,shannon!=0), aes(x=expinter,y=shannon))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.6,notch=F,aes(fill=type))+
  ggsignif::geom_signif(annotations =pvaluedata$anno1 ,y_position=c(4.9),tip_length = 0.001,
                        xmin =c(0.925:3.925),
                        xmax =c(1.225:4.225))+
  ggsignif::geom_signif(annotations =pvaluedata$anno2 ,y_position=c(4.5),tip_length = 0.001,
                        xmin = c(1.075:4.075),
                        xmax = c(1.225:4.225))+
  ggsignif::geom_signif(annotations =pvaluedata$anno3 ,y_position=c(4.2),tip_length = 0.001,
                        xmin = c(0.925:3.925),
                        xmax = c(1.075:4.075))+
  ggsignif::geom_signif(annotations =pvaluedata$anno4 ,y_position=c(5.2),tip_length = 0.001,
                        xmin = c(0.775:3.775),
                        xmax = c(0.925:3.925))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0,5.5))+
  scale_y_continuous(breaks = c(0,1,2,3,4))+
  scale_x_discrete(labels = c("[Min, 25%)","[25%, 50%)","[50%, 75%)","[75%, Max]"))+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "top",legend.direction = "vertical") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  labs(y = "Shannon-Wiener index",
       x =NULL,fill = NULL,color = NULL)+
  theme(legend.direction = "horizontal")
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.3gene_RBPdiversity.expinter.pdf",
       phs,width = 5, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 严格阈值对RBP多样性影响
```r
#函数根据RBP结合比例矩阵文件，统计RBP多样性diversity
#a输入基因&基因间区RBP结合比例矩阵文件，b输入物种名，c输入基因分类文件,d输入基因间区文件
getRBP_diversity <- function(a,b,c,d){
  #导入RBP结合比例矩阵
  gene_RBP_matrix <- read.delim(a, row.names=1)
  #计算每个基因的RBP多样性
  gene_RBP_diversity <- data.frame("geneid"=rownames(gene_RBP_matrix),
                                   "simpson"=vegan::diversity(gene_RBP_matrix, index = "simpson"),
                                   "shannon"=vegan::diversity(gene_RBP_matrix, index = "shannon"))
  #导入基因分类文件
  geneid_class <- read.delim(c)%>%select(2,4,7)%>%
    dplyr::rename(geneid=lncgid)
  #合并分类信息
  allclass <- geneid_class
  #基因分类以及RBP多样性信息合并
  geneid_class_RBP <- left_join(allclass,gene_RBP_diversity,by="geneid")%>%
    mutate(sp=b)
  geneid_class_RBP$simpson[is.na(geneid_class_RBP$simpson)==T] <- 0
  geneid_class_RBP$shannon[is.na(geneid_class_RBP$shannon)==T] <- 0
  return(geneid_class_RBP)
}
#人
hs_RBP_diversity <- getRBP_diversity(a="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.matrix.txt",
                                     b="Human",
                                     c="/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt",
                                     d="~/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.phastCons.txt")
hs_RBP_diversity=filter(hs_RBP_diversity,shannon!=0)
#
pas <- rbind(mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>200),cuttype="200"))
paa <- rbind(mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(hs_RBP_diversity,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>200),cuttype="200"))
#合并
alldata <- rbind(pas,paa)%>%
  mutate(type=case_when(type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA"))
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","50","100","150","200"))
#统计储存
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(shannon),median=median(shannon))
#PAS
p1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=shannon))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),
                        y_position=c(4.4),tip_length = 0.01,size=0.5,textsize=4)+
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0, 4.6))+
  labs(y = "Shannon-Wiener index",
       x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.paslnc_RBPdiversity.shannon.pdf",
       p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=shannon))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),
                        y_position=c(4.4),tip_length = 0.01,size=0.5,textsize=4)+
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0, 4.6))+
  labs(y = "Shannon-Wiener index",
       x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.paalnc_RBPdiversity.shannon.pdf",
       p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


##### PAS lncRNA基因与假基因的母基因之间的RBP相似性
```r
#R函数给指定物种每个pglnc基因分配随机的蛋白编码基因
#最终生成至少有一种RBP结合的pglnc-母基因，pglnc-随机蛋白编码基因，随机lnc-母基因的基因id对应文件
#（a输入物种每个基因id的RBP丰富度文件,b输入pglnc和其对应的母基因文件，c输出全部需要计算表达相关性的三类基因对）
pgrandompcg <- function(a,b,c,d){
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
  #人类
  hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                        b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.spexp_geneid_class.txt",
                        c = "Human")
  tj4exp <- group_by(hs_TPM,sp)%>%
    summarise(fen25=quantile(maxTPM,0.25),fen50=quantile(maxTPM,0.50),
              fen75=quantile(maxTPM,0.75))
  #增加表达区间
  hs_TPM <- mutate(hs_TPM,expinter=case_when(sp=="Human" & maxTPM < tj4exp$fen25[1] ~ "[Min,25%)",
                                             sp=="Human" & maxTPM < tj4exp$fen50[1] ~ "[25%,50%)",
                                             sp=="Human" & maxTPM < tj4exp$fen75[1] ~ "[50%,75%)",
                                             sp=="Human" & maxTPM >= tj4exp$fen75[1] ~ "[75%,Max]"
  ))%>%
    select(1,5)
  #导入每个基因id的RBP丰富度文件，筛选出与至少一种RBP结合的基因
  geneid_class <- read.delim(a)%>%
    filter(sp=="Human" & RBPnum>0)%>%
    select(1,2)%>%
    left_join(hs_TPM,by="geneid")
  geneid_class <- geneid_class[is.na(geneid_class$expinter)==F,]
  #获取pas lnc基因id
  paslncid <- filter(geneid_class,type=="PAS lncRNA")%>%
    select(1)
  #获取paa lnc基因id
  paalncid <- filter(geneid_class,type=="PAA lncRNA")%>%
    select(1)
  #获取蛋白编码基因id
  codingid <- filter(geneid_class,type=="Protein-coding")%>%
    select(1)
  #获取npalnc基因id
  npalncid <- filter(geneid_class,type=="NPA lncRNA")%>%
    select(1)
  #导入palnc和其对应的母基因文件,仅保留pglnc基因和母基因都与至少一种RBP结合的基因对
  #pas
  pgdlncRNA_transtarget <- read.delim(b)%>%
    filter(type=="Pseudogene-associated sense lncRNA")%>%
    select(2,10)
  colnames(pgdlncRNA_transtarget)[1] <- c("geneid")
  pgdlncRNA_transtarget <- inner_join(pgdlncRNA_transtarget,paslncid,by="geneid")
  colnames(pgdlncRNA_transtarget) <- c("pglncid","geneid")
  pgdlncRNA_transtarget <- inner_join(pgdlncRNA_transtarget,codingid,by="geneid")
  colnames(pgdlncRNA_transtarget) <- c("geneid","target")
  pgdlncRNA_transtarget <- mutate(pgdlncRNA_transtarget,type="PAS lncRNA gene vs Parent gene of pseudogene")
  pas_parent=pgdlncRNA_transtarget
  #paa
  pgdlncRNA_transtarget <- read.delim(b)%>%
    filter(type=="Pseudogene-associated antisense lncRNA")%>%
    select(2,10)
  colnames(pgdlncRNA_transtarget)[1] <- c("geneid")
  pgdlncRNA_transtarget <- inner_join(pgdlncRNA_transtarget,paalncid,by="geneid")
  colnames(pgdlncRNA_transtarget) <- c("pglncid","geneid")
  pgdlncRNA_transtarget <- inner_join(pgdlncRNA_transtarget,codingid,by="geneid")
  colnames(pgdlncRNA_transtarget) <- c("geneid","target")
  pgdlncRNA_transtarget <- mutate(pgdlncRNA_transtarget,type="PAA lncRNA gene vs Parent gene of pseudogene")
  paa_parent=pgdlncRNA_transtarget
  #获取需要分配随机蛋白基因的paslnc基因，产生相应基因对
  pglncnum <- nrow(pas_parent)
  pgdlncRNA_random <- data.frame()
  for (i in 1:pglncnum) {
    forpglnc <- pas_parent$geneid[i]
    forparent <- pas_parent$target[i]
    forexclude <- filter(pas_parent,geneid==forpglnc)$target
    forexpinter <- filter(geneid_class,geneid==forparent)$expinter
    forselect <- setdiff(filter(geneid_class,type=="Protein-coding"&expinter==forexpinter)$geneid,forexclude)
    randomparent <- sample_n(data.frame(forselect),100)
    forpgdlncRNA_random <- cbind(data.frame(forpglnc),randomparent)
    pgdlncRNA_random <- rbind(pgdlncRNA_random,forpgdlncRNA_random)
  }
  paslncRNA_random <- mutate(pgdlncRNA_random,type="PAS lncRNA gene vs Random protein-coding gene")
  colnames(paslncRNA_random) <- colnames(pas_parent)
  #获取需要分配随机蛋白基因的paalnc基因，产生相应基因对
  pglncnum <- nrow(paa_parent)
  pgdlncRNA_random <- data.frame()
  for (i in 1:pglncnum) {
    forpglnc <- paa_parent$geneid[i]
    forparent <- paa_parent$target[i]
    forexclude <- filter(paa_parent,geneid==forpglnc)$target
    forexpinter <- filter(geneid_class,geneid==forparent)$expinter
    forselect <- setdiff(filter(geneid_class,type=="Protein-coding"&expinter==forexpinter)$geneid,forexclude)
    randomparent <- sample_n(data.frame(forselect),100)
    forpgdlncRNA_random <- cbind(data.frame(forpglnc),randomparent)
    pgdlncRNA_random <- rbind(pgdlncRNA_random,forpgdlncRNA_random)
  }
  paalncRNA_random <- mutate(pgdlncRNA_random,type="PAA lncRNA gene vs Random protein-coding gene")
  colnames(paalncRNA_random) <- colnames(paa_parent)
  #获取需要分配随机lnc基因的pas母基因,产生相应基因对
  parentnum <- nrow(pas_parent)
  random_parent <- data.frame()
  for (i in 1:parentnum) {
    forpglnc <- pas_parent$geneid[i]
    forparent <- pas_parent$target[i]
    forexclude <- filter(pas_parent,target==forparent)$geneid
    forexpinter <- filter(geneid_class,geneid==forpglnc)$expinter
    forselect <- setdiff(filter(geneid_class,type!="Protein-coding"&expinter==forexpinter)$geneid,forexclude)
    randomlnc <- sample_n(data.frame(forselect),100)
    forrandom_parent <- cbind(randomlnc,data.frame(forparent))
    random_parent <- rbind(random_parent,forrandom_parent)
  }
  random_pasparent <- mutate(random_parent,type="Random lncRNA gene vs Parent gene of pseudogene")
  colnames(random_pasparent) <- colnames(pas_parent)
  #获取需要分配随机lnc基因的paa母基因,产生相应基因对
  parentnum <- nrow(paa_parent)
  random_parent <- data.frame()
  for (i in 1:parentnum) {
    forpglnc <- paa_parent$geneid[i]
    forparent <- paa_parent$target[i]
    forexclude <- filter(paa_parent,target==forparent)$geneid
    forexpinter <- filter(geneid_class,geneid==forpglnc)$expinter
    forselect <- setdiff(filter(geneid_class,type!="Protein-coding"&expinter==forexpinter)$geneid,forexclude)
    randomlnc <- sample_n(data.frame(forselect),100)
    forrandom_parent <- cbind(randomlnc,data.frame(forparent))
    random_parent <- rbind(random_parent,forrandom_parent)
  }
  random_paaparent <- mutate(random_parent,type="Random lncRNA gene vs Parent gene of pseudogene")
  colnames(random_paaparent) <- colnames(paa_parent)
  #合并全部需要计算表达相关性的基因对
  allgenepair <- rbind(pas_parent,paslncRNA_random,random_pasparent)
  allgenepair2 <- rbind(paa_parent,paalncRNA_random,random_paaparent)
  #储存
  data.table::fwrite(allgenepair,file =c,sep = '\t',row.names = F,quote = F,col.names = T)
  data.table::fwrite(allgenepair2,file =d,sep = '\t',row.names = F,quote = F,col.names = T)
}
#人类
pgrandompcg(a = "/home/zhzhang/PG/RBP/SPhs_3gene_RBPrichness.txt",
            b = "/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt",
            c = "/home/zhzhang/PG/RBP/Homo_sapiens.three_genepair_needRBPdistance.pas.1.txt",
            d = "/home/zhzhang/PG/RBP/Homo_sapiens.three_genepair_needRBPdistance.paa.1.txt")



#函数根据RBP结合比例矩阵文件，计算三类基因对之间的RBP Bray-Curtis similarity index
#(a输入RBP结合比例矩阵文件，b输入基因对文件，c输入物种名
getpgpcgRBPdistance <- function(a,b,c){
  #导入RBP结合比例矩阵文件
  gene_RBP_matrix <- read.delim(a)
  #导入基因对
  three_genepair_needdis <- read.delim(b)
  #循环计算每一对基因对的RBP相似性
  genepair_disdata <- data.frame()
  for (i in c(1:nrow(three_genepair_needdis))) {
    #提取两个基因ID
    geneid <- three_genepair_needdis[i,1]#
    targetid <- three_genepair_needdis[i,2]#
    #提取基因对RBP矩阵
    genepair_RBP <- filter(gene_RBP_matrix,X==geneid | X==targetid)%>%
      column_to_rownames("X")
    #计算基因对的RBP相似性
    if (nrow(genepair_RBP)==2) {
      vegdist <- vegan::vegdist(genepair_RBP, binary=F,method="bray")
      #每对基因对的RBP相似性存入新的dataframe
      genepair_disdata[i,1] <- geneid
      genepair_disdata[i,2] <- targetid
      genepair_disdata[i,3] <- three_genepair_needdis[i,3]
      genepair_disdata[i,4] <- 1-vegdist
    }
    #输出完成量
    print(i/nrow(three_genepair_needdis)*100)
  }
  #基因对的RBP相似性dataframe清除无法计算的对(dis=NA)
  genepair_disdata <- genepair_disdata[is.na(genepair_disdata[,4])==F,]
  colnames(genepair_disdata) <- c("geneid","target","type","similarity")
  #增加物种信息
  genepair_disdata <- mutate(genepair_disdata,sp=c)
}
#人
hs_pgpcgRBPsimi <- getpgpcgRBPdistance(a="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.matrix.txt",
                                       b="/home/zhzhang/PG/RBP/Homo_sapiens.three_genepair_needRBPdistance.pas.1.txt",
                                       c="Human")%>%
  mutate(fil="PAS")
hs_pgpcgRBPsimi2 <- getpgpcgRBPdistance(a="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.matrix.txt",
                                       b="/home/zhzhang/PG/RBP/Homo_sapiens.three_genepair_needRBPdistance.paa.1.txt",
                                       c="Human")%>%
  mutate(fil="PAA")
hs_pgpcgRBPsimi=rbind(hs_pgpcgRBPsimi,hs_pgpcgRBPsimi2)
#物种RBP相似性数据储存
hs_pgpcgRBPsimi$type <- factor(hs_pgpcgRBPsimi$type,
                               levels = c("PAS lncRNA gene vs Parent gene of pseudogene",
                                          "PAA lncRNA gene vs Parent gene of pseudogene",
                                          "PAS lncRNA gene vs Random protein-coding gene",
                                          "PAA lncRNA gene vs Random protein-coding gene",
                                          "Random lncRNA gene vs Parent gene of pseudogene"))
data.table::fwrite(hs_pgpcgRBPsimi,file ="/home/zhzhang/PG/RBP/Homo_sapiens.three_genepair_RBPsimidata.PASaPAA.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#统计储存
tj <- group_by(hs_pgpcgRBPsimi,fil,type)%>%
  summarise(mean=mean(similarity),median=median(similarity),sd=sd(similarity))
data.table::fwrite(tj,file ="/home/zhzhang/PG/RBP/Homo_sapiens.three_genepair_RBPsimidata.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#人
phs <- ggplot(data = filter(hs_pgpcgRBPsimi,fil=="PAS"),aes(x=type,y=similarity))+
  geom_boxplot(notch =T,fatten = 3,width=0.5,outlier.alpha = 0,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position = c(0.51,0.55),tip_length = 0.01,
                        comparisons=list(c("PAS lncRNA gene vs Parent gene of pseudogene",
                                           "PAS lncRNA gene vs Random protein-coding gene"),
                                         c("PAS lncRNA gene vs Parent gene of pseudogene",
                                           "Random lncRNA gene vs Parent gene of pseudogene")
                        ))+
  scale_fill_manual(values = c("#B20339","#4E6AB5","#614499"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0, 0.6))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6),labels = c("0","0.2","0.4","0.6"))+
  scale_x_discrete(labels = c("","",""))+
  labs(x = NULL, y ="Bray-Curtis\nsimilarity index of RBP",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 13)) +
  theme(axis.ticks.x = element_line(size = 0)) +
  theme(legend.text = element_text(size = 14))
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.PASlnc_parent_RBPsimi.pdf", 
       phs,width = 7.8, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")

phs <- ggplot(data = filter(hs_pgpcgRBPsimi,fil=="PAA"),aes(x=type,y=similarity))+
  geom_boxplot(notch =T,fatten = 3,width=0.5,outlier.alpha = 0,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position = c(0.51,0.545)-0.24,tip_length = 0.01,
                        comparisons=list(c("PAA lncRNA gene vs Parent gene of pseudogene",
                                           "PAA lncRNA gene vs Random protein-coding gene"),
                                         c("PAA lncRNA gene vs Parent gene of pseudogene",
                                           "Random lncRNA gene vs Parent gene of pseudogene")
                        ))+
  scale_fill_manual(values = c("#B20339","#4E6AB5","#614499"))+
  cowplot::theme_half_open()+
  coord_cartesian(ylim = c(0, 0.35))+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3),labels = c("0","0.1","0.2","0.3"))+
  scale_x_discrete(labels = c("","",""))+
  labs(x = NULL, y ="Bray-Curtis\nsimilarity index of RBP",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 13)) +
  theme(axis.ticks.x = element_line(size = 0)) +
  theme(legend.text = element_text(size = 14))
ggsave("/home/zhzhang/PG/RBP/Homo_sapiens.PAAlnc_parent_RBPsimi.pdf", 
       phs,width = 7.8, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### RBP富集分析
```r
#导入每个基因id的分类 RBP丰富度文件
geneid_class <- read.delim("/home/zhzhang/PG/RBP/SPhs_3gene_RBPrichness.txt")
#富集注释
allgeneRBP <- read.delim("~/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.txt")%>%
  select(-3)
anno <- mutate(allgeneRBP,GO_Description=RBP)
colnames(anno)[2] <- "GO"
go2gene <- anno[, c(2, 1)]
go2name <- anno[, c(2, 3)]
#富集结果
set.seed(10)
npalnc <- filter(geneid_class,type=="NPA lncRNA")%>%
  select(1)%>%
  dplyr::sample_n(2000)
paslnc <- filter(geneid_class,type=="PAS lncRNA")%>%
  select(1)
paalnc <- filter(geneid_class,type=="PAA lncRNA")%>%
  select(1)
ego_npg <- clusterProfiler::enricher(npalnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                     maxGSSize = 40000)@result%>%
  mutate(type="NPA lncRNA")
ego_pas <- clusterProfiler::enricher(paslnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                    maxGSSize = 40000)@result%>%
  mutate(type="PAS lncRNA")
ego_paa <- clusterProfiler::enricher(paalnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                     maxGSSize = 40000)@result%>%
  mutate(type="PAA lncRNA")
#合并
ego <- rbind(ego_pas,ego_paa,ego_npg)%>%
  select(10,1:9)
#储存
data.table::fwrite(ego,file ="/home/zhzhang/PG/RBP/Human.lncRNA.type.RBP_enrich.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)

#plot
ego <- read.delim("/home/zhzhang/PG/RBP/Human.lncRNA.type.RBP_enrich.txt")

ego <- separate(ego,GeneRatio,c("GeneRatio1","GeneRatio2"))%>%
  separate(BgRatio,c("BgRatio1","BgRatio2"))%>%
  mutate(GeneRatio1=as.numeric(GeneRatio1),GeneRatio2=as.numeric(GeneRatio2),
         BgRatio1=as.numeric(BgRatio1),BgRatio2=as.numeric(BgRatio2))%>%
  mutate(ES=(GeneRatio1/GeneRatio2)/(BgRatio1/BgRatio2))
ego$type <- factor(ego$type,levels=c("NPA lncRNA","PAA lncRNA",
                                     "PAS lncRNA"))
top10paslncRBP <- filter(ego,type=="PAS lncRNA")%>%
  top_n(10,desc(p.adjust))%>%
  arrange(p.adjust)
top10paalncRBP <- filter(ego,type=="PAA lncRNA")%>%
  top_n(10,desc(p.adjust))%>%
  arrange(p.adjust)
#富集图
egoplot <- filter(ego,ID %in% c(top10paslncRBP$ID,top10paalncRBP$ID))
egoplot$ID <- factor(egoplot$ID,levels = c(setdiff(top10paslncRBP$ID,c("UCHL5","DDX3X")),
                                           c("UCHL5","DDX3X"),
                                           setdiff(top10paalncRBP$ID,c("UCHL5","DDX3X"))))
pp1 <- ggplot(data = egoplot,aes(x=ID,y=type))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_vline(xintercept=c(1.5:17.5),color="gray")+
  geom_hline(yintercept=c(1.5:2.5),color="gray")+
  scale_fill_gradient(low="white",high="#C12039")+
  scale_size(range = c(1,15))+
  scale_y_discrete(labels=c("NPA lncRNA","PAA lncRNA",
                            "PAS lncRNA"))+
  cowplot::theme_half_open()+
  labs(x = NULL, y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",linewidth=1))
ggsave("/home/zhzhang/PG/RBP/Human.lnc.enrichRBP.pdf", 
       pp1,width = 12, height = 3,dpi=1200, units = "in", device='pdf',bg = "transparent")


```








##### 3.in human通过共表达基因推断假基因相关lncRNA基因的功能
```r
#人类数据合并，去除睾丸样本
allsample_TPM <- read.delim("~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt", row.names=1)
allsample_TPM <- select(allsample_TPM,-grep("testis",colnames(allsample_TPM)))
dd_TPM <- read.delim("~/PG/RNAseq/Homo_sapiens/ddsample_TPM.txt", row.names=1)
dd_TPM <- select(dd_TPM,-grep("Testis",colnames(dd_TPM)))
he <- cbind(allsample_TPM,dd_TPM)
data.table::fwrite(he,file ="/home/zhzhang/PG/RNAseq/Homo_sapiens/cbindsample_TPM.rmtestis.txt",
                   sep = '\t',row.names = T,quote = F,col.names = T)
#小鼠数据合并，去除睾丸样本
allsample_TPM <- read.delim("/home/zhzhang/PG/RNAseq/Mus_musculus/allsample_TPM.txt", row.names=1)
allsample_TPM <- select(allsample_TPM,-grep("gonad_male",colnames(allsample_TPM)))
dd_TPM <- read.delim("~/PG/RNAseq/Mus_musculus/ddsample_TPM.txt", row.names=1)
dd_TPM <- select(dd_TPM,-grep("Testis",colnames(dd_TPM)))
he <- cbind(allsample_TPM,dd_TPM)
data.table::fwrite(he,file ="/home/zhzhang/PG/RNAseq/Mus_musculus/cbindsample_TPM.rmtestis.txt",
                   sep = '\t',row.names = T,quote = F,col.names = T)



#获取合并数据中表达的基因
#（a输入基因表达TPM矩阵，b输入基因分类文件，c输入物种名，d输出在物种中表达的基因id和类型文件）
three_gene_exppercent <- function(a,b,c,d){
  #导入基因表达TPM矩阵
  allsample_TPM <- read.delim(a, row.names=1)
  #导入基因ID分类
  geneid_class <- read.delim(b)
  #循环计算不同阈值下，三类基因中表达的基因
  gene_data <- data.frame()
  for (i in c(1)) {
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
  data.table::fwrite(geneforsave,file =d,sep = '\t',row.names = F,quote = F,col.names = F)
}
#输出在物种中表达的基因id和类型文件
#小鼠
three_gene_exppercent(a = "~/PG/RNAseq/Mus_musculus/cbindsample_TPM.rmtestis.txt",
                                  b = "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                                  c = "Mouse",
                                  d = "~/PG/RNAseq/Mus_musculus/Mus_musculus.cbindsample.rmtestis.spexp_geneid_class.txt")
#人
three_gene_exppercent(a = "~/PG/RNAseq/Homo_sapiens/cbindsample_TPM.rmtestis.txt",
                                  b = "~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                                  c = "Human",
                                  d = "~/PG/RNAseq/Homo_sapiens/Homo_sapiens.cbindsample.rmtestis.spexp_geneid_class.txt")
#输出均传至/home/zhzhang/PG/RNAseqdata/物种/featureCounts_out/




#导入人矩阵
allsample_TPM <- read.delim("/home/zhzhang/PG/RNAseq/Homo_sapiens/cbindsample_TPM.rmtestis.txt", row.names=1)
#导入表达的geneid
expgene <- read.delim("~/PG/RNAseq/Homo_sapiens/Homo_sapiens.cbindsample.rmtestis.spexp_geneid_class.txt",header = F)%>%
  select(1)
#保留矩阵中表达的基因
allsample_TPM <- allsample_TPM[rownames(allsample_TPM)%in%expgene$V1,]
data.table::fwrite(allsample_TPM,file ="/home/zhzhang/PG/RNAseq/Homo_sapiens/cbindsample_TPM.rmtestis.expgene.txt",
                   sep = '\t',row.names = T,quote = F,col.names = T)
#导入小鼠矩阵
allsample_TPM <- read.delim("/home/zhzhang/PG/RNAseq/Mus_musculus/cbindsample_TPM.rmtestis.txt", row.names=1)
#导入表达的geneid
expgene <- read.delim("~/PG/RNAseq/Mus_musculus/Mus_musculus.cbindsample.rmtestis.spexp_geneid_class.txt",header = F)%>%
  select(1)
#保留矩阵中表达的基因
allsample_TPM <- allsample_TPM[rownames(allsample_TPM)%in%expgene$V1,]
data.table::fwrite(allsample_TPM,file ="/home/zhzhang/PG/RNAseq/Mus_musculus/cbindsample_TPM.rmtestis.expgene.txt",
                   sep = '\t',row.names = T,quote = F,col.names = T)
#输出均传至生科院/home/zhzhang/PG/RNAseqdata/物种/featureCounts_out/


#人
#分裂表达的lncid（1000 line per file）
grep -w "Pseudogene-derived lncRNA" "/home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/Homo_sapiens.cbindsample.rmtestis.spexp_geneid_class.txt"|awk '{print $1}' > /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/Homo_sapiens.spexp_pgdlncgene.id.txt
grep -w "Non-pseudogene-derived lncRNA" "/home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/Homo_sapiens.cbindsample.rmtestis.spexp_geneid_class.txt"|awk '{print $1}' > /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/Homo_sapiens.spexp_npgdlncgene.id.txt
cat /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/Homo_sapiens.spexp_pgdlncgene.id.txt /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/Homo_sapiens.spexp_npgdlncgene.id.txt |split -l 1000 -d - /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/lncid/lncid
ls /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/lncid/ > /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/split.lncid.txt
#分裂pcgid（10 line per file）
grep -w "Protein-coding" "/home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/Homo_sapiens.cbindsample.rmtestis.spexp_geneid_class.txt"|awk '{print $1}' > /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/Homo_sapiens.spexp_pcgene.id.txt
split -l 10 -d /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/Homo_sapiens.spexp_pcgene.id.txt /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/pcgid/pcgid
ls /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/pcgid/ > /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/split.pcgid.txt
#小鼠
#分裂表达的lncid（1000 line per file）
grep -w "Pseudogene-derived lncRNA" "/home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/Mus_musculus.cbindsample.rmtestis.spexp_geneid_class.txt"|awk '{print $1}' > /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/Mus_musculus.spexp_pgdlncgene.id.txt
grep -w "Non-pseudogene-derived lncRNA" "/home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/Mus_musculus.cbindsample.rmtestis.spexp_geneid_class.txt"|awk '{print $1}' > /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/Mus_musculus.spexp_npgdlncgene.id.txt
cat /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/Mus_musculus.spexp_pgdlncgene.id.txt /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/Mus_musculus.spexp_npgdlncgene.id.txt |split -l 1000 -d - /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/lncid/lncid
ls /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/lncid/ > /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/split.lncid.txt
#分裂pcgid（10 line per file）
grep -w "Protein-coding" "/home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/Mus_musculus.cbindsample.rmtestis.spexp_geneid_class.txt"|awk '{print $1}' > /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/Mus_musculus.spexp_pcgene.id.txt
split -l 10 -d /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/Mus_musculus.spexp_pcgene.id.txt /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/pcgid/pcgid
ls /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/pcgid/ > /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/split.pcgid.txt



#找到全部共表达的lnc-pcg基因对
#R
library(magrittr)
#生成参数列表
option_list <- list(optparse::make_option(
  opt_str ="--lncid",
  type = "character",
  default = NULL,
  help="splited spexp protein-coding geneid txt file"),
  optparse::make_option(
    opt_str ="--pcgid",
    type = "character",
    default = NULL,
    help="splited spexp protein-coding geneid txt file"),
  optparse::make_option(
    opt_str ="--outputone",
    type = "character",
    default = NULL,
    help="splited out"
  ),
  optparse::make_option(
    opt_str ="--outputtwo",
    type = "character",
    default = NULL,
    help="splited and filtered out"
  ))
#将参数列表传递给一个变量
args <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
getgenepaircor <- function(a,b,c,d,e){
  #导入表达的lnc基因list（1000）
  lncid <- read.delim(b)
  colnames(lncid) <- "geneid"
  #导入表达的蛋白编码基因list（10）
  pcgid <- read.delim(c)
  colnames(pcgid) <- "geneid"
  #导入基因表达量TPM矩阵
  allsample_TPM <- data.table::fread(a)%>%
    data.frame()
  colnames(allsample_TPM)[1] <- "geneid"
  #lnc基因表达矩阵
  lncid_allsample_TPM <- dplyr::inner_join(allsample_TPM,lncid,by="geneid")%>%
    tibble::column_to_rownames("geneid")%>%
    t()%>%
    data.frame()
  #蛋白基因表达矩阵
  pcgid_allsample_TPM <- dplyr::inner_join(allsample_TPM,pcgid,by="geneid")%>%
    tibble::column_to_rownames("geneid")%>%
    t()%>%
    data.frame()
      #计算基因对的表达相关性
      cor <- psych::corr.test(x=lncid_allsample_TPM,
                              y=pcgid_allsample_TPM,
                              method="spearman",adjust="holm")
      perrson <- cor[["r"]]%>%
        data.frame()%>%
        tibble::rownames_to_column("geneid")%>%
        tidyr::gather("pcgid","cor",2:(nrow(pcgid)+1))
      pvalue <- cor[["p.adj"]]%>%
        data.frame()%>%
        tibble::rownames_to_column("geneid")%>%
        tidyr::gather("pcgid","pvalue",2:(nrow(pcgid)+1))%>%
        dplyr::select(pvalue)
      corandpadj <- cbind(perrson,pvalue)
  #保存
  data.table::fwrite(corandpadj,file =d,
                     sep = '\t',row.names = F,quote = F,col.names = F)
  #筛选
  filtercorlong <- dplyr::filter(corandpadj,pvalue<0.001 & abs(cor)>0.75)
  data.table::fwrite(filtercorlong,file =e,
                     sep = '\t',row.names = F,quote = F,col.names = F)
}
#人
getgenepaircor(a="/home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/cbindsample_TPM.rmtestis.expgene.txt",
               b=args$lncid,
               c=args$pcgid,
               d=args$outputone,
               e=args$outputtwo)
#小鼠
getgenepaircor(a="/home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/cbindsample_TPM.rmtestis.expgene.txt",
               b=args$lncid,
               c=args$pcgid,
               d=args$outputone,
               e=args$outputtwo)


#parallel 并行执行R脚本
#人
rm /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/out2/*
rm /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/out1/*
cat "/home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/split.lncid.txt"|while read i
do
echo "${i}"
cat "/home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/split.pcgid.txt" | parallel -j 120 --verbose Rscript "/home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/sh/hs.R" --lncid /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/lncid/${i} --pcgid /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/pcgid/{} --outputone /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/out1/${i}.{}.out --outputtwo /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/out2/${i}.{}.out
done
#小鼠
rm /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/out2/*
rm /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/out1/*
cat "/home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/split.lncid.txt"|while read i
do
echo "${i}"
cat "/home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/split.pcgid.txt" | parallel -j 28 --verbose Rscript "/home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/sh/mm.R" --lncid /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/lncid/${i} --pcgid /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/pcgid/{} --outputone /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/out1/${i}.{}.out --outputtwo /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/out2/${i}.{}.out
done




#合并cor并行结果
#人
cat /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/out2/* > /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/merge/Homo_sapiens.spexp.lnc_pcgene.cordata.filter.txt
#获取全部lnc的coexp蛋白编码基因id（作为功能富集的bgd）
awk '{print $2}' "/home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/merge/Homo_sapiens.spexp.lnc_pcgene.cordata.filter.txt"|sort|uniq > /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/coexp_new/merge/Homo_sapiens.spexp.lnc_coexp_pcg.bgd.id.txt
#结果传至/home/zhzhang/PG/RNAseq/物种/addtarget_new



###pglnc富集的coexp基因功能富集
#导入基因分类
geneid_class <- read.delim("/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.cbindsample.rmtestis.spexp_geneid_class.txt", header=FALSE)
colnames(geneid_class) <- c("geneid","type")
#coexp富集注释
pglnc_pcgene <- read.delim("/home/zhzhang/PG/RNAseq/Homo_sapiens/addtarget_new/Homo_sapiens.spexp.lnc_pcgene.cordata.filter.txt", header=FALSE)%>%
  select(-3,-4)
anno <- mutate(pglnc_pcgene,GO_Description=V2)
colnames(anno)[1] <- "geneid"
colnames(anno)[2] <- "GO"
go2gene <- anno[, c(2, 1)]
go2name <- anno[, c(2, 3)]
#coexp富集结果
set.seed(10)
npalnc <- filter(geneid_class,type=="Non-pseudogene-associated lncRNA")%>%
  select(1)%>%
  dplyr::sample_n(2000)
paslnc <- filter(geneid_class,type=="Pseudogene-associated sense lncRNA")%>%
  select(1)
paalnc <- filter(geneid_class,type=="Pseudogene-associated antisense lncRNA")%>%
  select(1)
ego_npg <- clusterProfiler::enricher(npalnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "none",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                     maxGSSize = 40000)@result%>%
  mutate(type="NPA lncRNA")
ego_pas <- clusterProfiler::enricher(paslnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "none",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                    maxGSSize = 40000)@result%>%
  mutate(type="PAS lncRNA")
ego_paa <- clusterProfiler::enricher(paalnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "none",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                     maxGSSize = 40000)@result%>%
  mutate(type="PAA lncRNA")
#储存
data.table::fwrite(rbind(ego_pas,ego_paa,ego_npg),file ="/home/zhzhang/PG/RNAseq/Homo_sapiens/addtarget_new/Homo_sapiens.spexp.lnc.coexpgene_enrich.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#富集的coexp基因的p值分布
pvalued <- rbind(ego_pas,ego_paa,ego_npg)%>%
  filter(pvalue < 0.05)%>%
  mutate(ptype=case_when(pvalue < 0.00001 ~ "[Min, 1e-5)",
                         pvalue < 0.0001 ~ "[1e-5, 1e-4)",
                         pvalue < 0.001 ~ "[1e-4, 1e-3)",
                         pvalue < 0.01 ~ "[1e-3, 1e-2)",
                         pvalue < 0.05 ~ "[1e-2, 5e-2)"
  ))
pvaluedtj <- group_by(pvalued,type,ptype)%>%
  summarise(num=n())
data.table::fwrite(pvaluedtj,file ="/home/zhzhang/PG/RNAseq/Homo_sapiens/addtarget_new/Homo_sapiens.spexp.lnc.coexpgene_enrichnum.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
pvaluedtj$ptype <- factor(pvaluedtj$ptype,levels = c( "[Min, 1e-5)",
                                                      "[1e-5, 1e-4)",
                                                      "[1e-4, 1e-3)",
                                                      "[1e-3, 1e-2)","[1e-2, 5e-2)"))
pvaluedtj$type <- factor(pvaluedtj$type,levels = c("PAS lncRNA","PAA lncRNA",
                                                   "NPA lncRNA"))
pp <- ggplot(data = pvaluedtj,aes(x=ptype,y=num))+
  geom_col(aes(fill=type,color=type),width = 0.5)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                             "NPA lncRNA","Random intergenic"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA","Random intergenic"))+
  theme_half_open()+
  labs(x = "P-values (Fisher's exact test)", 
       y ="Number of co-expressed\ngenes enriched for lncRNA genes",
       fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  facet_wrap(~type,nrow=1,strip.position="top",as.table=F,scales="free_y")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))
ggsave("/home/zhzhang/PG/RNAseq/Homo_sapiens/addtarget_new/Homo_sapiens.spexp.lnc.coexpgene_enrichnum.pdf", 
       pp,width = 10, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")


#功能富集
#导入背景基因集
bgd <- read.table("/home/zhzhang/PG/RNAseq/Homo_sapiens/addtarget_new/Homo_sapiens.spexp.lnc_coexp_pcg.bgd.id.txt",
                  quote="\"", comment.char="")
colnames(bgd) <- "geneid"
#npglnc富集的coexp coding gene功能富集
npgcoexp_GO <- clusterProfiler::enrichGO(gene = filter(pvalued,type=="NPA lncRNA")$ID,universe=bgd$geneid,OrgDb = 'org.Hs.eg.db', keyType = "ENSEMBL", ont = 'BP', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 1)
npgcoexp_GO_result <- npgcoexp_GO@result%>%
  filter(p.adjust<0.05)
#paalnc富集的coexp coding gene功能富集
paacoexp_GO <- clusterProfiler::enrichGO(gene = filter(pvalued,type=="PAA lncRNA")$ID,universe=bgd$geneid,OrgDb = 'org.Hs.eg.db', keyType = "ENSEMBL", ont = 'BP', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 1)
paacoexp_GO_result <- paacoexp_GO@result%>%
  filter(p.adjust<0.05)%>%
  mutate(type="PAA lncRNA")
#pas
pascoexp_GO <- clusterProfiler::enrichGO(gene = filter(pvalued,type=="PAS lncRNA")$ID,universe=bgd$geneid,OrgDb = 'org.Hs.eg.db', keyType = "ENSEMBL", ont = 'BP', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 1)
pascoexp_GO_result <- pascoexp_GO@result%>%
  filter(p.adjust<0.05)%>%
  mutate(type="PAS lncRNA")
#储存
data.table::fwrite(rbind(pascoexp_GO_result,paacoexp_GO_result),file ="/home/zhzhang/PG/RNAseq/Homo_sapiens/addtarget_new/Homo_sapiens.spexp.palnc.enrichcoexpgene.GOBPenrich.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot_paslnc富集的coexp coding gene功能富集
res=c("mRNA metabolic process",
      "histone modification",
      "chromatin organization",
      "RNA splicing",
      "chromatin remodeling",
      "RNA catabolic process",
      "DNA repair",
      "double-strand break repair",
      "regulation of response to DNA damage stimulus",
      "regulation of RNA stability")
egoplot <- filter(pascoexp_GO_result,Description %in% res)
egoplot$Description[egoplot$Description=="regulation of response to DNA damage stimulus"] <- "regulation of response\nto DNA damage stimulus"
egoplot$Description <- factor(egoplot$Description,levels = egoplot$Description)
pp1 <- ggplot(data = egoplot,aes(x=Description,y=-log10(p.adjust)))+
  geom_col(color="gray",width = 0.001)+
  geom_point(aes(fill=-log10(p.adjust),size=Count),shape=21,color="black")+
  scale_fill_gradient(low="#E7EDF1",high="#7197AD",limits=c(1,25))+
  scale_size(range = c(5,20),breaks = c(100,200,300))+
  cowplot::theme_half_open()+
  labs(x = NULL, y ="-log10(p.adjust)",size="Count")+
  coord_cartesian(ylim = c(1, 23))+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 18),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(plot.margin = margin(t = 1,  # 顶部边缘距离
                             r = 1,  # 右边边缘距离
                             b = 0.5,  # 底部边缘距离
                             l = 1,  # 左边边缘距离
                             unit = "cm"))
ggsave("/home/zhzhang/PG/RNAseq/Homo_sapiens/addtarget_new/Homo_sapiens.spexp.paslnc.enrichcoexpgene.less1e-5.top10enrich.pdf", 
       pp1,width = 12.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#plot_paalnc富集的coexp coding gene功能富集
res=c("cilium organization",
      "cilium assembly",
      "cell-cell adhesion via plasma-membrane adhesion molecules",
      "cell projection assembly",
      "synapse organization",
      "regulation of trans-synaptic signaling",
      "synapse assembly",
      "chemical synaptic transmission",
      "cell junction organization",
      "vesicle-mediated transport in synapse")
egoplot <- filter(paacoexp_GO_result,Description %in% res)
egoplot$Description[egoplot$Description=="cell-cell adhesion via plasma-membrane adhesion molecules"] <- "cell-cell adhesion via plasma-\nmembrane adhesion molecules"
egoplot$Description[egoplot$Description=="regulation of trans-synaptic signaling"] <- "regulation of\ntrans-synaptic signaling"
egoplot$Description[egoplot$Description=="vesicle-mediated transport in synapse"] <- "vesicle-mediated\ntransport in synapse"
egoplot$Description <- factor(egoplot$Description,levels = egoplot$Description)
pp2 <- ggplot(data = egoplot,aes(x=Description,y=-log10(p.adjust)))+
  geom_col(color="gray",width = 0.001)+
  geom_point(aes(fill=-log10(p.adjust),size=Count),shape=21,color="black")+
  scale_fill_gradient(low="white",high="#5B794F",limits=c(1,6))+
  scale_size(range = c(5,15),breaks = c(40,60,80))+
  theme_half_open()+
  labs(x = NULL, y ="-log10(p.adjust)",size="Count")+
  coord_cartesian(ylim = c(0.3, 6.5))+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(plot.margin = margin(t = 1,  # 顶部边缘距离
                             r = 1,  # 右边边缘距离
                             b = 0.25,  # 底部边缘距离
                             l = 1,  # 左边边缘距离
                             unit = "cm"))
ggsave("/home/zhzhang/PG/RNAseq/Homo_sapiens/addtarget_new/Homo_sapiens.spexp.paalnc.enrichcoexpgene.less1e-5.top10enrich.pdf", 
       pp2,width = 12, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 4.in小鼠通过共表达基因推断假基因相关lncRNA基因的功能
```r
#合并cor并行结果
cat /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/out2/* > /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/merge/Mus_musculus.spexp.lnc_pcgene.cordata.filter.txt
#获取全部lnc的coexp蛋白编码基因id（富集的bgd）
awk '{print $2}' "/home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/merge/Mus_musculus.spexp.lnc_pcgene.cordata.filter.txt"|sort|uniq > /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/coexp_new/merge/Mus_musculus.spexp.lnc_coexp_pcg.bgd.id.txt
#结果传至/home/zhzhang/PG/RNAseq/物种/addtarget_new



###pglnc富集的coexp基因功能富集
#导入基因分类
geneid_class <- read.delim("/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.cbindsample.rmtestis.spexp_geneid_class.txt", header=FALSE)
colnames(geneid_class) <- c("geneid","type")
#coexp富集注释
pglnc_pcgene <- read.delim("/home/zhzhang/PG/RNAseq/Mus_musculus/addtarget_new/Mus_musculus.spexp.lnc_pcgene.cordata.filter.txt", header=FALSE)%>%
  select(-3,-4)
anno <- mutate(pglnc_pcgene,GO_Description=V2)
colnames(anno)[1] <- "geneid"
colnames(anno)[2] <- "GO"
go2gene <- anno[, c(2, 1)]
go2name <- anno[, c(2, 3)]
#coexp富集结果
set.seed(10)
npalnc <- filter(geneid_class,type=="Non-pseudogene-associated lncRNA")%>%
  select(1)%>%
  dplyr::sample_n(2000)
paslnc <- filter(geneid_class,type=="Pseudogene-associated sense lncRNA")%>%
  select(1)
paalnc <- filter(geneid_class,type=="Pseudogene-associated antisense lncRNA")%>%
  select(1)
ego_npg <- clusterProfiler::enricher(npalnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "none",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                     maxGSSize = 40000)@result%>%
  mutate(type="NPA lncRNA")
ego_pas <- clusterProfiler::enricher(paslnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "none",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                    maxGSSize = 40000)@result%>%
  mutate(type="PAS lncRNA")
ego_paa <- clusterProfiler::enricher(paalnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "none",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                     maxGSSize = 40000)@result%>%
  mutate(type="PAA lncRNA")
#储存
data.table::fwrite(rbind(ego_pas,ego_paa,ego_npg),file ="/home/zhzhang/PG/RNAseq/Mus_musculus/addtarget_new/Mus_musculus.spexp.lnc.coexpgene_enrich.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#富集的coexp基因的p值分布
pvalued <- rbind(ego_pas,ego_paa,ego_npg)%>%
  filter(pvalue < 0.05)%>%
  mutate(ptype=case_when(pvalue < 0.00001 ~ "[Min, 1e-5)",
                         pvalue < 0.0001 ~ "[1e-5, 1e-4)",
                         pvalue < 0.001 ~ "[1e-4, 1e-3)",
                         pvalue < 0.01 ~ "[1e-3, 1e-2)",
                         pvalue < 0.05 ~ "[1e-2, 5e-2)"
  ))
pvaluedtj <- group_by(pvalued,type,ptype)%>%
  summarise(num=n())
data.table::fwrite(pvaluedtj,file ="/home/zhzhang/PG/RNAseq/Mus_musculus/addtarget_new/Mus_musculus.spexp.lnc.coexpgene_enrichnum.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
pvaluedtj$ptype <- factor(pvaluedtj$ptype,levels = c( "[Min, 1e-5)",
                                                      "[1e-5, 1e-4)",
                                                      "[1e-4, 1e-3)",
                                                      "[1e-3, 1e-2)","[1e-2, 5e-2)"))
pvaluedtj$type <- factor(pvaluedtj$type,levels = c("PAS lncRNA","PAA lncRNA",
                                                   "NPA lncRNA"))
pp <- ggplot(data = pvaluedtj,aes(x=ptype,y=num))+
  geom_col(aes(fill=type,color=type),width = 0.5)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                             "NPA lncRNA","Random intergenic"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA","Random intergenic"))+
  theme_half_open()+
  labs(x = "P-values (Fisher's exact test)", 
       y ="Number of co-expressed\ngenes enriched for lncRNA genes",
       fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "none", legend.direction = "horizontal")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  facet_wrap(~type,nrow=1,strip.position="top",as.table=F,scales="free_y")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))
ggsave("/home/zhzhang/PG/RNAseq/Mus_musculus/addtarget_new/Mus_musculus.spexp.lnc.coexpgene_enrichnum.pdf", 
       pp,width = 10, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")


#功能富集
#导入背景基因集
bgd <- read.table("/home/zhzhang/PG/RNAseq/Mus_musculus/addtarget_new/Mus_musculus.spexp.lnc_coexp_pcg.bgd.id.txt",
                  quote="\"", comment.char="")
colnames(bgd) <- "geneid"
#npglnc富集的coexp coding gene功能富集
npgcoexp_GO <- clusterProfiler::enrichGO(gene = filter(pvalued,type=="NPA lncRNA")$ID,universe=bgd$geneid,OrgDb = 'org.Mm.eg.db', keyType = "ENSEMBL", ont = 'BP', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 1)
npgcoexp_GO_result <- npgcoexp_GO@result%>%
  filter(p.adjust<0.05)
#paalnc富集的coexp coding gene功能富集
paacoexp_GO <- clusterProfiler::enrichGO(gene = filter(pvalued,type=="PAA lncRNA")$ID,universe=bgd$geneid,OrgDb = 'org.Mm.eg.db', keyType = "ENSEMBL", ont = 'BP', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 1)
paacoexp_GO_result <- paacoexp_GO@result%>%
  filter(p.adjust<0.05)%>%
  mutate(type="PAA lncRNA")
#pas
pascoexp_GO <- clusterProfiler::enrichGO(gene = filter(pvalued,type=="PAS lncRNA")$ID,universe=bgd$geneid,OrgDb = 'org.Mm.eg.db', keyType = "ENSEMBL", ont = 'BP', pAdjustMethod = 'fdr', pvalueCutoff = 0.05, qvalueCutoff = 1)
pascoexp_GO_result <- pascoexp_GO@result%>%
  filter(p.adjust<0.05)%>%
  mutate(type="PAS lncRNA")
#储存
data.table::fwrite(rbind(pascoexp_GO_result,paacoexp_GO_result),file ="/home/zhzhang/PG/RNAseq/Mus_musculus/addtarget_new/Mus_musculus.spexp.palnc.enrichcoexpgene.GOBPenrich.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot_paslnc富集的coexp coding gene功能富集
res=c("regulation of mRNA metabolic process",
      "histone modification",
      "chromatin organization",
      "RNA splicing",
      "chromatin remodeling",
      "RNA catabolic process",
      "DNA repair",
      "double-strand break repair",
      "regulation of response to DNA damage stimulus",
      "regulation of RNA stability")
egoplot <- filter(pascoexp_GO_result,Description %in% res)
egoplot$Description[egoplot$Description=="regulation of mRNA metabolic process"] <- "regulation of\nmRNA metabolic process"
egoplot$Description[egoplot$Description=="regulation of response to DNA damage stimulus"] <- "regulation of response\nto DNA damage stimulus"
egoplot$Description <- factor(egoplot$Description,levels = egoplot$Description)
pp1 <- ggplot(data = egoplot,aes(x=Description,y=-log10(p.adjust)))+
  geom_col(color="gray",width = 0.001)+
  geom_point(aes(fill=-log10(p.adjust),size=Count),shape=21,color="black")+
  scale_fill_gradient(low="#E7EDF1",high="#7197AD",limits=c(1,15))+
  scale_size(range = c(5,15),breaks = c(40,60,80))+
  theme_half_open()+
  labs(x = NULL, y ="-log10(p.adjust)",size="Count")+
  coord_cartesian(ylim = c(1, 17))+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(plot.margin = margin(t = 1,  # 顶部边缘距离
                             r = 1,  # 右边边缘距离
                             b = 0.5,  # 底部边缘距离
                             l = 1,  # 左边边缘距离
                             unit = "cm"))
ggsave("/home/zhzhang/PG/RNAseq/Mus_musculus/addtarget_new/Mus_musculus.spexp.paslnc.enrichcoexpgene.less1e-5.top10enrich.pdf", 
       pp1,width = 12, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#plot_paalnc富集的coexp coding gene功能富集
res=c("regulation of mRNA metabolic process",
      "histone modification",
      "chromatin organization",
      "RNA splicing",
      "chromatin remodeling",
      "RNA catabolic process",
      "DNA repair",
      "double-strand break repair",
      "regulation of response to DNA damage stimulus",
      "regulation of RNA stability")
egoplot <- filter(paacoexp_GO_result,Description %in% res)
egoplot$Description[egoplot$Description=="regulation of mRNA metabolic process"] <- "regulation of\nmRNA metabolic process"
egoplot$Description[egoplot$Description=="regulation of response to DNA damage stimulus"] <- "regulation of response\nto DNA damage stimulus"
egoplot$Description <- factor(egoplot$Description,levels = egoplot$Description)
pp2 <- ggplot(data = egoplot,aes(x=Description,y=-log10(p.adjust)))+
  geom_col(color="gray",width = 0.001)+
  geom_point(aes(fill=-log10(p.adjust),size=Count),shape=21,color="black")+
  scale_fill_gradient(low="white",high="#5B794F",limits=c(1,20))+
  scale_size(range = c(5,15),breaks = c(60,100,140))+
  theme_half_open()+
  labs(x = NULL, y ="-log10(p.adjust)",size="Count")+
  coord_cartesian(ylim = c(1, 21))+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(plot.margin = margin(t = 1,  # 顶部边缘距离
                             r = 1,  # 右边边缘距离
                             b = 0.25,  # 底部边缘距离
                             l = 1,  # 左边边缘距离
                             unit = "cm"))
ggsave("/home/zhzhang/PG/RNAseq/Mus_musculus/addtarget_new/Mus_musculus.spexp.paalnc.enrichcoexpgene.less1e-5.top10enrich.pdf", 
       pp2,width = 12, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#plot_npalnc富集的coexp coding gene功能富集
res=c("chromosome segregation",
      "nuclear chromosome segregation",
      "nuclear division",
      "mitotic spindle organization",
      "DNA replication",
      "mitotic sister chromatid segregation",
      "meiotic cell cycle",
      "chromosome separation",
      "cell division",
      "cell cycle phase transition")
egoplot <- filter(npgcoexp_GO_result,Description %in% res)
egoplot$Description <- factor(egoplot$Description,levels = egoplot$Description)
pp3 <- ggplot(data = egoplot,aes(x=Description,y=-log10(p.adjust)))+
  geom_col(color="gray",width = 0.001)+
  geom_point(aes(fill=-log10(p.adjust),size=Count),shape=21,color="black")+
  scale_fill_gradient(low="white",high="#F8A364",limits=c(1,16))+
  scale_size(range = c(5,15),breaks = c(25,50))+
  theme_half_open()+
  labs(x = NULL, y ="-log10(p.adjust)",size="Count")+
  coord_cartesian(ylim = c(4, 17))+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(plot.margin = margin(t = 1,  # 顶部边缘距离
                             r = 1,  # 右边边缘距离
                             b = 0.25,  # 底部边缘距离
                             l = 1,  # 左边边缘距离
                             unit = "cm"))
ggsave("/home/zhzhang/PG/RNAseq/Mus_musculus/addtarget_new/Mus_musculus.spexp.npalnc.enrichcoexpgene.less1e-5.top10enrich.pdf", 
       pp3,width = 12, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```





