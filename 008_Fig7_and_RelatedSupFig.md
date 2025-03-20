#### 八. 小鼠多组织转录组数据进行衰老相关分析
##### 1.小鼠衰老相关转录组数据每个组织差异表达分析
```r
a="/home/zhzhang/PG/RNAseq/Mus_musculus/agingsample_TPM.txt"
b="/home/zhzhang/PG/RNAseq/Mus_musculus/agingsample_readscount.txt"
tissue=c("Aorta","Brain","Heart","Kidney","Liver","Lung","Muscle","Skin")
#导入基因表达TPM矩阵
allsample_TPM <- read.delim(a, row.names=1)
#导入基因表达readcount矩阵
allsample_RC <- read.delim(b, row.names=1)
#
allDE=data.frame()
for (i in 1:8) {
  fort=tissue[i]
  #根据TPM>=1提取指定组织表达的基因id
  fort_TPM=allsample_TPM[,grep(fort,colnames(allsample_TPM))]
  gene_exp <- data.frame(max=apply(fort_TPM,1,max))%>%
    rownames_to_column("geneid")%>%
    filter(max>=1)%>%
    select(geneid)
  rm(fort_TPM)
  gc()
  #提取指定组织RC矩阵,保留表达的基因
  for_RC=allsample_RC[,grep(fort,colnames(allsample_RC))]
  for_RC=for_RC[rownames(for_RC) %in% gene_exp$geneid,]
  #生成样本信息表
  sample_info=data.frame(row.names = colnames(for_RC),row=colnames(for_RC))%>%
    separate(row,c("one","age","three"),remove = T)%>%
    select(age)
  sample_info$age=factor(sample_info$age,levels = c("6mo","24mo","30mo"))
  #构建DEseq矩阵
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = for_RC, colData = sample_info, design= ~age)
  #DE分析
  dds <- DESeq2::DESeq(dds)
  #根据指定的样品组别导出DE分析结果
  res <- DESeq2::results(dds,contrast=c("age","24mo","6mo"),tidy=T,pAdjustMethod="fdr")%>%
    dplyr::rename(geneid=row,log2FC_24m=log2FoldChange,padj_24m=padj)%>%
    select(geneid,log2FC_24m,padj_24m)
  res$padj_24m[is.na(res$padj_24m)==T]=1
  res2 <- DESeq2::results(dds,contrast=c("age","30mo","6mo"),tidy=T,pAdjustMethod="fdr")%>%
    dplyr::rename(geneid=row,log2FC_30m=log2FoldChange,padj_30m=padj)%>%
    select(geneid,log2FC_30m,padj_30m)%>%
    mutate(tissue=fort)
  res2$padj_30m[is.na(res2$padj_30m)==T]=1
  reshe <- left_join(res,res2,by="geneid")
  rm(dds)
  rm(for_RC)
  rm(res)
  rm(res2)
  gc()
  #
  allDE=rbind(allDE,reshe)
  rm(reshe)
  gc()
}
#
data.table::fwrite(allDE,
                   file ="/home/zhzhang/PG/aging_mouse/Alltissue_expgene_OYDEA.aging.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)


```
##### 2.组织水平（each tissue），比较三类基因的DE比例
```r
tissue=c("Aorta","Brain","Heart","Kidney","Liver","Lung","Muscle","Skin")
#导入基因类型id
geneid_class=read.delim("/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt")
#导入DE信息
DE <- read.delim("/home/zhzhang/PG/aging_mouse/Alltissue_expgene_OYDEA.aging.txt",)
#添加基因类型
DE=left_join(DE,geneid_class,by="geneid")%>%
  mutate(de=case_when(padj_24m<0.05 & abs(log2FC_24m)>0.25 ~ 1,
                      padj_30m<0.05 & abs(log2FC_30m)>0.25 ~ 1,
                      T~0))
DE=DE[is.na(DE$type)==F,]

#
#统计组织水平上，3类基因中de基因和表达基因
zongti <- group_by(DE,tissue,type)%>%
  summarise(allnum=n(),denum=sum(de))%>%
  mutate(nodenum=allnum-denum,deratio=denum*100/allnum)%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
zongti$type <- factor(zongti$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
zongti$tissue <- factor(zongti$tissue,levels =rev(c("Brain","Skin","Aorta","Heart","Kidney","Muscle","Liver","Lung")))
zongti <- arrange(zongti,tissue,type)
data.table::fwrite(zongti,
                   file ="/home/zhzhang/PG/aging_mouse/mouse_alltissue.deratio.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#PLOT
pm=ggplot(data =zongti,aes(x=tissue,y=deratio))+
  geom_point(aes(color=type),size=2)+
  geom_line(aes(group=type,color=type))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                     limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  scale_y_continuous(limits = c(0,27))+
  labs(x =NULL, y ="Proportions (%)",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 18),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")+
  coord_flip()
ggsave("/home/zhzhang/PG/aging_mouse/mouse_alltissue.deratio.pdf", 
       pm,width = 5, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#odd ratio plot
#fishertest【pas vs npa】
tizt=zongti
tissuelist=tissue
tiztpv <- data.frame()
for (i in 1:length(tissuelist)) {
  fisherre <- fisher.test(filter(tizt,tissue==stringr::str_to_title(tissuelist[i]) & type!="Protein-coding" & type!="PAA lncRNA")[,4:5])
  pv <- fisherre[["p.value"]]
  or <- fisherre[["estimate"]][["odds ratio"]]
  orlow <- fisherre[["conf.int"]][1]
  orhigh <- fisherre[["conf.int"]][2]
  fpv <- data.frame(tissue=stringr::str_to_title(tissuelist[i]),pvalue=pv,OR=or,
                    ORlow=orlow,ORhigh=orhigh)
  tiztpv <- rbind(tiztpv,fpv)
}
tiztpv <- mutate(tiztpv,pstr=case_when(pvalue<0.001~"***",
                                       pvalue<0.01~"**",
                                       pvalue<0.05~"*",
                                       T~"N.S."))
data.table::fwrite(tiztpv,
                   file ="/home/zhzhang/PG/aging_mouse/mouse_alltissue.deratio.fisherOR.pasVSnpa.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#
tiztpv$tissue <- factor(tiztpv$tissue,levels =rev(c("Brain","Skin","Aorta","Heart","Kidney","Muscle","Liver","Lung")))
ph1 <- ggplot(data =tiztpv,aes(x=tissue,y=OR))+
  geom_hline(yintercept=1,linetype=2,color="#FAA465")+
  geom_errorbar(aes(ymin=ORlow,ymax=ORhigh),width=0.1,position = position_nudge(x=-0.1))+
  geom_point(color="#7197AD",size=5,shape=18,position = position_nudge(x=-0.1))+
  geom_text(aes(label=pstr),position = position_nudge(x=0.1))+
  cowplot::theme_half_open()+
  scale_y_continuous(limits = c(0.15,max(tiztpv$ORhigh+0.2)))+
  labs(x =NULL, y ="Odds ratio",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 18),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")+
  coord_flip()
ggsave("/home/zhzhang/PG/aging_mouse/mouse_alltissue.deratio.fisherOR.pasVSnpa.pdf", 
       ph1,width = 4, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#fishertest【paa vs npa】
tiztpv <- data.frame()
for (i in 1:length(tissuelist)) {
  fisherre <- fisher.test(filter(tizt,tissue==stringr::str_to_title(tissuelist[i]) & type!="Protein-coding" & type!="PAS lncRNA")[,4:5])
  pv <- fisherre[["p.value"]]
  or <- fisherre[["estimate"]][["odds ratio"]]
  orlow <- fisherre[["conf.int"]][1]
  orhigh <- fisherre[["conf.int"]][2]
  fpv <- data.frame(tissue=stringr::str_to_title(tissuelist[i]),pvalue=pv,OR=or,
                    ORlow=orlow,ORhigh=orhigh)
  tiztpv <- rbind(tiztpv,fpv)
}
tiztpv <- mutate(tiztpv,pstr=case_when(pvalue<0.001~"***",
                                       pvalue<0.01~"**",
                                       pvalue<0.05~"*",
                                       T~"N.S."))
data.table::fwrite(tiztpv,
                   file ="/home/zhzhang/PG/aging_mouse/mouse_alltissue.deratio.fisherOR.paaVSnpa.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#
tiztpv$tissue <- factor(tiztpv$tissue,levels =rev(c("Brain","Skin","Aorta","Heart","Kidney","Muscle","Liver","Lung")))
ph1 <- ggplot(data =tiztpv,aes(x=tissue,y=OR))+
  geom_hline(yintercept=1,linetype=2,color="#FAA465")+
  geom_errorbar(aes(ymin=ORlow,ymax=ORhigh),width=0.1,position = position_nudge(x=-0.1))+
  geom_point(color="#628255",size=5,shape=18,position = position_nudge(x=-0.1))+
  geom_text(aes(label=pstr),position = position_nudge(x=0.1))+
  cowplot::theme_half_open()+
  scale_y_continuous(limits = c(0.15,max(tiztpv$ORhigh+0.1)))+
  labs(x =NULL, y ="Odds ratio",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 18),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")+
  coord_flip()
ggsave("/home/zhzhang/PG/aging_mouse/mouse_alltissue.deratio.fisherOR.paaVSnpa.pdf", 
       ph1,width = 4, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#fishertest【pas vs paa】
tiztpv <- data.frame()
for (i in 1:length(tissuelist)) {
  fisherre <- fisher.test(filter(tizt,tissue==stringr::str_to_title(tissuelist[i]) & type!="Protein-coding" & type!="NPA lncRNA")[,4:5])
  pv <- fisherre[["p.value"]]
  or <- fisherre[["estimate"]][["odds ratio"]]
  orlow <- fisherre[["conf.int"]][1]
  orhigh <- fisherre[["conf.int"]][2]
  fpv <- data.frame(tissue=stringr::str_to_title(tissuelist[i]),pvalue=pv,OR=or,
                    ORlow=orlow,ORhigh=orhigh)
  tiztpv <- rbind(tiztpv,fpv)
}
tiztpv <- mutate(tiztpv,pstr=case_when(pvalue<0.001~"***",
                                       pvalue<0.01~"**",
                                       pvalue<0.05~"*",
                                       T~"N.S."))
data.table::fwrite(tiztpv,
                   file ="/home/zhzhang/PG/aging_mouse/mouse_alltissue.deratio.fisherOR.pasVSpaa.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#
tiztpv$tissue <- factor(tiztpv$tissue,levels =rev(c("Brain","Skin","Aorta","Heart","Kidney","Muscle","Liver","Lung")))
ph1 <- ggplot(data =tiztpv,aes(x=tissue,y=OR))+
  geom_hline(yintercept=1,linetype=2,color="#628255")+
  geom_errorbar(aes(ymin=ORlow,ymax=ORhigh),width=0.1,position = position_nudge(x=-0.1))+
  geom_point(color="#7197AD",size=5,shape=18,position = position_nudge(x=-0.1))+
  geom_text(aes(label=pstr),position = position_nudge(x=0.1))+
  cowplot::theme_half_open()+
  scale_y_continuous(limits = c(0.1,max(tiztpv$ORhigh+0.2)),breaks = c(0,5,10) )+
  labs(x =NULL, y ="Odds ratio",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 18),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")+
  coord_flip()
ggsave("/home/zhzhang/PG/aging_mouse/mouse_alltissue.deratio.fisherOR.pasVSpaa.pdf", 
       ph1,width = 4, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 3.DEG中，三类基因中组织共享性agingDEG的比例
```r
#观察DE的三类基因中，tissue-shared衰老DEG比例
tissuelist <- c("Brain","Skin","Aorta","Heart","Kidney","Muscle","Liver","Lung")
#导入基因类型id
geneid_class=read.delim("/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt")
#导入DE信息
DE <- read.delim("/home/zhzhang/PG/aging_mouse/Alltissue_expgene_OYDEA.aging.txt",)
#添加基因类型
DE=left_join(DE,geneid_class,by="geneid")%>%
  mutate(de=case_when(padj_24m<0.05 & abs(log2FC_24m)>0.25 ~ 1,
                      padj_30m<0.05 & abs(log2FC_30m)>0.25 ~ 1,
                      T~0))
DE=DE[is.na(DE$type)==F,]
DE=mutate(DE,type=case_when(type=="Protein-coding" ~ "Protein-coding",
                         type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                         type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                         type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
#统计在至少一个组织中表达的基因，表达的组织数以及差异表达的组织数
tjalltissue <- group_by(DE,geneid,type)%>%
  summarise(expnum=n(),denum=sum(de))%>%
  data.frame()%>%
  mutate(exptype=case_when(expnum>1~"Non-tissue-specific",
                           T~"Tissue-specific"),
         detype=case_when(denum>1~"Tissue-shared",
                          T~"Non-tissue-shared"))
#筛选在至少一个组织中表现出衰老差异表达的基因，并根据denum分配类型
fil_tjalltissue <- filter(tjalltissue,denum>0)
#统计DE的三类基因中，tissue-shared衰老DEG比例
tj11 <- group_by(fil_tjalltissue,type,detype)%>%
  summarise(num=n())%>%
  data.frame()
tj12 <- group_by(fil_tjalltissue,type)%>%
  summarise(allnum=n())%>%
  data.frame()
tj1 <- left_join(tj11,tj12,by="type")%>%
  filter(detype=="Tissue-shared")%>%
  mutate(nonum=allnum-num,ratio=num*100/allnum)
tj1$type <- factor(tj1$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
data.table::fwrite(tj1,
                   file ="/home/zhzhang/PG/aging_mouse/mouse_alltissue.tissueshared_degratio.inallDEG.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot tissuesharedDE基因比例差异
cedif <- ggplot(data = tj1,aes(x=type,y=ratio))+
  geom_col(aes(fill=type),width = 0.5)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  scale_y_continuous(breaks = c(0,20,40,60,80))+
  theme_half_open()+
  labs(x = NULL, y ="Tissue-shared\naging-related DEGs (%)",fill = NULL)+
  theme(axis.title = element_text(size = 17),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 0),
        legend.position = "none")+
  theme(axis.text.x = element_text(angle =30)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  theme(axis.ticks.x = element_line(colour = NA))
ggsave("/home/zhzhang/PG/aging_mouse/mouse_alltissue.tissueshared_degratio.inallDEG.pdf", 
       cedif,width = 4.8,height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#统计de基因中，三类基因的组织特异性表达情况
tj21 <- group_by(fil_tjalltissue,type,exptype)%>%
  summarise(num=n())%>%
  data.frame()
tj2 <- left_join(tj21,tj12,by="type")%>%
  filter(exptype=="Tissue-specific")%>%
  mutate(nonum=allnum-num,ratio=num*100/allnum)
tj2$type <- factor(tj2$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
data.table::fwrite(tj2,
                   file ="/home/zhzhang/PG/aging_mouse/mouse_alltissue.tissuespecificgene_ratio.inallDEG.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot exptissuespecific基因inDEG比例差异
cedif2 <- ggplot(data = tj2,aes(x=type,y=ratio))+
  geom_col(aes(fill=type),width = 0.5)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  scale_y_continuous(breaks = c(0,10,20,30,40))+
  theme_half_open()+
  labs(x = NULL, y ="Tissue-specific genes (%)",fill = NULL)+
  theme(axis.title = element_text(size = 17),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        legend.position = "none")+
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave("/home/zhzhang/PG/aging_mouse/mouse_alltissue.tissuespecificgene_ratio.inallDEG.pdf", 
       cedif2,width = 4, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#统计非组织特异性且DE的基因中，在几个组织中出现aging-related DE
tjzz1 <- group_by(fil_tjalltissue,expnum,type)%>%
  summarise(mediandenum=median(denum),meandenum=mean(denum),num=n())
tjzz1$type <- factor(tjzz1$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
tjzz1 <- arrange(tjzz1,expnum,type)
data.table::fwrite(tjzz1,
                   file ="/home/zhzhang/PG/aging_mouse/mouse_alltissue.detissuenum.allDEG.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#wilcoxtest
tiztpv1 <- data.frame()
for (i in 2:8) {
  fisherre <- wilcox.test(filter(fil_tjalltissue,expnum==i & type=="PAS lncRNA")$denum,
                          filter(fil_tjalltissue,expnum==i & type=="NPA lncRNA")$denum)
  pv <- fisherre[["p.value"]]
  fpv <- data.frame(expnum=i,pvalue=pv)
  tiztpv1 <- rbind(tiztpv1,fpv)
}
tiztpv1 <- mutate(tiztpv1,pstr=case_when(pvalue<0.001~"***",
                                         pvalue<0.01~"**",
                                         pvalue<0.05~"*",
                                         T~"N.S."))
tiztpv2 <- data.frame()
for (i in 2:8) {
  fisherre <- wilcox.test(filter(fil_tjalltissue,expnum==i & type=="PAS lncRNA")$denum,
                          filter(fil_tjalltissue,expnum==i & type=="Protein-coding")$denum)
  pv <- fisherre[["p.value"]]
  fpv <- data.frame(expnum=i,pvalue=pv)
  tiztpv2 <- rbind(tiztpv2,fpv)
}
tiztpv2 <- mutate(tiztpv2,pstr=case_when(pvalue<0.001~"***",
                                         pvalue<0.01~"**",
                                         pvalue<0.05~"*",
                                         T~"N.S."))
tiztpv3 <- data.frame()
for (i in c(2:5,7:8)) {
  fisherre <- wilcox.test(filter(fil_tjalltissue,expnum==i & type=="PAA lncRNA")$denum,
                          filter(fil_tjalltissue,expnum==i & type=="NPA lncRNA")$denum)
  pv <- fisherre[["p.value"]]
  fpv <- data.frame(expnum=i,pvalue=pv)
  tiztpv3 <- rbind(tiztpv3,fpv)
}
tiztpv3 <- mutate(tiztpv3,pstr=case_when(pvalue<0.001~"***",
                                         pvalue<0.01~"**",
                                         pvalue<0.05~"*",
                                         T~"N.S."))
#PLOT
fil_tjalltissue$expnum <- factor(fil_tjalltissue$expnum ,
                                 levels = c(1:8))
fil_tjalltissue$type <- factor(fil_tjalltissue$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
tjzz1$expnum <- factor(tjzz1$expnum ,
                       levels = c(1:8))
ph3 <- ggplot(data =filter(fil_tjalltissue,expnum!=1),aes(x=expnum,y=denum))+
  geom_boxplot(width=0.6,outlier.alpha = 0,aes(color=type))+
  geom_point(data =filter(tjzz1,expnum!=1&type=="Protein-coding"),aes(x=expnum,y=meandenum),shape=23,color="#BB5A5D",position = position_nudge(x=-0.225))+
  geom_point(data =filter(tjzz1,expnum!=1&type=="PAS lncRNA"),aes(x=expnum,y=meandenum),shape=23,color="#7197AD",position = position_nudge(x=-0.075))+
  geom_point(data =filter(tjzz1,expnum!=1&type=="PAA lncRNA"),aes(x=expnum,y=meandenum),shape=23,color="#628255",position = position_nudge(x=0.075))+
  geom_point(data =filter(tjzz1,expnum!=1&type=="NPA lncRNA"),aes(x=expnum,y=meandenum),shape=23,color="#FAA465",position = position_nudge(x=0.225))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                     limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6))+
  coord_cartesian(ylim = c(1,7))+
  ggsignif::geom_signif(annotations=tiztpv1$pstr,textsize=4,
                        y_position=c(2.2,2.2,2.3,2.3,3.2,2.3,6.2),tip_length = 0,size=0.4,
                        xmin = c(0.925:3.925,5,5.925,6.925),xmax = c(1.225:4.225,5.2,6.225,7.225))+
  ggsignif::geom_signif(annotations=tiztpv2$pstr,textsize=4,
                        y_position=c(3,3,3.1,3.1,3.8,6.2,6.5),tip_length = 0,size=0.4,
                        xmin = c(0.775:3.775,4.8,5.775,6.775),xmax = c(0.925:3.925,5,5.925,6.925))+
  ggsignif::geom_signif(annotations=tiztpv3$pstr,textsize=4,
                        y_position=c(1.4,1.5,1.5,3.1,5.2,5.2),tip_length = 0,size=0.4,
                        xmin = c(1.075:4.075,6.075,7.075),xmax = c(1.225:4.225,6.225,7.225))+
  labs(x ="Number of tissues with expression", y ="Number of tissues with\naging-related differential expression",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")
ggsave("/home/zhzhang/PG/aging_mouse/mouse_alltissue.detissuenum.allDEG.pdf", 
       ph3,width = 12, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
##### 4.筛选验证
```r
#mouse
#观察DE的三类基因中，tissue-shared衰老DEG比例
tissuelist <- c("Brain","Skin","Aorta","Heart","Kidney","Muscle","Liver","Lung")
#导入基因类型id
geneid_class=read.delim("/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt")
#导入DE信息
DE <- read.delim("/home/zhzhang/PG/aging_mouse/Alltissue_expgene_OYDEA.aging.txt")
#添加基因类型
DE=left_join(DE,geneid_class,by="geneid")%>%
  mutate(de=case_when(padj_24m<0.05 & abs(log2FC_24m)>0.25 ~ 1,
                      padj_30m<0.05 & abs(log2FC_30m)>0.25 ~ 1,
                      T~0))
DE=DE[is.na(DE$type)==F,]
DE=mutate(DE,type=case_when(type=="Protein-coding" ~ "Protein-coding",
                            type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                            type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                            type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
#统计在至少一个组织中表达的基因，表达的组织数以及差异表达的组织数
tjalltissue <- group_by(DE,geneid,type)%>%
  summarise(expnum=n(),denum=sum(de))%>%
  data.frame()%>%
  mutate(exptype=case_when(expnum>1~"Non-tissue-specific",
                           T~"Tissue-specific"),
         detype=case_when(denum>1~"Tissue-shared",
                          T~"Non-tissue-shared"))
#筛选mouse中DE的PAS,PAA
mousePAS=filter(tjalltissue,type=="PAS lncRNA" & denum>0)
mousePAA=filter(tjalltissue,type=="PAA lncRNA" & denum>0)
top_n(mousePAS,5,denum)
top_n(mousePAA,5,denum)
grep("PB",mousePAS$geneid,value = T)
#human
#观察DE的三类基因中，tissue-shared衰老DEG比例
tissuelist <- c("brain","skin","heart","muscle","bonemarrow","liver","lung","ovary","testis")
alltissue <- data.frame()
for (i in 1:length(tissuelist)) {
  #导入
  OYDEA <- read.delim(paste("/home/zhzhang/PG/aging/",tissuelist[i],"/seurat/Celltype_expgene_OYDEA.txt",sep=""))
  #统计基因在特定的组织水平是否de
  tjzt <- group_by(OYDEA,genename,type)%>%
    summarise(denum=sum(de))%>%
    mutate(de=case_when(denum>=1~1,denum==0~0))%>%
    data.frame()%>%
    select(-denum)
  alltissue <- rbind(alltissue,tjzt)
}
#统计在至少一个组织中表达的基因，表达的组织数以及差异表达的组织数
hstjalltissue <- group_by(alltissue,genename,type)%>%
  summarise(expnum=n(),denum=sum(de))%>%
  data.frame()%>%
  mutate(exptype=case_when(expnum>1~"Non-tissue-specific",
                           T~"Tissue-specific"),
         detype=case_when(denum>1~"Tissue-shared",
                          T~"Non-tissue-shared"))
geneIDname <- read.delim("~/PG/scRNAseq/gene.ID.name.human.txt", header=FALSE)
colnames(geneIDname) <- c("geneid","genename")
geneIDname=mutate(geneIDname,genename=stringr::str_replace_all(genename,"_","-"))
hstjalltissue=left_join(hstjalltissue,geneIDname,by="genename")%>%
  select("geneid" , "type"  ,  "expnum" , "denum"  , "exptype" ,"detype")
#筛选human中DE的PAS,PAA
humanPAS=filter(hstjalltissue,type=="PAS lncRNA" & denum>0)
humanPAA=filter(hstjalltissue,type=="PAA lncRNA" & denum>0)
#同源信息
hsmm=read.delim("/home/zhzhang/PG/aging_mouse/homo/hs.Mus_musculus.11homo_lncRNApair.pro.txt")
colnames(hsmm)=c("human","mouse")
mmhs=read.delim("/home/zhzhang/PG/aging_mouse/homo/mm.Homo_sapiens.11homo_lncRNApair.pro.txt")
colnames(mmhs)=c("mouse","human")
he=rbind(hsmm,select(mmhs,human,mouse))%>%
  distinct(human,mouse)
#DEG同源数量确定
humanPAS=left_join(dplyr::rename(humanPAS,human=geneid),he,by="human")%>%
  left_join(dplyr::rename(mousePAS,mouse=geneid),by="mouse")
humanPAA=left_join(dplyr::rename(humanPAA,human=geneid),he,by="human")%>%
  left_join(dplyr::rename(mousePAA,mouse=geneid),by="mouse")
humanPAS[is.na(humanPAS$denum.y)==F,]
humanPAA[is.na(humanPAA$denum.y)==F,]

#导入DE信息
DE <- read.delim("/home/zhzhang/PG/aging_mouse/Alltissue_expgene_OYDEA.aging.txt")%>%
  filter(geneid %in% qPCRlist$geneid)


qPCRlist <- read.table("~/PG/Z_DF/QPCR/qPCRlist.tsv", quote="\"", comment.char="")%>%
  filter(V1!="ENSMUSG00000121478"&V1!="ENSMUSG00000121485"&V1!="Mou_XLOC_018103"&V1!="Mou_XLOC_037407")
colnames(qPCRlist)="geneid"
#筛选包含更多lnc的组织
group_by(DE,tissue)%>%summarise(num=n())


```
```r
"/home/zhzhang/PG/Z_DF/qPCRlist.tsv"
#基因对应的gtf
grep -w -Ff "/share/home/zhzhang24/PG/Z_DF/QPCR/qPCRlist.tsv" "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" > /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.gtf
#gtf提取基因转录本序列 fasta
gffread /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.gtf -g "/share/home/zhzhang24/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa" -w /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.gtf.fa
#gtf链信息反转
cat /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.gtf|sed -e 's/+/@/g; s/-/+/g; s/@/-/g' > /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.temp.gtf
gffread /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.temp.gtf -g "/share/home/zhzhang24/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa" -w /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.temp.gtf.fa

#提取基因外显子bed（merge并计算共享数）
cat "/share/home/zhzhang24/PG/Z_DF/QPCR/qPCRlist.tsv"|while read gene
do
strand=`grep -w "$gene" "/share/home/zhzhang24/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed"|awk '{print $6}'|uniq`
grep -w "$gene" "/share/home/zhzhang24/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed"|awk '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6}'|bedtools genomecov -i - -g /share/home/zhzhang24/PG/refgenome/Mus_musculus.GRCm39.dna.chrsize.txt -bg|awk -v gene=$gene -v strand=$strand '{print $1"\t"$2"\t"$3"\t"gene"_"NR"_"$4"\texon\t"strand}' >> /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.exon.bed
done
#bed去除与假基因重叠区域，避免同源区域误差
bedtools subtract -a /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.exon.bed -b "/share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg.bed" > /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.exon.rmpg.bed
#转换模板连
cat /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.exon.rmpg.bed|sed -e 's/+/@/g; s/-/+/g; s/@/-/g' > /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.exon.rmpg.temp.bed
bedtools getfasta -fi "/share/home/zhzhang24/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa" -bed /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.exon.rmpg.temp.bed -s -name > /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.exon.rmpg.temp.bed.fa

#提取基因外显子bed（不merge）
cat "/share/home/zhzhang24/PG/Z_DF/QPCR/qPCRlist.tsv"|while read gene
do
grep -w "$gene" "/share/home/zhzhang24/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed"|awk '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6}' >> /share/home/zhzhang24/PG/Z_DF/QPCR/PAlncRNAGene_qPCRlist.refexon.bed
done


```
```r
grep -w -v "Protein-coding" /share/home/zhzhang24/PG/RNAseqdata/newGTF/Homo_sapiens.geneid_class.txt|grep ENSG|grep -v chr|grep -v -F "." > "/share/home/zhzhang24/crisp/lncRNAlist.txt"
awk '{print $1}' "/share/home/zhzhang24/crisp/lncRNAlist.txt" > "/share/home/zhzhang24/crisp/lncRNAlist.tiqu.txt"
#基因对应的gtf
grep -w -Ff "/share/home/zhzhang24/crisp/lncRNAlist.tiqu.txt" "/share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/Homo_sapiens.GRCh38.108.chr.gtf" > /share/home/zhzhang24/crisp/lncRNAGene.gtf
#gtf提取基因转录本序列 fasta
gffread /share/home/zhzhang24/crisp/lncRNAGene.gtf -g "/share/home/zhzhang24/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa" -w /share/home/zhzhang24/crisp/lncRNAGene.gtf.fa


```
##### 5.进行RT-qPCR的基因的表达改变热图in RNA-seq
```r
tissuel=c("Brain","Skin","Aorta","Heart","Kidney","Muscle","Liver","Lung")
genel=c("ENSMUSG00000121500",
        "ENSMUSG00000120992",
        "ENSMUSG00000121378",
        "Mou_XLOC_017668",
        "ENSMUSG00000114608",
        "Mou_XLOC_000586",
        
        "PB.7089",
        
        "Mou_XLOC_028665",
        "Mou_XLOC_013922",
        
        "ENSMUSG00000097039",
        "Mou_XLOC_000120")
#导入差异表达情况
OYDEA <- read.delim("~/PG/aging_mouse/Alltissue_expgene_OYDEA.aging.txt")
#提取11 gene
predata=data.frame(he=paste(rep(genel,times=1,each=8),
                            rep(tissuel,times=8,each=1),
                            sep = "-"))%>%
  left_join(mutate(OYDEA,he=paste(geneid,tissue,sep = "-"))%>%select(-geneid,-tissue),
            by="he")
predata$log2FC_24m[is.na(predata$log2FC_24m)==T]=0
predata$log2FC_30m[is.na(predata$log2FC_30m)==T]=0
predata$padj_24m[is.na(predata$padj_24m)==T]=1
predata$padj_30m[is.na(predata$padj_30m)==T]=1
#提取log2FC
FCdata=select(predata,1,2,4)%>%
  gather(age,FC,2:3)%>%
  mutate(age=case_when(age=="log2FC_24m" ~ "24mo",
                       age=="log2FC_30m" ~ "30mo"))%>%
  separate(he,c("geneid","tissue"),sep = "-")
FCdata$geneid=factor(FCdata$geneid,levels = genel)
FCdata$tissue=factor(FCdata$tissue,levels = tissuel)
FCdata=spread(FCdata,tissue,FC)%>%
  arrange(geneid,age)%>%
  mutate(he=paste(age,geneid,sep="-"))%>%
  select(-age,-geneid)%>%
  column_to_rownames("he")
#提取pvalue
pvdata=select(predata,1,3,5)%>%
  gather(age,pv,2:3)%>%
  mutate(age=case_when(age=="padj_24m" ~ "24mo",
                       age=="padj_30m" ~ "30mo"))%>%
  separate(he,c("geneid","tissue"),sep = "-")
pvdata$geneid=factor(pvdata$geneid,levels = genel)
pvdata$tissue=factor(pvdata$tissue,levels = tissuel)
pvdata=spread(pvdata,tissue,pv)%>%
  arrange(geneid,age)%>%
  mutate(he=paste(age,geneid,sep="-"))%>%
  select(-age,-geneid)%>%
  column_to_rownames("he")
#pvalue矩阵中不显著对应的fc矩阵的位置，FC变为0
FCdata2=FCdata
FCdata2[pvdata>0.05]=0
#pvalue矩阵变为*
pvdata2=pvdata
pvdata2[pvdata>0.05]=""
pvdata2[pvdata<0.05]="*"
pvdata2[pvdata<0.01]="**"
pvdata2[pvdata<0.001]="***"
#top3热图
pFCdata=FCdata[1:12,]
ppvdata=pvdata2[1:12,]
geneln=c("ENSMUSG00000121500",
        "D17H6S56E-5",
        "ENSMUSG00000121378",
        "Mou_XLOC_017668",
        "Gm36161",
        "Mou_XLOC_000586")
rownames(pFCdata)=paste(rep(c("24mo","30mo"),times=6,each=1),rep(geneln,times=1,each=2),sep = "-")
rownames(ppvdata)=paste(rep(c("24mo","30mo"),times=6,each=1),rep(geneln,times=1,each=2),sep = "-")
#样本信息矩阵
samplem=data.frame(row=paste(rep(c("24mo","30mo"),times=6,each=1),rep(geneln,times=1,each=2),sep = "-"),
                   Age=rep(c("24mo vs. 6mo","30mo vs. 6mo"),times=6,each=1))%>%
  column_to_rownames("row")
#plot
bk <- c(seq(min(pFCdata),0,by=0.01),seq(0.1,max(pFCdata),by=0.01))
ph <- pheatmap::pheatmap(pFCdata,
                         scale="none",cluster_cols = F, cluster_rows=F,
                         color = c(colorRampPalette(colors = c("#6082B0","white"))(length(bk)*abs(min(pFCdata))/(abs(min(pFCdata))+max(pFCdata))),colorRampPalette(colors = c("white","#CE6463"))(length(bk)*max(pFCdata)/(abs(min(pFCdata))+max(pFCdata)))),
                         breaks=bk,
                         border_color=NA,angle_col=45,
                         fontsize=15,fontsize_col=12,fontsize_row=12,
                         display_numbers=ppvdata,
                         gaps_row = c(2,4,6,8,10),
                         annotation_row=select(samplem,Age),
                         annotation_colors=list("Age"=c('24mo vs. 6mo'='#FDDF8A','30mo vs. 6mo'='#FCAD60')))
ggsave("/home/zhzhang/PG/Z_DF/QPCR/TOP3.PAlnc.log2FC.pdf", 
       ph,width = 11, height = 8,dpi=1200, units = "in", device='pdf',bg = "transparent")



#4homo热图
pFCdata=FCdata[15:22,]
ppvdata=pvdata2[15:22,]
geneln=c("Mou_XLOC_028665",
         "Mou_XLOC_013922",
         "Pvt1",
         "Mou_XLOC_000120")
rownames(pFCdata)=paste(rep(c("24mo","30mo"),times=4,each=1),rep(geneln,times=1,each=2),sep = "-")
rownames(ppvdata)=paste(rep(c("24mo","30mo"),times=4,each=1),rep(geneln,times=1,each=2),sep = "-")
#样本信息矩阵
samplem=data.frame(row=paste(rep(c("24mo","30mo"),times=4,each=1),rep(geneln,times=1,each=2),sep = "-"),
                   Age=rep(c("24mo vs. 6mo","30mo vs. 6mo"),times=4,each=1))%>%
  column_to_rownames("row")
#plot
bk <- c(seq(min(pFCdata),-0.01,by=0.01),seq(0,max(pFCdata),by=0.01))
ph <- pheatmap::pheatmap(pFCdata,
                         scale="none",cluster_cols = F, cluster_rows=F,
                         color = c(colorRampPalette(colors = c("#6082B0","white"))(length(bk)*abs(min(pFCdata))/(abs(min(pFCdata))+max(pFCdata))),colorRampPalette(colors = c("white","#CE6463"))(length(bk)*max(pFCdata)/(abs(min(pFCdata))+max(pFCdata)))),
                         breaks=bk,
                         border_color=NA,angle_col=45,
                         fontsize=15,fontsize_col=15,fontsize_row=15,
                         display_numbers=ppvdata,
                         gaps_row = c(2,4,6,8),
                         annotation_row=select(samplem,Age),
                         annotation_colors=list("Age"=c('24mo vs. 6mo'='#FDDF8A','30mo vs. 6mo'='#FCAD60')))
ggsave("/home/zhzhang/PG/Z_DF/QPCR/HOMO4.PAlnc.log2FC.pdf", 
       ph,width = 11, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


##### 6.sup table
```r
#小鼠中PAlnc的差异表达
#导入基因类型id
geneid_class=read.delim("/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt")
#导入DE信息
DE <- read.delim("/home/zhzhang/PG/aging_mouse/Alltissue_expgene_OYDEA.aging.txt")
#添加基因类型
DE=left_join(DE,geneid_class,by="geneid")%>%
  mutate(de=case_when(padj_24m<0.05 & abs(log2FC_24m)>0.25 ~ 1,
                      padj_30m<0.05 & abs(log2FC_30m)>0.25 ~ 1,
                      T~0))
DE=DE[is.na(DE$type)==F,]
DE$type=factor(DE$type,levels =c("Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA") )
DE=filter(DE,type=="Pseudogene-associated sense lncRNA" | type=="Pseudogene-associated antisense lncRNA")%>%
  filter(de==1)%>%
  select(-de)%>%
  arrange(type,geneid)%>%
  select(1,7,6,2:5)
colnames(DE)=c("Gene ID","Type","Tissue","Log2FC (24mo vs. 6mo)","Adjusted p-value (24mo vs. 6mo)","Log2FC (30mo vs. 6mo)","Adjusted p-value (30mo vs. 6mo)")
data.table::fwrite(DE,file ="/home/zhzhang/PG/Z_DF/septable/PAlnc_DEG_mouse.tsv",
                   sep = '\t',row.names = F,quote = F,col.names = T)


#人类中PAlnc的差异表达
#观察DE的三类基因中，tissue-shared衰老DEG比例
tissuelist <- c("brain","skin","heart","muscle","bonemarrow","liver","lung","ovary","testis")
alltissue <- data.frame()
for (i in 1:length(tissuelist)) {
  #导入
  OYDEA <- read.delim(paste("/home/zhzhang/PG/aging/",tissuelist[i],"/seurat/Celltype_expgene_OYDEA.txt",sep=""))%>%
    mutate(Tissue=tissuelist[i])
  alltissue <- rbind(alltissue,OYDEA)
}
hstjalltissue=filter(alltissue,de==1)%>%
  filter(type=="PAS lncRNA" | type=="PAA lncRNA")%>%
  mutate(type=case_when(type=="PAS lncRNA" ~ "Pseudogene-associated sense lncRNA",
                        type=="PAA lncRNA" ~ "Pseudogene-associated antisense lncRNA"))
geneIDname <- read.delim("~/PG/scRNAseq/gene.ID.name.human.txt", header=FALSE)
colnames(geneIDname) <- c("geneid","genename")
geneIDname=mutate(geneIDname,genename=stringr::str_replace_all(genename,"_","-"))
hstjalltissue=left_join(hstjalltissue,geneIDname,by="genename")%>%
  select(13,8,12,7,3,6,4,5,11)%>%
  mutate(pct.1=pct.1*100,pct.2=pct.2*100)
rank=distinct(arrange(hstjalltissue,desc(abs(avg_log2FC))),geneid)$geneid
hstjalltissue$geneid=factor(hstjalltissue$geneid,levels =rank )
hstjalltissue$type=factor(hstjalltissue$type,levels =c("Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA") )
hstjalltissue=arrange(hstjalltissue,type,geneid,Tissue)
hstjalltissue$Tissue[hstjalltissue$Tissue=="bonemarrow"]="bone marrow"
colnames(hstjalltissue)=c("Gene ID","Type","Tissue","Cell type","Log2FC (old vs. young)","Adjusted p-value (old vs. young)",
                          "Percentage of cells expressing the gene (old)","Percentage of cells expressing the gene (young)",
                          "Type of differential expression")
data.table::fwrite(hstjalltissue,file ="/home/zhzhang/PG/Z_DF/septable/PAlnc_DEG_human.tsv",
                   sep = '\t',row.names = F,quote = F,col.names = T)


#人类其他总表
lncRNA_class <- read.delim("~/PG/lncRNA_class/new_r1/Homo_sapiens.lncRNA_class.txt")
lncRNAage <- read.delim("~/PG/age/Human.lncRNA.age.txt")
a <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paalnc_lncRNA_exon.merge.phastCons.txt"
paslnc_lncRNA_exon <- rbind(read.delim(a, header=FALSE),
                            read.delim(b, header=FALSE))%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("pgid"),sep="____")%>%
  group_by(pgid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,4)%>%
  dplyr::rename(geneid=pgid)
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
  geneclass_tranlen <- left_join(geneid_class,alltran_len,by="geneid")
}
tlen <- gettlen(a="~/PG/Evolution/tran_len/Homo_sapiens.allgeneexon.bed",
                  b="~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                  c="Human")%>%
  select(1,3)
gettexonnum <- function(a,b,c){
  #导入外显子bed文件
  allgeneexon <- read.delim(a, header=FALSE)%>%
    select(-1,-6)
  colnames(allgeneexon) <- c("start","end","geneid","transcriptid")
  #导入基因分类文件
  geneid_class <- read.delim(b)
  #计算每个基因每个转录本长度
  alltran_exonnum <- mutate(allgeneexon,num=1)%>%
    group_by(geneid,transcriptid)%>%
    summarise(exonnum=sum(num))%>%
    data.frame()%>%
    group_by(geneid)%>%
    summarise(maxexonnum=max(exonnum))%>%
    data.frame()
  #合并分类和转录本长度
  geneclass_tranlen <- left_join(geneid_class,alltran_exonnum,by="geneid")
}
exon <- gettexonnum(a="~/PG/Evolution/tran_len/Homo_sapiens.allgeneexon.bed",
                      b="~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                      c="Human")%>%
  select(1,3)
tgenemaxTPM <- function(a,b,c){
  #导入基因表达TPM矩阵，计算出每个基因在所有样本中的最大表达量
  allsample_TPM <- read.delim(a, row.names=1)
  gene_max_TPM <- apply(allsample_TPM,1,max)%>%
    data.frame()
  colnames(gene_max_TPM) <- "maxTPM"
  gene_max_TPM <- rownames_to_column(gene_max_TPM,"geneid")
  #基因分类和maxtpm合并
  geneid_class <- read.delim(b)%>%
    left_join(gene_max_TPM,by="geneid")
  return(geneid_class)
}
TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                      c = "Human")%>%
  select(1,3)
gene_tissuenum <- read.delim("~/PG/RNAseq/Homo_sapiens/Homo_sapiens.gene_tissuenum.txt")%>%
  select(1,3)
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
  geneid_class <- read.delim(c)
  allgeneexon_m6a_len <- left_join(geneid_class,allgeneexon_m6a_len,by="geneid")%>%
    filter(type!="Interference lncRNA")
  #合并
  return(allgeneexon_m6a_len)
}
m6ad <- geneexonM6Ad(a="/home/zhzhang/PG/m6a/Homo_sapiens.allgeneexon_intersect_m6A.txt",
                       c="/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                       e="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexon.bed")%>%
  select(1,5)
getRBP_NUM <- function(a,b,c){
  #导入基因具有的RBP种类及其结合比例信息文件,计算RBP种类数
  allgene_RBP <- read.delim(a)%>%
    mutate(ratio=1)%>%
    group_by(geneid)%>%
    summarise(RBPnum=sum(ratio))
  #导入基因分类
  geneid_class <- read.delim(c)
  #合并
  allclass <- geneid_class
  geneid_class_RBP <- left_join(allclass,allgene_RBP,by="geneid")%>%
    filter(type!="Interference lncRNA")%>%
    mutate(sp=b)
  geneid_class_RBP$RBPnum[is.na(geneid_class_RBP$RBPnum)==T] <- 0
  return(geneid_class_RBP)
}
RBP_NUM <- getRBP_NUM(a="/home/zhzhang/PG/RBP/Homo_sapiens.allgeneexonAintergenic.RBP.txt",
                         b="Human",
                         c="~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt"
                         )%>%
  select(1,3)
RBPrich <- read.delim("~/PG/RBP/SPhs_3gene_RBPrichness.txt")%>%
  select(1,3)
RBPdiversity <- read.delim("~/PG/RBP/SPhs_3gene_RBPdiversity.txt")%>%
  select(1,4)
he=filter(lncRNA_class,type=="Pseudogene-associated sense lncRNA" | type=="Pseudogene-associated antisense lncRNA")%>%
  left_join(lncRNAage,by="geneid")%>%
  left_join(paslnc_lncRNA_exon,by="geneid")%>%
  left_join(tlen,by="geneid")%>%
  left_join(exon,by="geneid")%>%
  left_join(TPM,by="geneid")%>%
  left_join(gene_tissuenum,by="geneid")%>%
  left_join(m6ad,by="geneid")%>%
  left_join(RBPrich,by="geneid")%>%
  left_join(RBPdiversity,by="geneid")
he$type=factor(he$type,levels = c("Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA"))
he1=arrange(he,type,geneid)
colnames(he1)=c("Gene ID","Type","Evolutionary age","PhastCons score",
                "Maximal transcript length (nt)","Maximal exon number",
                "Maximal expression level (TPM)","Number of tissues expressing the gene",
                "Density of m6A sites (pcs/kb)","RBP richness","RBP diversity (Shannon-Wiener index)")
data.table::fwrite(he1,file ="/home/zhzhang/PG/Z_DF/septable/PAlnc_zong_human.tsv",
                   sep = '\t',row.names = F,quote = F,col.names = T)



#mouse
lncRNA_class <- read.delim("~/PG/lncRNA_class/new_r1/Mus_musculus.lncRNA_class.txt")
lncRNAage <- read.delim("~/PG/age/Mouse.lncRNA.age.txt")
a <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paalnc_lncRNA_exon.merge.phastCons.txt"
paslnc_lncRNA_exon <- rbind(read.delim(a, header=FALSE),
                            read.delim(b, header=FALSE))%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("pgid"),sep="____")%>%
  group_by(pgid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,4)%>%
  dplyr::rename(geneid=pgid)
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
  geneclass_tranlen <- left_join(geneid_class,alltran_len,by="geneid")
}
tlen <- gettlen("~/PG/Evolution/tran_len/Mus_musculus.allgeneexon.bed",
                "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                "Mouse")%>%
  select(1,3)
gettexonnum <- function(a,b,c){
  #导入外显子bed文件
  allgeneexon <- read.delim(a, header=FALSE)%>%
    select(-1,-6)
  colnames(allgeneexon) <- c("start","end","geneid","transcriptid")
  #导入基因分类文件
  geneid_class <- read.delim(b)
  #计算每个基因每个转录本长度
  alltran_exonnum <- mutate(allgeneexon,num=1)%>%
    group_by(geneid,transcriptid)%>%
    summarise(exonnum=sum(num))%>%
    data.frame()%>%
    group_by(geneid)%>%
    summarise(maxexonnum=max(exonnum))%>%
    data.frame()
  #合并分类和转录本长度
  geneclass_tranlen <- left_join(geneid_class,alltran_exonnum,by="geneid")
}
exon <- gettexonnum("~/PG/Evolution/tran_len/Mus_musculus.allgeneexon.bed",
                    "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                    "Mouse")%>%
  select(1,3)
tgenemaxTPM <- function(a,b,c){
  #导入基因表达TPM矩阵，计算出每个基因在所有样本中的最大表达量
  allsample_TPM <- read.delim(a, row.names=1)
  gene_max_TPM <- apply(allsample_TPM,1,max)%>%
    data.frame()
  colnames(gene_max_TPM) <- "maxTPM"
  gene_max_TPM <- rownames_to_column(gene_max_TPM,"geneid")
  #基因分类和maxtpm合并
  geneid_class <- read.delim(b)%>%
    left_join(gene_max_TPM,by="geneid")
  return(geneid_class)
}
TPM <- tgenemaxTPM(a = "/home/zhzhang/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                   b = "/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                   c = "Mouse")%>%
  select(1,3)
gene_tissuenum <- read.delim("~/PG/RNAseq/Mus_musculus/Mus_musculus.gene_type_tissueexp.txt")%>%
  group_by(geneid,type)%>%
  summarise(tnum=sum(exp))%>%
  data.frame()%>%select(1,3)
he=filter(lncRNA_class,type=="Pseudogene-associated sense lncRNA" | type=="Pseudogene-associated antisense lncRNA")%>%
  left_join(lncRNAage,by="geneid")%>%
  left_join(paslnc_lncRNA_exon,by="geneid")%>%
  left_join(tlen,by="geneid")%>%
  left_join(exon,by="geneid")%>%
  left_join(TPM,by="geneid")%>%
  left_join(gene_tissuenum,by="geneid")
he$type=factor(he$type,levels = c("Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA"))
he1=arrange(he,type,geneid)
colnames(he1)=c("Gene ID","Type","Evolutionary age","PhastCons score",
                "Maximal transcript length (nt)","Maximal exon number",
                "Maximal expression level (TPM)","Number of tissues expressing the gene"
)
data.table::fwrite(he1,file ="/home/zhzhang/PG/Z_DF/septable/PAlnc_zong_mouse.tsv",
                   sep = '\t',row.names = F,quote = F,col.names = T)




```





