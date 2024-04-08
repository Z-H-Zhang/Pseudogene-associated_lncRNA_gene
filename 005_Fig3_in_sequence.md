### FIG3
### in sequence
##### 1.两类lncRNA外显子序列保守性对比
```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的保守性
###lncrna外显子 保守性打分计算
#人类
#提取pglnc外显子bed,重叠部分merge
cat /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.lncRNA_class.txt|awk -F "\t" '$2=="Pseudogene-derived lncRNA"{print $1}' > /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.pgdlnc.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.all_lncRNA_exon.bed|grep -w -Ff /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.pgdlnc.geneid.txt |awk '{print "chr"$1"\t"$2"\t"$3}'|bedtools sort|bedtools merge|awk '{print $0"\tpglnc_"NR}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.pgdlnc_lncRNA_mergeexon.ucsc.bed
#提取npglnc外显子bed，重叠部分merge
cat /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.lncRNA_class.txt|awk -F "\t" '$2=="Non-pseudogene-derived lncRNA"{print $1}' > /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.npgdlnc.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.all_lncRNA_exon.bed|grep -w -Ff /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.npgdlnc.geneid.txt |awk '{print "chr"$1"\t"$2"\t"$3}'|bedtools sort|bedtools merge|awk '{print $0"\tnpglnc_"NR}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.npgdlnc_lncRNA_mergeexon.ucsc.bed
#分别保守性打分
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/hg38.phastCons100way.bw /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.pgdlnc_lncRNA_mergeexon.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.pgdlnc_lncRNA_exon.phastCons.txt
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/hg38.phastCons100way.bw /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.npgdlnc_lncRNA_mergeexon.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.npgdlnc_lncRNA_exon.phastCons.txt
#小鼠
#提取pglnc外显子bed，重叠部分merge
cat /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.lncRNA_class.txt|awk -F "\t" '$2=="Pseudogene-derived lncRNA"{print $1}' > /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.pgdlnc.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed|grep -w -Ff /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.pgdlnc.geneid.txt |awk '{print "chr"$1"\t"$2"\t"$3}'|bedtools sort|bedtools merge|awk '{print $0"\tpglnc_"NR}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.pgdlnc_lncRNA_mergeexon.ucsc.bed
#提取npglnc外显子bed，重叠部分merge
cat /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.lncRNA_class.txt|awk -F "\t" '$2=="Non-pseudogene-derived lncRNA"{print $1}' > /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.npgdlnc.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed|grep -w -Ff /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.npgdlnc.geneid.txt |awk '{print "chr"$1"\t"$2"\t"$3}'|bedtools sort|bedtools merge|awk '{print $0"\tnpglnc_"NR}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.npgdlnc_lncRNA_mergeexon.ucsc.bed
#分别保守性打分
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/mm39.phastCons35way.bw /home/zhzhang/PG/Evolution/conserve/Mus_musculus.pgdlnc_lncRNA_mergeexon.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.pgdlnc_lncRNA_exon.phastCons.txt
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/mm39.phastCons35way.bw /home/zhzhang/PG/Evolution/conserve/Mus_musculus.npgdlnc_lncRNA_mergeexon.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.npgdlnc_lncRNA_exon.phastCons.txt



###基因间区 保守性打分计算
##获取基因间区bed文件
#获取全部基因bed文件
#人类
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.pre_allgene.bed
#小鼠
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.pre_allgene.bed
#人
a <- "~/PG/Evolution/conserve/Homo_sapiens.pre_allgene.bed"
b <- "~/PG/Evolution/conserve/Homo_sapiens.allgene.bed"
#小鼠
a <- "~/PG/Evolution/conserve/Mus_musculus.pre_allgene.bed"
b <- "~/PG/Evolution/conserve/Mus_musculus.allgene.bed"
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
#获取全部基因启动子bed文件
cat "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgene.bed"|awk '$6=="+" {print $1"\t"$2-2000"\t"$2"\t"$4"\tpromoter\t"$6} $6=="-" {print $1"\t"$3"\t"$3+2000"\t"$4"\tpromoter\t"$6}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgenepromoter.bed
cat "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgene.bed"|awk '$6=="+" {print $1"\t"$2-2000"\t"$2"\t"$4"\tpromoter\t"$6} $6=="-" {print $1"\t"$3"\t"$3+2000"\t"$4"\tpromoter\t"$6}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgenepromoter.bed
#合并基因区和启动子区域bed文件
cp /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgene.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgeneApromoter.bed
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgenepromoter.bed >> /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgeneApromoter.bed
cp /home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgene.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgeneApromoter.bed
cat /home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgenepromoter.bed >> /home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgeneApromoter.bed
bedtools sort -i /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgeneApromoter.bed|bedtools merge > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgeneApromoter.sort.bed
bedtools sort -i /home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgeneApromoter.bed|bedtools merge > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgeneApromoter.sort.bed
rm /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgeneApromoter.bed
rm /home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgeneApromoter.bed
#补集获得基因间区bed文件
cat "/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa.fai" | awk '{print $1"\t"$2}' > /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt
cat "/home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa.fai" | awk '{print $1"\t"$2}' > /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chrsize.txt
bedtools complement -i /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgeneApromoter.sort.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.intergenic.bed
bedtools complement -i /home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgeneApromoter.sort.bed -g /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chrsize.txt > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.intergenic.bed
#生成随机的20000个3kb区域bed文件
bedtools random -n 20000 -l 3000 -seed 1024 -g "/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt" > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_region.bed
bedtools random -n 20000 -l 3000 -seed 1024 -g "/home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chrsize.txt" > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.20000random_3kb_region.bed
##生成随机的20000个3kb基因间区bed文件,并打分计算
#人类
bedtools shuffle -i /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_region.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -incl /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.intergenic.bed -noOverlapping|bedtools sort|awk '{print $1"\t"$2"\t"$3"\tintergenic_"$4"\tintergenic\t"$6}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed|awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.ucsc.bed
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/hg38.phastCons100way.bw /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.phastCons.txt
#小鼠
bedtools shuffle -i /home/zhzhang/PG/Evolution/conserve/Mus_musculus.20000random_3kb_region.bed -g /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chrsize.txt -incl /home/zhzhang/PG/Evolution/conserve/Mus_musculus.intergenic.bed -noOverlapping|bedtools sort|awk '{print $1"\t"$2"\t"$3"\tintergenic_"$4"\tintergenic\t"$6}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.20000random_3kb_intergenic.bed
cat /home/zhzhang/PG/Evolution/conserve/Mus_musculus.20000random_3kb_intergenic.bed|awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.20000random_3kb_intergenic.ucsc.bed
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/mm39.phastCons35way.bw /home/zhzhang/PG/Evolution/conserve/Mus_musculus.20000random_3kb_intergenic.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.20000random_3kb_intergenic.phastCons.txt




###蛋白编码基因外显子（CDS,5UTR,3UTR）的保守性打分
##获取蛋白编码基因外显子CDS的bed（CDS所属基因类型在22/24列）,重叠部分merge
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="CDS" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="CDS" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.bed
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.bed|bedtools sort|bedtools merge|awk '{print "chr"$1"\t"$2"\t"$3"\tCDS_"NR}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.ucsc.bed
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="CDS" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="CDS" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.bed
cat /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.bed|bedtools sort|bedtools merge|awk '{print "chr"$1"\t"$2"\t"$3"\tCDS_"NR}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.ucsc.bed
#保守性打分
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/hg38.phastCons100way.bw /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.phastCons.txt
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/mm39.phastCons35way.bw /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.phastCons.txt



##获取蛋白编码基因外显子5UTR的bed（5UTR所属基因类型在20/22/24列）,重叠部分merge
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="five_prime_utr" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="five_prime_utr" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="five_prime_utr" && $20=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.bed
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.bed|bedtools sort|bedtools merge|awk '{print "chr"$1"\t"$2"\t"$3"\t5UTR_"NR}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.ucsc.bed
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="five_prime_utr" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="five_prime_utr" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="five_prime_utr" && $20=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.bed
cat /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.bed|bedtools sort|bedtools merge|awk '{print "chr"$1"\t"$2"\t"$3"\t5UTR_"NR}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.ucsc.bed
#保守性打分
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/hg38.phastCons100way.bw /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.phastCons.txt
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/mm39.phastCons35way.bw /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.phastCons.txt



##获取蛋白编码基因外显子3UTR的bed（3UTR所属基因类型在20/22/24列）,重叠部分merge
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="three_prime_utr" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="three_prime_utr" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="three_prime_utr" && $20=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.bed
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.bed|bedtools sort|bedtools merge|awk '{print "chr"$1"\t"$2"\t"$3"\t3UTR_"NR}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.ucsc.bed
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="three_prime_utr" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="three_prime_utr" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="three_prime_utr" && $20=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.bed
cat /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.bed|bedtools sort|bedtools merge|awk '{print "chr"$1"\t"$2"\t"$3"\t3UTR_"NR}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.ucsc.bed
#保守性打分
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/hg38.phastCons100way.bw /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.phastCons.txt
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/mm39.phastCons35way.bw /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.phastCons.txt

#传至实验室服务器/home/zhzhang/PG/Evolution/conserve
```


```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的保守性
#人类
a <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.pgdlnc_lncRNA_exon.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.npgdlnc_lncRNA_exon.phastCons.txt"
c <- "~/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.phastCons.txt"
d <- "~/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.phastCons.txt"
e <- "~/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.phastCons.txt"
f <- "~/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.lncRNA_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.lncRNA_PhastCons.png"
#小鼠
a <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.pgdlnc_lncRNA_exon.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.npgdlnc_lncRNA_exon.phastCons.txt"
c <- "~/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.phastCons.txt"
d <- "~/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.phastCons.txt"
e <- "~/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.phastCons.txt"
f <- "~/PG/Evolution/conserve/Mus_musculus.20000random_3kb_intergenic.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.lncRNA_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.lncRNA_PhastCons.png"
#导入pgdlnc外显子的保守性打分
pgdlnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  select(6)%>%
  mutate(type="Pseudogene-derived lncRNA")
colnames(pgdlnc_lncRNA_exon)[1] <- "score"
#导入npgdlnc外显子的保守性打分
npgdlnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  select(6)%>%
  mutate(type="Non-pseudogene-derived lncRNA")
colnames(npgdlnc_lncRNA_exon)[1] <- "score"
#导入蛋白编码基因CDS,5UTR,3UTR保守性打分
proteingeneCDS <- read.delim(c, header=FALSE)%>%
  select(6)%>%
  mutate(type="CDS")
colnames(proteingeneCDS) <- colnames(pgdlnc_lncRNA_exon)
proteingene5UTR <- read.delim(d, header=FALSE)%>%
  select(6)%>%
  mutate(type="5'UTR")
colnames(proteingene5UTR) <- colnames(pgdlnc_lncRNA_exon)
proteingene3UTR <- read.delim(e, header=FALSE)%>%
  select(6)%>%
  mutate(type="3'UTR")
colnames(proteingene3UTR) <- colnames(pgdlnc_lncRNA_exon)
#导入基因间区保守性打分
intergenic <- read.delim(f, header=FALSE)%>%
  select(6)%>%
  mutate(type="Random intergenic")
colnames(intergenic) <- colnames(pgdlnc_lncRNA_exon)
#合并
alldata <- rbind(pgdlnc_lncRNA_exon,npgdlnc_lncRNA_exon,proteingeneCDS,proteingene5UTR,proteingene3UTR,intergenic)
alldata$type <- factor(alldata$type,levels = c("Random intergenic","Non-pseudogene-derived lncRNA",
                                       "Pseudogene-derived lncRNA","CDS","5'UTR","3'UTR"))
#统计
tj <- group_by(alldata,type)%>%
  summarise(mean=mean(score),median=median(score))
data.table::fwrite(tj,file =g,sep = '\t',row.names = F,quote = F,col.names = T)
#plot
p1 <- ggplot(data=alldata, aes(x=type,y=score))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  geom_signif(map_signif_level=T,y_position=c(0.4,1,1.1,1.2),
              comparisons = list(c("Random intergenic","Non-pseudogene-derived lncRNA"),
                                 c("Non-pseudogene-derived lncRNA","Pseudogene-derived lncRNA"),
                                 c("Pseudogene-derived lncRNA","CDS"),
                                 c("Pseudogene-derived lncRNA","5'UTR")))+
  scale_fill_manual(values=c("#8491B4","#FAA465","#7197AD","#BB5A5D","#BB5A5D","#BB5A5D"),
                    limits=c("Random intergenic","Non-pseudogene-derived lncRNA",
                             "Pseudogene-derived lncRNA","CDS","5'UTR","3'UTR"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  scale_x_discrete(labels = c("Random\nintergenic","Non-pseudogene-\nderived lncRNA",
                              "Pseudogene-\nderived lncRNA","CDS","5'UTR","3'UTR"))+
  labs(y = "PhastCons score", x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave(h,p1,width = 6, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")



```





##### 3.两类lncRNA基因转录本长度
```r
#提取外显子bed，为求转录本长度做准备【chr start end geneid transcriptid strand】
#人类
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" |grep -w exon|awk '$13=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $11=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/Evolution/tran_len/Homo_sapiens.allgeneexon.bed
#小鼠
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" |grep -w exon|awk '$13=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $11=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/Evolution/tran_len/Mus_musculus.allgeneexon.bed
#鸡
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Gallus_gallus.GRCg7b.108.chr.rmpg.novellncRNA.gtf" |grep -w exon|awk '$13=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $11=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/Evolution/tran_len/Gallus_gallus.allgeneexon.bed


```
```r
#对比两类lnc转录本长度
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
    summarise(len=sum(len))
  #合并分类和转录本长度
  geneclass_tranlen <- left_join(geneid_class,alltran_len,by="geneid")%>%
    filter(type!="Interference lncRNA")%>%
    mutate(sp=c)
}
#人
hstlen <- gettlen("~/PG/Evolution/tran_len/Homo_sapiens.allgeneexon.bed",
                  "~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                  "Human")
#小鼠
mmtlen <- gettlen("~/PG/Evolution/tran_len/Mus_musculus.allgeneexon.bed",
                  "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                  "Mouse")
#合并不同物种信息
allsp_tlen <- rbind(hstlen,mmtlen)
allsp_tlen$type <- factor(allsp_tlen$type,
                          levels = c("Protein-coding",
                                     "Pseudogene-derived lncRNA",
                                     "Non-pseudogene-derived lncRNA"))
#统计储存
tj <- group_by(allsp_tlen,sp,type)%>%
  summarise(mean=mean(len),median=median(len))
data.table::fwrite(tj,file ="/home/zhzhang/PG/Evolution/tran_len/ALLSP_3gene_transcript_len.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
pmm <- ggplot(data = filter(allsp_tlen,sp=="Mouse"),aes(x=type,y=log10(len)))+
  geom_violin(width=0.9,aes(fill=type,color=type))+
  geom_boxplot(fatten = 3,width=0.1,outlier.alpha = 0)+
  geom_signif(map_signif_level=T,y_position = c(5.5,5),
              comparisons=list(c("Protein-coding","Pseudogene-derived lncRNA"),
                               c("Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0.5, 6))+
  scale_y_continuous(breaks = c(1,2,3,4,5),labels = c("1","2","3","4","5"))+
  scale_x_discrete(labels = c("Protein-coding","Pseudogene-derived\nlncRNA","Non-pseudogene-\nderived lncRNA"))+
  labs(x = NULL, y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[10]*"("*"b"*"p"*")"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =35))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/Evolution/tran_len/Mus_musculus.lnc_tranlen.png", 
       pmm,width = 4, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")
phs <- ggplot(data = filter(allsp_tlen,sp=="Human"),aes(x=type,y=log10(len)))+
  geom_violin(width=0.9,aes(fill=type,color=type))+
  geom_boxplot(fatten = 3,width=0.1,outlier.alpha = 0)+
  geom_signif(map_signif_level=T,y_position = c(5.5,5),
              comparisons=list(c("Protein-coding","Pseudogene-derived lncRNA"),
                               c("Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  coord_cartesian(ylim = c(1, 6))+
  scale_y_continuous(breaks = c(1,2,3,4,5),labels = c("1","2","3","4","5"))+
  scale_x_discrete(labels = c("Protein-coding","Pseudogene-derived\nlncRNA","Non-pseudogene-\nderived lncRNA"))+
  labs(x = NULL, y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[10]*"("*"b"*"p"*")"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =35))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/Evolution/tran_len/Homo_sapiens.lnc_tranlen.png", 
       phs,width = 4, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")


```
##### 4.两类lncRNA基因转录本外显子数量
```r
#对比两类lnc转录本外显子数量
#函数计算三类基因每个转录本的外显子数量
#a输入全部基因外显子bed，b输入基因分类文件,c输入物种名
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
    summarise(exonnum=sum(num))
  #合并分类和转录本长度
  geneclass_tranlen <- left_join(geneid_class,alltran_exonnum,by="geneid")%>%
    filter(type!="Interference lncRNA")%>%
    mutate(sp=c)
}
#人
hstlen <- gettexonnum(a="~/PG/Evolution/tran_len/Homo_sapiens.allgeneexon.bed",
                      b="~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                      c="Human")
#小鼠
mmtlen <- gettexonnum("~/PG/Evolution/tran_len/Mus_musculus.allgeneexon.bed",
                      "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                      "Mouse")
#合并不同物种信息
allsp_ten <- rbind(hstlen,mmtlen)
allsp_ten$type <- factor(allsp_ten$type,
                         levels = c("Protein-coding",
                                    "Pseudogene-derived lncRNA",
                                    "Non-pseudogene-derived lncRNA"))
#统计储存
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(1024)
tj <- group_by(allsp_ten,sp,type)%>%
  summarise(mean=mean(exonnum),median=median(exonnum),
            confmin=boot::boot.ci(boot::boot(exonnum, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(exonnum, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][5])
data.table::fwrite(tj,file ="/home/zhzhang/PG/Evolution/tran_en/ALLSP_3gene_transcript_exonnum.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
tj$type <- factor(tj$type,levels = c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))
t.test(filter(allsp_ten,sp=="Mouse" & type=="Protein-coding")$exonnum,
       filter(allsp_ten,sp=="Mouse" & type=="Pseudogene-derived lncRNA")$exonnum)
t.test(filter(allsp_ten,sp=="Mouse" & type=="Pseudogene-derived lncRNA")$exonnum,
       filter(allsp_ten,sp=="Mouse" & type=="Non-pseudogene-derived lncRNA")$exonnum)
pmm <- ggplot(data = filter(tj,sp=="Mouse"),aes(x=type,y=mean))+
  geom_point(size=2,aes(fill=type,color=type))+
  geom_errorbar(aes(ymin=confmin,ymax=confmax,color=type),width=0.15)+
  geom_signif(annotations=c("***","***"),y_position=c(7.8,4.5),
              xmin = c(1,2),xmax = c(2,3),tip_length = 0.01)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  scale_y_continuous(limits = c(2.5,8),breaks = c(0,2.5,5,7.5))+
  scale_x_discrete(labels = c("Protein-coding","Pseudogene-derived\nlncRNA","Non-pseudogene-\nderived lncRNA"))+
  labs(x = NULL, y ="Transcript exon number",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =35))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/Evolution/tran_en/Mus_musculus.lnc_tranexonn.png", 
       pmm,width = 4, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")
t.test(filter(allsp_ten,sp=="Human" & type=="Protein-coding")$exonnum,
       filter(allsp_ten,sp=="Human" & type=="Pseudogene-derived lncRNA")$exonnum)
t.test(filter(allsp_ten,sp=="Human" & type=="Pseudogene-derived lncRNA")$exonnum,
       filter(allsp_ten,sp=="Human" & type=="Non-pseudogene-derived lncRNA")$exonnum)
phs <- ggplot(data = filter(tj,sp=="Human"),aes(x=type,y=mean))+
  geom_point(size=2,aes(fill=type,color=type))+
  geom_errorbar(aes(ymin=confmin,ymax=confmax,color=type),width=0.15)+
  geom_signif(annotations=c("***","***"),y_position=c(8.5,5.7),
              xmin = c(1,2),xmax = c(2,3),tip_length = 0.01)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  scale_y_continuous(limits = c(2.5,9),breaks = c(0,3,6,9))+
  scale_x_discrete(labels = c("Protein-coding","Pseudogene-derived\nlncRNA","Non-pseudogene-\nderived lncRNA"))+
  labs(x = NULL, y ="Transcript exon number",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =35))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/Evolution/tran_en/Homo_sapiens.lnc_tranexonn.png", 
       phs,width = 4, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")




```


##### 5.两类lncRNA基因外显子 SNP密度
```r
#dsSNP数据库下载数据
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz
gzip -d "/home/zhzhang/PG/Evolution/snp/00-common_all.vcf.gz"
#1-based vcf转为0-based bed
tail -n +58 "/home/zhzhang/PG/Evolution/snp/00-common_all.vcf"|awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > /home/zhzhang/PG/Evolution/snp/common_allsnp.bed


wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20230702.vcf.gz
gzip -d "/home/zhzhang/PG/Evolution/snp/clinvar_20230702.vcf.gz"
tail -n +28 "/home/zhzhang/PG/Evolution/snp/clinvar_20230702.vcf"|awk '{print $1"\t"$2"\t"$2"\t"$3}' > /home/zhzhang/PG/Evolution/snp/clinvar_20230702snp.bed


```
```r
#pgdlnc_lncRNA_merge exon bed文件
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.pgdlnc_lncRNA_mergeexon.ucsc.bed|sed "s/^chr//g" > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.pgdlnc_lncRNA_mergeexon.bed
#npgdlnc_lncRNA_merge exon bed文件
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.npgdlnc_lncRNA_mergeexon.ucsc.bed|sed "s/^chr//g" > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.npgdlnc_lncRNA_mergeexon.bed
#基因间区 bed文件
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed|awk '{print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.bed
#proteingeneCDS bed文件
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.ucsc.bed|sed "s/^chr//g" > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.bed
#proteingene5UTR bed文件
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.ucsc.bed|sed "s/^chr//g" > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.bed
#proteingene3UTR bed文件
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.ucsc.bed|sed "s/^chr//g" > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.bed


#六种区域去除cpg岛区域
bedtools subtract -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.pgdlnc_lncRNA_mergeexon.bed -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.pgdlnc_lncRNA_mergeexon.nocpg.bed
bedtools subtract -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.npgdlnc_lncRNA_mergeexon.bed -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.npgdlnc_lncRNA_mergeexon.nocpg.bed
bedtools subtract -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.bed -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.nocpg.bed
bedtools subtract -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.bed -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.nocpg.bed
bedtools subtract -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.bed -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.nocpg.bed
bedtools subtract -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.bed -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.nocpg.bed




#获取六种区域中，SNP的数量
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.pgdlnc_lncRNA_mergeexon.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.pgdlnc_lncRNA_mergeexon.SNPnum.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.npgdlnc_lncRNA_mergeexon.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.npgdlnc_lncRNA_mergeexon.SNPnum.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.SNPnum.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.SNPnum.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.SNPnum.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.SNPnum.txt
#传至实验室服务器/home/zhzhang/PG/Evolution/snp


```
```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的SNP密度
#人类commonSNP
a <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.pgdlnc_lncRNA_mergeexon.SNPnum.txt"
b <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.npgdlnc_lncRNA_mergeexon.SNPnum.txt"
c <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.SNPnum.txt"
d <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.SNPnum.txt"
e <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.SNPnum.txt"
f <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.SNPnum.txt"
g <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.lncRNA_SNPdensity.tj.txt"
h <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.lncRNA_SNPdensity.png"

#导入pgdlnc外显子的SNP数量
pgdlnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="Pseudogene-derived lncRNA")
colnames(pgdlnc_lncRNA_exon)[c(1,2)] <- c("SNPnum","SNPdensity")
#导入npgdlnc外显子的SNP数量
npgdlnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="Non-pseudogene-derived lncRNA")
colnames(npgdlnc_lncRNA_exon)[c(1,2)] <- c("SNPnum","SNPdensity")
#导入蛋白编码基因CDS,5UTR,3UTR的SNP数量
proteingeneCDS <- read.delim(c, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="CDS")
colnames(proteingeneCDS) <- colnames(pgdlnc_lncRNA_exon)
proteingene5UTR <- read.delim(d, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="5'UTR")
colnames(proteingene5UTR) <- colnames(pgdlnc_lncRNA_exon)
proteingene3UTR <- read.delim(e, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="3'UTR")
colnames(proteingene3UTR) <- colnames(pgdlnc_lncRNA_exon)
#导入基因间区的SNP数量
intergenic <- read.delim(f, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="Random intergenic")
colnames(intergenic) <- colnames(pgdlnc_lncRNA_exon)
#合并
alldata <- rbind(pgdlnc_lncRNA_exon,npgdlnc_lncRNA_exon,proteingeneCDS,proteingene5UTR,proteingene3UTR,intergenic)
alldata$type <- factor(alldata$type,levels = c("Random intergenic","Non-pseudogene-derived lncRNA",
                                               "Pseudogene-derived lncRNA","CDS","5'UTR","3'UTR"))
#统计
tj <- group_by(alldata,type)%>%
  summarise(mean=mean(SNPdensity),median=median(SNPdensity),
            meannum=mean(SNPnum),mediannum=median(SNPnum))
data.table::fwrite(tj,file =g,sep = '\t',row.names = F,quote = F,col.names = T)
#plot
p1 <- ggplot(data=alldata, aes(x=type,y=log10(SNPdensity+1)))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  geom_signif(map_signif_level=T,y_position=c(2,2.4,2.8,3.2),
              comparisons = list(c("Random intergenic","Non-pseudogene-derived lncRNA"),
                                 c("Non-pseudogene-derived lncRNA","Pseudogene-derived lncRNA"),
                                 c("Pseudogene-derived lncRNA","CDS"),
                                 c("Pseudogene-derived lncRNA","3'UTR")))+
  scale_fill_manual(values=c("#8491B4","#FAA465","#7197AD","#BB5A5D","#BB5A5D","#BB5A5D"),
                    limits=c("Random intergenic","Non-pseudogene-derived lncRNA",
                             "Pseudogene-derived lncRNA","CDS","5'UTR","3'UTR"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,1,2,3),labels = c("0","1","2","3"))+
  scale_x_discrete(labels = c("Random\nintergenic","Non-pseudogene-\nderived lncRNA",
                              "Pseudogene-\nderived lncRNA","CDS","5'UTR","3'UTR"))+
  labs(y =expression("l"*"o"*"g"[10]*"("*"S"*"N"*"P"~"d"*"e"*"n"*"s"*"i"*"t"*"y"*")"),
       x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave(h,p1,width = 6, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")




```




##### 7.CADD评分 SNP
```r
#CADD数据库下载SNP的CADD分数https://cadd.gs.washington.edu/download
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.snv.tsv.gz -P /home/zhzhang/PG/Evolution/CADD/
#产生positionID（chr_position）
tail -n +3 "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.r3.0.snv.tsv"|awk '{print $0"\t"$1"_"$2}' > /home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt

#产生common SNP的 SNPID和positionID对照表
cat "/home/zhzhang/PG/Evolution/snp/common_allsnp.bed"|awk '{print $4"\t"$1"_"$3}' > /home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt

#获取六种区域中，SNP的ID
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.pgdlnc_lncRNA_mergeexon.bed -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.pgdlnc_lncRNA_mergeexon.common_SNPID.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.npgdlnc_lncRNA_mergeexon.bed -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.npgdlnc_lncRNA_mergeexon.common_SNPID.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.bed -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.20000random_3kb_intergenic.common_SNPID.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.bed -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingeneCDS.common_SNPID.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.bed -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingene5UTR.common_SNPID.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.bed -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingene3UTR.common_SNPID.txt


#根据六种区域内SNPID提取对应positionID
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.pgdlnc_lncRNA_mergeexon.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.pgdlnc_lncRNA_mergeexon.common_SNP.positionID.txt
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.npgdlnc_lncRNA_mergeexon.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.npgdlnc_lncRNA_mergeexon.common_SNP.positionID.txt
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.20000random_3kb_intergenic.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.20000random_3kb_intergenic.common_SNP.positionID.txt
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingeneCDS.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingeneCDS.common_SNP.positionID.txt
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingene5UTR.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene5UTR.common_SNP.positionID.txt
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingene3UTR.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene3UTR.common_SNP.positionID.txt


#根据六种区域内common SNP的positionID提取CADD分数
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.pgdlnc_lncRNA_mergeexon.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.pgdlnc_lncRNA_mergeexon.common_SNP.CADD.txt
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.npgdlnc_lncRNA_mergeexon.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.npgdlnc_lncRNA_mergeexon.common_SNP.CADD.txt
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.20000random_3kb_intergenic.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.20000random_3kb_intergenic.common_SNP.CADD.txt
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingeneCDS.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingeneCDS.common_SNP.CADD.txt
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene5UTR.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene5UTR.common_SNP.CADD.txt
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene3UTR.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene3UTR.common_SNP.CADD.txt

```
```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的SNP密度
#人类commonSNP
a <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.pgdlnc_lncRNA_mergeexon.common_SNP.CADD.txt"
b <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.npgdlnc_lncRNA_mergeexon.common_SNP.CADD.txt"
c <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingeneCDS.common_SNP.CADD.txt"
d <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene5UTR.common_SNP.CADD.txt"
e <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene3UTR.common_SNP.CADD.txt"
f <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.20000random_3kb_intergenic.common_SNP.CADD.txt"
g <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.SNPCADD.tj.txt"
h <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.SNPCADD.png"
#导入pgdlnc外显子的SNP的CADD分数
pgdlnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  group_by(V7)%>%
  top_n(1,V6)%>%
  select(6)%>%
  mutate(type="Pseudogene-derived lncRNA")
colnames(pgdlnc_lncRNA_exon) <- c("pID","CADD","type")
#导入npgdlnc外显子的SNP
npgdlnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  group_by(V7)%>%
  top_n(1,V6)%>%
  select(6)%>%
  mutate(type="Non-pseudogene-derived lncRNA")
colnames(npgdlnc_lncRNA_exon) <- colnames(pgdlnc_lncRNA_exon)
#导入蛋白编码基因CDS,5UTR,3UTR的SNP
proteingeneCDS <- read.delim(c, header=FALSE)%>%
  group_by(V7)%>%
  top_n(1,V6)%>%
  select(6)%>%
  mutate(type="CDS")
colnames(proteingeneCDS) <- colnames(pgdlnc_lncRNA_exon)
proteingene5UTR <- read.delim(d, header=FALSE)%>%
  group_by(V7)%>%
  top_n(1,V6)%>%
  select(6)%>%
  mutate(type="5'UTR")
colnames(proteingene5UTR) <- colnames(pgdlnc_lncRNA_exon)
proteingene3UTR <- read.delim(e, header=FALSE)%>%
  group_by(V7)%>%
  top_n(1,V6)%>%
  select(6)%>%
  mutate(type="3'UTR")
colnames(proteingene3UTR) <- colnames(pgdlnc_lncRNA_exon)
#导入基因间区的SNP数量
intergenic <- read.delim(f, header=FALSE)%>%
  group_by(V7)%>%
  top_n(1,V6)%>%
  select(6)%>%
  mutate(type="Random intergenic")
colnames(intergenic) <- colnames(pgdlnc_lncRNA_exon)
#合并
alldata <- rbind(pgdlnc_lncRNA_exon,npgdlnc_lncRNA_exon,proteingeneCDS,proteingene5UTR,proteingene3UTR,intergenic)
alldata$type <- factor(alldata$type,levels = c("Random intergenic","Non-pseudogene-derived lncRNA",
                                               "Pseudogene-derived lncRNA","CDS","5'UTR","3'UTR"))
#统计
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(1024)
tj <- group_by(alldata,type)%>%
  summarise(num=n(),median=median(CADD),mean=mean(CADD),
            confmin=boot::boot.ci(boot::boot(CADD, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(CADD, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][5])
data.table::fwrite(tj,file =g,sep = '\t',row.names = F,quote = F,col.names = T)
#plot log转换的CADD根据排名的标准化分数
p1 <- ggplot(data=alldata, aes(x=type,y=log10(CADD+1)))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  geom_signif(map_signif_level=T,y_position=c(1.6,1.8,2,2.3),
              comparisons = list(c("Random intergenic","Non-pseudogene-derived lncRNA"),
                                 c("Non-pseudogene-derived lncRNA","Pseudogene-derived lncRNA"),
                                 c("Pseudogene-derived lncRNA","CDS"),
                                 c("Pseudogene-derived lncRNA","3'UTR")))+
  scale_fill_manual(values=c("#8491B4","#FAA465","#7197AD","#BB5A5D","#BB5A5D","#BB5A5D"),
                    limits=c("Random intergenic","Non-pseudogene-derived lncRNA",
                             "Pseudogene-derived lncRNA","CDS","5'UTR","3'UTR"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 2.5))+
  scale_x_discrete(labels = c("Random\nintergenic","Non-pseudogene-\nderived lncRNA",
                              "Pseudogene-\nderived lncRNA","CDS","5'UTR","3'UTR"))+
  labs(y =expression("S"*"N"*"P"~"l"*"o"*"g"[10]*"("*"C"*"A"*"D"*"D"~"s"*"c"*"o"*"r"*"e"*")"),x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave(h,p1,width = 6, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")



```


##### 6.SNP DAF
```r
#从ensemble数据库下载千人基因组项目phase 3的SNP结果
wget https://ftp.ensembl.org/pub/release-110/variation/gvf/homo_sapiens/1000GENOMES-phase_3.gvf.gz -P /home/zhzhang/PG/Evolution/snp/
gzip -d /home/zhzhang/PG/Evolution/snp/1000GENOMES-phase_3.gvf.gz

#获取common SNP ID
awk '{print $4}' /home/zhzhang/PG/Evolution/snp/common_allsnp.bed > /home/zhzhang/PG/Evolution/snp/common_allsnp.ID.txt
#提取千人基因组项目phase 3的SNP结果中的common SNP
grep -w -Ff /home/zhzhang/PG/Evolution/snp/common_allsnp.ID.txt "/home/zhzhang/PG/Evolution/snp/1000GENOMES-phase_3.gvf" > /home/zhzhang/PG/Evolution/snp/1000GENOMES-phase_3.common.gvf

#提取原文件每行的第1列，第四列，和第五列。以及最后一列中ID，AMR，SAS，EUR，EAS，AFR。值得注意的是，每行中AMR，SAS，EUR，EAS，AFR的顺序并不固定，并且有的行没有这些属性中的一个或几个。此外，AMR，SAS，EUR，EAS，AFR这几个属性，有的值包含多个由英文逗号分隔的值，请将这些由英文逗号分隔的值加和。并且，不具有这些属性中的一个或几个的行，对应的属性输出为NA。最后添加一列，计算AMR，SAS，EUR，EAS，AFR这五个属性的平均值，NA值不参与计算
awk 'BEGIN {OFS="\t"} 
     {
         id="NA"; amr="NA"; sas="NA"; eur="NA"; eas="NA"; afr="NA";
         sum=0; count=0;
         n=split($9, arr, ";");
         for (i=1; i<=n; i++) {
             split(arr[i], temp, "=");
             if (temp[1] == "ID") {
                 id = temp[2];
             } else if (temp[1] == "AMR") {
                 split(temp[2], values, ",");
                 for (j=1; j<=length(values); j++) {
                     amr += values[j];
                     if (values[j] != "NA") {
                         sum += values[j]
                     }
                 }
                 if (amr != "NA") {
                         count++
                     }
             } else if (temp[1] == "SAS") {
                 split(temp[2], values, ",");
                 for (j=1; j<=length(values); j++) {
                     sas += values[j];
                     if (values[j] != "NA") {
                         sum += values[j]
                     }
                 }
                 if (sas != "NA") {
                         count++
                     }
             } else if (temp[1] == "EUR") {
                 split(temp[2], values, ",");
                 for (j=1; j<=length(values); j++) {
                     eur += values[j];
                     if (values[j] != "NA") {
                         sum += values[j]
                     }
                 }
                 if (eur != "NA") {
                         count++
                     }
             } else if (temp[1] == "EAS") {
                 split(temp[2], values, ",");
                 for (j=1; j<=length(values); j++) {
                     eas += values[j];
                     if (values[j] != "NA") {
                         sum += values[j]
                     }
                 }
                 if (eas != "NA") {
                         count++
                     }
             } else if (temp[1] == "AFR") {
                 split(temp[2], values, ",");
                 for (j=1; j<=length(values); j++) {
                     afr += values[j];
                     if (values[j] != "NA") {
                         sum += values[j]
                     }
                 }
                 if (afr != "NA") {
                         count++
                     }
             }
         }
         if (count > 0) {
             avg = sum / count;
         } else {
             avg = "NA";
         }
         print $1"\t"$4"\t"$5"\t"id"\t"amr"\t"sas"\t"eur"\t"eas"\t"afr"\t"avg;
     }' /home/zhzhang/PG/Evolution/snp/1000GENOMES-phase_3.common.gvf > "/home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.DAF.bed"





#脚本
# 打开原文件
with open('/home/zhzhang/PG/Evolution/snp/1000GENOMES-phase_3.common.gvf', 'r') as file:
    for line in file:
        # 分割每行数据并提取需要的信息
        columns = line.strip().split('\t')
        col1 = columns[0]
        col4 = columns[3]
        col5 = columns[4]

        # 处理最后一列中ID和AFR属性
        last_column = columns[-1]
        id_value = ''
        afr_value = 'NA'

        if 'ID=' in last_column:
            id_index = last_column.index('ID=') + 3
            id_value = last_column[id_index:last_column.index(';', id_index)]

        if 'AFR=' in last_column:
            afr_index = last_column.index('AFR=') + 4
            afr_values = [float(val) for val in last_column[afr_index:].split(';')[0].split(',') if val]  # 提取AFR属性的值并转换为浮点数列表
            afr_sum = sum(afr_values)  # 对AFR属性的值求和
            afr_value = str(afr_sum)

        print(col1, col4, col5, id_value, afr_value, sep='\t')

#执行脚本
python "/home/zhzhang/PG/Evolution/snp/DAF/chuli.py" > /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed


```
```r
#获取六种区域中的SNP
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.pgdlnc_lncRNA_mergeexon.bed -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.pgdlnc_lncRNA_mergeexon.SNPDAF.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.npgdlnc_lncRNA_mergeexon.bed -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.npgdlnc_lncRNA_mergeexon.SNPDAF.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.bed -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.20000random_3kb_intergenic.SNPDAF.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.bed -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingeneCDS.SNPDAF.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.bed -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene5UTR.SNPDAF.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.bed -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene3UTR.SNPDAF.txt



#从ucsc获取CPG岛区域（防止CPG岛区域高突变率的影响）
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz -P /home/zhzhang/PG/Evolution/snp/
gzip -d /home/zhzhang/PG/Evolution/snp/cpgIslandExt.txt.gz
awk '$2!~/_/{print $2"\t"$3"\t"$4}' /home/zhzhang/PG/Evolution/snp/cpgIslandExt.txt |sed 's/chr//g' > /home/zhzhang/PG/Evolution/snp/cpgIsland.bed

#去除上述数据中的CPG岛区域SNP
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.pgdlnc_lncRNA_mergeexon.SNPDAF.txt -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed -v > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.pgdlnc_lncRNA_mergeexon.SNPDAF.noncpg.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.npgdlnc_lncRNA_mergeexon.SNPDAF.txt -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed -v > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.npgdlnc_lncRNA_mergeexon.SNPDAF.noncpg.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.20000random_3kb_intergenic.SNPDAF.txt -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed -v > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.20000random_3kb_intergenic.SNPDAF.noncpg.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingeneCDS.SNPDAF.txt -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed -v > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingeneCDS.SNPDAF.noncpg.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene5UTR.SNPDAF.txt -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed -v > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene5UTR.SNPDAF.noncpg.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene3UTR.SNPDAF.txt -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed -v > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene3UTR.SNPDAF.noncpg.txt



```
```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的SNP DAF
#人类commonSNP
a <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.pgdlnc_lncRNA_mergeexon.SNPDAF.noncpg.txt"
b <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.npgdlnc_lncRNA_mergeexon.SNPDAF.noncpg.txt"
c <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingeneCDS.SNPDAF.noncpg.txt"
d <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene5UTR.SNPDAF.noncpg.txt"
e <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene3UTR.SNPDAF.noncpg.txt"
f <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.20000random_3kb_intergenic.SNPDAF.noncpg.txt"
g <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.lncRNA_SNPDAF.tj.txt"
h <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.lncRNA_SNPDAF.png"
#导入pgdlnc外显子的SNP,去除averageDAF为NA 的SNP
pgdlnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  select(5)%>%
  filter(V5<=1)%>%
  mutate(type="Pseudogene-derived lncRNA")
colnames(pgdlnc_lncRNA_exon)[c(1,2)] <- c("DAF","type")
#导入npgdlnc外显子的SNP
npgdlnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  select(5)%>%
  filter(V5<=1)%>%
  mutate(type="Non-pseudogene-derived lncRNA")
colnames(npgdlnc_lncRNA_exon)[c(1,2)] <- colnames(pgdlnc_lncRNA_exon)
#导入蛋白编码基因CDS,5UTR,3UTR的SNP
proteingeneCDS <- read.delim(c, header=FALSE)%>%
  select(5)%>%
  filter(V5<=1)%>%
  mutate(type="CDS")
colnames(proteingeneCDS) <- colnames(pgdlnc_lncRNA_exon)
proteingene5UTR <- read.delim(d, header=FALSE)%>%
  select(5)%>%
  filter(V5<=1)%>%
  mutate(type="5'UTR")
colnames(proteingene5UTR) <- colnames(pgdlnc_lncRNA_exon)
proteingene3UTR <- read.delim(e, header=FALSE)%>%
  select(5)%>%
  filter(V5<=1)%>%
  mutate(type="3'UTR")
colnames(proteingene3UTR) <- colnames(pgdlnc_lncRNA_exon)
#导入基因间区的SNP数量
intergenic <- read.delim(f, header=FALSE)%>%
  select(5)%>%
  filter(V5<=1)%>%
  mutate(type="Random intergenic")
colnames(intergenic) <- colnames(pgdlnc_lncRNA_exon)
#合并
alldata <- rbind(pgdlnc_lncRNA_exon,npgdlnc_lncRNA_exon,proteingeneCDS,proteingene5UTR,proteingene3UTR,intergenic)
alldata$type <- factor(alldata$type,levels = c("Random intergenic","Non-pseudogene-derived lncRNA",
                                               "Pseudogene-derived lncRNA","CDS","5'UTR","3'UTR"))
#统计
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(1024)
tj <- group_by(alldata,type)%>%
  summarise(num=n(),median=median(DAF),mean=mean(DAF),
            confmin=boot::boot.ci(boot::boot(DAF, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(DAF, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][5])

data.table::fwrite(tj,file =g,sep = '\t',row.names = F,quote = F,col.names = T)
#plot
t.test(pgdlnc_lncRNA_exon$DAF,intergenic$DAF) #p-value = 0.0000000477
t.test(npgdlnc_lncRNA_exon$DAF,intergenic$DAF) #p-value < 2.2e-16
t.test(pgdlnc_lncRNA_exon$DAF,npgdlnc_lncRNA_exon$DAF) #p-value = 0.8547
t.test(pgdlnc_lncRNA_exon$DAF,proteingeneCDS$DAF) #p-value < 2.2e-16
t.test(pgdlnc_lncRNA_exon$DAF,proteingene3UTR$DAF) #p-value < 2.2e-16

phs <- ggplot(data = tj,aes(x=type,y=mean))+
  geom_point(size=2,aes(fill=type,color=type))+
  geom_errorbar(aes(ymin=confmin,ymax=confmax,color=type),width=0.15)+
  geom_signif(annotations=c("***","***","N.S.","***","***"),
              y_position=c(0.121,0.125,0.116,0.119,0.123),
              xmin = c(1,1,2,3,3),xmax = c(2,3,3,4,6))+
  scale_fill_manual(values=c("#8491B4","#FAA465","#7197AD","#BB5A5D","#BB5A5D","#BB5A5D"),
                    limits=c("Random intergenic","Non-pseudogene-derived lncRNA",
                             "Pseudogene-derived lncRNA","CDS","5'UTR","3'UTR"))+
  scale_color_manual(values=c("#8491B4","#FAA465","#7197AD","#BB5A5D","#BB5A5D","#BB5A5D"),
                    limits=c("Random intergenic","Non-pseudogene-derived lncRNA",
                             "Pseudogene-derived lncRNA","CDS","5'UTR","3'UTR"))+
  theme_half_open()+
  scale_x_discrete(labels = c("Random\nintergenic","Non-pseudogene-\nderived lncRNA",
                              "Pseudogene-\nderived lncRNA","CDS","5'UTR","3'UTR"))+
  labs(y = "SNP derived allele frequency", x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))

ggsave(h,phs,width = 6, height = 5.6,dpi=1200, units = "in", device='png',bg = "transparent")


```




##### 7.lncRNA演化年龄
```r
### download UCSC chain files
# Human
#全部lncRNA的外显子bed
/home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.all_lncRNA_exon.bed
#######人类
#主物种信息
vipsp="Homo_sapiens"
gn="hg38"
vipspdir="hs"
vipgtf="/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf"
vipfa="/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa"
#
cd /home/zhzhang/PG/age/chain/${vipspdir}/
grep -w -v "${vipsp}" "/home/zhzhang/PG/age/chain/sh/sp.txt"|while read i
do
source /home/zhzhang/miniconda3/bin/activate SEQ
arr=($i)
quan=${arr[0]}
jian=${arr[1]}
chain1=`ls ${jian}/${gn}To*.over.chain.gz`
chrsize=`ls /home/zhzhang/PG/refgenome/${quan}.*.dna.chrsize.txt`
gtf=`ls /home/zhzhang/PG/RNAseqdata/newGTF/${quan}.*.chr.rmpg.novellncRNA.gtf`
fa=`ls /home/zhzhang/PG/refgenome/${quan}.*.dna.chr.fa`
#start
echo "# == 1.align human lncRNA to ${jian} =="
awk '{print "chr"$0}' /home/zhzhang/PG/lncRNA_class/${vipsp}/${vipsp}.all_lncRNA_exon.bed | liftOver stdin ${chain1} stdout /dev/null -minMatch=0.1 |sed 's/^chr//g' > /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.bed
if [ "$jian" == "Chicken" ]
then
awk '{print "chr"$0}' /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.bed| liftOver stdin "/home/zhzhang/PG/age/chain/hs/Chicken/galGal6ToGCF_016699485.2.over.chain.gz" stdout /dev/null -minMatch=0.1 |sed "s/NC_052532.1/1/g;s/NC_052533.1/2/g;s/NC_052534.1/3/g;s/NC_052535.1/4/g;s/NC_052536.1/5/g;s/NC_052537.1/6/g;s/NC_052538.1/7/g;s/NC_052539.1/8/g;s/NC_052540.1/9/g;s/NC_052541.1/10/g;s/NC_052542.1/11/g;s/NC_052543.1/12/g;s/NC_052544.1/13/g;s/NC_052545.1/14/g;s/NC_052546.1/15/g;s/NC_052547.1/16/g;s/NC_052548.1/17/g;s/NC_052549.1/18/g;s/NC_052550.1/19/g;s/NC_052551.1/20/g;s/NC_052552.1/21/g;s/NC_052553.1/22/g;s/NC_052554.1/23/g;s/NC_052555.1/24/g;s/NC_052556.1/25/g;s/NC_052557.1/26/g;s/NC_052558.1/27/g;s/NC_052559.1/28/g;s/NC_052560.1/29/g;s/NC_052561.1/30/g;s/NC_052562.1/31/g;s/NC_052563.1/32/g;s/NC_052564.1/33/g;s/NC_052565.1/34/g;s/NC_052566.1/35/g;s/NC_052567.1/36/g;s/NC_052568.1/37/g;s/NC_052569.1/38/g;s/NC_052570.1/39/g;s/NC_052571.1/W/g;s/NC_052572.1/Z/g" > /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.1.bed
rm /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.bed
mv /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.1.bed /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.bed
fi
awk '$1~/^[1234567890]$/ || $1~/^[1234567890][1234567890]$/ || $1~/^[XYWZ]$/ {print $0}' /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.bed |bedtools slop -g ${chrsize} -b 10000 |bedtools sort > /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.sort.bed
echo "# == 2.intersect human lncRNA with ${jian} lncRNA =="
bedtools intersect -a /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.sort.bed -b /home/zhzhang/PG/lncRNA_class/${quan}/${quan}.all_lncRNA_exon.bed -wo -s > /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.intersect.${quan}.all_lncRNA_exon.bed
echo "# == 3.extract human lncRNAs which have ${jian} syntenic-lncRNA =="
source /home/zhzhang/miniconda3/bin/activate rbase
Rscript /home/zhzhang/PG/age/chain/sh/step3.R --input /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipsp}.${jian}.all_lncRNA_exon.intersect.${quan}.all_lncRNA_exon.bed --sp ${jian} --outputone /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.${quan}.homo_lncRNApair.txt --outputtwo /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.own_${quan}_homo.lncRNAid.txt
done


#######小鼠
vipsp="Mus_musculus"
gn="mm39"
vipspdir="mm"
vipgtf="/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf"
vipfa="/home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa"


################# R脚本Step3 ################
library(dplyr)
library(tidyverse)
#生成参数列表
option_list <- list(optparse::make_option(
  opt_str ="--input",
  type = "character",
  default = NULL,
  help=""),
  optparse::make_option(
    opt_str ="--sp",
    type = "character",
    default = NULL,
    help=""
  ),
  optparse::make_option(
    opt_str ="--outputone",
    type = "character",
    default = NULL,
    help=""
  ),
  optparse::make_option(
    opt_str ="--outputtwo",
    type = "character",
    default = NULL,
    help=""
  ))
#将参数列表传递给一个变量
args <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
#导入
intersect <- read.delim(args$input, header=FALSE)
colnames(intersect)[c(4,5,10,11,13)] <- c("geneid_1","transcriptid_1","geneid_2","transcriptid_2","len")
#统计每个基因每个转录本和另一物种每个基因每个转录本之间的交集
tj <- group_by(intersect,geneid_1,transcriptid_1,geneid_2,transcriptid_2)%>%
  summarise(len=sum(len))
#去除原分组
newtj <- data.frame(tj$geneid_1,tj$transcriptid_1,tj$geneid_2,tj$transcriptid_2,tj$len)
colnames(newtj) <- colnames(tj)
#每种syntenic基因对仅保留交集长度最长的转录本id，输出
genepair <- mutate(newtj,gp=paste(geneid_1,geneid_2,sep="----"))%>%
  group_by(gp)%>%
  top_n(1,len)%>%
  distinct(gp,.keep_all=T)
data.table::fwrite(genepair,file =args$outputone,sep = '\t',row.names = F,quote = F,col.names = T)
#仅保留在另一物种具有syntenic lnc的主物种lnc基因id
ownhomolnc <- distinct(newtj,geneid_1)%>%
  mutate(homo=1)
colnames(ownhomolnc) <- c("geneid",args$sp)
data.table::fwrite(ownhomolnc,file =args$outputtwo,sep = '\t',row.names = F,quote = F,col.names = T)
################################################



#
echo "# == 4.lastz between human lncRNAs and ${jian} syntenic-lncRNA =="
tail -n +2 /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.${quan}.homo_lncRNApair.txt | while read ii
do
arrarr=($ii)
lnc1=${arrarr[0]}
trans1=${arrarr[1]}
lnc2=${arrarr[2]}
trans2=${arrarr[3]}
grep -w "${trans1}" "${vipgtf}" > /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/trans1.gtf
grep -w "${trans2}" "${gtf}" > /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/trans2.gtf
gffread /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/trans1.gtf -g "${vipfa}" -w "/home/zhzhang/PG/age/chain/${vipspdir}/${jian}/trans1.fa"
gffread /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/trans2.gtf -g "${fa}" -w "/home/zhzhang/PG/age/chain/${vipspdir}/${jian}/trans2.fa"
alscore=`lastz /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/trans1.fa /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/trans2.fa --hspthresh=0 --format=MAF|grep "a score="|sed 's/a score=//g'|awk 'BEGIN {max=-1} {if ($1 > max) max=$1} END {print max}'`
echo -e "${lnc1}\t${lnc2}\t${alscore}" >> /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.${quan}.homo_lncRNApair.alscore.txt
done



################### R脚本Step5_pro ######################
library(dplyr)
library(tidyverse)
#生成参数列表
option_list <- list(optparse::make_option(
  opt_str ="--input",
  type = "character",
  default = NULL,
  help=""),
  optparse::make_option(
    opt_str ="--sp",
    type = "character",
    default = NULL,
    help=""
  ),
  optparse::make_option(
    opt_str ="--outputone",
    type = "character",
    default = NULL,
    help=""
  ),
  optparse::make_option(
    opt_str ="--outputtwo",
    type = "character",
    default = NULL,
    help=""
  ))
#将参数列表传递给一个变量
args <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
#导入同源lnc对比对得分
alscore <- read.delim(args$input, header=FALSE)
colnames(alscore) <- c("geneid_1","geneid_2","score")
#筛选出互惠最佳1-1同源lnc对
aaa <- 9999999
homo1v1pair <- data.frame()
while (aaa != 0) {
  if (aaa == 9999999) {
    whilealscore <- alscore
  }else {
    whilealscore <- whilealscore[whilealscore$geneid_1 %in% homo1v1pairwhile$geneid_1 ==F,]
    whilealscore <- whilealscore[whilealscore$geneid_2 %in% homo1v1pairwhile$geneid_2 ==F,]
  }
  filtergene1 <- group_by(whilealscore,geneid_1)%>%
    top_n(1,score)%>%
    mutate(lncpair=paste(geneid_1,geneid_2,sep = "----"))
  filtergene2 <- group_by(whilealscore,geneid_2)%>%
    top_n(1,score)%>%
    mutate(lncpair=paste(geneid_1,geneid_2,sep = "----"))
  homo1v1pairwhile <- data.frame("lncpair"=intersect(filtergene1$lncpair,filtergene2$lncpair))%>%
    separate(lncpair,c("geneid_1","geneid_2"),sep="----")
  aaa <- nrow(homo1v1pairwhile)
  homo1v1pair <- rbind(homo1v1pair,homo1v1pairwhile)
}
#储存1-1同源对
data.table::fwrite(homo1v1pair,
                   file =args$outputone,
                   sep = '\t',row.names = F,quote = F,col.names = T)
#储存具有1-1同源lnc的主物种lncid
ownhomo <- data.frame(homo1v1pair$geneid_1)%>%
  mutate(homo=1)
colnames(ownhomo) <- c("geneid",args$sp)
data.table::fwrite(ownhomo,
                   file =args$outputtwo,
                   sep = '\t',row.names = F,quote = F,col.names = T)
#############################################################






echo "# == 5.extract human lncRNAs which have 1-1 ${jian} homo-lncRNA =="
source /home/zhzhang/miniconda3/bin/activate rbase
Rscript /home/zhzhang/PG/age/chain/sh/step5.pro.R --input /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.${quan}.homo_lncRNApair.alscore.txt --sp ${jian} --outputone /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.${quan}.11homo_lncRNApair.pro.txt --outputtwo /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.own_${quan}_11homo.lncRNAid.pro.txt





```
```r
#Human
############
#导入全部lncRNA id和分类
lncRNA_class <- read.delim("~/PG/lncRNA_class/Homo_sapiens.lncRNA_class.txt")
#导入第一个年龄段
Macaque_homo <- read.delim("~/PG/age/chain/hs/Macaque/hs.own_Macaca_mulatta_11homo.lncRNAid.pro.txt")%>%
  mutate(Macaque=1)%>%
  distinct(geneid,.keep_all = T)
#导入第二个年龄段
Mouse_homo <- read.delim("~/PG/age/chain/hs/Mouse/hs.own_Mus_musculus_11homo.lncRNAid.pro.txt")%>%
  mutate(Mouse=2)%>%
  distinct(geneid,.keep_all = T)
Rat_homo <- read.delim("~/PG/age/chain/hs/Rat/hs.own_Rattus_norvegicus_11homo.lncRNAid.pro.txt")%>%
  mutate(Rat=2)%>%
  distinct(geneid,.keep_all = T)
Rabbit_homo <- read.delim("~/PG/age/chain/hs/Rabbit/hs.own_Oryctolagus_cuniculus_11homo.lncRNAid.pro.txt")%>%
  mutate(Rabbit=2)%>%
  distinct(geneid,.keep_all = T)
#导入第三个年龄段
Opossum_homo <- read.delim("~/PG/age/chain/hs/Opossum/hs.own_Monodelphis_domestica_11homo.lncRNAid.pro.txt")%>%
  mutate(Opossum=3)%>%
  distinct(geneid,.keep_all = T)
#导入第四个年龄段
Chicken_homo <- read.delim("~/PG/age/chain/hs/Chicken/hs.own_Gallus_gallus_11homo.lncRNAid.pro.txt")%>%
  mutate(Chicken=4)%>%
  distinct(geneid,.keep_all = T)
#导入第五个年龄段
Zebrafish_homo <- read.delim("~/PG/age/chain/hs/Zebrafish/hs.own_Danio_rerio_11homo.lncRNAid.pro.txt")%>%
  mutate(Zebrafish=5)%>%
  distinct(geneid,.keep_all = T)
#形成矩阵
agematrix <- left_join(select(lncRNA_class,1),Macaque_homo,by="geneid")%>%
  left_join(Mouse_homo,by="geneid")%>%
  left_join(Rat_homo,by="geneid")%>%
  left_join(Rabbit_homo,by="geneid")%>%
  left_join(Opossum_homo,by="geneid")%>%
  left_join(Chicken_homo,by="geneid")%>%
  left_join(Zebrafish_homo,by="geneid")%>%
  column_to_rownames("geneid")
agematrix$Macaque[is.na(agematrix$Macaque)==T] <- 0
agematrix$Mouse[is.na(agematrix$Mouse)==T] <- 0
agematrix$Rat[is.na(agematrix$Rat)==T] <- 0
agematrix$Rabbit[is.na(agematrix$Rabbit)==T] <- 0
agematrix$Opossum[is.na(agematrix$Opossum)==T] <- 0
agematrix$Chicken[is.na(agematrix$Chicken)==T] <- 0
agematrix$Zebrafish[is.na(agematrix$Zebrafish)==T] <- 0
#根据lnc最古老同源物确定年龄
lncage <- apply(agematrix,1,max)%>%
  data.frame()
colnames(lncage) <- "max"
lncage$max <- as.character(lncage$max)
lncage$max[lncage$max=="1"] <- "29 Ma"
lncage$max[lncage$max=="2"] <- "87 Ma"
lncage$max[lncage$max=="3"] <- "160 Ma"
lncage$max[lncage$max=="4"] <- "319 Ma"
lncage$max[lncage$max=="5"] <- "429 Ma"
lncage$max[lncage$max=="0"] <- "Human"
lncage <- rownames_to_column(lncage,"geneid")
#导入具有潜在同源对的主物种lncid
Macaque_phomo <- read.delim("~/PG/age/chain/hs/Macaque/hs.own_Macaca_mulatta_homo.lncRNAid.txt")%>%
  select(1)
Mouse_phomo <- read.delim("~/PG/age/chain/hs/Mouse/hs.own_Mus_musculus_homo.lncRNAid.txt")%>%
  select(1)
Rat_phomo <- read.delim("~/PG/age/chain/hs/Rat/hs.own_Rattus_norvegicus_homo.lncRNAid.txt")%>%
  select(1)
Rabbit_phomo <- read.delim("~/PG/age/chain/hs/Rabbit/hs.own_Oryctolagus_cuniculus_homo.lncRNAid.txt")%>%
  select(1)
Opossum_phomo <- read.delim("~/PG/age/chain/hs/Opossum/hs.own_Monodelphis_domestica_homo.lncRNAid.txt")%>%
  select(1)
Chicken_phomo <- read.delim("~/PG/age/chain/hs/Chicken/hs.own_Gallus_gallus_homo.lncRNAid.txt")%>%
  select(1)
Zebrafish_phomo <- read.delim("~/PG/age/chain/hs/Zebrafish/hs.own_Danio_rerio_homo.lncRNAid.txt")%>%
  select(1)
#确定至少一次拥有潜在homo的lnc
lnc_phomo <- distinct(rbind(Macaque_phomo,Mouse_phomo,Rat_phomo,Rabbit_phomo,Opossum_phomo,
                            Chicken_phomo,Zebrafish_phomo),geneid)%>%
  mutate(phomo=1)
#没有1-1同源lnc但有潜在homo的lnc(年龄模糊)
lncage <- left_join(lncage,lnc_phomo,by="geneid")%>%
  mutate(ambiguous=0)
lncage$phomo[is.na(lncage$phomo)==T] <- 0
lncage$ambiguous[lncage$phomo==1 & lncage$max=="Human"] <- 1
lncage$max[lncage$ambiguous==1] <- "Ambiguous"
lncage <- select(lncage,1,2)
#储存
data.table::fwrite(lncage,file ="/home/zhzhang/PG/age/Human.lncRNA.age.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#富集
anno <- mutate(lncage,GO_Description=max)
colnames(anno)[2] <- "GO"
go2gene <- anno[, c(2, 1)]
go2name <- anno[, c(2, 3)]
npglnc <- filter(lncRNA_class,type=="Non-pseudogene-derived lncRNA")%>%
  select(1)
pglnc <- filter(lncRNA_class,type=="Pseudogene-derived lncRNA")%>%
  select(1)
ego_npg <- clusterProfiler::enricher(npglnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                    maxGSSize = 40000)@result%>%
  mutate(type="Non-pseudogene-derived lncRNA")
ego_pg <- clusterProfiler::enricher(pglnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                    maxGSSize = 40000)@result%>%
  mutate(type="Pseudogene-derived lncRNA")
#合并
ego <- rbind(ego_pg,ego_npg)%>%
  select(10,1:9)
#储存
data.table::fwrite(ego,file ="/home/zhzhang/PG/age/Human.lncRNA.type.age_enrich.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)

#plot
ego <- read.delim("~/PG/age/Human.lncRNA.type.age_enrich.txt")

ego <- filter(ego,ID!="Ambiguous")%>%
  separate(GeneRatio,c("GeneRatio1","GeneRatio2"))%>%
  separate(BgRatio,c("BgRatio1","BgRatio2"))%>%
  mutate(GeneRatio1=as.numeric(GeneRatio1),GeneRatio2=as.numeric(GeneRatio2),
         BgRatio1=as.numeric(BgRatio1),BgRatio2=as.numeric(BgRatio2))%>%
  mutate(ES=(GeneRatio1/GeneRatio2)/(BgRatio1/BgRatio2))
ego$type <- factor(ego$type,levels=c("Non-pseudogene-derived lncRNA",
                                    "Pseudogene-derived lncRNA"))
ego$ID <- factor(ego$ID,levels=c("429 Ma","319 Ma","160 Ma",
                                    "87 Ma","29 Ma","Human"))
#fishertest检验两类lnc指定年龄基因比例的差异p值
pdata <- distinct(ego,ID)
for (i in 1:nrow(pdata)) {
  forid <- pdata$ID[i]
  spglncnum <- filter(ego,type=="Pseudogene-derived lncRNA" & ID==forid)$GeneRatio1
  snpglncnum <- filter(ego,type=="Non-pseudogene-derived lncRNA" & ID==forid)$GeneRatio1
  pdata[i,2] <- fisher.test(data.frame(c(spglncnum,snpglncnum),c(1213-spglncnum,37003-snpglncnum)))[["p.value"]]
}
#富集dotplot
pp1 <- ggplot(data = ego,aes(x=type,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:5.5),color="gray")+
  geom_vline(xintercept=c(1.5),color="gray")+
  scale_fill_gradient(low="white",high="#C12039",breaks=c(0,10,20,30))+
  scale_size(breaks = c(0.5,1,1.5,2,2.5),range = c(1,15))+
  scale_x_discrete(labels=c("Non-pseudogene-\nderived lncRNA",
                            "Pseudogene-derived\nlncRNA"))+
  theme_half_open()+
  labs(x = NULL, y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave("/home/zhzhang/PG/age/Human.lnc.enrich.pdf", 
       pp1,width = 4.8, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")










#Mouse
############
#导入全部lncRNA id和分类
lncRNA_class <- read.delim("~/PG/lncRNA_class/Mus_musculus.lncRNA_class.txt")
#导入第一个年龄段
Rat_homo <- read.delim("~/PG/age/chain/mm/Rat/mm.own_Rattus_norvegicus_11homo.lncRNAid.pro.txt")%>%
  mutate(Rat=1)%>%
  distinct(geneid,.keep_all = T)
#导入第二个年龄段
Rabbit_homo <- read.delim("~/PG/age/chain/mm/Rabbit/mm.own_Oryctolagus_cuniculus_11homo.lncRNAid.pro.txt")%>%
  mutate(Rabbit=2)%>%
  distinct(geneid,.keep_all = T)
#导入第三个年龄段
Human_homo <- read.delim("~/PG/age/chain/mm/Human/mm.own_Homo_sapiens_11homo.lncRNAid.pro.txt")%>%
  mutate(Human=3)%>%
  distinct(geneid,.keep_all = T)
Macaque_homo <- read.delim("~/PG/age/chain/mm/Macaque/mm.own_Macaca_mulatta_11homo.lncRNAid.pro.txt")%>%
  mutate(Macaque=3)%>%
  distinct(geneid,.keep_all = T)
#导入第四个年龄段
Opossum_homo <- read.delim("~/PG/age/chain/mm/Opossum/mm.own_Monodelphis_domestica_11homo.lncRNAid.pro.txt")%>%
  mutate(Opossum=4)%>%
  distinct(geneid,.keep_all = T)
#导入第五个年龄段
Chicken_homo <- read.delim("~/PG/age/chain/mm/Chicken/mm.own_Gallus_gallus_11homo.lncRNAid.pro.txt")%>%
  mutate(Chicken=5)%>%
  distinct(geneid,.keep_all = T)
#导入第六个年龄段
Zebrafish_homo <- read.delim("~/PG/age/chain/mm/Zebrafish/mm.own_Danio_rerio_11homo.lncRNAid.pro.txt")%>%
  mutate(Zebrafish=6)%>%
  distinct(geneid,.keep_all = T)
#形成矩阵
agematrix <- left_join(select(lncRNA_class,1),Macaque_homo,by="geneid")%>%
  left_join(Human_homo,by="geneid")%>%
  left_join(Rat_homo,by="geneid")%>%
  left_join(Rabbit_homo,by="geneid")%>%
  left_join(Opossum_homo,by="geneid")%>%
  left_join(Chicken_homo,by="geneid")%>%
  left_join(Zebrafish_homo,by="geneid")%>%
  column_to_rownames("geneid")
agematrix$Macaque[is.na(agematrix$Macaque)==T] <- 0
agematrix$Human[is.na(agematrix$Human)==T] <- 0
agematrix$Rat[is.na(agematrix$Rat)==T] <- 0
agematrix$Rabbit[is.na(agematrix$Rabbit)==T] <- 0
agematrix$Opossum[is.na(agematrix$Opossum)==T] <- 0
agematrix$Chicken[is.na(agematrix$Chicken)==T] <- 0
agematrix$Zebrafish[is.na(agematrix$Zebrafish)==T] <- 0
#根据lnc最古老同源物确定年龄
lncage <- apply(agematrix,1,max)%>%
  data.frame()
colnames(lncage) <- "max"
lncage$max <- as.character(lncage$max)
lncage$max[lncage$max=="1"] <- "13 Ma"
lncage$max[lncage$max=="2"] <- "79 Ma"
lncage$max[lncage$max=="3"] <- "87 Ma"
lncage$max[lncage$max=="4"] <- "160 Ma"
lncage$max[lncage$max=="5"] <- "319 Ma"
lncage$max[lncage$max=="6"] <- "429 Ma"
lncage$max[lncage$max=="0"] <- "Mouse"
lncage <- rownames_to_column(lncage,"geneid")
#导入具有潜在同源对的主物种lncid
Macaque_phomo <- read.delim("~/PG/age/chain/mm/Macaque/mm.own_Macaca_mulatta_homo.lncRNAid.txt")%>%
  select(1)
Human_phomo <- read.delim("~/PG/age/chain/mm/Human/mm.own_Homo_sapiens_homo.lncRNAid.txt")%>%
  select(1)
Rat_phomo <- read.delim("~/PG/age/chain/mm/Rat/mm.own_Rattus_norvegicus_homo.lncRNAid.txt")%>%
  select(1)
Rabbit_phomo <- read.delim("~/PG/age/chain/mm/Rabbit/mm.own_Oryctolagus_cuniculus_homo.lncRNAid.txt")%>%
  select(1)
Opossum_phomo <- read.delim("~/PG/age/chain/mm/Opossum/mm.own_Monodelphis_domestica_homo.lncRNAid.txt")%>%
  select(1)
Chicken_phomo <- read.delim("~/PG/age/chain/mm/Chicken/mm.own_Gallus_gallus_homo.lncRNAid.txt")%>%
  select(1)
Zebrafish_phomo <- read.delim("~/PG/age/chain/mm/Zebrafish/mm.own_Danio_rerio_homo.lncRNAid.txt")%>%
  select(1)
#确定至少一次拥有潜在homo的lnc
lnc_phomo <- distinct(rbind(Macaque_phomo,Human_phomo,Rat_phomo,Rabbit_phomo,Opossum_phomo,
                            Chicken_phomo,Zebrafish_phomo),geneid)%>%
  mutate(phomo=1)
#去除没有1-1同源lnc但有潜在homo的lnc(年龄模糊)
lncage <- left_join(lncage,lnc_phomo,by="geneid")%>%
  mutate(ambiguous=0)
lncage$phomo[is.na(lncage$phomo)==T] <- 0
lncage$ambiguous[lncage$phomo==1 & lncage$max=="Mouse"] <- 1
lncage$max[lncage$ambiguous==1] <- "Ambiguous"
lncage <- select(lncage,1,2)
#储存
data.table::fwrite(lncage,file ="/home/zhzhang/PG/age/Mouse.lncRNA.age.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)

#富集
anno <- mutate(lncage,GO_Description=max)
colnames(anno)[2] <- "GO"
go2gene <- anno[, c(2, 1)]
go2name <- anno[, c(2, 3)]
npglnc <- filter(lncRNA_class,type=="Non-pseudogene-derived lncRNA")%>%
  select(1)
pglnc <- filter(lncRNA_class,type=="Pseudogene-derived lncRNA")%>%
  select(1)
ego_npg <- clusterProfiler::enricher(npglnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                     maxGSSize = 40000)@result%>%
  mutate(type="Non-pseudogene-derived lncRNA")
ego_pg <- clusterProfiler::enricher(pglnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                    maxGSSize = 40000)@result%>%
  mutate(type="Pseudogene-derived lncRNA")
#合并
ego <- rbind(ego_pg,ego_npg)%>%
  select(10,1:9)
#储存
data.table::fwrite(ego,file ="/home/zhzhang/PG/age/Mouse.lncRNA.type.age_enrich.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)

#plot
ego <- read.delim("~/PG/age/Mouse.lncRNA.type.age_enrich.txt")

ego <- filter(ego,ID!="Ambiguous")%>%
  separate(GeneRatio,c("GeneRatio1","GeneRatio2"))%>%
  separate(BgRatio,c("BgRatio1","BgRatio2"))%>%
  mutate(GeneRatio1=as.numeric(GeneRatio1),GeneRatio2=as.numeric(GeneRatio2),
         BgRatio1=as.numeric(BgRatio1),BgRatio2=as.numeric(BgRatio2))%>%
  mutate(ES=(GeneRatio1/GeneRatio2)/(BgRatio1/BgRatio2))
ego$type <- factor(ego$type,levels=c("Non-pseudogene-derived lncRNA",
                                     "Pseudogene-derived lncRNA"))
ego$ID <- factor(ego$ID,levels=c("429 Ma","319 Ma","160 Ma",
                                 "87 Ma","79 Ma","13 Ma","Mouse"))
#fishertest检验两类lnc指定年龄基因比例的差异p值
pdata <- distinct(ego,ID)
for (i in 1:nrow(pdata)) {
  forid <- pdata$ID[i]
  spglncnum <- filter(ego,type=="Pseudogene-derived lncRNA" & ID==forid)$GeneRatio1
  snpglncnum <- filter(ego,type=="Non-pseudogene-derived lncRNA" & ID==forid)$GeneRatio1
  pdata[i,2] <- fisher.test(data.frame(c(spglncnum,snpglncnum),c(714-spglncnum,28844-snpglncnum)))[["p.value"]]
}
#富集图
pp1 <- ggplot(data = ego,aes(x=type,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:6.5),color="gray")+
  geom_vline(xintercept=c(1.5),color="gray")+
  scale_fill_gradient(low="white",high="#C12039")+
  scale_size(breaks = c(0.5,1,1.5,2,2.5),range = c(1,13))+
  scale_x_discrete(labels=c("Non-pseudogene-\nderived lncRNA",
                            "Pseudogene-derived\nlncRNA"))+
  theme_half_open()+
  labs(x = NULL, y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave("/home/zhzhang/PG/age/Mouse.lnc.enrich.pdf", 
       pp1,width = 4.8, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")




```
