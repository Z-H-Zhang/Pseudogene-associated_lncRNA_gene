### 三.演化相关分析
##### 1.两类lncRNA外显子序列保守性对比
```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的保守性
###lncrna外显子 保守性打分计算
#人类小鼠
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
#提取pas lnc外显子bed,重叠部分merge
cat "/share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.lncRNA_class.txt"|awk -F "\t" '$2=="Pseudogene-associated sense lncRNA"{print $1}' > /share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.paslnc.geneid.txt
cat "/share/home/zhzhang24/PG/lncRNA_class/${i}/${i}.all_lncRNA_exon.merge.bed"|grep -w -Ff /share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.paslnc.geneid.txt |awk '{print $1"\t"$2"\t"$3"\t"$4}' > /share/home/zhzhang24/PG/Evolution/conserve/${i}.paslnc_lncRNA_exon.merge.bed
awk '{print "chr"$0"____"NR}' /share/home/zhzhang24/PG/Evolution/conserve/${i}.paslnc_lncRNA_exon.merge.bed > /share/home/zhzhang24/PG/Evolution/conserve/${i}.paslnc_lncRNA_exon.merge.ucsc.bed
#提取paa lnc外显子bed,重叠部分merge
cat "/share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.lncRNA_class.txt"|awk -F "\t" '$2=="Pseudogene-associated antisense lncRNA"{print $1}' > /share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.paalnc.geneid.txt
cat "/share/home/zhzhang24/PG/lncRNA_class/${i}/${i}.all_lncRNA_exon.merge.bed"|grep -w -Ff /share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.paalnc.geneid.txt |awk '{print $1"\t"$2"\t"$3"\t"$4}' > /share/home/zhzhang24/PG/Evolution/conserve/${i}.paalnc_lncRNA_exon.merge.bed
awk '{print "chr"$0"____"NR}' /share/home/zhzhang24/PG/Evolution/conserve/${i}.paalnc_lncRNA_exon.merge.bed > /share/home/zhzhang24/PG/Evolution/conserve/${i}.paalnc_lncRNA_exon.merge.ucsc.bed
#提取npa lnc外显子bed，重叠部分merge
cat "/share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.lncRNA_class.txt"|awk -F "\t" '$2=="Non-pseudogene-associated lncRNA"{print $1}' > /share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.npalnc.geneid.txt
cat "/share/home/zhzhang24/PG/lncRNA_class/${i}/${i}.all_lncRNA_exon.merge.bed"|grep -w -Ff /share/home/zhzhang24/PG/lncRNA_class/new_r1/${i}.npalnc.geneid.txt |awk '{print $1"\t"$2"\t"$3"\t"$4}' > /share/home/zhzhang24/PG/Evolution/conserve/${i}.npalnc_lncRNA_exon.merge.bed
awk '{print "chr"$0"____"NR}' /share/home/zhzhang24/PG/Evolution/conserve/${i}.npalnc_lncRNA_exon.merge.bed > /share/home/zhzhang24/PG/Evolution/conserve/${i}.npalnc_lncRNA_exon.merge.ucsc.bed
#分别保守性打分
micromamba run -n SEQ bigWigAverageOverBed ${bw} /share/home/zhzhang24/PG/Evolution/conserve/${i}.paslnc_lncRNA_exon.merge.ucsc.bed /share/home/zhzhang24/PG/Evolution/conserve/${i}.paslnc_lncRNA_exon.merge.phastCons.txt
micromamba run -n SEQ bigWigAverageOverBed ${bw} /share/home/zhzhang24/PG/Evolution/conserve/${i}.paalnc_lncRNA_exon.merge.ucsc.bed /share/home/zhzhang24/PG/Evolution/conserve/${i}.paalnc_lncRNA_exon.merge.phastCons.txt
micromamba run -n SEQ bigWigAverageOverBed ${bw} /share/home/zhzhang24/PG/Evolution/conserve/${i}.npalnc_lncRNA_exon.merge.ucsc.bed /share/home/zhzhang24/PG/Evolution/conserve/${i}.npalnc_lncRNA_exon.merge.phastCons.txt
done
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/Evolution/conserve/*_lncRNA_exon.merge.phastCons.txt zhzhang@122.205.95.67:/home/zhzhang/PG/Evolution/conserve/


#假基因区域 保守性打分计算
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
#bed文件
cat "/share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg.bed"|awk '{print "chr"$1"\t"$2"\t"$3"\t"$7}' > /share/home/zhzhang24/PG/Evolution/conserve/${i}.hpg.bed
#保守性打分
micromamba run -n SEQ bigWigAverageOverBed ${bw} /share/home/zhzhang24/PG/Evolution/conserve/${i}.hpg.bed /share/home/zhzhang24/PG/Evolution/conserve/${i}.hpg.phastCons.txt
done
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/Evolution/conserve/*.hpg.phastCons.txt zhzhang@122.205.95.67:/home/zhzhang/PG/Evolution/conserve/


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
#human
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="CDS" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="CDS" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.bed
awk '{print $4}' /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.bed|sort -u|while read ii
do
grep -w "${ii}" /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.bed|bedtools sort|bedtools merge|awk -v name="${ii}" '{print $1"\t"$2"\t"$3"\t"name}' >> /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDSmerge.bed
done
rm /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.bed
awk '{print "chr"$0"____"NR}' /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDSmerge.bed > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDSmerge.ucsc.bed
#mouse
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="CDS" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="CDS" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.bed
awk '{print $4}' /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.bed|sort -u|while read ii
do
grep -w "${ii}" /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.bed|bedtools sort|bedtools merge|awk -v name="${ii}" '{print $1"\t"$2"\t"$3"\t"name}' >> /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDSmerge.bed
done
rm /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.bed
awk '{print "chr"$0"____"NR}' /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDSmerge.bed > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDSmerge.ucsc.bed
#保守性打分
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/hg38.phastCons100way.bw /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDSmerge.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.phastCons.txt
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/mm39.phastCons35way.bw /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDSmerge.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.phastCons.txt



##获取蛋白编码基因外显子5UTR的bed（5UTR所属基因类型在20/22/24列）,重叠部分merge
#human
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="five_prime_utr" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="five_prime_utr" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="five_prime_utr" && $20=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.bed
awk '{print $4}' /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.bed|sort -u|while read ii
do
grep -w "${ii}" /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.bed|bedtools sort|bedtools merge|awk -v name="${ii}" '{print $1"\t"$2"\t"$3"\t"name}' >> /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTRmerge.bed
done
rm /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.bed
awk '{print "chr"$0"____"NR}' /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTRmerge.bed > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTRmerge.ucsc.bed
#mouse
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="five_prime_utr" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="five_prime_utr" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="five_prime_utr" && $20=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.bed
awk '{print $4}' /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.bed|sort -u|while read ii
do
grep -w "${ii}" /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.bed|bedtools sort|bedtools merge|awk -v name="${ii}" '{print $1"\t"$2"\t"$3"\t"name}' >> /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTRmerge.bed
done
rm /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.bed
awk '{print "chr"$0"____"NR}' /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTRmerge.bed > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTRmerge.ucsc.bed
#保守性打分
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/hg38.phastCons100way.bw /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTRmerge.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.phastCons.txt
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/mm39.phastCons35way.bw /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTRmerge.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.phastCons.txt



##获取蛋白编码基因外显子3UTR的bed（3UTR所属基因类型在20/22/24列）,重叠部分merge
#human
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="three_prime_utr" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="three_prime_utr" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="three_prime_utr" && $20=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.bed
awk '{print $4}' /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.bed|sort -u|while read ii
do
grep -w "${ii}" /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.bed|bedtools sort|bedtools merge|awk -v name="${ii}" '{print $1"\t"$2"\t"$3"\t"name}' >> /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTRmerge.bed
done
rm /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.bed
awk '{print "chr"$0"____"NR}' /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTRmerge.bed > /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTRmerge.ucsc.bed
#mouse
tail -n +6 "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" |sed "s/\"//g;s/\;//g"|awk '$3=="three_prime_utr" && $24=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="three_prime_utr" && $22=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7} $3=="three_prime_utr" && $20=="protein_coding" {print $1"\t"$4"\t"$5"\t"$10"\t"$3"\t"$7}' > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.bed
awk '{print $4}' /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.bed|sort -u|while read ii
do
grep -w "${ii}" /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.bed|bedtools sort|bedtools merge|awk -v name="${ii}" '{print $1"\t"$2"\t"$3"\t"name}' >> /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTRmerge.bed
done
rm /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.bed
awk '{print "chr"$0"____"NR}' /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTRmerge.bed > /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTRmerge.ucsc.bed
#保守性打分
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/hg38.phastCons100way.bw /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTRmerge.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.phastCons.txt
bigWigAverageOverBed /home/zhzhang/PG/Evolution/conserve/mm39.phastCons35way.bw /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTRmerge.ucsc.bed /home/zhzhang/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.phastCons.txt


#传至实验室服务器/home/zhzhang/PG/Evolution/conserve
```


```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的保守性
#人类
a <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paalnc_lncRNA_exon.merge.phastCons.txt"
ab <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.npalnc_lncRNA_exon.merge.phastCons.txt"
c <- "~/PG/Evolution/conserve/Homo_sapiens.proteingeneCDS.phastCons.txt"
d <- "~/PG/Evolution/conserve/Homo_sapiens.proteingene5UTR.phastCons.txt"
e <- "~/PG/Evolution/conserve/Homo_sapiens.proteingene3UTR.phastCons.txt"
f <- "~/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.lncRNA_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.lncRNA_PhastCons.pdf"
#小鼠
a <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paalnc_lncRNA_exon.merge.phastCons.txt"
ab <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.npalnc_lncRNA_exon.merge.phastCons.txt"
c <- "~/PG/Evolution/conserve/Mus_musculus.proteingeneCDS.phastCons.txt"
d <- "~/PG/Evolution/conserve/Mus_musculus.proteingene5UTR.phastCons.txt"
e <- "~/PG/Evolution/conserve/Mus_musculus.proteingene3UTR.phastCons.txt"
f <- "~/PG/Evolution/conserve/Mus_musculus.20000random_3kb_intergenic.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.lncRNA_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.lncRNA_PhastCons.pdf"
#导入paslnc外显子的保守性打分
paslnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("pgid"),sep="____")%>%
  group_by(pgid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(4)%>%
  mutate(type="Pseudogene-associated sense lncRNA")
#导入paalnc外显子的保守性打分
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("pgid"),sep="____")%>%
  group_by(pgid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(4)%>%
  mutate(type="Pseudogene-associated antisense lncRNA")
#导入npa lnc外显子的保守性打分
npalnc_lncRNA_exon <- read.delim(ab, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("pgid"),sep="____")%>%
  group_by(pgid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(4)%>%
  mutate(type="Non-pseudogene-associated lncRNA")
#导入蛋白编码基因CDS,5UTR,3UTR保守性打分
proteingeneCDS <- read.delim(c, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("pgid"),sep="____")%>%
  group_by(pgid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(4)%>%
  mutate(type="CDS")
proteingene5UTR <- read.delim(d, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("pgid"),sep="____")%>%
  group_by(pgid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(4)%>%
  mutate(type="5'UTR")
proteingene3UTR <- read.delim(e, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("pgid"),sep="____")%>%
  group_by(pgid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(4)%>%
  mutate(type="3'UTR")
#导入基因间区保守性打分
intergenic <- read.delim(f, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("pgid"),sep="____")%>%
  group_by(pgid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(4)%>%
  mutate(type="Random intergenic")
#合并
alldata <- rbind(paslnc_lncRNA_exon,paalnc_lncRNA_exon,npalnc_lncRNA_exon,proteingeneCDS,proteingene5UTR,proteingene3UTR,intergenic)
alldata$type <- factor(alldata$type,levels = c("CDS","5'UTR","3'UTR",
                                               "Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                                               "Non-pseudogene-associated lncRNA","Random intergenic"))
#统计
tj <- group_by(alldata,type)%>%
  summarise(mean=mean(score),median=median(score))
data.table::fwrite(tj,file =g,sep = '\t',row.names = F,quote = F,col.names = T)
#PLOT
#HUMAN:
  ggsignif::geom_signif(map_signif_level=T,y_position=c(0.4,0.9,0.75,1.1,1,1.2),
                      comparisons = list(c("Random intergenic","Non-pseudogene-associated lncRNA"),
                                         c("Non-pseudogene-associated lncRNA","Pseudogene-associated sense lncRNA"),
                                         c("Non-pseudogene-associated lncRNA","Pseudogene-associated antisense lncRNA"),
                                         c("Pseudogene-associated sense lncRNA","5'UTR"),
                                         c("Pseudogene-associated sense lncRNA","3'UTR"),
                                         c("Pseudogene-associated antisense lncRNA","5'UTR")
                      ))+
#MOUSE:
ggsignif::geom_signif(map_signif_level=T,y_position=c(0.45,1.1,1,1.1,1,1.2,1.3),
                        comparisons = list(c("Random intergenic","Non-pseudogene-associated lncRNA"),
                                           c("Non-pseudogene-associated lncRNA","Pseudogene-associated sense lncRNA"),
                                           c("Non-pseudogene-associated lncRNA","Pseudogene-associated antisense lncRNA"),
                                           c("Pseudogene-associated sense lncRNA","5'UTR"),
                                           c("Pseudogene-associated sense lncRNA","3'UTR"),
                                           c("Pseudogene-associated antisense lncRNA","5'UTR"),
                                           c("Pseudogene-associated antisense lncRNA","3'UTR")
                        ))+
#plot
  p1 <- ggplot(data=alldata, aes(x=type,y=score))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position=c(0.45,1.1,1,1.1,1,1.2,1.3),
                        comparisons = list(c("Random intergenic","Non-pseudogene-associated lncRNA"),
                                           c("Non-pseudogene-associated lncRNA","Pseudogene-associated sense lncRNA"),
                                           c("Non-pseudogene-associated lncRNA","Pseudogene-associated antisense lncRNA"),
                                           c("Pseudogene-associated sense lncRNA","5'UTR"),
                                           c("Pseudogene-associated sense lncRNA","3'UTR"),
                                           c("Pseudogene-associated antisense lncRNA","5'UTR"),
                                           c("Pseudogene-associated antisense lncRNA","3'UTR")
                        ))+
  scale_fill_manual(values=c("#8491B4","#FAA465","#7197AD","#628255","#BB5A5D","#BB5A5D","#BB5A5D"),
                    limits=c("Random intergenic","Non-pseudogene-associated lncRNA",
                             "Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA","CDS","5'UTR","3'UTR"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  scale_x_discrete(labels = c("CDS","5'UTR","3'UTR",
                              "PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA","Random\nintergenic"))+
  labs(y = "PhastCons score", x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave(h,p1,width = 6, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


##### 2.对比不同年龄组中，两类lnc的序列保守性
```r
#对比不同年龄组中两类lncRNA外显子的保守性
#人类
aaa <- "~/PG/age/Human.lncRNA.age.txt"
agelist <- c("429 Ma","319 Ma","160 Ma",
             "87 Ma","29 Ma","Human")
a <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paalnc_lncRNA_exon.merge.phastCons.txt"
ab <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.npalnc_lncRNA_exon.merge.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.lncRNA_age_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.lncRNA_age_PhastCons.pdf"
#小鼠
aaa <- "~/PG/age/Mouse.lncRNA.age.txt"
agelist <- c("429 Ma","319 Ma","160 Ma",
             "87 Ma","79 Ma","13 Ma","Mouse")
a <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paalnc_lncRNA_exon.merge.phastCons.txt"
ab <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.npalnc_lncRNA_exon.merge.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.lncRNA_age_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.lncRNA_age_PhastCons.pdf"
#导入lnc年龄
age <- read.delim(aaa)
#导入paslnc外显子的保守性打分
paslnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,4)%>%
  mutate(type="PAS lncRNA")%>%
  left_join(age,by="geneid")%>%
  filter(max %in% agelist)
#导入paalnc外显子的保守性打分
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,4)%>%
  mutate(type="PAA lncRNA")%>%
  left_join(age,by="geneid")%>%
  filter(max %in% agelist)
#导入npalnc外显子的保守性打分
npalnc_lncRNA_exon <- read.delim(ab, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,4)%>%
  mutate(type="NPA lncRNA")%>%
  left_join(age,by="geneid")%>%
  filter(max %in% agelist)
#合并
alldata <- rbind(paslnc_lncRNA_exon,paalnc_lncRNA_exon,npalnc_lncRNA_exon)
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA","NPA lncRNA"))
alldata$max <- factor(alldata$max,levels=agelist)
#统计
tj <- group_by(alldata,type,max)%>%
  summarise(mean=mean(score),median=median(score))
data.table::fwrite(tj,file =g,sep = '\t',row.names = F,quote = F,col.names = T)
#plot
#查看检验结果
ggplot(data=alldata, aes(x=type,y=score))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,
                        comparisons = list(c("PAS lncRNA","NPA lncRNA"),
                                           c("PAS lncRNA","PAA lncRNA"),
                                           c("PAA lncRNA","NPA lncRNA"))
                        )+
  scale_fill_manual(values=c("#7197AD","#628255","#FAA465"),
                    limits=c("PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  labs(y = "PhastCons score", x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "top") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  facet_wrap(~max,nrow=1)
#human
ggsignif::geom_signif(annotations=c("N.S.","N.S.","N.S.","N.S.","N.S.","***"),
                      y_position=c(0.75,0.93,0.68,0.88,0.83,0.97),tip_length = 0,size=0.5,textsize=4,
                      xmin = c(0.8:5.8),
                      xmax = c(1:6))+
  ggsignif::geom_signif(annotations=c("N.S.","***","***","***","***","***"),
                        y_position=c(0.95,0.85,0.6,0.8,0.5,0.53),tip_length = 0,size=0.5,textsize=4,
                        xmin = c(1:6),
                        xmax = c(1.2:6.2))+
  ggsignif::geom_signif(annotations=c("N.S.","***","***","***","***","***"),
                        y_position=c(1.05,1.05,0.78,0.98,0.93,1.05),tip_length = 0,size=0.5,textsize=4,
                        xmin = c(0.8:5.8),
                        xmax = c(1.2:6.2))+
#mouse
  ggsignif::geom_signif(annotations=c("N.S.","N.S.","N.S.","***","N.S.","N.S.","N.S."),
                        y_position=c(1.13),tip_length = 0,size=0.5,textsize=4,
                        xmin = c(0.8:6.8),
                        xmax = c(1:7))+
  ggsignif::geom_signif(annotations=c("N.S.","***","***","***","***","***","***"),
                        y_position=c(1.03),tip_length = 0,size=0.5,textsize=4,
                        xmin = c(1:7),
                        xmax = c(1.2:7.2))+
  ggsignif::geom_signif(annotations=c("N.S.","***","**","***","***","***","***"),
                        y_position=c(1.23),tip_length = 0,size=0.5,textsize=4,
                        xmin = c(0.8:6.8),
                        xmax = c(1.2:7.2))+
#plot
  p1 <- ggplot(data=alldata, aes(x=max,y=score))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.6,notch=T,aes(fill=type))+
  ggsignif::geom_signif(annotations=c("N.S.","N.S.","N.S.","***","N.S.","N.S.","N.S."),
                        y_position=c(1.13),tip_length = 0,size=0.5,textsize=4,
                        xmin = c(0.8:6.8),
                        xmax = c(1:7))+
  ggsignif::geom_signif(annotations=c("N.S.","***","***","***","***","***","***"),
                        y_position=c(1.03),tip_length = 0,size=0.5,textsize=4,
                        xmin = c(1:7),
                        xmax = c(1.2:7.2))+
  ggsignif::geom_signif(annotations=c("N.S.","***","**","***","***","***","***"),
                        y_position=c(1.23),tip_length = 0,size=0.5,textsize=4,
                        xmin = c(0.8:6.8),
                        xmax = c(1.2:7.2))+
  scale_fill_manual(values=c("#7197AD","#628255","#FAA465"),
                    limits=c("PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  labs(y = "PhastCons score", x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "top") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave(h,p1,width = 6, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 2.5.更严格的阈值对两类palnc的序列保守性的影响
```r
#对比不同组中两类lncRNA外显子的保守性
#人类
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt"
a <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paalnc_lncRNA_exon.merge.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.palncRNA_cutofftype_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paslncRNA_cutofftype_PhastCons.pdf"
i <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paalncRNA_cutofftype_PhastCons.pdf"
#小鼠
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.palncRNA_message.txt"
a <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paalnc_lncRNA_exon.merge.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.palncRNA_cutofftype_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paslncRNA_cutofftype_PhastCons.pdf"
i <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paalncRNA_cutofftype_PhastCons.pdf"
#导入lnc分组
age <- read.delim(aaa)%>%select(2,4)
colnames(age)[1]="geneid"
#导入paslnc外显子的保守性打分
paslnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,4)%>%
  mutate(type="PAS lncRNA")%>%
  left_join(age,by="geneid")
#
pas <- rbind(mutate(filter(paslnc_lncRNA_exon,total_overlap_len>0),cuttype="0"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>50),cuttype="50"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>100),cuttype="100"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>150),cuttype="150"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>200),cuttype="200"))
#导入paalnc外显子的保守性打分
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,4)%>%
  mutate(type="PAA lncRNA")%>%
  left_join(age,by="geneid")
#
paa <- rbind(mutate(filter(paalnc_lncRNA_exon,total_overlap_len>0),cuttype="0"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>50),cuttype="50"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>100),cuttype="100"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>150),cuttype="150"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>200),cuttype="200"))
#合并
alldata <- rbind(pas,paa)
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","50","100","150","200"))
#统计
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(score),median=median(score))
data.table::fwrite(tj,file =g,sep = '\t',row.names = F,quote = F,col.names = T)
#plot
#PAS
p1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=score))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,test.args = c("greater"),
                        comparisons = list(c("200","0")),
                        y_position=c(1),tip_length = 0.01,size=0.5,textsize=4)+#human:y_position=0.9
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  labs(y = "PhastCons score", x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
    theme(axis.title.y = element_text(size = 15),
          axis.title.x= element_text(size = 14),
          axis.text.y  = element_text(size = 13),
          axis.text.x = element_text(size = 12),legend.position = "none")+
    theme(axis.text.x = element_text(angle =45))+
    theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(h,p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=score))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,test.args = c("greater"),
                        comparisons = list(c("200","0")),
                        y_position=c(1),tip_length = 0.01,size=0.5,textsize=4)+#human:y_position=0.9
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  labs(y = "PhastCons score", x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(i,p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")






```
##### 2.6更严格的表达阈值对两类palnc的序列保守性的影响
```r
#计算基因最大表达量
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
                      b = "/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                      c = "Mouse")
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                      c = "Human")
ALLsp_TPM <- rbind(mm_TPM,hs_TPM)%>%select(1,3)
#对比不同组中两类lncRNA外显子的保守性
#人类
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt"
a <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paalnc_lncRNA_exon.merge.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.palncRNA_expcutofftype_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paslncRNA.expcutoff.PhastCons.pdf"
i <- "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.paalncRNA.expcutoff.PhastCons.pdf"
#小鼠
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.palncRNA_message.txt"
a <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paslnc_lncRNA_exon.merge.phastCons.txt"
b <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paalnc_lncRNA_exon.merge.phastCons.txt"
g <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.palncRNA_expcutofftype_PhastCons.tj.txt"
h <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paslncRNA.expcutoff.PhastCons.pdf"
i <- "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.paalncRNA.expcutoff.PhastCons.pdf"
#导入lnc分组
age <- read.delim(aaa)%>%
  select(2)%>%
  dplyr::rename(geneid=lncgid)%>%
  left_join(ALLsp_TPM,by="geneid")
#导入paslnc外显子的保守性打分
paslnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,4)%>%
  mutate(type="PAS lncRNA")%>%
  left_join(age,by="geneid")
#
pas <- rbind(mutate(filter(paslnc_lncRNA_exon,maxTPM>0),cuttype="0"),
             mutate(filter(paslnc_lncRNA_exon,maxTPM>0.5),cuttype="0.5"),
             mutate(filter(paslnc_lncRNA_exon,maxTPM>1),cuttype="1"),
             mutate(filter(paslnc_lncRNA_exon,maxTPM>1.5),cuttype="1.5"),
             mutate(filter(paslnc_lncRNA_exon,maxTPM>2),cuttype="2"))
#导入paalnc外显子的保守性打分
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  filter(V3!=0)%>%
  select(1,3,4)%>%
  separate(V1,c("geneid"),sep="____")%>%
  group_by(geneid)%>%
  summarise(sumscore=sum(V4),sumlen=sum(V3))%>%
  mutate(score=sumscore/sumlen)%>%
  select(1,4)%>%
  mutate(type="PAA lncRNA")%>%
  left_join(age,by="geneid")
#
paa <- rbind(mutate(filter(paalnc_lncRNA_exon,maxTPM>0),cuttype="0"),
             mutate(filter(paalnc_lncRNA_exon,maxTPM>0.5),cuttype="0.5"),
             mutate(filter(paalnc_lncRNA_exon,maxTPM>1),cuttype="1"),
             mutate(filter(paalnc_lncRNA_exon,maxTPM>1.5),cuttype="1.5"),
             mutate(filter(paalnc_lncRNA_exon,maxTPM>2),cuttype="2"))
#合并
alldata <- rbind(pas,paa)
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","0.5","1","1.5","2"))
#统计
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(score),median=median(score))
data.table::fwrite(tj,file =g,sep = '\t',row.names = F,quote = F,col.names = T)
#plot
#PAS
p1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=score))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("2","0")),
                        y_position=c(1),tip_length = 0.01,size=0.5,textsize=4)+#human:y_position=0.9
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  labs(y = "PhastCons score", x ="Cutoff of\nexpression (TPM)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(h,p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=score))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("2","0")),
                        y_position=c(1),tip_length = 0.01,size=0.5,textsize=4)+#human:y_position=0.9
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = c("0","0.25","0.5","0.75","1"))+
  labs(y = "PhastCons score", x ="Cutoff of\nexpression (TPM)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(i,p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


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
#斑马鱼
tail -n +6 "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Danio_rerio.GRCz11.108.chr.rmpg.novellncRNA.gtf" |grep -w exon|awk '$13=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $11=="transcript_id" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /share/home/zhzhang24/PG/Evolution/tran_len/Danio_rerio.allgeneexon.bed


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
                  b="~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                  c="Human")
#小鼠
mmtlen <- gettlen("~/PG/Evolution/tran_len/Mus_musculus.allgeneexon.bed",
                  "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                  "Mouse")
#合并不同物种信息
allsp_tlen <- rbind(hstlen,mmtlen)
allsp_tlen$type <- factor(allsp_tlen$type,
                          levels = c("Protein-coding",
                                     "Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                                     "Non-pseudogene-associated lncRNA"))
#统计储存
tj <- group_by(allsp_tlen,sp,type)%>%
  summarise(mean=mean(len),median=median(len))
data.table::fwrite(tj,file ="/home/zhzhang/PG/Evolution/tran_len/ALLSP_3gene_transcript_len.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
pmm <- ggplot(data = filter(allsp_tlen,sp=="Mouse"),aes(x=type,y=log2(len)))+
  geom_violin(width=0.9,aes(fill=type,color=type))+
  geom_boxplot(fatten = 3,width=0.1,outlier.alpha = 0)+
  ggsignif::geom_signif(map_signif_level=T,y_position = c(19,18,17),
                        comparisons=list(c("Protein-coding","Pseudogene-associated sense lncRNA"),
                                         c("Pseudogene-associated sense lncRNA","Non-pseudogene-associated lncRNA"),
                                         c("Pseudogene-associated antisense lncRNA","Non-pseudogene-associated lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA","Random intergenic"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-associated lncRNA","Random intergenic"))+
  theme_half_open()+
  #coord_cartesian(ylim = c(0.5, 6))+
  #scale_y_continuous(breaks = c(1,2,3,4,5),labels = c("1","2","3","4","5"))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  labs(x = NULL, y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[2]*"("*"b"*"p"*")"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/Evolution/tran_len/Mus_musculus.lnc_tranlen.pdf", 
       pmm,width = 4, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
phs <- ggplot(data = filter(allsp_tlen,sp=="Human"),aes(x=type,y=log2(len)))+
  geom_violin(width=0.9,aes(fill=type,color=type))+
  geom_boxplot(fatten = 1,width=0.1,outlier.alpha = 0)+
  ggsignif::geom_signif(map_signif_level=T,y_position = c(17,19.5,18.5),
                        comparisons=list(c("Protein-coding","Pseudogene-associated sense lncRNA"),
                                         c("Pseudogene-associated sense lncRNA","Non-pseudogene-associated lncRNA"),
                                         c("Pseudogene-associated antisense lncRNA","Non-pseudogene-associated lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA","Random intergenic"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-associated lncRNA","Random intergenic"))+
  theme_half_open()+
  #coord_cartesian(ylim = c(1, 6))+
  #scale_y_continuous(breaks = c(1,2,3,4,5),labels = c("1","2","3","4","5"))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  labs(x = NULL, y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[2]*"("*"b"*"p"*")"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/Evolution/tran_len/Homo_sapiens.lnc_tranlen.pdf", 
       phs,width = 4, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


##### 3.5 更严格阈值对palnc转录本长度的影响
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
  geneid_class <- read.delim(b)%>%select(2,4,7)%>%
    dplyr::rename(geneid=lncgid)
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
#人
tlen <- gettlen(a="~/PG/Evolution/tran_len/Homo_sapiens.allgeneexon.bed",
                  b="/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt",
                  c="Human")
h="/home/zhzhang/PG/Evolution/tran_len/Homo_sapiens.paslnc_tranlen.pdf"
i="/home/zhzhang/PG/Evolution/tran_len/Homo_sapiens.paalnc_tranlen.pdf"
#小鼠
tlen <- gettlen("~/PG/Evolution/tran_len/Mus_musculus.allgeneexon.bed",
                "/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.palncRNA_message.txt",
                  "Mouse")
h="/home/zhzhang/PG/Evolution/tran_len/Mus_musculus.paslnc_tranlen.pdf"
i="/home/zhzhang/PG/Evolution/tran_len/Mus_musculus.paalnc_tranlen.pdf"
#不同物种信息
#
pas <- rbind(mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>200),cuttype="200"))
paa <- rbind(mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>200),cuttype="200"))
#合并
alldata <- rbind(pas,paa)%>%
  mutate(type=case_when(type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA"))
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","50","100","150","200"))
#统计储存
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(len),median=median(len))
#PAS
p1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=log2(len)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=log2(mean)),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=log2(mean),group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),
                        y_position=c(14.5),tip_length = 0.01,size=0.5,textsize=4)+
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  theme_half_open()+
  coord_cartesian(ylim = c(7.5, 15))+#human:ylim = c(7, 15)
  labs(y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[2]*"("*"b"*"p"*")"*")"),
       x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(h,p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=log2(len)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=log2(mean)),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=log2(mean),group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),
                        y_position=c(14.5),tip_length = 0.01,size=0.5,textsize=4)+#human:y_position=0.9
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  theme_half_open()+
  coord_cartesian(ylim = c(7.5, 15))+ #human c(6.6,15)
  labs(y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[2]*"("*"b"*"p"*")"*")"),
       x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(i,p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
##### 3.6更严格表达阈值对palnc转录本长度的影响
```r
#计算基因最大表达量
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
                      b = "/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                      c = "Mouse")
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                      c = "Human")
ALLsp_TPM <- rbind(mm_TPM,hs_TPM)%>%select(1,3)
#对比两类lnc转录本长度
#函数计算三类基因每个转录本的长度
#a输入全部基因外显子bed，b输入基因分类文件,c输入物种名
gettlen <- function(a,b,c){
  #导入外显子bed文件
  allgeneexon <- read.delim(a, header=FALSE)%>%
    select(-1,-6)
  colnames(allgeneexon) <- c("start","end","geneid","transcriptid")
  #导入基因分类文件
  geneid_class <- read.delim(b)%>%select(2,4,7)%>%
    dplyr::rename(geneid=lncgid)
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
#人
tlen <- gettlen(a="~/PG/Evolution/tran_len/Homo_sapiens.allgeneexon.bed",
                b="/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt",
                c="Human")%>%
  left_join(ALLsp_TPM,by="geneid")
h="/home/zhzhang/PG/Evolution/tran_len/Homo_sapiens.paslnc.expcutoff.tranlen.pdf"
i="/home/zhzhang/PG/Evolution/tran_len/Homo_sapiens.paalnc.expcutoff.tranlen.pdf"
#小鼠
tlen <- gettlen("~/PG/Evolution/tran_len/Mus_musculus.allgeneexon.bed",
                "/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.palncRNA_message.txt",
                "Mouse")%>%
  left_join(ALLsp_TPM,by="geneid")
h="/home/zhzhang/PG/Evolution/tran_len/Mus_musculus.paslnc.expcutoff.tranlen.pdf"
i="/home/zhzhang/PG/Evolution/tran_len/Mus_musculus.paalnc.expcutoff.tranlen.pdf"
#不同物种信息
#
pas <- rbind(mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&maxTPM>0),cuttype="0"),
             mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&maxTPM>0.5),cuttype="0.5"),
             mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&maxTPM>1),cuttype="1"),
             mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&maxTPM>1.5),cuttype="1.5"),
             mutate(filter(tlen,type=="Pseudogene-associated sense lncRNA"&maxTPM>2),cuttype="2"))
paa <- rbind(mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&maxTPM>0),cuttype="0"),
             mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&maxTPM>0.5),cuttype="0.5"),
             mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&maxTPM>1),cuttype="1"),
             mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&maxTPM>1.5),cuttype="1.5"),
             mutate(filter(tlen,type=="Pseudogene-associated antisense lncRNA"&maxTPM>2),cuttype="2"))
#合并
alldata <- rbind(pas,paa)%>%
  mutate(type=case_when(type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA"))
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","0.5","1","1.5","2"))
#统计储存
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(len),median=median(len))
#PAS
p1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=log2(len)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=log2(mean)),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=log2(mean),group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("2","0")),
                        y_position=c(14.5),tip_length = 0.01,size=0.5,textsize=4)+
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  theme_half_open()+
  coord_cartesian(ylim = c(7.5, 15))+#human:ylim = c(7, 15)
  labs(y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[2]*"("*"b"*"p"*")"*")"),
       x ="Cutoff of\nexpression (TPM)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(h,p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=log2(len)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=log2(mean)),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=log2(mean),group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("2","0")),
                        y_position=c(14.4),tip_length = 0.01,size=0.5,textsize=4)+#human:y_position=0.9
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  theme_half_open()+
  coord_cartesian(ylim = c(7.5, 15))+ #human c(6.6,15)
  labs(y =expression("T"*"r"*"a"*"n"*"s"*"c"*"r"*"i"*"p"*"t"~"l"*"e"*"n"*"g"*"t"*"h"~"("*"l"*"o"*"g"[2]*"("*"b"*"p"*")"*")"),
       x ="Cutoff of\nexpression (TPM)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(i,p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



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
                                     "Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                                     "Non-pseudogene-associated lncRNA"))
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
tj$type <- factor(tj$type,levels = c("Protein-coding",
                                     "Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                                     "Non-pseudogene-associated lncRNA"))
t.test(filter(allsp_ten,sp=="Mouse" & type=="Protein-coding")$exonnum,
       filter(allsp_ten,sp=="Mouse" & type=="Pseudogene-associated sense lncRNA")$exonnum)
t.test(filter(allsp_ten,sp=="Mouse" & type=="Pseudogene-associated sense lncRNA")$exonnum,
       filter(allsp_ten,sp=="Mouse" & type=="Non-pseudogene-associated lncRNA")$exonnum)
t.test(filter(allsp_ten,sp=="Mouse" & type=="Pseudogene-associated antisense lncRNA")$exonnum,
       filter(allsp_ten,sp=="Mouse" & type=="Non-pseudogene-associated lncRNA")$exonnum)
pmm <- ggplot(data = filter(tj,sp=="Mouse"),aes(x=type,y=mean))+
  geom_point(size=2,aes(fill=type,color=type))+
  geom_errorbar(aes(ymin=confmin,ymax=confmax,color=type),width=0.15)+
  ggsignif::geom_signif(annotations=c("***","***","***"),y_position=c(7.8,4.5,4),
              xmin = c(1,2,3),xmax = c(2,4,4),tip_length = 0.01)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA","Random intergenic"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-associated lncRNA","Random intergenic"))+
  theme_half_open()+
  scale_y_continuous(limits = c(2.5,8),breaks = c(0,2.5,5,7.5))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  labs(x = NULL, y ="Transcript exon number",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/Evolution/tran_en/Mus_musculus.lnc_tranexonn.pdf", 
       pmm,width = 4, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
t.test(filter(allsp_ten,sp=="Human" & type=="Protein-coding")$exonnum,
       filter(allsp_ten,sp=="Human" & type=="Pseudogene-associated sense lncRNA")$exonnum)
t.test(filter(allsp_ten,sp=="Human" & type=="Pseudogene-associated sense lncRNA")$exonnum,
       filter(allsp_ten,sp=="Human" & type=="Non-pseudogene-associated lncRNA")$exonnum)
t.test(filter(allsp_ten,sp=="Human" & type=="Pseudogene-associated antisense lncRNA")$exonnum,
       filter(allsp_ten,sp=="Human" & type=="Non-pseudogene-associated lncRNA")$exonnum)
phs <- ggplot(data = filter(tj,sp=="Human"),aes(x=type,y=mean))+
  geom_point(size=2,aes(fill=type,color=type))+
  geom_errorbar(aes(ymin=confmin,ymax=confmax,color=type),width=0.15)+
  ggsignif::geom_signif(annotations=c("***","***","***"),y_position=c(8.5,5.7,4.5),
              xmin = c(1,2,3),xmax = c(2,4,4),tip_length = 0.01)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA","Random intergenic"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-associated lncRNA","Random intergenic"))+
  theme_half_open()+
  scale_y_continuous(limits = c(2.5,9),breaks = c(0,3,6,9))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  labs(x = NULL, y ="Transcript exon number",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =35))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/Evolution/tran_en/Homo_sapiens.lnc_tranexonn.pdf", 
       phs,width = 4, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


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
#从ucsc获取CPG岛区域（防止CPG岛区域高突变率的影响）
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz -P /home/zhzhang/PG/Evolution/snp/
gzip -d /home/zhzhang/PG/Evolution/snp/cpgIslandExt.txt.gz
awk '$2!~/_/{print $2"\t"$3"\t"$4}' /home/zhzhang/PG/Evolution/snp/cpgIslandExt.txt |sed 's/chr//g' > /home/zhzhang/PG/Evolution/snp/cpgIsland.bed

#六种区域去除cpg岛区域
bedtools subtract -a /share/home/zhzhang24/PG/Evolution/conserve/Homo_sapiens.paslnc_lncRNA_exon.merge.bed -b /share/home/zhzhang24/PG/Evolution/snp/cpgIsland.bed > /share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paslnc_lncRNA_mergeexon.nocpg.bed
bedtools subtract -a /share/home/zhzhang24/PG/Evolution/conserve/Homo_sapiens.paalnc_lncRNA_exon.merge.bed -b /share/home/zhzhang24/PG/Evolution/snp/cpgIsland.bed > /share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paalnc_lncRNA_mergeexon.nocpg.bed
bedtools subtract -a /share/home/zhzhang24/PG/Evolution/conserve/Homo_sapiens.npalnc_lncRNA_exon.merge.bed -b /share/home/zhzhang24/PG/Evolution/snp/cpgIsland.bed > /share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.npalnc_lncRNA_mergeexon.nocpg.bed
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed|awk '{print $1"\t"$2"\t"$3"\t"$4}'|bedtools subtract -a - -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.nocpg.bed
bedtools subtract -a /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDSmerge.bed -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.nocpg.bed
bedtools subtract -a /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTRmerge.bed -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.nocpg.bed
bedtools subtract -a /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTRmerge.bed -b /home/zhzhang/PG/Evolution/snp/cpgIsland.bed > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.nocpg.bed




#获取六种区域中，SNP的数量
bedtools intersect -a /share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paslnc_lncRNA_mergeexon.nocpg.bed -b /share/home/zhzhang24/PG/Evolution/snp/common_allsnp.bed -c > /share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paslnc_lncRNA_mergeexon.SNPnum.txt
bedtools intersect -a /share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paalnc_lncRNA_mergeexon.nocpg.bed -b /share/home/zhzhang24/PG/Evolution/snp/common_allsnp.bed -c > /share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paalnc_lncRNA_mergeexon.SNPnum.txt
bedtools intersect -a /share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.npalnc_lncRNA_mergeexon.nocpg.bed -b /share/home/zhzhang24/PG/Evolution/snp/common_allsnp.bed -c > /share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.npalnc_lncRNA_mergeexon.SNPnum.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.SNPnum.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.SNPnum.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.SNPnum.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.nocpg.bed -b /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -c > /home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.SNPnum.txt

#传至实验室服务器/home/zhzhang/PG/Evolution/snp


```
```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的SNP密度
#人类commonSNP
a <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paslnc_lncRNA_mergeexon.SNPnum.txt"
b <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paalnc_lncRNA_mergeexon.SNPnum.txt"
ab <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.npalnc_lncRNA_mergeexon.SNPnum.txt"
c <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.SNPnum.txt"
d <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.SNPnum.txt"
e <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.SNPnum.txt"
f <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.SNPnum.txt"
g <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.lncRNA_SNPdensity.tj.txt"
h <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.lncRNA_SNPdensity.pdf"

#导入paslnc外显子的SNP数量
paslnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="PAS lncRNA")
colnames(paslnc_lncRNA_exon)[c(1,2)] <- c("SNPnum","SNPdensity")
#导入paalnc外显子的SNP数量
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="PAA lncRNA")
colnames(paalnc_lncRNA_exon)[c(1,2)] <- c("SNPnum","SNPdensity")
#导入npgdlnc外显子的SNP数量
npalnc_lncRNA_exon <- read.delim(ab, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="NPA lncRNA")
colnames(npalnc_lncRNA_exon)[c(1,2)] <- c("SNPnum","SNPdensity")
#导入蛋白编码基因CDS,5UTR,3UTR的SNP数量
proteingeneCDS <- read.delim(c, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="CDS")
colnames(proteingeneCDS) <- colnames(paslnc_lncRNA_exon)
proteingene5UTR <- read.delim(d, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="5'UTR")
colnames(proteingene5UTR) <- colnames(paslnc_lncRNA_exon)
proteingene3UTR <- read.delim(e, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="3'UTR")
colnames(proteingene3UTR) <- colnames(paslnc_lncRNA_exon)
#导入基因间区的SNP数量
intergenic <- read.delim(f, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(5,6)%>%
  mutate(type="Random intergenic")
colnames(intergenic) <- colnames(paslnc_lncRNA_exon)
#合并
alldata <- rbind(paslnc_lncRNA_exon,paalnc_lncRNA_exon,npalnc_lncRNA_exon,proteingeneCDS,proteingene5UTR,proteingene3UTR,intergenic)
alldata$type <- factor(alldata$type,levels = c("CDS","5'UTR","3'UTR",
                                               "PAS lncRNA","PAA lncRNA",
                                               "NPA lncRNA","Random intergenic"))
#统计
tj <- group_by(alldata,type)%>%
  summarise(mean=mean(SNPdensity),median=median(SNPdensity),
            meannum=mean(SNPnum),mediannum=median(SNPnum))
data.table::fwrite(tj,file =g,sep = '\t',row.names = F,quote = F,col.names = T)
#plot
p1 <- ggplot(data=alldata, aes(x=type,y=log2(SNPdensity+1)))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position=c(6.7,8.5,7.7,9.5,10.2),
              comparisons = list(c("Random intergenic","NPA lncRNA"),
                                 c("NPA lncRNA","PAA lncRNA"),
                                 c("PAS lncRNA","PAA lncRNA"),
                                 c("PAS lncRNA","NPA lncRNA"),
                                 c("PAS lncRNA","3'UTR")))+
  scale_fill_manual(values=c("#8491B4","#FAA465","#7197AD","#628255","#BB5A5D","#BB5A5D","#BB5A5D"),
                    limits=c("Random intergenic","NPA lncRNA",
                             "PAS lncRNA","PAA lncRNA","CDS","5'UTR","3'UTR"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,3,6,9),labels = c("0","3","6","9"))+
  scale_x_discrete(labels = c("CDS","5'UTR","3'UTR","PAS lncRNA","PAA lncRNA","NPA lncRNA","Random\nintergenic"))+
  labs(y =expression("l"*"o"*"g"[2]*"("*"S"*"N"*"P"~"d"*"e"*"n"*"s"*"i"*"t"*"y"*")"),
       x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave(h,p1,width = 6, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 5.5 更严格阈值对SNP密度影响
```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的SNP密度
#人类commonSNP
a <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paslnc_lncRNA_mergeexon.SNPnum.txt"
b <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paalnc_lncRNA_mergeexon.SNPnum.txt"
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt"
h <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paslncRNA_SNPdensity.pdf"
i <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paalncRNA_SNPdensity.pdf"
#导入lnc分组
age <- read.delim(aaa)%>%select(2,4)
colnames(age)[1]="geneid"
#导入paslnc外显子的SNP数量
paslnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(4,6)%>%
  mutate(type="PAS lncRNA")%>%
  dplyr::rename(geneid=V4)%>%
  left_join(age,by="geneid")
#导入paalnc外显子的SNP数量
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(4,6)%>%
  mutate(type="PAA lncRNA")%>%
  dplyr::rename(geneid=V4)%>%
  left_join(age,by="geneid")
#
pas <- rbind(mutate(filter(paslnc_lncRNA_exon,total_overlap_len>0),cuttype="0"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>50),cuttype="50"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>100),cuttype="100"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>150),cuttype="150"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>200),cuttype="200"))
#
paa <- rbind(mutate(filter(paalnc_lncRNA_exon,total_overlap_len>0),cuttype="0"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>50),cuttype="50"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>100),cuttype="100"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>150),cuttype="150"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>200),cuttype="200"))
#合并
alldata <- rbind(pas,paa)
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","50","100","150","200"))
#统计
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(38)
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(SNPd),median=median(SNPd),sd=sd(SNPd),
            confmin=boot::boot.ci(boot::boot(SNPd, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(SNPd, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][5])
#plot
#PAS
p1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=log2(SNPd+1)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),
                        y_position=c(8),tip_length = 0.01,size=0.5,textsize=4)+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=log2(mean),fill=cuttype),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=log2(mean),group=type),linetype="dashed",size=0.3)+
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  theme_half_open()+
  labs(y =expression("l"*"o"*"g"[2]*"("*"S"*"N"*"P"~"d"*"e"*"n"*"s"*"i"*"t"*"y"*")"), x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(h,p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=log2(SNPd+1)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("200","0")),
                        y_position=c(8),tip_length = 0.01,size=0.5,textsize=4)+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=log2(mean),fill=cuttype),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=log2(mean),group=type),linetype="dashed",size=0.3)+
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  theme_half_open()+
  labs(y =expression("l"*"o"*"g"[2]*"("*"S"*"N"*"P"~"d"*"e"*"n"*"s"*"i"*"t"*"y"*")"), x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(i,p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")




#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的SNP密度
#人类commonSNP
a <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paslnc_lncRNA_mergeexon.SNPnum.txt"
b <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paalnc_lncRNA_mergeexon.SNPnum.txt"
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt"
h <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paslncRNA_SNPdensity.pdf"
i <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paalncRNA_SNPdensity.pdf"
#导入lnc分组
age <- read.delim(aaa)%>%select(2,4)
colnames(age)[1]="geneid"
#导入paslnc外显子的SNP数量
paslnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(4,6)%>%
  mutate(type="PAS lncRNA")%>%
  dplyr::rename(geneid=V4)%>%
  left_join(age,by="geneid")
#导入paalnc外显子的SNP数量
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(4,6)%>%
  mutate(type="PAA lncRNA")%>%
  dplyr::rename(geneid=V4)%>%
  left_join(age,by="geneid")
#
pas <- rbind(mutate(filter(paslnc_lncRNA_exon,total_overlap_len>0),cuttype="0"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>50),cuttype="50"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>100),cuttype="100"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>150),cuttype="150"),
             mutate(filter(paslnc_lncRNA_exon,total_overlap_len>200),cuttype="200"))
#
paa <- rbind(mutate(filter(paalnc_lncRNA_exon,total_overlap_len>0),cuttype="0"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>50),cuttype="50"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>100),cuttype="100"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>150),cuttype="150"),
             mutate(filter(paalnc_lncRNA_exon,total_overlap_len>200),cuttype="200"))
#合并
alldata <- rbind(pas,paa)
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","50","100","150","200"))
#统计
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(38)
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(SNPd),median=median(SNPd),sd=sd(SNPd),
            confmin=boot::boot.ci(boot::boot(SNPd, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(SNPd, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][5])
#plot
#PAS
p1 <- ggplot(data=filter(tj,type=="PAS lncRNA"), aes(x=cuttype,y=mean))+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean,fill=cuttype),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  geom_errorbar(aes(ymin=confmin,ymax=confmax,color=cuttype),width=0.15)+
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  scale_color_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  theme_half_open()+
  labs(y = "SNP density", x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(h,p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(tj,type=="PAA lncRNA"), aes(x=cuttype,y=mean))+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean,fill=cuttype),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  geom_errorbar(aes(ymin=confmin,ymax=confmax,color=cuttype),width=0.15)+
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  scale_color_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  theme_half_open()+
  labs(y = "SNP density", x ="Cutoff of\noverlap length (bp)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(i,p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 5.6更严格表达阈值对SNP密度影响
```r
#计算基因最大表达量
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
                      b = "/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                      c = "Mouse")
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                      c = "Human")
ALLsp_TPM <- rbind(mm_TPM,hs_TPM)%>%select(1,3)
#人类commonSNP
a <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paslnc_lncRNA_mergeexon.SNPnum.txt"
b <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paalnc_lncRNA_mergeexon.SNPnum.txt"
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt"
h <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paslncRNA.expcutoff.SNPdensity.pdf"
i <- "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.paalncRNA.expcutoff.SNPdensity.pdf"
#导入lnc分组
age <- read.delim(aaa)%>%
  select(2)%>%
  dplyr::rename(geneid=lncgid)%>%
  left_join(ALLsp_TPM,by="geneid")
#导入paslnc外显子的SNP数量
paslnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(4,6)%>%
  mutate(type="PAS lncRNA")%>%
  dplyr::rename(geneid=V4)%>%
  left_join(age,by="geneid")
#导入paalnc外显子的SNP数量
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  mutate(SNPd=V5/((V3-V2+1)/1000))%>%
  select(4,6)%>%
  mutate(type="PAA lncRNA")%>%
  dplyr::rename(geneid=V4)%>%
  left_join(age,by="geneid")
#
pas <- rbind(mutate(filter(paslnc_lncRNA_exon,maxTPM>0),cuttype="0"),
             mutate(filter(paslnc_lncRNA_exon,maxTPM>0.5),cuttype="0.5"),
             mutate(filter(paslnc_lncRNA_exon,maxTPM>1),cuttype="1"),
             mutate(filter(paslnc_lncRNA_exon,maxTPM>1.5),cuttype="1.5"),
             mutate(filter(paslnc_lncRNA_exon,maxTPM>2),cuttype="2"))
#
paa <- rbind(mutate(filter(paalnc_lncRNA_exon,maxTPM>0),cuttype="0"),
             mutate(filter(paalnc_lncRNA_exon,maxTPM>0.5),cuttype="0.5"),
             mutate(filter(paalnc_lncRNA_exon,maxTPM>1),cuttype="1"),
             mutate(filter(paalnc_lncRNA_exon,maxTPM>1.5),cuttype="1.5"),
             mutate(filter(paalnc_lncRNA_exon,maxTPM>2),cuttype="2"))
#合并
alldata <- rbind(pas,paa)
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","0.5","1","1.5","2"))
#统计
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(38)
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(SNPd),median=median(SNPd),sd=sd(SNPd),
            confmin=boot::boot.ci(boot::boot(SNPd, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(SNPd, mean_forboot, R = 100),conf=0.95,type=c('perc'))[["percent"]][5])
#plot
#PAS
p1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=log2(SNPd+1)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("2","0")),
                        y_position=c(8),tip_length = 0.01,size=0.5,textsize=4)+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=log2(mean),fill=cuttype),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=log2(mean),group=type),linetype="dashed",size=0.3)+
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  theme_half_open()+
  labs(y =expression("l"*"o"*"g"[2]*"("*"S"*"N"*"P"~"d"*"e"*"n"*"s"*"i"*"t"*"y"*")"), x ="Cutoff of\nexpression (TPM)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(h,p1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#PAA
p2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=log2(SNPd+1)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  ggsignif::geom_signif(map_signif_level = T,
                        comparisons = list(c("2","0")),
                        y_position=c(8),tip_length = 0.01,size=0.5,textsize=4)+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=log2(mean),fill=cuttype),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=log2(mean),group=type),linetype="dashed",size=0.3)+
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  theme_half_open()+
  labs(y =expression("l"*"o"*"g"[2]*"("*"S"*"N"*"P"~"d"*"e"*"n"*"s"*"i"*"t"*"y"*")"), x ="Cutoff of\nexpression (TPM)",fill = NULL,color = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(i,p2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```




##### 6.CADD评分 SNP
```r
#CADD数据库下载SNP的CADD分数https://cadd.gs.washington.edu/download
wget https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.snv.tsv.gz -P /home/zhzhang/PG/Evolution/CADD/
#产生positionID（chr_position）
tail -n +3 "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.r3.0.snv.tsv"|awk '{print $0"\t"$1"_"$2}' > /home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt

#产生common SNP的 SNPID和positionID对照表
cat "/home/zhzhang/PG/Evolution/snp/common_allsnp.bed"|awk '{print $4"\t"$1"_"$3}' > /home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt

#获取六种区域中，SNP的ID
bedtools intersect -a /share/home/zhzhang24/PG/Evolution/snp/common_allsnp.bed -b "/share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paslnc_lncRNA_mergeexon.nocpg.bed" -wa |awk '{print $4}'|sort|uniq > /share/home/zhzhang24/PG/Evolution/GERP/Homo_sapiens.paslnc_lncRNA_mergeexon.common_SNPID.txt
bedtools intersect -a /share/home/zhzhang24/PG/Evolution/snp/common_allsnp.bed -b "/share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paalnc_lncRNA_mergeexon.nocpg.bed" -wa |awk '{print $4}'|sort|uniq > /share/home/zhzhang24/PG/Evolution/GERP/Homo_sapiens.paalnc_lncRNA_mergeexon.common_SNPID.txt
bedtools intersect -a /share/home/zhzhang24/PG/Evolution/snp/common_allsnp.bed -b "/share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.npalnc_lncRNA_mergeexon.nocpg.bed" -wa |awk '{print $4}'|sort|uniq > /share/home/zhzhang24/PG/Evolution/GERP/Homo_sapiens.npalnc_lncRNA_mergeexon.common_SNPID.txt
cat /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed|awk '{print $1"\t"$2"\t"$3"\t"$4}'|bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b - -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.20000random_3kb_intergenic.common_SNPID.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingeneCDSmerge.bed" -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingeneCDS.common_SNPID.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene5UTRmerge.bed" -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingene5UTR.common_SNPID.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/common_allsnp.bed -b "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.proteingene3UTRmerge.bed" -wa |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingene3UTR.common_SNPID.txt


#根据六种区域内SNPID提取对应positionID
grep -w -Ff /share/home/zhzhang24/PG/Evolution/GERP/Homo_sapiens.paslnc_lncRNA_mergeexon.common_SNPID.txt "/share/home/zhzhang24/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /share/home/zhzhang24/PG/Evolution/CADD/Homo_sapiens.paslnc_lncRNA_mergeexon.common_SNP.positionID.txt
grep -w -Ff /share/home/zhzhang24/PG/Evolution/GERP/Homo_sapiens.paalnc_lncRNA_mergeexon.common_SNPID.txt "/share/home/zhzhang24/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /share/home/zhzhang24/PG/Evolution/CADD/Homo_sapiens.paalnc_lncRNA_mergeexon.common_SNP.positionID.txt
grep -w -Ff /share/home/zhzhang24/PG/Evolution/GERP/Homo_sapiens.npalnc_lncRNA_mergeexon.common_SNPID.txt "/share/home/zhzhang24/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /share/home/zhzhang24/PG/Evolution/CADD/Homo_sapiens.npalnc_lncRNA_mergeexon.common_SNP.positionID.txt
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.20000random_3kb_intergenic.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.20000random_3kb_intergenic.common_SNP.positionID.txt
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingeneCDS.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingeneCDS.common_SNP.positionID.txt
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingene5UTR.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene5UTR.common_SNP.positionID.txt
grep -w -Ff /home/zhzhang/PG/Evolution/GERP/Homo_sapiens.proteingene3UTR.common_SNPID.txt "/home/zhzhang/PG/Evolution/CADD/common_allsnp.SNPID.positionID.txt"|awk '{print $2}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene3UTR.common_SNP.positionID.txt


#根据六种区域内common SNP的positionID提取CADD分数
grep -w -Ff /share/home/zhzhang24/PG/Evolution/CADD/Homo_sapiens.paslnc_lncRNA_mergeexon.common_SNP.positionID.txt "/share/home/zhzhang24/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /share/home/zhzhang24/PG/Evolution/CADD/Homo_sapiens.paslnc_lncRNA_mergeexon.common_SNP.CADD.txt
grep -w -Ff /share/home/zhzhang24/PG/Evolution/CADD/Homo_sapiens.paalnc_lncRNA_mergeexon.common_SNP.positionID.txt "/share/home/zhzhang24/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /share/home/zhzhang24/PG/Evolution/CADD/Homo_sapiens.paalnc_lncRNA_mergeexon.common_SNP.CADD.txt
grep -w -Ff /share/home/zhzhang24/PG/Evolution/CADD/Homo_sapiens.npalnc_lncRNA_mergeexon.common_SNP.positionID.txt "/share/home/zhzhang24/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /share/home/zhzhang24/PG/Evolution/CADD/Homo_sapiens.npalnc_lncRNA_mergeexon.common_SNP.CADD.txt
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.20000random_3kb_intergenic.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.20000random_3kb_intergenic.common_SNP.CADD.txt
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingeneCDS.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingeneCDS.common_SNP.CADD.txt
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene5UTR.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene5UTR.common_SNP.CADD.txt
grep -w -Ff /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene3UTR.common_SNP.positionID.txt "/home/zhzhang/PG/Evolution/CADD/gnomad.genomes.snv.positionID.txt"|awk '{print $0}' > /home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene3UTR.common_SNP.CADD.txt

```
```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的SNP CADD
#人类commonSNP
a <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.paslnc_lncRNA_mergeexon.common_SNP.CADD.txt"
b <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.paalnc_lncRNA_mergeexon.common_SNP.CADD.txt"
ab <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.npalnc_lncRNA_mergeexon.common_SNP.CADD.txt"
c <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingeneCDS.common_SNP.CADD.txt"
d <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene5UTR.common_SNP.CADD.txt"
e <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.proteingene3UTR.common_SNP.CADD.txt"
f <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.20000random_3kb_intergenic.common_SNP.CADD.txt"
g <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.SNPCADD.tj.txt"
h <- "/home/zhzhang/PG/Evolution/CADD/Homo_sapiens.SNPCADD.pdf"
#导入pgdlnc外显子的SNP的CADD分数
pgdlnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  group_by(V7)%>%
  top_n(1,V6)%>%
  select(6)%>%
  mutate(type="PAS lncRNA")
colnames(pgdlnc_lncRNA_exon) <- c("pID","CADD","type")
#导入pgdlnc外显子的SNP的CADD分数
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  group_by(V7)%>%
  top_n(1,V6)%>%
  select(6)%>%
  mutate(type="PAA lncRNA")
colnames(paalnc_lncRNA_exon) <- c("pID","CADD","type")
#导入npgdlnc外显子的SNP
npgdlnc_lncRNA_exon <- read.delim(ab, header=FALSE)%>%
  group_by(V7)%>%
  top_n(1,V6)%>%
  select(6)%>%
  mutate(type="NPA lncRNA")
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
alldata <- rbind(pgdlnc_lncRNA_exon,paalnc_lncRNA_exon,npgdlnc_lncRNA_exon,proteingeneCDS,proteingene5UTR,proteingene3UTR,intergenic)
alldata$type <- factor(alldata$type,levels = c("CDS","5'UTR","3'UTR",
                                               "PAS lncRNA","PAA lncRNA",
                                               "NPA lncRNA","Random intergenic"))
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
  ggsignif::geom_signif(map_signif_level=T,y_position=c(1.6,1.8,2,1.6,1.8),
              comparisons = list(c("Random intergenic","NPA lncRNA"),
                                 c("NPA lncRNA","PAA lncRNA"),
                                 c("NPA lncRNA","PAS lncRNA"),
                                 c("PAA lncRNA","PAS lncRNA"),
                                 c("PAS lncRNA","3'UTR")))+
  scale_fill_manual(values=c("#8491B4","#FAA465","#7197AD","#628255","#BB5A5D","#BB5A5D","#BB5A5D"),
                    limits=c("Random intergenic","NPA lncRNA",
                             "PAS lncRNA","PAA lncRNA","CDS","5'UTR","3'UTR"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 2.2))+
  scale_x_discrete(labels = c("CDS","5'UTR","3'UTR","PAS lncRNA","PAA lncRNA","NPA lncRNA","Random\nintergenic"))+
  labs(y =expression("S"*"N"*"P"~"l"*"o"*"g"[10]*"("*"C"*"A"*"D"*"D"~"s"*"c"*"o"*"r"*"e"*")"),x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave(h,p1,width = 6, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```


##### 6.5 SNP DAF
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
bedtools intersect -a /share/home/zhzhang24/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b "/share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paslnc_lncRNA_mergeexon.nocpg.bed" -wa > /share/home/zhzhang24/PG/Evolution/snp/DAF/Homo_sapiens.paslnc_lncRNA_mergeexon.SNPDAF.txt
bedtools intersect -a /share/home/zhzhang24/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b "/share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.paalnc_lncRNA_mergeexon.nocpg.bed" -wa > /share/home/zhzhang24/PG/Evolution/snp/DAF/Homo_sapiens.paalnc_lncRNA_mergeexon.SNPDAF.txt
bedtools intersect -a /share/home/zhzhang24/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b "/share/home/zhzhang24/PG/Evolution/snp/Homo_sapiens.npalnc_lncRNA_mergeexon.nocpg.bed" -wa > /share/home/zhzhang24/PG/Evolution/snp/DAF/Homo_sapiens.npalnc_lncRNA_mergeexon.SNPDAF.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.20000random_3kb_intergenic.nocpg.bed" -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.20000random_3kb_intergenic.SNPDAF.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingeneCDS.nocpg.bed" -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingeneCDS.SNPDAF.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene5UTR.nocpg.bed" -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene5UTR.SNPDAF.txt
bedtools intersect -a /home/zhzhang/PG/Evolution/snp/DAF/1000GENOMES_phase_3.common.AFRDAF.bed -b "/home/zhzhang/PG/Evolution/snp/Homo_sapiens.proteingene3UTR.nocpg.bed" -wa > /home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene3UTR.SNPDAF.txt


```
```r
#对比两类lncRNA外显子，随机基因间区（阴性对照），蛋白编码基因外显子（CDS,5UTR,3UTR）（阳性对照）的SNP DAF
#人类commonSNP
a <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.paslnc_lncRNA_mergeexon.SNPDAF.txt"
b <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.paalnc_lncRNA_mergeexon.SNPDAF.txt"
ab <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.npalnc_lncRNA_mergeexon.SNPDAF.txt"
c <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingeneCDS.SNPDAF.txt"
d <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene5UTR.SNPDAF.txt"
e <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.proteingene3UTR.SNPDAF.txt"
f <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.20000random_3kb_intergenic.SNPDAF.txt"
g <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.lncRNA_SNPDAF.tj.txt"
h <- "/home/zhzhang/PG/Evolution/snp/DAF/Homo_sapiens.lncRNA_SNPDAF.pdf"
#导入pgdlnc外显子的SNP,去除averageDAF为NA 的SNP
pgdlnc_lncRNA_exon <- read.delim(a, header=FALSE)%>%
  select(5)%>%
  filter(V5<=1)%>%
  mutate(type="PAS lncRNA")
colnames(pgdlnc_lncRNA_exon)[c(1,2)] <- c("DAF","type")
#导入paalnc外显子的SNP,去除averageDAF为NA 的SNP
paalnc_lncRNA_exon <- read.delim(b, header=FALSE)%>%
  select(5)%>%
  filter(V5<=1)%>%
  mutate(type="PAA lncRNA")
colnames(paalnc_lncRNA_exon)[c(1,2)] <- c("DAF","type")
#导入npgdlnc外显子的SNP
npgdlnc_lncRNA_exon <- read.delim(ab, header=FALSE)%>%
  select(5)%>%
  filter(V5<=1)%>%
  mutate(type="NPA lncRNA")
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
alldata <- rbind(pgdlnc_lncRNA_exon,paalnc_lncRNA_exon,npgdlnc_lncRNA_exon,proteingeneCDS,proteingene5UTR,proteingene3UTR,intergenic)
alldata$type <- factor(alldata$type,levels = c("CDS","5'UTR","3'UTR",
                                               "PAS lncRNA","PAA lncRNA",
                                               "NPA lncRNA","Random intergenic"))
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
t.test(pgdlnc_lncRNA_exon$DAF,paalnc_lncRNA_exon$DAF) #p-value = 0.5939
t.test(npgdlnc_lncRNA_exon$DAF,intergenic$DAF) #p-value < 2.2e-16
t.test(paalnc_lncRNA_exon$DAF,npgdlnc_lncRNA_exon$DAF) #p-value = 0.6866
t.test(pgdlnc_lncRNA_exon$DAF,proteingeneCDS$DAF) #p-value < 2.2e-16
t.test(pgdlnc_lncRNA_exon$DAF,proteingene3UTR$DAF) #p-value < 2.2e-16

phs <- ggplot(data = tj,aes(x=type,y=mean))+
  geom_point(size=2,aes(fill=type,color=type))+
  geom_errorbar(aes(ymin=confmin,ymax=confmax,color=type),width=0.15)+
  ggsignif::geom_signif(annotations=c("***","N.S.","N.S.","***"),
              y_position=c(0.115,0.117,0.115,0.12),
              xmin = c(3,4,5,6),xmax = c(4,5,6,7))+
  scale_fill_manual(values=c("#8491B4","#FAA465","#7197AD","#628255","#BB5A5D","#BB5A5D","#BB5A5D"),
                    limits=c("Random intergenic","NPA lncRNA",
                             "PAS lncRNA","PAA lncRNA","CDS","5'UTR","3'UTR"))+
  scale_color_manual(values=c("#8491B4","#FAA465","#7197AD","#628255","#BB5A5D","#BB5A5D","#BB5A5D"),
                     limits=c("Random intergenic","NPA lncRNA",
                              "PAS lncRNA","PAA lncRNA","CDS","5'UTR","3'UTR"))+
  theme_half_open()+
  scale_x_discrete(labels = c("CDS","5'UTR","3'UTR","PAS lncRNA","PAA lncRNA","NPA lncRNA","Random\nintergenic"))+
  labs(y = "SNP derived allele frequency", x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =45)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave(h,phs,width = 6, height = 5.6,dpi=1200, units = "in", device='pdf',bg = "transparent")



```




##### 7.lncRNA基因年龄的鉴定和比较
```r
### download reciprocal-best互惠最佳 netted chain files
#for human
cd /home/zhzhang/PG/age/chain/hs/
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/vsHg38/reciprocalBest/galGal6.hg38.rbest.chain.gz &
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/vsHg38/reciprocalBest/hg38.galGal6.rbest.chain.gz &
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/vsGCF_016699485.2/reciprocalBest/GCF_016699485.2.galGal6.rbest.chain.gz &
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/vsGCF_016699485.2/reciprocalBest/galGal6.GCF_016699485.2.rbest.chain.gz &
https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsHg38/reciprocalBest/rheMac10.hg38.rbest.chain.gz
https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsHg38/reciprocalBest/hg38.rheMac10.rbest.chain.gz
nohup wget -c -P Chimp/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/reciprocalBest/panTro6.hg38.rbest.chain.gz &
nohup wget -c -P Chimp/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/reciprocalBest/hg38.panTro6.rbest.chain.gz &
nohup wget -c -P Dog/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsCanFam6/reciprocalBest/canFam6.hg38.rbest.chain.gz  &
nohup wget -c -P Dog/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsCanFam6/reciprocalBest/hg38.canFam6.rbest.chain.gz &
nohup wget -c -P Mouse/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm39/reciprocalBest/hg38.mm39.rbest.chain.gz &
nohup wget -c -P Mouse/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMm39/reciprocalBest/mm39.hg38.rbest.chain.gz &
nohup wget -c -P Opossum/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMonDom5/reciprocalBest/hg38.monDom5.rbest.chain.gz &
nohup wget -c -P Opossum/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsMonDom5/reciprocalBest/monDom5.hg38.rbest.chain.gz &
nohup wget -c -P Platypus/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsOrnAna2/reciprocalBest/hg38.ornAna2.rbest.chain.gz &
nohup wget -c -P Platypus/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsOrnAna2/reciprocalBest/ornAna2.hg38.rbest.chain.gz &
nohup wget -c -P Rabbit/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsOryCun2/reciprocalBest/hg38.oryCun2.rbest.chain.gz &
nohup wget -c -P Rabbit/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsOryCun2/reciprocalBest/oryCun2.hg38.rbest.chain.gz &
nohup wget -c -P Rat/ https://hgdownload.soe.ucsc.edu/goldenPath/rn7/vsHg38/reciprocalBest/hg38.rn7.rbest.chain.gz &
nohup wget -c -P Rat/ https://hgdownload.soe.ucsc.edu/goldenPath/rn7/vsHg38/reciprocalBest/rn7.hg38.rbest.chain.gz &
nohup wget -c -P Rhesus/ https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsHg38/reciprocalBest/hg38.rheMac10.rbest.chain.gz &
nohup wget -c -P Rhesus/ https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsHg38/reciprocalBest/rheMac10.hg38.rbest.chain.gz &
nohup wget -c -P XTropicalis https://hgdownload.soe.ucsc.edu/goldenPath/xenTro10/vsHg38/reciprocalBest/hg38.xenTro10.rbest.chain.gz &
nohup wget -c -P XTropicalis https://hgdownload.soe.ucsc.edu/goldenPath/xenTro10/vsHg38/reciprocalBest/xenTro10.hg38.rbest.chain.gz &
nohup wget -c -P Zebrafish https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsDanRer11/reciprocalBest/danRer11.hg38.rbest.chain.gz &
nohup wget -c -P Zebrafish https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsDanRer11/reciprocalBest/hg38.danRer11.rbest.chain.gz &


#for mouse
cd /home/zhzhang/PG/age/chain/mm/
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsGalGal6/reciprocalBest/galGal6.mm39.rbest.chain.gz &
nohup wget -c -P Chicken/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsGalGal6/reciprocalBest/mm39.galGal6.rbest.chain.gz &
https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsMm39/reciprocalBest/mm39.rheMac10.rbest.chain.gz
nohup wget -c -P Chimp/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsPanTro6/reciprocalBest/mm39.panTro6.rbest.chain.gz &
nohup wget -c -P Chimp/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsPanTro6/reciprocalBest/panTro6.mm39.rbest.chain.gz &
nohup wget -c -P Dog/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsCanFam6/reciprocalBest/canFam6.mm39.rbest.chain.gz & 
nohup wget -c -P Dog/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsCanFam6/reciprocalBest/mm39.canFam6.rbest.chain.gz &
nohup wget -c -P Opossum/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsMonDom5/reciprocalBest/mm39.monDom5.rbest.chain.gz &
nohup wget -c -P Opossum/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsMonDom5/reciprocalBest/monDom5.mm39.rbest.chain.gz &
nohup wget -c -P platypus/ http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsOrnAna1/reciprocalBest/mm10.ornAna1.rbest.chain.gz &
nohup wget -c -P platypus/ http://hgdownload.soe.ucsc.edu/goldenPath/mm10/vsOrnAna1/reciprocalBest/ornAna1.mm10.rbest.chain.gz &
nohup wget -c -P Rabbit/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsOryCun2/reciprocalBest/mm39.oryCun2.rbest.chain.gz &
nohup wget -c -P Rabbit/ https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsOryCun2/reciprocalBest/oryCun2.mm39.rbest.chain.gz &
nohup wget -c -P Rat/ https://hgdownload.soe.ucsc.edu/goldenPath/rn7/vsMm39/reciprocalBest/mm39.rn7.rbest.chain.gz &
nohup wget -c -P Rat/ https://hgdownload.soe.ucsc.edu/goldenPath/rn7/vsMm39/reciprocalBest/rn7.mm39.rbest.chain.gz &
nohup wget -c -P Rhesus/ https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsMm39/reciprocalBest/mm39.rheMac10.rbest.chain.gz &
nohup wget -c -P Rhesus/ https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/vsMm39/reciprocalBest/rheMac10.mm39.rbest.chain.gz &
nohup wget -c -P XTropicalis https://hgdownload.soe.ucsc.edu/goldenPath/xenTro10/vsMm39/reciprocalBest/mm39.xenTro10.rbest.chain.gz &
nohup wget -c -P XTropicalis https://hgdownload.soe.ucsc.edu/goldenPath/xenTro10/vsMm39/reciprocalBest/xenTro10.mm39.rbest.chain.gz &
nohup wget -c -P Zebrafish https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsDanRer11/reciprocalBest/danRer11.mm39.rbest.chain.gz &
nohup wget -c -P Zebrafish https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsDanRer11/reciprocalBest/mm39.danRer11.rbest.chain.gz &
nohup wget -c -P Human https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsHg38/reciprocalBest/hg38.mm39.rbest.chain.gz &
nohup wget -c -P Human https://hgdownload.soe.ucsc.edu/goldenPath/mm39/vsHg38/reciprocalBest/mm39.hg38.rbest.chain.gz &





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
echo "# == 4.lastz between human lncRNAs and ${jian} homo-lncRNA =="
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




################# R脚本Step5 ################
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
filtergene1 <- group_by(alscore,geneid_1)%>%
  top_n(1,score)%>%
  mutate(lncpair=paste(geneid_1,geneid_2,sep = "----"))
filtergene2 <- group_by(alscore,geneid_2)%>%
  top_n(1,score)%>%
  mutate(lncpair=paste(geneid_1,geneid_2,sep = "----"))
homo1v1pair <- data.frame("lncpair"=intersect(filtergene1$lncpair,filtergene2$lncpair))%>%
  separate(lncpair,c("geneid_1","geneid_2"),sep="----")
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
############################################




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
#Rscript /home/zhzhang/PG/age/chain/sh/step5.R --input /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.${quan}.homo_lncRNApair.alscore.txt --sp ${jian} --outputone /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.${quan}.11homo_lncRNApair.txt --outputtwo /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.own_${quan}_11homo.lncRNAid.txt
Rscript /home/zhzhang/PG/age/chain/sh/step5.pro.R --input /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.${quan}.homo_lncRNApair.alscore.txt --sp ${jian} --outputone /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.${quan}.11homo_lncRNApair.pro.txt --outputtwo /home/zhzhang/PG/age/chain/${vipspdir}/${jian}/${vipspdir}.own_${quan}_11homo.lncRNAid.pro.txt





#结果传至实验室服务器
scp -r -P 22 zhzhang@211.69.141.147:"/home/zhzhang/PG/age/" /home/zhzhang/PG/


```
```r
#Human
############
#导入全部lncRNA id和分类
lncRNA_class <- read.delim("~/PG/lncRNA_class/new/Homo_sapiens.lncRNA_class.txt")
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
#去除没有1-1同源lnc但有潜在homo的lnc(年龄模糊)
lncage <- left_join(lncage,lnc_phomo,by="geneid")%>%
  mutate(ambiguous=0)
lncage$phomo[is.na(lncage$phomo)==T] <- 0
lncage$ambiguous[lncage$phomo==1 & lncage$max=="Human"] <- 1
lncage$max[lncage$ambiguous==1] <- "Ambiguous"
lncage <- select(lncage,1,2)
#储存
data.table::fwrite(lncage,file ="/home/zhzhang/PG/age/Human.lncRNA.age.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)





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

```
```r
#human
a="/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.lncRNA_class.txt"
b="/home/zhzhang/PG/age/Human.lncRNA.age.txt"
c="/home/zhzhang/PG/age/Human.lncRNA.type.age_enrich.txt"
d=c("429 Ma","319 Ma","160 Ma",
    "87 Ma","29 Ma","Human")
e="/home/zhzhang/PG/age/Human.lnc.enrich.pdf"
#mouse
a="/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.lncRNA_class.txt"
b="/home/zhzhang/PG/age/Mouse.lncRNA.age.txt"
c="/home/zhzhang/PG/age/Mouse.lncRNA.type.age_enrich.txt"
d=c("429 Ma","319 Ma","160 Ma",
    "87 Ma","79 Ma","13 Ma","Mouse")
e="/home/zhzhang/PG/age/Mouse.lnc.enrich.pdf"
#导入全部lncRNA id和分类,age
lncRNA_class <- read.delim(a)
lncage <- read.delim(b)
#富集
anno <- mutate(lncage,GO_Description=max)
colnames(anno)[2] <- "GO"
go2gene <- anno[, c(2, 1)]
go2name <- anno[, c(2, 3)]
paslnc <- filter(lncRNA_class,type=="Pseudogene-associated sense lncRNA")%>%
  select(1)
paalnc <- filter(lncRNA_class,type=="Pseudogene-associated antisense lncRNA")%>%
  select(1)
set.seed(10)
npalnc <- filter(lncRNA_class,type=="Non-pseudogene-associated lncRNA")%>%
  select(1)%>%
  dplyr::sample_n(2000)
ego_npa <- clusterProfiler::enricher(npalnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                     maxGSSize = 40000)@result%>%
  mutate(type="NPA lncRNA")
ego_pas <- clusterProfiler::enricher(paslnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                    maxGSSize = 40000)@result%>%
  mutate(type="PAS lncRNA")
ego_paa <- clusterProfiler::enricher(paalnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                     maxGSSize = 40000)@result%>%
  mutate(type="PAA lncRNA")
#合并
ego <- rbind(ego_pas,ego_paa,ego_npa)%>%
  select(10,1:9)
#储存
data.table::fwrite(ego,file =c,
                   sep = '\t',row.names = F,quote = F,col.names = T)

#plot
ego <- read.delim(c)
ego <- filter(ego,ID!="Ambiguous")%>%
  separate(GeneRatio,c("GeneRatio1","GeneRatio2"))%>%
  separate(BgRatio,c("BgRatio1","BgRatio2"))%>%
  mutate(GeneRatio1=as.numeric(GeneRatio1),GeneRatio2=as.numeric(GeneRatio2),
         BgRatio1=as.numeric(BgRatio1),BgRatio2=as.numeric(BgRatio2))%>%
  mutate(ES=(GeneRatio1/GeneRatio2)/(BgRatio1/BgRatio2))
ego$type <- factor(ego$type,levels=c("PAS lncRNA","PAA lncRNA","NPA lncRNA"))
ego$ID <- factor(ego$ID,levels=d)
#富集dotplot
#human
pp1 <- ggplot(data = ego,aes(x=type,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:5.5),color="gray")+
  geom_vline(xintercept=c(1.5:2.5),color="gray")+
  scale_fill_gradient(low="white",high="#C12039",breaks=c(0,10,20,30))+
  scale_size(breaks = c(0.5,1,1.5,2,2.5),range = c(3,15))+
  theme_half_open()+
  labs(x = NULL, y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(e, 
       pp1,width = 5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#or mouse
pp1 <- ggplot(data = ego,aes(x=type,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:6.5),color="gray")+
  geom_vline(xintercept=c(1.5:2.5),color="gray")+
  scale_fill_gradient(low="white",high="#C12039",breaks=c(0,5,10,15))+
  scale_size(breaks = c(1,2,3),range = c(3,15))+
  theme_half_open()+
  labs(x = NULL, y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(e, 
       pp1,width = 4.8, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#fishertest检验三类lnc指定年龄基因比例的差异p值
npalnc2 <- filter(lncRNA_class,type=="Non-pseudogene-associated lncRNA")%>%
  select(1)
ego_npa2 <- clusterProfiler::enricher(npalnc2$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                      maxGSSize = 40000)@result%>%
  mutate(type="NPA lncRNA")
ego2 <- rbind(ego_pas,ego_paa,ego_npa2)%>%
  select(10,1:9)%>%
  filter(ID!="Ambiguous")%>%
  separate(GeneRatio,c("GeneRatio1","GeneRatio2"))%>%
  separate(BgRatio,c("BgRatio1","BgRatio2"))%>%
  mutate(GeneRatio1=as.numeric(GeneRatio1),GeneRatio2=as.numeric(GeneRatio2),GeneRatio=GeneRatio1*100/GeneRatio2,
         BgRatio1=as.numeric(BgRatio1),BgRatio2=as.numeric(BgRatio2))%>%
  mutate(ES=(GeneRatio1/GeneRatio2)/(BgRatio1/BgRatio2))
ego2$type <- factor(ego2$type,levels=c("PAS lncRNA","PAA lncRNA","NPA lncRNA"))
ego2$ID <- factor(ego2$ID,levels=d)
pdata <- distinct(ego2,ID)
for (i in 1:nrow(pdata)) {
  forid <- pdata$ID[i]
  spglncnum <- filter(ego2,type=="PAS lncRNA" & ID==forid)$GeneRatio1
  spglncnum2 <- filter(ego2,type=="PAS lncRNA" & ID==forid)$GeneRatio2
  snpglncnum <- filter(ego2,type=="NPA lncRNA" & ID==forid)$GeneRatio1
  snpglncnum2 <- filter(ego2,type=="NPA lncRNA" & ID==forid)$GeneRatio2
  paalncnum <- filter(ego2,type=="PAA lncRNA" & ID==forid)$GeneRatio1
  paalncnum2 <- filter(ego2,type=="PAA lncRNA" & ID==forid)$GeneRatio2
  pdata[i,2] <- fisher.test(data.frame(c(spglncnum,snpglncnum),c(spglncnum2-spglncnum,snpglncnum2-snpglncnum)))[["p.value"]]
  pdata[i,3] <- fisher.test(data.frame(c(paalncnum,snpglncnum),c(paalncnum2-paalncnum,snpglncnum2-snpglncnum)))[["p.value"]]
  pdata[i,4] <- fisher.test(data.frame(c(spglncnum,paalncnum),c(spglncnum2-spglncnum,paalncnum2-paalncnum)))[["p.value"]]
  }
colnames(pdata)[2:4]=c("PASvsNPA","PAAvsNPA","PASvsPAA")



```
#### 7.5 严格阈值对palnc年龄的影响
```r
#human
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt"
a="/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.lncRNA_class.txt"
b="/home/zhzhang/PG/age/Human.lncRNA.age.txt"
d=c("429 Ma","319 Ma","160 Ma",
    "87 Ma","29 Ma","Human")
e="/home/zhzhang/PG/age/Human.paslnc.enrich.pdf"
f="/home/zhzhang/PG/age/Human.paalnc.enrich.pdf"
#mouse
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.palncRNA_message.txt"
a="/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.lncRNA_class.txt"
b="/home/zhzhang/PG/age/Mouse.lncRNA.age.txt"
d=c("429 Ma","319 Ma","160 Ma",
    "87 Ma","79 Ma","13 Ma","Mouse")
e="/home/zhzhang/PG/age/Mouse.paslnc.enrich.pdf"
f="/home/zhzhang/PG/age/Mouse.paalnc.enrich.pdf"
#导入lnc分组
lncRNA_class <- read.delim(aaa)%>%select(2,4,7)
colnames(lncRNA_class)[1]="geneid"
#导入全部lncRNA id和分类,age
lncage <- read.delim(b)
#富集
anno <- mutate(lncage,GO_Description=max)
colnames(anno)[2] <- "GO"
go2gene <- anno[, c(2, 1)]
go2name <- anno[, c(2, 3)]
ego=data.frame()
for (i in c(0,50,100,150,200)) {
  paslnc <- filter(lncRNA_class,type=="Pseudogene-associated sense lncRNA" & total_overlap_len>i)%>%
    select(1)
  paalnc <- filter(lncRNA_class,type=="Pseudogene-associated antisense lncRNA" & total_overlap_len>i)%>%
    select(1)
  ego_pas <- clusterProfiler::enricher(paslnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                       maxGSSize = 40000)@result%>%
    mutate(type="PAS lncRNA",cuttype=as.character(i))
  ego_paa <- clusterProfiler::enricher(paalnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                       maxGSSize = 40000)@result%>%
    mutate(type="PAA lncRNA",cuttype=as.character(i))
  #合并
  ego <- rbind(ego,ego_pas,ego_paa)
}
ego <- arrange(ego,type,cuttype)%>%
  select(10:11,1:9)
#
ego <- filter(ego,ID!="Ambiguous")%>%
  separate(GeneRatio,c("GeneRatio1","GeneRatio2"))%>%
  separate(BgRatio,c("BgRatio1","BgRatio2"))%>%
  mutate(GeneRatio1=as.numeric(GeneRatio1),GeneRatio2=as.numeric(GeneRatio2),
         BgRatio1=as.numeric(BgRatio1),BgRatio2=as.numeric(BgRatio2))%>%
  mutate(ES=(GeneRatio1/GeneRatio2)/(BgRatio1/BgRatio2))
ego$ID <- factor(ego$ID,levels = d)
ego$cuttype <- factor(ego$cuttype,levels = c("0","50","100","150","200"))
ego <- arrange(ego,type,cuttype)
#富集dotplot
#human
pp1 <- ggplot(data = filter(ego,type=="PAS lncRNA"),aes(x=cuttype,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:5.5),color="gray")+
  geom_vline(xintercept=c(1.5:4.5),color="gray")+
  scale_fill_gradient(low="white",high="#7197AD",breaks=c(0,5,10))+
  scale_size(breaks = c(1,1.5,2,2.5),range = c(3,15))+
  theme_half_open()+
  labs(x = "Cutoff of overlap length (bp)", y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(e,pp1,width = 5.5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
pp2 <- ggplot(data = filter(ego,type=="PAA lncRNA"),aes(x=cuttype,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:5.5),color="gray")+
  geom_vline(xintercept=c(1.5:4.5),color="gray")+
  scale_fill_gradient(low="white",high="#628255",breaks=c(0,10,20,30))+
  scale_size(breaks = c(1,1.5,2,2.5),range = c(3,15))+
  theme_half_open()+
  labs(x = "Cutoff of overlap length (bp)", y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(f,pp2,width = 5.5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#or mouse
pp1 <- ggplot(data = filter(ego,type=="PAS lncRNA"),aes(x=cuttype,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:6.5),color="gray")+
  geom_vline(xintercept=c(1.5:4.5),color="gray")+
  scale_fill_gradient(low="white",high="#7197AD",breaks=c(0,2.5,5,10))+
  scale_size(breaks = c(1,1.5,2,2.5),range = c(3,15))+
  theme_half_open()+
  labs(x = "Cutoff of overlap length (bp)", y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(e,pp1,width = 5.5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
pp2 <- ggplot(data = filter(ego,type=="PAA lncRNA"),aes(x=cuttype,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:6.5),color="gray")+
  geom_vline(xintercept=c(1.5:4.5),color="gray")+
  scale_fill_gradient(low="white",high="#628255",breaks=c(0,5,10,15))+
  scale_size(breaks = c(1,1.5,2,2.5),range = c(3,15))+
  theme_half_open()+
  labs(x = "Cutoff of overlap length (bp)", y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(f,pp2,width = 5.5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
#### 7.6严格表达阈值对palnc年龄的影响
```r
#计算基因最大表达量
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
                      b = "/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                      c = "Mouse")
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                      c = "Human")
ALLsp_TPM <- rbind(mm_TPM,hs_TPM)%>%select(1,3)
#human
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt"
a="/home/zhzhang/PG/lncRNA_class/new_r1/Homo_sapiens.lncRNA_class.txt"
b="/home/zhzhang/PG/age/Human.lncRNA.age.txt"
d=c("429 Ma","319 Ma","160 Ma",
    "87 Ma","29 Ma","Human")
e="/home/zhzhang/PG/age/Human.paslnc.expcutoff.enrich.pdf"
f="/home/zhzhang/PG/age/Human.paalnc.expcutoff.enrich.pdf"
#mouse
aaa <- "/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.palncRNA_message.txt"
a="/home/zhzhang/PG/lncRNA_class/new_r1/Mus_musculus.lncRNA_class.txt"
b="/home/zhzhang/PG/age/Mouse.lncRNA.age.txt"
d=c("429 Ma","319 Ma","160 Ma",
    "87 Ma","79 Ma","13 Ma","Mouse")
e="/home/zhzhang/PG/age/Mouse.paslnc.expcutoff.enrich.pdf"
f="/home/zhzhang/PG/age/Mouse.paalnc.expcutoff.enrich.pdf"
#导入lnc分组
lncRNA_class <- read.delim(aaa)%>%
  select(2,7)%>%
  dplyr::rename(geneid=lncgid)%>%
  left_join(ALLsp_TPM,by="geneid")
#导入全部lncRNA id和分类,age
lncage <- read.delim(b)
#富集
anno <- mutate(lncage,GO_Description=max)
colnames(anno)[2] <- "GO"
go2gene <- anno[, c(2, 1)]
go2name <- anno[, c(2, 3)]
ego=data.frame()
for (i in c(0,0.5,1,1.5,2)) {
  paslnc <- filter(lncRNA_class,type=="Pseudogene-associated sense lncRNA" & maxTPM>i)%>%
    select(1)
  paalnc <- filter(lncRNA_class,type=="Pseudogene-associated antisense lncRNA" & maxTPM>i)%>%
    select(1)
  ego_pas <- clusterProfiler::enricher(paslnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                       maxGSSize = 40000)@result%>%
    mutate(type="PAS lncRNA",cuttype=as.character(i))
  ego_paa <- clusterProfiler::enricher(paalnc$geneid, TERM2GENE = go2gene, TERM2NAME = go2name, pAdjustMethod = "fdr",pvalueCutoff  = 0, qvalueCutoff  = 0,minGSSize = 0,
                                       maxGSSize = 40000)@result%>%
    mutate(type="PAA lncRNA",cuttype=as.character(i))
  #合并
  ego <- rbind(ego,ego_pas,ego_paa)
}
ego <- arrange(ego,type,cuttype)%>%
  select(10:11,1:9)
#
ego <- filter(ego,ID!="Ambiguous")%>%
  separate(GeneRatio,c("GeneRatio1","GeneRatio2"))%>%
  separate(BgRatio,c("BgRatio1","BgRatio2"))%>%
  mutate(GeneRatio1=as.numeric(GeneRatio1),GeneRatio2=as.numeric(GeneRatio2),
         BgRatio1=as.numeric(BgRatio1),BgRatio2=as.numeric(BgRatio2))%>%
  mutate(ES=(GeneRatio1/GeneRatio2)/(BgRatio1/BgRatio2))
ego$ID <- factor(ego$ID,levels = d)
ego$cuttype <- factor(ego$cuttype,levels = c("0","0.5","1","1.5","2"))
ego <- arrange(ego,type,cuttype)
#富集dotplot
#human
pp1 <- ggplot(data = filter(ego,type=="PAS lncRNA"),aes(x=cuttype,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:5.5),color="gray")+
  geom_vline(xintercept=c(1.5:4.5),color="gray")+
  scale_fill_gradient(low="white",high="#7197AD",breaks=c(0,5,10))+
  scale_size(breaks = c(1,1.5,2,2.5),range = c(3,15))+
  cowplot::theme_half_open()+
  labs(x = "Cutoff of expression (TPM)", y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(e,pp1,width = 5.5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
pp2 <- ggplot(data = filter(ego,type=="PAA lncRNA"),aes(x=cuttype,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:5.5),color="gray")+
  geom_vline(xintercept=c(1.5:4.5),color="gray")+
  scale_fill_gradient(low="white",high="#628255",breaks=c(0,10,20,30))+
  scale_size(breaks = c(1,1.5,2,2.5),range = c(3,15))+
  cowplot::theme_half_open()+
  labs(x = "Cutoff of expression (TPM)", y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(f,pp2,width = 5.5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#or mouse
pp1 <- ggplot(data = filter(ego,type=="PAS lncRNA"),aes(x=cuttype,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:6.5),color="gray")+
  geom_vline(xintercept=c(1.5:4.5),color="gray")+
  scale_fill_gradient(low="white",high="#7197AD",breaks=c(0,2.5,5,10))+
  scale_size(breaks = c(1,1.5,2,2.5),range = c(3,15))+
  theme_half_open()+
  labs(x = "Cutoff of expression (TPM)", y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(e,pp1,width = 5.5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
pp2 <- ggplot(data = filter(ego,type=="PAA lncRNA"),aes(x=cuttype,y=ID))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:6.5),color="gray")+
  geom_vline(xintercept=c(1.5:4.5),color="gray")+
  scale_fill_gradient(low="white",high="#628255",breaks=c(0,5,10,15))+
  scale_size(breaks = c(1,1.5,2,2.5),range = c(3,15))+
  theme_half_open()+
  labs(x = "Cutoff of expression (TPM)", y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave(f,pp2,width = 5.5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```

### 四.基因表达谱分析
##### 1.RNA-seq数据
```r
#RNA-SEQ
scp -P 5614 zhzhang@47.102.45.227:/NAS/chenl/lncRNA/mouse/Illumina/*.fq.gz /home/zhzhang/PG/RNAseqdata/Mus_musculus/
scp -P 5614 zhzhang@47.102.45.227:/NAS/chenl/lncRNA/mouse/Illumina/other/68_clean* /home/zhzhang/PG/RNAseqdata/Mus_musculus/

scp -P 5614 zhzhang@47.102.45.227:/NAS/chenl/lnc_zebrafish/data/Zebrafish_*.fq.gz /home/zhzhang/PG/RNAseqdata/Danio_rerio/

scp -P 5614 zhzhang@47.102.45.227:/NAS/chenl/lncRNA/Illumina_data/chicken/*.fq.gz /home/zhzhang/PG/RNAseqdata/Gallus_gallus/

scp -P 5614 zhzhang@47.102.45.227:/NAS/chenl/lncRNA/Illumina_data/medaka/*.fq.gz /home/zhzhang/PG/RNAseqdata/Oryzias_latipes/
```
```r
#人类数据（https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1733/sdrf: https://www.ebi.ac.uk/ena/browser/view/PRJEB4337 and https://www.ebi.ac.uk/ena/browser/view/PRJEB6971）
for i in $(seq 325 495)
do
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/ERR315${i}/ERR315${i}_1.fastq.gz
done


for i in $(seq 325 495)
do
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR315/ERR315${i}/ERR315${i}_2.fastq.gz
done


for i in $(seq 122 155)
do
if [ "${i}" != "126" -a "${i}" != "137" -a "${i}" != "144" -a "${i}" != "145" -a "${i}" != "154" ]
then
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR579/ERR579${i}/ERR579${i}_1.fastq.gz
fi
done


for i in $(seq 122 155)
do
if [ "${i}" != "126" -a "${i}" != "137" -a "${i}" != "144" -a "${i}" != "145" -a "${i}" != "154" ]
then
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR579/ERR579${i}/ERR579${i}_2.fastq.gz
fi
done



#批量改名
#R
E_MTAB_2836 <- read.delim("~/E-MTAB-2836.tsv")%>%
  select(1,6,31)
colnames(E_MTAB_2836) <- c("sample","sex","ERR")
E_MTAB_2836_r <-   distinct(E_MTAB_2836,ERR,.keep_all = T)%>%
  separate(sample,c("t","r"),sep = "_")

E_MTAB_2836_r[E_MTAB_2836_r[,1]=="endometrium",3] <- "female"
E_MTAB_2836_r[E_MTAB_2836_r[,1]=="fallopiantube",3] <- "female"
E_MTAB_2836_r[E_MTAB_2836_r[,1]=="ovary",3] <- "female"
E_MTAB_2836_r[E_MTAB_2836_r[,1]=="prostate",3] <- "male"
E_MTAB_2836_r[E_MTAB_2836_r[,1]=="testis",3] <- "male"
E_MTAB_2836_r[E_MTAB_2836_r[,3]=="  ",3] <- "no"

E_MTAB_2836_r <- arrange(E_MTAB_2836_r,t)%>%
  mutate(id=c(1:200),sp="Human")%>%
  mutate(yuan1=paste(ERR,"1.fastq.gz",sep="_"),hou1=paste(sp,t,sex,id,"R1.fq.gz",sep="_"))%>%
  mutate(yuan2=paste(ERR,"2.fastq.gz",sep="_"),hou2=paste(sp,t,sex,id,"R2.fq.gz",sep="_"))

E_MTAB_2836_r1 <- select(E_MTAB_2836_r,yuan1,hou1)
data.table::fwrite(E_MTAB_2836_r1,
                   file ="~/E_MTAB_2836_r1.txt",sep = '\t',row.names = F,quote = F,col.names = F)
E_MTAB_2836_r2 <- select(E_MTAB_2836_r,yuan2,hou2)
data.table::fwrite(E_MTAB_2836_r2,
                   file ="~/E_MTAB_2836_r2.txt",sep = '\t',row.names = F,quote = F,col.names = F)

#linux
cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens/
for i in $(seq 1 200)
do
 oldname=`cat "/home/zhzhang/PG/RNAseqdata/E_MTAB_2836_r1.txt" |awk '{print $1}' |sed -n "${i}p"`
 newfile=`cat "/home/zhzhang/PG/RNAseqdata/E_MTAB_2836_r1.txt" |awk '{print $2}' |sed -n "${i}p"`
 mv ${oldname} ${newfile}
done

for i in $(seq 1 200)
do
 oldname=`cat "/home/zhzhang/PG/RNAseqdata/E_MTAB_2836_r2.txt" |awk '{print $1}' |sed -n "${i}p"`
 newfile=`cat "/home/zhzhang/PG/RNAseqdata/E_MTAB_2836_r2.txt" |awk '{print $2}' |sed -n "${i}p"`
 mv ${oldname} ${newfile}
done



```
```r
#人DD数据（https://www.ebi.ac.uk/ena/browser/view/PRJEB26969
"/home/zhzhang/PG/RNAseqdata/humandd.sh"

cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/
for i in $(seq 1 313)
do
 oldname=`awk '{print $1}' "/home/zhzhang/PG/RNAseqdata/filereport_read_run_PRJEB26969_tsv.txt" |sed -n "${i}p"`
 newfile=`awk '{print $2}' "/home/zhzhang/PG/RNAseqdata/filereport_read_run_PRJEB26969_tsv.txt" |sed -n "${i}p"|awk -F '.' '{print $2"_"$3"_"$5"_"$4"_"$1".fq.gz"}'`
 mv "${oldname}.fastq.gz" "${newfile}"
done


#小鼠DD数据（https://www.ebi.ac.uk/ena/browser/view/PRJEB26869
"/home/zhzhang/PG/RNAseqdata/mousedd.sh"
sed -n '1,40p' "/home/zhzhang/PG/RNAseqdata/mousedd.sh" > "/home/zhzhang/PG/RNAseqdata/mousedd1.sh"
sed -n '41,80p' "/home/zhzhang/PG/RNAseqdata/mousedd.sh" > "/home/zhzhang/PG/RNAseqdata/mousedd2.sh"
sed -n '81,120p' "/home/zhzhang/PG/RNAseqdata/mousedd.sh" > "/home/zhzhang/PG/RNAseqdata/mousedd3.sh"
sed -n '121,160p' "/home/zhzhang/PG/RNAseqdata/mousedd.sh" > "/home/zhzhang/PG/RNAseqdata/mousedd4.sh"
sed -n '161,200p' "/home/zhzhang/PG/RNAseqdata/mousedd.sh" > "/home/zhzhang/PG/RNAseqdata/mousedd5.sh"
sed -n '201,240p' "/home/zhzhang/PG/RNAseqdata/mousedd.sh" > "/home/zhzhang/PG/RNAseqdata/mousedd6.sh"
sed -n '241,280p' "/home/zhzhang/PG/RNAseqdata/mousedd.sh" > "/home/zhzhang/PG/RNAseqdata/mousedd7.sh"
sed -n '281,317p' "/home/zhzhang/PG/RNAseqdata/mousedd.sh" > "/home/zhzhang/PG/RNAseqdata/mousedd8.sh"
cd /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/
for i in $(seq 1 317)
do
 oldname=`awk '{print $1}' "/home/zhzhang/PG/RNAseqdata/filereport_read_run_PRJEB26869_tsv.txt" |sed -n "${i}p"`
 newfile=`awk '{print $2}' "/home/zhzhang/PG/RNAseqdata/filereport_read_run_PRJEB26869_tsv.txt" |sed -n "${i}p"|awk -F '.' '{print $2"_"$3"_"$5"_"$4"_"$1".fq.gz"}'`
 mv "${oldname}.fastq.gz" "${newfile}"
done



```
```r
#小鼠衰老
tail -n +2 "/share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/filereport_read_run_PRJNA936435_tsv.txt"|awk '{print $2"\t"$3}'|sed 's/srr/fastq/g'|grep -v 15mo  > /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/sample.tsv
#下载R1
cat "/share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/sample.tsv"|while read aaa
do
sample=`echo -e "${aaa}"|awk '{print $1}'`
srr=`echo -e "${aaa}"|awk '{print $2}'`
wget -c ${srr}_1.fastq.gz -O /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/fq/${sample}_R1.fq.gz
echo -e "${sample}" >> /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/R1_D.log
done
#下载R2
cat "/share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/sample.tsv"|while read aaa
do
sample=`echo -e "${aaa}"|awk '{print $1}'`
srr=`echo -e "${aaa}"|awk '{print $2}'`
wget -c ${srr}_2.fastq.gz -O /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/fq/${sample}_R2.fq.gz
echo -e "${sample}" >> /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/R2_D.log
done




```


##### 2.基因组索引构建
```r
#斑马鱼
for i in $(seq 1 25|sort)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome.${i}.fa >> /home/zhzhang/PG/refgenome/Danio_rerio.GRCz11.dna.chr.fa
done



#鸡
for i in $(seq 1 39|sort)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.${i}.fa >> /home/zhzhang/PG/refgenome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.chr.fa
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.W.fa >> /home/zhzhang/PG/refgenome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.chr.fa
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.Z.fa >> /home/zhzhang/PG/refgenome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.chr.fa


#人
for i in $(seq 1 22|sort)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${i}.fa >> /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa
done
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa" >> /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa" >> /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa


#小鼠
for i in $(seq 1 19|sort)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.${i}.fa >> /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa
done
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.X.fa" >> /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.Y.fa" >> /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa



#猕猴
for i in $(seq 1 20|sort)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.${i}.fa >> /home/zhzhang/PG/refgenome/Macaca_mulatta.Mmul_10.dna.chr.fa
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.X.fa >> /home/zhzhang/PG/refgenome/Macaca_mulatta.Mmul_10.dna.chr.fa
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.Y.fa >> /home/zhzhang/PG/refgenome/Macaca_mulatta.Mmul_10.dna.chr.fa


#大鼠
for i in $(seq 1 20|sort)
do
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.${i}.fa" >> /home/zhzhang/PG/refgenome/Rattus_norvegicus.mRatBN7.2.dna.chr.fa
done
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.X.fa" >> /home/zhzhang/PG/refgenome/Rattus_norvegicus.mRatBN7.2.dna.chr.fa
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.Y.fa" >> /home/zhzhang/PG/refgenome/Rattus_norvegicus.mRatBN7.2.dna.chr.fa


#兔子
for i in $(seq 1 21|sort)
do
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.${i}.fa" >> /home/zhzhang/PG/refgenome/Oryctolagus_cuniculus.OryCun2.0.dna.chr.fa
done
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.X.fa" >> /home/zhzhang/PG/refgenome/Oryctolagus_cuniculus.OryCun2.0.dna.chr.fa

#负鼠
for i in $(seq 1 8|sort)
do
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna.primary_assembly.${i}.fa" >> /home/zhzhang/PG/refgenome/Monodelphis_domestica.ASM229v1.dna.chr.fa
done
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna.primary_assembly.X.fa" >> /home/zhzhang/PG/refgenome/Monodelphis_domestica.ASM229v1.dna.chr.fa







#STAR基因组索引
STAR --runThreadN 28 --runMode genomeGenerate --genomeDir /home/zhzhang/PG/refgenome/STAR_Danio_rerio_GRCz11/ --genomeFastaFiles /home/zhzhang/PG/refgenome/Danio_rerio.GRCz11.dna.chr.fa --sjdbGTFfile /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/Danio_rerio.GRCz11.108.chr.gtf --sjdbOverhang 149
STAR --runThreadN 28 --runMode genomeGenerate --genomeDir /home/zhzhang/PG/refgenome/STAR_Gallus_gallus_GRCg7b/ --genomeFastaFiles /home/zhzhang/PG/refgenome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.chr.fa --sjdbGTFfile /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.chr.gtf --sjdbOverhang 149 --genomeSAindexNbases 13
STAR --runThreadN 28 --runMode genomeGenerate --genomeDir /home/zhzhang/PG/refgenome/STAR_Homo_sapiens_GRCh38/ --genomeFastaFiles /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa --sjdbGTFfile /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/Homo_sapiens.GRCh38.108.chr.gtf --sjdbOverhang 149
STAR --runThreadN 28 --runMode genomeGenerate --genomeDir /home/zhzhang/PG/refgenome/STAR_Mus_musculus_GRCm39/ --genomeFastaFiles /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa --sjdbGTFfile /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.chr.gtf --sjdbOverhang 149
STAR --runThreadN 28 --runMode genomeGenerate --genomeDir /home/zhzhang/PG/refgenome/STAR_Oryzias_latipes_ASM223467v1/ --genomeFastaFiles /home/zhzhang/PG/refgenome/Oryzias_latipes.ASM223467v1.dna.chr.fa --sjdbGTFfile /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryzias_latipes/mysql/Oryzias_latipes.ASM223467v1.108.chr.gtf --sjdbOverhang 149 --genomeSAindexNbases 13

STAR --runThreadN 28 --runMode genomeGenerate --genomeDir /share/home/zhzhang24/PG/refgenome/STAR_Mus_musculus_GRCm39_default/ --genomeFastaFiles /share/home/zhzhang24/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa --sjdbGTFfile /share/home/zhzhang24/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.chr.gtf

```
##### 3.Rawdata质控
```r
#斑马鱼
cd /home/zhzhang/PG/RNAseqdata/Danio_rerio/
for i in Zebrafish_*_R1.fq.gz
do
sample=${i%_R1.*}
trim_galore -j 28 --paired -q 20 --phred33 --length 20 -e 0.1 --stringency 3 -o /home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/ /home/zhzhang/PG/RNAseqdata/Danio_rerio/${sample}_R1.fq.gz /home/zhzhang/PG/RNAseqdata/Danio_rerio/${sample}_R2.fq.gz
done






#小鼠
cd /home/zhzhang/PG/RNAseqdata/Mus_musculus/
for i in Mouse_*_R1.fq.gz
do
sample=${i%_R1.*}
trim_galore -j 28 --paired -q 20 --phred33 --length 20 -e 0.1 --stringency 3 -o /home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/ /home/zhzhang/PG/RNAseqdata/Mus_musculus/${sample}_R1.fq.gz /home/zhzhang/PG/RNAseqdata/Mus_musculus/${sample}_R2.fq.gz
done



#青鳉
cd /home/zhzhang/PG/RNAseqdata/Oryzias_latipes/
for i in Medaka_*_R1.fq.gz
do
sample=${i%_R1.*}
trim_galore -j 28 --paired -q 20 --phred33 --length 20 -e 0.1 --stringency 3 -o /home/zhzhang/PG/RNAseqdata/Oryzias_latipes/qc_fq/ /home/zhzhang/PG/RNAseqdata/Oryzias_latipes/${sample}_R1.fq.gz /home/zhzhang/PG/RNAseqdata/Oryzias_latipes/${sample}_R2.fq.gz
done





#鸡
cd /home/zhzhang/PG/RNAseqdata/Gallus_gallus/
for i in Chicken_*_R1.fq.gz
do
sample=${i%_R1.*}
trim_galore -j 28 --paired -q 20 --phred33 --length 20 -e 0.1 --stringency 3 -o /home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/ /home/zhzhang/PG/RNAseqdata/Gallus_gallus/${sample}_R1.fq.gz /home/zhzhang/PG/RNAseqdata/Gallus_gallus/${sample}_R2.fq.gz
done



#人
cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens/
for i in $(seq 1 200)
do
q=`ls Human_*_R1.fq.gz|sed -n "${i}p"`
sample=${q%_R1.*}
trim_galore -j 28 --paired -q 20 --phred33 --length 20 -e 0.1 --stringency 3 -o /home/zhzhang/PG/RNAseqdata/Homo_sapiens/qc_fq/ /home/zhzhang/PG/RNAseqdata/Homo_sapiens/${sample}_R1.fq.gz /home/zhzhang/PG/RNAseqdata/Homo_sapiens/${sample}_R2.fq.gz
done

#人dd数据
cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/
for i in $(seq 1 313)
do
q=`ls Human_*.fq.gz|sed -n "${i}p"`
sample=${q%.fq.gz}
trim_galore -j 28 -q 20 --phred33 --length 20 -e 0.1 --stringency 3 -o /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/qc_fq/ /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/${sample}.fq.gz
done

#小鼠dd数据
cd /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/
for i in $(seq 1 317)
do
q=`ls Mouse_*.fq.gz|sed -n "${i}p"`
sample=${q%.fq.gz}
trim_galore -j 28 -q 20 --phred33 --length 20 -e 0.1 --stringency 3 -o /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/qc_fq/ /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/${sample}.fq.gz
done



```
```r
#小鼠衰老数据
#质控
cat "/share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/sample.tsv"|while read aaa
do
i=`echo -e "${aaa}"|awk '{print $1}'`
micromamba run -n SEQ trim_galore -j 64 --paired -q 20 --phred33 --length 20 -e 0.1 --stringency 3 -o /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/fq_qc/ /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/fq/${i}_R1.fq.gz /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/fq/${i}_R2.fq.gz
done


#比对
module load arm/star/2.7.11b
cat "/share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/sample.tsv"|while read aaa
do
i=`echo -e "${aaa}"|awk '{print $1}'`
STAR --runMode alignReads --runThreadN 64 --outSAMtype BAM SortedByCoordinate --genomeDir /share/home/zhzhang24/PG/refgenome/STAR_Mus_musculus_GRCm39 --readFilesIn /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/fq_qc/${i}_R1_val_1.fq.gz /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/fq_qc/${i}_R2_val_2.fq.gz --outFileNamePrefix /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/STAR_output/${i} --readFilesCommand zcat --twopassMode Basic --outBAMsortingThreadN 28 --outSAMmultNmax 1 --outSJfilterReads Unique --outSAMattrIHstart 0 --outSAMstrandField intronMotif
done




#比对结果质控（根据MAPQ过滤）
module load arm/samtools/1.21
cat "/share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/sample.tsv"|while read aaa
do
i=`echo -e "${aaa}"|awk '{print $1}'`
samtools view -@ 28 -b /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/STAR_output/${i}Aligned.sortedByCoord.out.bam -q 20 -o /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/bam/${i}.q20.bam
done


#筛选测序深度过低需要排除的样本
cat "/share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/sample.tsv"|while read aaa
do
i=`echo -e "${aaa}"|awk '{print $1}'`
grep "Uniquely mapped reads number" /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/STAR_output/${i}Log.final.out|awk -F '\t' -v sample="${i}" '{print sample"\t"$2}' >> /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/sample_UniqueReadNum.tsv
done
#汇总低质量样本
cat /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/sample_UniqueReadNum.tsv |awk '$2<4000000{print $1}' > /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/lowQ_sample.tsv
#移除低质量样本
cat "/share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/lowQ_sample.tsv"|while read i
do
mv "/share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/bam/${i}.q20.bam" /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/rm_bam/
done

#表达定量
cd /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/bam/
micromamba run -n SEQ featureCounts -T 64 -p --countReadPairs -t exon -g gene_id -a "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" -o /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/featureCounts_out/fc_aging_readscount.txt *.bam
#now
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/RNAseqdata/Mus_musculus_aging/featureCounts_out/fc_aging_readscount.txt zhzhang@122.205.95.67:/home/zhzhang/PG/RNAseq/Mus_musculus/

```


##### 4.比对
```r
#斑马鱼
cd /home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/
for i in Zebrafish_*_R1_val_1.fq.gz
do
sample=${i%_R1*}
STAR --runMode alignReads --runThreadN 28 --outSAMtype BAM SortedByCoordinate --genomeDir /home/zhzhang/PG/refgenome/STAR_Danio_rerio_GRCz11 --readFilesIn /home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/${sample}_R1_val_1.fq.gz /home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/${sample}_R2_val_2.fq.gz --outFileNamePrefix /home/zhzhang/PG/RNAseqdata/Danio_rerio/STAR_output/${sample} --readFilesCommand zcat --sjdbOverhang 149 --twopassMode Basic --outBAMsortingThreadN 20 --outSAMmultNmax 1 --outSJfilterReads Unique --outSAMattrIHstart 0 --outSAMstrandField intronMotif
done


#小鼠
cd /home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/
for i in Mouse_*_R1_val_1.fq.gz
do
sample=${i%_R1*}
STAR --runMode alignReads --runThreadN 28 --outSAMtype BAM SortedByCoordinate --genomeDir /home/zhzhang/PG/refgenome/STAR_Mus_musculus_GRCm39 --readFilesIn /home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/${sample}_R1_val_1.fq.gz /home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/${sample}_R2_val_2.fq.gz --outFileNamePrefix /home/zhzhang/PG/RNAseqdata/Mus_musculus/STAR_output/${sample} --readFilesCommand zcat --sjdbOverhang 149 --twopassMode Basic --outBAMsortingThreadN 20 --outSAMmultNmax 1 --outSJfilterReads Unique --outSAMattrIHstart 0 --outSAMstrandField intronMotif
done



#鸡
cd /home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/
for i in Chicken_*_R1_val_1.fq.gz
do
sample=${i%_R1*}
STAR --runMode alignReads --runThreadN 28 --outSAMtype BAM SortedByCoordinate --genomeDir /home/zhzhang/PG/refgenome/STAR_Gallus_gallus_GRCg7b --readFilesIn /home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/${sample}_R1_val_1.fq.gz /home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/${sample}_R2_val_2.fq.gz --outFileNamePrefix /home/zhzhang/PG/RNAseqdata/Gallus_gallus/STAR_output/${sample} --readFilesCommand zcat --sjdbOverhang 149 --twopassMode Basic --outBAMsortingThreadN 20 --outSAMmultNmax 1 --outSJfilterReads Unique --outSAMattrIHstart 0 --outSAMstrandField intronMotif
done


#人
cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens/qc_fq/
for i in $(seq 1 200)
do
q=`ls Human_*_R1_val_1.fq.gz|sed -n "${i}p"`
sample=${q%_R1*}
STAR --runMode alignReads --runThreadN 28 --outSAMtype BAM SortedByCoordinate --genomeDir /home/zhzhang/PG/refgenome/STAR_Homo_sapiens_GRCh38/ --readFilesIn /home/zhzhang/PG/RNAseqdata/Homo_sapiens/qc_fq/${sample}_R1_val_1.fq.gz /home/zhzhang/PG/RNAseqdata/Homo_sapiens/qc_fq/${sample}_R2_val_2.fq.gz --outFileNamePrefix /home/zhzhang/PG/RNAseqdata/Homo_sapiens/STAR_output/${sample} --readFilesCommand zcat --sjdbOverhang 149 --twopassMode Basic --outBAMsortingThreadN 20 --outSAMmultNmax 1 --outSJfilterReads Unique --outSAMattrIHstart 0 --outSAMstrandField intronMotif
done


#人dd
cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/qc_fq/
for i in $(seq 1 313)
do
q=`ls Human_*_trimmed.fq.gz|sed -n "${i}p"`
sample=${q%_trimmed.fq.gz}
STAR --runMode alignReads --runThreadN 28 --outSAMtype BAM SortedByCoordinate --genomeDir /home/zhzhang/PG/refgenome/STAR_Homo_sapiens_GRCh38/ --readFilesIn /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/qc_fq/${sample}_trimmed.fq.gz --outFileNamePrefix /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/STAR_output/${sample} --readFilesCommand zcat --sjdbOverhang 149 --twopassMode Basic --outBAMsortingThreadN 20 --outSAMmultNmax 1 --outSJfilterReads Unique --outSAMattrIHstart 0 --outSAMstrandField intronMotif
done




#小鼠dd
cd /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/qc_fq/
for i in $(seq 1 317)
do
q=`ls Mouse_*_trimmed.fq.gz|sed -n "${i}p"`
sample=${q%_trimmed.fq.gz}
STAR --runMode alignReads --runThreadN 28 --outSAMtype BAM SortedByCoordinate --genomeDir /home/zhzhang/PG/refgenome/STAR_Mus_musculus_GRCm39/ --readFilesIn /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/qc_fq/${sample}_trimmed.fq.gz --outFileNamePrefix /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/STAR_output/${sample} --readFilesCommand zcat --sjdbOverhang 149 --twopassMode Basic --outBAMsortingThreadN 20 --outSAMmultNmax 1 --outSJfilterReads Unique --outSAMattrIHstart 0 --outSAMstrandField intronMotif
done




```
##### 5.Bam质控
```r
#斑马鱼
cd /home/zhzhang/PG/RNAseqdata/Danio_rerio/STAR_output/
for i in Zebrafish_*Aligned.sortedByCoord.out.bam
do
sample=${i%Aligned*}
samtools view -@ 28 -b /home/zhzhang/PG/RNAseqdata/Danio_rerio/STAR_output/${sample}Aligned.sortedByCoord.out.bam -q 20 -o /home/zhzhang/PG/RNAseqdata/Danio_rerio/q20bam/${sample}.q20.bam
done


#小鼠Mus_musculus
cd /home/zhzhang/PG/RNAseqdata/Mus_musculus/STAR_output/
for i in Mouse_*Aligned.sortedByCoord.out.bam
do
sample=${i%Aligned*}
samtools view -@ 28 -b /home/zhzhang/PG/RNAseqdata/Mus_musculus/STAR_output/${sample}Aligned.sortedByCoord.out.bam -q 20 -o /home/zhzhang/PG/RNAseqdata/Mus_musculus/q20bam/${sample}.q20.bam
done



#鸡Gallus_gallus
cd /home/zhzhang/PG/RNAseqdata/Gallus_gallus/STAR_output/
for i in Chicken_*Aligned.sortedByCoord.out.bam
do
sample=${i%Aligned*}
samtools view -@ 120 -b /home/zhzhang/PG/RNAseqdata/Gallus_gallus/STAR_output/${sample}Aligned.sortedByCoord.out.bam -q 20 -o /home/zhzhang/PG/RNAseqdata/Gallus_gallus/q20bam/${sample}.q20.bam
done


#人
cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens/STAR_output/
for i in $(seq 1 200)
do
q=`ls Human_*Aligned.sortedByCoord.out.bam|sed -n "${i}p"`
sample=${q%Aligned*}
samtools view -@ 120 -b /home/zhzhang/PG/RNAseqdata/Homo_sapiens/STAR_output/${sample}Aligned.sortedByCoord.out.bam -q 20 -o /home/zhzhang/PG/RNAseqdata/Homo_sapiens/q20bam/${sample}.q20.bam
done


#人dd
cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/STAR_output/
for i in $(seq 1 313)
do
q=`ls Human_*Aligned.sortedByCoord.out.bam|sed -n "${i}p"`
sample=${q%Aligned*}
samtools view -@ 28 -b /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/STAR_output/${sample}Aligned.sortedByCoord.out.bam -q 20 -o /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/q20bam/${sample}.q20.bam
done




#小鼠dd
cd /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/STAR_output/
for i in $(seq 1 317)
do
q=`ls Mouse_*Aligned.sortedByCoord.out.bam|sed -n "${i}p"`
sample=${q%Aligned*}
samtools view -@ 28 -b /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/STAR_output/${sample}Aligned.sortedByCoord.out.bam -q 20 -o /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/q20bam/${sample}.q20.bam
done



#斑马鱼鸡小鼠转bw
Zebrafish 1345101833
Chicken 1041122857
Mouse 2723414844
#banmayu
cd /home/zhzhang/PG/RNAseqdata/Danio_rerio/q20bam/
find *.bam|while read i
do
sample=${i%.q20.bam}
bamCoverage -b /home/zhzhang/PG/RNAseqdata/Danio_rerio/q20bam/${sample}.q20.bam -o /home/zhzhang/PG/RNAseqdata/Danio_rerio/bw/${sample}.bw --binSize 10 --normalizeUsing BPM -p 28 --effectiveGenomeSize 1345101833
done
#xiaoshu
cd /home/zhzhang/PG/RNAseqdata/Mus_musculus/q20bam/
find *.bam|while read i
do
sample=${i%.q20.bam}
bamCoverage -b /home/zhzhang/PG/RNAseqdata/Mus_musculus/q20bam/${sample}.q20.bam -o /home/zhzhang/PG/RNAseqdata/Mus_musculus/bw/${sample}.bw --binSize 10 --normalizeUsing BPM -p 28 --effectiveGenomeSize 2723414844
done
#ji
cd /home/zhzhang/PG/RNAseqdata/Gallus_gallus/q20bam/
find *.bam|while read i
do
sample=${i%.q20.bam}
bamCoverage -b /home/zhzhang/PG/RNAseqdata/Gallus_gallus/q20bam/${sample}.q20.bam -o /home/zhzhang/PG/RNAseqdata/Gallus_gallus/bw/${sample}.bw --binSize 10 --normalizeUsing BPM -p 28 --effectiveGenomeSize 1041122857
done


#UCSC小工具bigWigMerge可合并多个bw文件，合并输出为bdg格式，合并方法为每个bin的信号值在全部文件中sum(也可参数取max)
#
cd /home/zhzhang/PG/RNAseqdata/Danio_rerio/bw/
find *bw|awk -F "_" '{print $1"_"$2"_"$3}'|sort|uniq|while read i
do
bigWigMerge /home/zhzhang/PG/RNAseqdata/Danio_rerio/bw/${i}*.bw /home/zhzhang/PG/RNAseqdata/Danio_rerio/bwmerge/${i}.bedGraph
sort -k1,1 -k2,2n /home/zhzhang/PG/RNAseqdata/Danio_rerio/bwmerge/${i}.bedGraph > /home/zhzhang/PG/RNAseqdata/Danio_rerio/bwmerge/${i}.sort.bedGraph
bedGraphToBigWig /home/zhzhang/PG/RNAseqdata/Danio_rerio/bwmerge/${i}.sort.bedGraph "/home/zhzhang/PG/refgenome/Danio_rerio.GRCz11.dna.chrsize.txt" /home/zhzhang/PG/RNAseqdata/Danio_rerio/bwmerge/${i}.bw
done
#
cd /home/zhzhang/PG/RNAseqdata/Mus_musculus/bw/
find *bw|awk -F "_" '{print $1"_"$2"_"$3}'|sort|uniq|while read i
do
bigWigMerge /home/zhzhang/PG/RNAseqdata/Mus_musculus/bw/${i}*.bw /home/zhzhang/PG/RNAseqdata/Mus_musculus/bwmerge/${i}.bedGraph
sort -k1,1 -k2,2n /home/zhzhang/PG/RNAseqdata/Mus_musculus/bwmerge/${i}.bedGraph > /home/zhzhang/PG/RNAseqdata/Mus_musculus/bwmerge/${i}.sort.bedGraph
bedGraphToBigWig /home/zhzhang/PG/RNAseqdata/Mus_musculus/bwmerge/${i}.sort.bedGraph "/home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chrsize.txt" /home/zhzhang/PG/RNAseqdata/Mus_musculus/bwmerge/${i}.bw
done
#
cd /home/zhzhang/PG/RNAseqdata/Gallus_gallus/bw/
find *bw|awk -F "_" '{print $1"_"$2"_"$3}'|sort|uniq|while read i
do
bigWigMerge /home/zhzhang/PG/RNAseqdata/Gallus_gallus/bw/${i}*.bw /home/zhzhang/PG/RNAseqdata/Gallus_gallus/bwmerge/${i}.bedGraph
sort -k1,1 -k2,2n /home/zhzhang/PG/RNAseqdata/Gallus_gallus/bwmerge/${i}.bedGraph > /home/zhzhang/PG/RNAseqdata/Gallus_gallus/bwmerge/${i}.sort.bedGraph
bedGraphToBigWig /home/zhzhang/PG/RNAseqdata/Gallus_gallus/bwmerge/${i}.sort.bedGraph "/home/zhzhang/PG/refgenome/Gallus_gallus.GRCg7b.dna.chrsize.txt" /home/zhzhang/PG/RNAseqdata/Gallus_gallus/bwmerge/${i}.bw
done


```
##### 6.GTF文件准备
```r
#排除掉假基因的参考基因组注释，追加写入三代和二代novel lnc的外显子水平注释
#鸡
cp "/home/zhzhang/PG/HPG/GTF_PG/Gallus_gallus.GRCg7b.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Gallus_gallus.GRCg7b.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Gallus_gallus.GRCg7b.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.3Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Gallus_gallus.GRCg7b.108.chr.rmpg.novellncRNA.gtf
#人
cp "/home/zhzhang/PG/HPG/GTF_PG/Homo_sapiens.GRCh38.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.3Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf
#小鼠
cp "/home/zhzhang/PG/HPG/GTF_PG/Mus_musculus.GRCm39.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.3Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf

#猕猴
cp "/home/zhzhang/PG/HPG/GTF_PG/Macaca_mulatta.Mmul_10.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Macaca_mulatta.Mmul_10.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Macaca_mulatta.Mmul_10.108.chr.rmpg.novellncRNA.gtf
#大鼠
cp "/home/zhzhang/PG/HPG/GTF_PG/Rattus_norvegicus.mRatBN7.2.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Rattus_norvegicus.mRatBN7.2.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Rattus_norvegicus.mRatBN7.2.108.chr.rmpg.novellncRNA.gtf
#兔子
cp "/home/zhzhang/PG/HPG/GTF_PG/Oryctolagus_cuniculus.OryCun2.0.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Oryctolagus_cuniculus.OryCun2.0.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Oryctolagus_cuniculus.OryCun2.0.108.chr.rmpg.novellncRNA.gtf
#负鼠
cp "/home/zhzhang/PG/HPG/GTF_PG/Monodelphis_domestica.ASM229v1.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Monodelphis_domestica.ASM229v1.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Monodelphis_domestica.ASM229v1.108.chr.rmpg.novellncRNA.gtf
#斑马鱼
cp "/home/zhzhang/PG/HPG/GTF_PG/Danio_rerio.GRCz11.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Danio_rerio.GRCz11.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Danio_rerio.GRCz11.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.3Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Danio_rerio.GRCz11.108.chr.rmpg.novellncRNA.gtf


```
##### 7.表达定量+定量结果整理为表达矩阵+标准化TPM
```r
cd /home/zhzhang/PG/RNAseqdata/Gallus_gallus/q20bam
featureCounts -T 28 -s 2 -p --countReadPairs -t exon -g gene_id -a "/home/zhzhang/PG/RNAseqdata/newGTF/Gallus_gallus.GRCg7b.108.chr.rmpg.novellncRNA.gtf" -o /home/zhzhang/PG/RNAseqdata/Gallus_gallus/featureCounts_out/fc_allsample_readscount.txt *.bam

cd /share/home/zhzhang24/PG/RNAseqdata/Danio_rerio/q20bam
micromamba run -n SEQ featureCounts -T 64 -s 2 -p --countReadPairs -t exon -g gene_id -a "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Danio_rerio.GRCz11.108.chr.rmpg.novellncRNA.gtf" -o /share/home/zhzhang24/PG/RNAseqdata/Danio_rerio/featureCounts_out/fc_allsample_readscount.txt *.bam

cd /home/zhzhang/PG/RNAseqdata/Mus_musculus/q20bam
featureCounts -T 28 -s 2 -p --countReadPairs -t exon -g gene_id -a "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" -o /home/zhzhang/PG/RNAseqdata/Mus_musculus/featureCounts_out/fc_allsample_readscount.txt *.bam



cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens/q20bam
featureCounts -T 64 -p --countReadPairs -t exon -g gene_id -a "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" -o /home/zhzhang/PG/RNAseqdata/Homo_sapiens/featureCounts_out/fc_allsample_readscount.txt *.bam


cd /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/q20bam
featureCounts -T 64 -t exon -g gene_id -a "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf" -o /home/zhzhang/PG/RNAseqdata/Homo_sapiens_dd/featureCounts_out/fc_ddsample_readscount.txt *.bam



cd /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/q20bam
featureCounts -T 64 -t exon -g gene_id -a "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf" -o /home/zhzhang/PG/RNAseqdata/Mus_musculus_dd/featureCounts_out/fc_ddsample_readscount.txt *.bam


```
```r
#R:定量结果整理为表达矩阵+标准化TPM
#fc_out为featruecounts对物种全部样本基因定量输出，generc_out为基因表达量readscount矩阵，gene_TPM_out为基因表达量TPM矩阵
chrpcgnum <- function(fc_out,generc_out,gene_TPM_out){
  #导入featruecounts输出
  fc_allsample_readscount <- read.delim(fc_out, comment.char="#")
  #生成readscount表达矩阵
  gene_rc <- select(fc_allsample_readscount,-c(2:6))%>%
    column_to_rownames("Geneid")
  colnamedf <- data.frame(cn=colnames(gene_rc))%>%
    separate(cn,c("cn","no"),sep = "[.]")%>%
    select(1)
  colnames(gene_rc) <- colnamedf$cn
  #储存
  data.table::fwrite(gene_rc,file =generc_out,sep = '\t',row.names = T,quote = F,col.names = T)
  #生成归一化矩阵
  #归一化长度
  gene_TPM <- gene_rc/(fc_allsample_readscount$Length/1000)
  #归一化深度
  for (i in c(1:ncol(gene_TPM))) {
    suofang <- sum(gene_TPM[,i])/1000000
    gene_TPM[,i] <- gene_TPM[,i]/suofang
  }
  #储存
  data.table::fwrite(gene_TPM,file =gene_TPM_out,sep = '\t',row.names = T,quote = F,col.names = T)
}


#鸡
chrpcgnum("~/PG/RNAseq/Gallus_gallus/fc_allsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Gallus_gallus/allsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Gallus_gallus/allsample_TPM.txt")
#斑马鱼
chrpcgnum("/home/zhzhang/PG/RNAseq/Danio_rerio/fc_allsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Danio_rerio/allsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Danio_rerio/allsample_TPM.txt")
#小鼠
chrpcgnum("~/PG/RNAseq/Mus_musculus/fc_allsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Mus_musculus/allsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Mus_musculus/allsample_TPM.txt")
#人
chrpcgnum("~/PG/RNAseq/Homo_sapiens/fc_allsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Homo_sapiens/allsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Homo_sapiens/allsample_TPM.txt")
#人dd data
chrpcgnum("/home/zhzhang/PG/RNAseq/Homo_sapiens/fc_ddsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Homo_sapiens/ddsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Homo_sapiens/ddsample_TPM.txt")
#小鼠dd data
chrpcgnum("/home/zhzhang/PG/RNAseq/Mus_musculus/fc_ddsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Mus_musculus/ddsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Mus_musculus/ddsample_TPM.txt")
#小鼠aging
chrpcgnum("/home/zhzhang/PG/RNAseq/Mus_musculus/fc_aging_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Mus_musculus/agingsample_readscount.txt",
          "/home/zhzhang/PG/RNAseq/Mus_musculus/agingsample_TPM.txt")


```
##### 8.基因分类文件(蛋白编码/lnc的geneid-分类对照list文件)
```r
#lncrna基因部分
cp "/share/home/zhzhang24/PG/lncRNA_class/new_r1/Gallus_gallus.lncRNA_class.txt" /share/home/zhzhang24/PG/RNAseqdata/newGTF/Gallus_gallus.geneid_class.txt
cp "/share/home/zhzhang24/PG/lncRNA_class/new_r1/Homo_sapiens.lncRNA_class.txt" /share/home/zhzhang24/PG/RNAseqdata/newGTF/Homo_sapiens.geneid_class.txt
cp "/share/home/zhzhang24/PG/lncRNA_class/new_r1/Mus_musculus.lncRNA_class.txt" /share/home/zhzhang24/PG/RNAseqdata/newGTF/Mus_musculus.geneid_class.txt

#蛋白编码基因部分
cat "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Gallus_gallus.GRCg7b.108.chr.rmpg.novellncRNA.gtf"|awk '$3== "gene" && $17=="gene_biotype" {print $10"\t"$18} $3== "gene" && $15=="gene_biotype" {print $10"\t"$16}'|sed "s/\"//g;s/\;//g"|grep protein_coding|awk '{print $1"\t""Protein-coding"}' >> /share/home/zhzhang24/PG/RNAseqdata/newGTF/Gallus_gallus.geneid_class.txt
cat "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf"|awk '$3== "gene" && $17=="gene_biotype" {print $10"\t"$18} $3== "gene" && $15=="gene_biotype" {print $10"\t"$16}'|sed "s/\"//g;s/\;//g"|grep "protein_coding"|awk '{print $1"\t""Protein-coding"}' >> /share/home/zhzhang24/PG/RNAseqdata/newGTF/Homo_sapiens.geneid_class.txt
cat "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf"|awk '$3== "gene" && $17=="gene_biotype" {print $10"\t"$18} $3== "gene" && $15=="gene_biotype" {print $10"\t"$16}'|sed "s/\"//g;s/\;//g"|grep "protein_coding"|awk '{print $1"\t""Protein-coding"}' >> /share/home/zhzhang24/PG/RNAseqdata/newGTF/Mus_musculus.geneid_class.txt

传至实验室服务器/home/zhzhang/PG/RNAseq/物种/下
```
##### 9.两类lncRNA与蛋白编码基因的表达比例(基因在某样本中TPM>=1视为在该样本中表达)
```r
#计算两类lnc基因与蛋白编码基因的表达比例
#1.在物种中三类基因的表达比例，只要在至少一个样本中表达(TPM>=1)即在物种中表达。
#3.在每种组织中三类基因的平均表达比例。柱。
#1.函数输出：三类基因（蛋白，两类lncrna）在不同阈值下在物种中的表达比例
#（a输入基因表达TPM矩阵，b输入基因分类文件，c输入物种名，d输出在物种中表达的基因id和类型文件）
three_gene_exppercent <- function(a,b,c,d){
  #导入基因表达TPM矩阵
  allsample_TPM <- read.delim(a, row.names=1)
  #导入基因ID分类
  geneid_class <- read.delim(b)
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
#小鼠
mm_spexp <- three_gene_exppercent(a ="/home/zhzhang/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                                  b = "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                                  c = "Mouse",
                                  d="~/PG/RNAseq/Mus_musculus/Mus_musculus.spexp_geneid_class.txt")
#人
hs_spexp <- three_gene_exppercent(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                                  b = "~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                                  c = "Human",
                                  d="~/PG/RNAseq/Homo_sapiens/Homo_sapiens.spexp_geneid_class.txt")
#合并
allsp_spexp <- rbind(mm_spexp,hs_spexp)
allsp_spexp$type <- factor(allsp_spexp$type,levels = c("Protein-coding","Pseudogene-associated sense lncRNA",
                                                       "Pseudogene-associated antisense lncRNA",
                                                       "Non-pseudogene-associated lncRNA"))
#储存
data.table::fwrite(allsp_spexp,file ="/home/zhzhang/PG/RNAseq/ALLSP_3gene_spexppercent.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#画图
#人
ph <- ggplot(data = filter(allsp_spexp,sp=="Human"),aes(x=threshold,y=spexpratio*100))+
  geom_point(aes(color=type),size=2)+
  geom_line(aes(group=type,color=type))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-associated lncRNA"),
                     labels=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA"))+
  theme_half_open()+
  scale_y_continuous(limits = c(70,105))+
  labs(x = "Expression cutoff (TPM)", y ="Expressed genes (%)",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "top", legend.direction = "vertical")+
  theme(legend.position = c(0.1, 0.9))
ggsave("/home/zhzhang/PG/RNAseq/plot/Human.pglnc_sp_exp_percent.pdf", 
       ph,width = 4.5, height = 4.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#小鼠
pm <- ggplot(data = filter(allsp_spexp,sp=="Mouse"),aes(x=threshold,y=spexpratio*100))+
  geom_point(aes(color=type),size=2)+
  geom_line(aes(group=type,color=type))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-associated lncRNA"),
                     labels=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA"))+
  theme_half_open()+
  scale_y_continuous(limits = c(20,85))+
  labs(x = "Expression cutoff (TPM)", y ="Expressed genes (%)",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "top", legend.direction = "vertical")+
  theme(legend.position = c(0.1, 0.7))
ggsave("/home/zhzhang/PG/RNAseq/plot/Mouse.pglnc_sp_exp_percent.pdf", 
       pm,width = 4.5, height = 4.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#2.函数输出每个样本中三类基因的表达比例。
#（a输入基因表达TPM矩阵，b输入基因分类文件，c输入物种名）
three_gene_sampleexppercent <- function(a,b,c){
  #导入基因表达TPM矩阵
  allsample_TPM <- read.delim(a, row.names=1)
  #导入基因ID分类
  geneid_class <- read.delim(b)
  #基因在样本中TPM>=1视为表达，制作基因表达的定性矩阵，表达为1，不表达为0
  gene_exp <- allsample_TPM
  for (i in c(1:ncol(gene_exp))) {
    gene_exp[gene_exp[,i]>=1,i] <- 1
    gene_exp[gene_exp[,i]<1,i] <- 0
  }
  #合并基因类型和基因表达的定性矩阵
  gene_exp <- left_join(geneid_class,rownames_to_column(gene_exp,"geneid"),by="geneid")
  #计算每个样本中coding基因的表达比例
  coding_gene <- filter(gene_exp,type=="Protein-coding")%>%
    select(-1,-2)
  coding_gene_sampleexpratio <- data.frame(expratio=apply(coding_gene,2,sum)/nrow(coding_gene),
                                           type="Protein-coding")%>%
    rownames_to_column("sample")
  #计算每个样本中paslnc基因的表达比例
  pglnc_gene <- filter(gene_exp,type=="Pseudogene-associated sense lncRNA")%>%
    select(-1,-2)
  pglnc_gene_sampleexpratio <- data.frame(expratio=apply(pglnc_gene,2,sum)/nrow(pglnc_gene),
                                          type="Pseudogene-associated sense lncRNA")%>%
    rownames_to_column("sample")
  #计算每个样本中paalnc基因的表达比例
  pglnc_gene2 <- filter(gene_exp,type=="Pseudogene-associated antisense lncRNA")%>%
    select(-1,-2)
  pglnc_gene_sampleexpratio2 <- data.frame(expratio=apply(pglnc_gene2,2,sum)/nrow(pglnc_gene2),
                                           type="Pseudogene-associated antisense lncRNA")%>%
    rownames_to_column("sample")
  #计算每个样本中npalnc基因的表达比例
  npglnc_gene <- filter(gene_exp,type=="Non-pseudogene-associated lncRNA")%>%
    select(-1,-2)
  npglnc_gene_sampleexpratio <- data.frame(expratio=apply(npglnc_gene,2,sum)/nrow(npglnc_gene),
                                           type="Non-pseudogene-associated lncRNA")%>%
    rownames_to_column("sample")
  #合并
  gene_sampleexpratio <- rbind(coding_gene_sampleexpratio,pglnc_gene_sampleexpratio,pglnc_gene_sampleexpratio2,npglnc_gene_sampleexpratio)
  #添加物种信息
  gene_sampleexpratio <- mutate(gene_sampleexpratio,sp=c)
  return(gene_sampleexpratio)
}
#计算各物种三类基因表达比例
#小鼠
mm_sampleexp <- three_gene_sampleexppercent(a = "/home/zhzhang/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                                            b = "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                                            c = "Mouse")
#人
hs_sampleexp <- three_gene_sampleexppercent(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                                            b = "~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                                            c = "Human")
#合并
allsp_sampleexp <- rbind(mm_sampleexp,hs_sampleexp)
allsp_sampleexp$type <- factor(allsp_sampleexp$type,levels = c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                                                               "Non-pseudogene-associated lncRNA"))
#统计储存
allsp_sampleexp_tj <- group_by(allsp_sampleexp,sp,type)%>%
  summarise(exp_percent_min=min(expratio),exp_percent_max=max(expratio),
            exp_percent_mean=mean(expratio),exp_percent_median=median(expratio))
data.table::fwrite(allsp_sampleexp_tj,file ="/home/zhzhang/PG/RNAseq/ALLSP_3gene_sampleexppercent.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#3.在每种组织中三类基因的平均表达比例。柱。
#样本表达比例信息添加组织信息
allsp_tissueexp <- separate(allsp_sampleexp,sample,c("no","tissue","sex"),sep="_",remove = F)%>%
  select(-no)
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="gonad" & allsp_tissueexp$sex=="male"] <- "testis"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="gonad" & allsp_tissueexp$sex=="female"] <- "ovary"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="gut"] <- "colon"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="bonemarrow"] <- "bone marrow"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="fallopiantube"] <- "fallopian tube"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="gallbladder"] <- "gall bladder"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="lymphnode"] <- "lymph node"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="salivarygland"] <- "salivary gland"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="skeletalmuscle"] <- "skeletal muscle"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="smallintestine"] <- "small intestine"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="smoothmuscle"] <- "smooth muscle"
allsp_tissueexp$tissue[allsp_tissueexp$tissue=="urinarybladder"] <- "urinary bladder"
mean_forboot <- function(data, index) {
  return(mean(data[index]))
}
set.seed(1024)
allsp_tissueexp <- select(allsp_tissueexp,-sex)%>%
  group_by(sp,tissue,type)%>%
  summarise(mean=mean(expratio),confmin=boot::boot.ci(boot::boot(expratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][4],
            confmax=boot::boot.ci(boot::boot(expratio, mean_forboot, R = 1000),conf=0.95,type=c('perc'))[["percent"]][5])
data.table::fwrite(allsp_tissueexp,file ="/home/zhzhang/PG/RNAseq/ALLSP_3gene_tissueexppercent.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot展示每种组织三类基因的平均表达比例
pmouse <- ggplot(data = filter(allsp_tissueexp,sp=="Mouse"),aes(x=tissue,y=mean*100))+
  geom_col(width=0.6,position="dodge",aes(fill=type))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Mouse",type=="Protein-coding"),
                position=position_nudge(x=-0.225),aes(ymax=confmax*100,ymin=confmin*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Mouse",type=="Pseudogene-associated sense lncRNA"),
                position=position_nudge(x=-0.075),aes(ymax=confmax*100,ymin=confmin*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Mouse",type=="Pseudogene-associated antisense lncRNA"),
                position=position_nudge(x=0.075),aes(ymax=confmax*100,ymin=confmin*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Mouse",type=="Non-pseudogene-associated lncRNA"),
                position=position_nudge(x=0.225),aes(ymax=confmax*100,ymin=confmin*100))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,20,40,60,80))+
  scale_x_discrete(labels = Hmisc::capitalize(distinct(filter(data.frame(allsp_tissueexp),sp=="Mouse"),tissue)$tissue))+
  labs(x = NULL, y ="Expression proportion (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  theme(axis.text.x = element_text(angle =30))
ggsave("/home/zhzhang/PG/RNAseq/plot/Mus_musculus.tissuelnc_exp_percent.pdf", 
       pmouse,width = 8, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
phuman <- ggplot(data = filter(allsp_tissueexp,sp=="Human"),aes(x=tissue,y=mean*100))+
  geom_col(width=0.6,position="dodge",aes(fill=type))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Human",type=="Protein-coding"),
                position=position_nudge(x=-0.225),aes(ymax=confmax*100,ymin=confmin*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Human",type=="Pseudogene-associated sense lncRNA"),
                position=position_nudge(x=-0.075),aes(ymax=confmax*100,ymin=confmin*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Human",type=="Pseudogene-associated antisense lncRNA"),
                position=position_nudge(x=0.075),aes(ymax=confmax*100,ymin=confmin*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Human",type=="Non-pseudogene-associated lncRNA"),
                position=position_nudge(x=0.225),aes(ymax=confmax*100,ymin=confmin*100))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,25,50,75,100),labels = c("0","25","50","75","100"))+
  scale_x_discrete(labels = Hmisc::capitalize(distinct(filter(data.frame(allsp_tissueexp),sp=="Human"),tissue)$tissue))+
  labs(x = NULL, y ="Expression proportion (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  theme(axis.text.x = element_text(angle =30))
ggsave("/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.tissuelnc_exp_percent.pdf", 
       phuman,width = 18, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
##### 10.两类lncRNA与蛋白编码基因（表达的）在物种/各组织中的最大表达量
```r
##画图对比各物种两类lnc基因的最大表达量
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
                      b = "/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.spexp_geneid_class.txt",
                      c = "Mouse")
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.spexp_geneid_class.txt",
                      c = "Human")
#不同物种三类基因表达量数据合并
ALLsp_TPM <- rbind(mm_TPM,hs_TPM)
ALLsp_TPM$type <- factor(ALLsp_TPM$type,levels =c("Protein-coding","Pseudogene-associated sense lncRNA",
                                                  "Pseudogene-associated antisense lncRNA",
                                                  "Non-pseudogene-associated lncRNA"))
#统计储存
tj <- group_by(ALLsp_TPM,sp,type)%>%
  summarise(mean=mean(log2(maxTPM)),median=median(log2(maxTPM)))
data.table::fwrite(tj,file ="/home/zhzhang/PG/RNAseq/ALLSP_3gene_maxexplog2TPM.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
pmm <- ggplot(data = filter(ALLsp_TPM,sp=="Mouse"),aes(x=type,y=log2(maxTPM)))+
  geom_boxplot(fatten = 3,notch = T,width=0.5,outlier.alpha = 0,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position = c(10,6.5,5),tip_length = 0.01,
              comparisons=list(c("Protein-coding","Pseudogene-associated sense lncRNA"),
                               c("Pseudogene-associated sense lncRNA","Non-pseudogene-associated lncRNA"),
                               c("Pseudogene-associated antisense lncRNA","Non-pseudogene-associated lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,2,4,6,8,10))+
  coord_cartesian(ylim = c(0, 11))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  labs(x = NULL, y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =30))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/RNAseq/plot/Mus_musculus.lnc_maxexp.pdf", 
       pmm,width = 4.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
phs <- ggplot(data = filter(ALLsp_TPM,sp=="Human"),aes(x=type,y=log2(maxTPM)))+
  geom_boxplot(fatten = 3,notch = T,width=0.5,outlier.alpha = 0,aes(fill=type))+
  ggsignif::geom_signif(map_signif_level=T,y_position = c(11.5,10,8.5,7),tip_length = 0.01,
                        comparisons=list(c("Protein-coding","Pseudogene-associated sense lncRNA"),
                                         c("Pseudogene-associated sense lncRNA","Non-pseudogene-associated lncRNA"),
                                         c("Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA"),
                                         c("Pseudogene-associated antisense lncRNA","Non-pseudogene-associated lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 13))+
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10,12.5),labels = c("0","2.5","5","7.5","10","12.5"))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  labs(x = NULL, y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =30))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.lnc_maxexp.pdf", 
       phs,width = 4.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
##### 
```r
#函数输出每个组织中表达的基因ID及其最大表达量。（list：tissue-ID-maxexp-type）
tissuemaxexp <- function(a,b,c){
  #导入基因ID分类
  geneid_class <- read.delim(b)
  #导入基因表达TPM矩阵
  allsample_TPM <- read.delim(a)
  colnames(allsample_TPM)[1] <- "geneid"
  lallsample_TPM <- gather(allsample_TPM,sample,exp,2:ncol(allsample_TPM))%>%
    separate(sample,c("sp","tissue","sex"),sep="_")
  #
  lallsample_TPM$tissue[lallsample_TPM$tissue=="gonad" & lallsample_TPM$sex=="male"] <- "testis"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="gonad" & lallsample_TPM$sex=="female"] <- "ovary"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="gut"] <- "colon"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="bonemarrow"] <- "bone marrow"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="fallopiantube"] <- "fallopian tube"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="gallbladder"] <- "gall bladder"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="lymphnode"] <- "lymph node"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="salivarygland"] <- "salivary gland"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="skeletalmuscle"] <- "skeletal muscle"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="smallintestine"] <- "small intestine"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="smoothmuscle"] <- "smooth muscle"
  lallsample_TPM$tissue[lallsample_TPM$tissue=="urinarybladder"] <- "urinary bladder"
  #每个组织中，每个基因最大表达量
  lallsample_TPM <- select(lallsample_TPM,-sp,-sex)%>%
    group_by(tissue,geneid)%>%
    summarise(exp=max(exp))
  #去除组织中不表达的基因
  unglallsample_TPM <- data.frame(lallsample_TPM)%>%
    filter(exp>=1)%>%
    left_join(geneid_class,by="geneid")%>%
    mutate(sp=c)
}
#人
hs_tissumaxexp <- tissuemaxexp(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                               b = "~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                               c = "Human")
#小鼠
mm_tissumaxexp <- tissuemaxexp(a = "~/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                               b = "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                               c = "Mouse")
#合并
allsp_sampleexp <- rbind(hs_tissumaxexp,mm_tissumaxexp)
allsp_sampleexp$type <- factor(allsp_sampleexp$type,levels = c("Protein-coding","Pseudogene-associated sense lncRNA",
                                                               "Pseudogene-associated antisense lncRNA",
                                                               "Non-pseudogene-associated lncRNA"))
allsp_sampleexp <- allsp_sampleexp[is.na(allsp_sampleexp$type)==F,]
#统计储存
allsp_sampleexp_tj <- group_by(allsp_sampleexp,sp,tissue,type)%>%
  summarise(median=median(exp))
data.table::fwrite(allsp_sampleexp_tj,file ="/home/zhzhang/PG/RNAseq/ALLSP_alltissue_3gene_maxTPM.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#wilcox test p值计算(pgdlnc vs npgdlnc)
pdata <- distinct(allsp_sampleexp_tj,sp,tissue)
for (i in 1:nrow(pdata)) {
  forsp <- pdata$sp[i]
  fortissue <- pdata$tissue[i]
  pva <- wilcox.test(filter(allsp_sampleexp,sp==forsp & tissue==fortissue & type=="Pseudogene-associated sense lncRNA")$exp,
                     filter(allsp_sampleexp,sp==forsp & tissue==fortissue & type=="Non-pseudogene-associated lncRNA")$exp)[["p.value"]]
  if (pva > 0.05) {
    annom <- "N.S."
  }
  if (pva < 0.05) {
    annom <- "*"
  }
  if (pva < 0.01) {
    annom <- "**"
  }
  if (pva < 0.001) {
    annom <- "***"
  }
  pdata[i,3] <- pva
  pdata[i,4] <- annom
  pva2 <- wilcox.test(filter(allsp_sampleexp,sp==forsp & tissue==fortissue & type=="Pseudogene-associated antisense lncRNA")$exp,
                     filter(allsp_sampleexp,sp==forsp & tissue==fortissue & type=="Non-pseudogene-associated lncRNA")$exp)[["p.value"]]
  if (pva2 > 0.05) {
    annom2 <- "N.S."
  }
  if (pva2 < 0.05) {
    annom2 <- "*"
  }
  if (pva2 < 0.01) {
    annom2 <- "**"
  }
  if (pva2 < 0.001) {
    annom2 <- "***"
  }
  pdata[i,5] <- pva2
  pdata[i,6] <- annom2
}
colnames(pdata)[c(3:6)] <- c("PASvsNPApvalue","PASvsNPAanno","PAAvsNPApvalue","PAAvsNPAanno")
data.table::fwrite(pdata,file ="/home/zhzhang/PG/RNAseq/ALLSP_alltissue_3gene_maxTPM.diffpvalue.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot每种组织，三类基因最大表达量
pmouse <- ggplot(data = filter(allsp_sampleexp,sp=="Mouse"),aes(x=tissue,y=log2(exp)))+
  geom_boxplot(fatten = 3,notch = T,width=0.6,outlier.alpha = 0,aes(fill=type))+
  ggsignif::geom_signif(annotations = filter(pdata,sp=="Mouse")$PASvsNPAanno,y_position = c(11),
              xmin = c(0.925:5.925),xmax = c(1.225:6.225),tip_length = 0)+
  ggsignif::geom_signif(annotations = filter(pdata,sp=="Mouse")$PAAvsNPAanno,y_position = c(10),
                        xmin = c(1.075:6.075),xmax = c(1.225:6.225),tip_length = 0)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 11))+
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10,12.5),labels = c("0","2.5","5","7.5","10","12.5"))+
  scale_x_discrete(labels = str_to_title(distinct(filter(data.frame(allsp_sampleexp),sp=="Mouse"),tissue)$tissue))+
  labs(x = NULL, y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  theme(axis.text.x = element_text(angle =30))
ggsave("/home/zhzhang/PG/RNAseq/plot/Mus_musculus.tissue_3gene_maxexp.pdf", 
       pmouse,width = 8, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
phuman <- ggplot(data = filter(allsp_sampleexp,sp=="Human"),aes(x=tissue,y=log2(exp)))+
  geom_boxplot(fatten = 3,notch = T,width=0.6,outlier.alpha = 0,aes(fill=type))+
  ggsignif::geom_signif(annotations = filter(pdata,sp=="Human")$PASvsNPAanno,y_position = c(11),
                        xmin = c(0.925:31.925),xmax = c(1.225:32.225),tip_length = 0)+
  ggsignif::geom_signif(annotations = filter(pdata,sp=="Human")$PAAvsNPAanno,y_position = c(10),
                        xmin = c(1.075:32.075),xmax = c(1.225:32.225),tip_length = 0)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 11))+
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10,12.5),labels = c("0","2.5","5","7.5","10","12.5"))+
  scale_x_discrete(labels = Hmisc::capitalize(distinct(filter(data.frame(allsp_sampleexp),sp=="Human"),tissue)$tissue))+
  labs(x = NULL, y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  theme(axis.text.x = element_text(angle =30))
ggsave("/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.tissue_3gene_maxexp.pdf", 
       phuman,width = 18, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")




```




##### 10.5 严格阈值对表达水平影响
```r
##画图对比各物种两类lnc基因的最大表达量
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
ALLsp_TPM <- tgenemaxTPM(a = "/home/zhzhang/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.spexp_geneid_class.txt",
                      c = "Mouse")
d="~/PG/lncRNA_class/new_r1/Mus_musculus.palncRNA_message.txt"
e="/home/zhzhang/PG/RNAseq/plot/Mus_musculus.paslnc_maxexp.pdf"
f="/home/zhzhang/PG/RNAseq/plot/Mus_musculus.paalnc_maxexp.pdf"
#人类
ALLsp_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.spexp_geneid_class.txt",
                      c = "Human")
d="~/PG/lncRNA_class/new_r1/Homo_sapiens.palncRNA_message.txt"
e="/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.paslnc_maxexp.pdf"
f="/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.paalnc_maxexp.pdf"
#导入lncRNA信息
palncRNA_message <- read.delim(d)%>%
  select(2,4)%>%
  dplyr::rename(geneid=lncgid)
ALLsp_TPM <- left_join(palncRNA_message,ALLsp_TPM,by="geneid")
#
pas <- rbind(mutate(filter(ALLsp_TPM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(ALLsp_TPM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(ALLsp_TPM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(ALLsp_TPM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(ALLsp_TPM,type=="Pseudogene-associated sense lncRNA"&total_overlap_len>200),cuttype="200"))
paa <- rbind(mutate(filter(ALLsp_TPM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>0),cuttype="0"),
             mutate(filter(ALLsp_TPM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>50),cuttype="50"),
             mutate(filter(ALLsp_TPM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>100),cuttype="100"),
             mutate(filter(ALLsp_TPM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>150),cuttype="150"),
             mutate(filter(ALLsp_TPM,type=="Pseudogene-associated antisense lncRNA"&total_overlap_len>200),cuttype="200"))
#合并
alldata <- rbind(pas,paa)%>%
  mutate(type=case_when(type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA"))
alldata$type <- factor(alldata$type,levels = c("PAS lncRNA","PAA lncRNA"))
alldata$cuttype <- factor(alldata$cuttype,levels = c("0","50","100","150","200"))
#统计储存
tj <- group_by(alldata,type,cuttype)%>%
  summarise(mean=mean(log2(maxTPM)),median=median(log2(maxTPM)))
#plot
ph1 <- ggplot(data=filter(alldata,type=="PAS lncRNA"), aes(x=cuttype,y=log2(maxTPM)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAS lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,test.args = c("greater"),
                        comparisons = list(c("200","0")),
                        y_position=c(7),tip_length = 0.01,size=0.5,textsize=4)+#human:y_position=c(8)
  scale_fill_manual(values=c("#7197AD","#70A2BF","#6CADD3","#5FB3E5","#4DB9F8"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,2,4,6,8,10))+
  coord_cartesian(ylim = c(0, 8))+#human:ylim = c(0, 9)
  labs(x = "Cutoff of\noverlap length (bp)", y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(e, 
       ph1,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
ph2 <- ggplot(data=filter(alldata,type=="PAA lncRNA"), aes(x=cuttype,y=log2(maxTPM)))+
  geom_boxplot(fatten = 2,outlier.alpha = 0,width=0.4,notch=T,aes(fill=cuttype))+
  geom_point(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean),shape=21,size=0.8)+
  geom_line(data =filter(tj,type=="PAA lncRNA") ,aes(y=mean,group=type),linetype="dashed",size=0.3)+
  ggsignif::geom_signif(map_signif_level = T,test.args = c("greater"),
                        comparisons = list(c("200","0")),
                        y_position=c(6),tip_length = 0.01,size=0.5,textsize=4)+#human:y_position=c(7)
  scale_fill_manual(values=c("#628255","#6D9B5B","#77B65D","#7DD15B","#80EC54"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,2,4,6,8,10))+
  coord_cartesian(ylim = c(0, 7))+#human:ylim = c(0, 8)
  labs(x = "Cutoff of\noverlap length (bp)", y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x= element_text(size = 14),
        axis.text.y  = element_text(size = 13),
        axis.text.x = element_text(size = 12),legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(f, 
       ph2,width = 2.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```




##### 11.lnc组织特异性和组织分布情况
```r
#函数输出每个基因在每种组织中是否表达(在某组织至少一个样本中表达即为在该组织中表达)
#（a输入基因表达TPM矩阵，b输入基因分类文件，c输出路径）
three_gene_sampleexppercent <- function(a,b,c){
  #导入基因表达TPM矩阵
  allsample_TPM <- read.delim(a, row.names=1)
  #导入基因ID分类
  geneid_class <- read.delim(b)
  #基因在样本中TPM>=1视为表达，制作基因表达的定性矩阵，表达为1，不表达为0
  gene_exp <- allsample_TPM
  for (i in c(1:ncol(gene_exp))) {
    gene_exp[gene_exp[,i]>=1,i] <- 1
    gene_exp[gene_exp[,i]<1,i] <- 0
  }
  #合并基因类型和基因表达的定性矩阵,并转换为基因在不同组织中的表达情况
  gene_matrix <- left_join(geneid_class,rownames_to_column(gene_exp,"geneid"),by="geneid")
  gene_tissue_exp <- gather(gene_matrix,sample,exp,3:ncol(gene_matrix))%>%
    separate(sample,c("sp","tissue","sex"),sep="_")
  gene_tissue_exp$tissue[gene_tissue_exp$tissue=="gonad" & gene_tissue_exp$sex=="male"] <- "testis"
  gene_tissue_exp$tissue[gene_tissue_exp$tissue=="gonad" & gene_tissue_exp$sex=="female"] <- "ovary"
  gene_tissue_exp$tissue[gene_tissue_exp$tissue=="gut"] <- "colon"
  gene_tissue_exp <- group_by(gene_tissue_exp,geneid,type,tissue)%>%
    summarise(exp=max(exp))
  #储存
  data.table::fwrite(gene_tissue_exp,file =c,sep = '\t',row.names = F,quote = F,col.names = T)
}
#输出每个基因在每种组织中是否表达
#小鼠
mm_sampleexp <- three_gene_sampleexppercent(a = "~/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                                            b = "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                                            c = "/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.gene_type_tissueexp.txt")
#人
hs_sampleexp <- three_gene_sampleexppercent(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                                            b = "~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                                            c = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.gene_type_tissueexp.txt")
#对比三类基因表达的广泛程度(人)
hs_gene_type_tissueexp <- read.delim("~/PG/RNAseq/Homo_sapiens/Homo_sapiens.gene_type_tissueexp.txt")%>%
  group_by(geneid,type)%>%
  summarise(tnum=sum(exp))
#储存
data.table::fwrite(hs_gene_type_tissueexp,file ="/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.gene_tissuenum.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
hs_gene_type_tissueexp <- data.frame(hs_gene_type_tissueexp)%>%
  filter(tnum!=0)
hs_gene_type_tissueexp$type <- factor(hs_gene_type_tissueexp$type,
                                      levels = rev(c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                                                 "Non-pseudogene-associated lncRNA")))
library(gghalves)
wp1 <- signif(wilcox.test(filter(hs_gene_type_tissueexp,type=="Pseudogene-associated sense lncRNA")$tnum,
                          filter(hs_gene_type_tissueexp,type=="Non-pseudogene-associated lncRNA")$tnum)[["p.value"]],
              2)
wp12 <- signif(wilcox.test(filter(hs_gene_type_tissueexp,type=="Pseudogene-associated antisense lncRNA")$tnum,
                          filter(hs_gene_type_tissueexp,type=="Non-pseudogene-associated lncRNA")$tnum)[["p.value"]],
              2)
wp2 <- signif(wilcox.test(filter(hs_gene_type_tissueexp,type=="Pseudogene-associated sense lncRNA")$tnum,
                          filter(hs_gene_type_tissueexp,type=="Protein-coding")$tnum)[["p.value"]],
              2)
hs <- ggplot(data = hs_gene_type_tissueexp,aes(x=type,y=tnum))+
  geom_half_violin(width=1,side = "r",position = position_nudge(x=0.2),aes(color=type,fill=type))+
  geom_point(data = filter(hs_gene_type_tissueexp,type=="Protein-coding"),stroke=0,
             alpha=0.8/(18264/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(hs_gene_type_tissueexp,type=="Non-pseudogene-associated lncRNA"),stroke=0,
             alpha=0.8/(25582/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(hs_gene_type_tissueexp,type=="Pseudogene-associated sense lncRNA"),stroke=0,
             alpha=0.8,size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(hs_gene_type_tissueexp,type=="Pseudogene-associated antisense lncRNA"),stroke=0,
             alpha=0.8,size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_boxplot(fatten=3,width=0.05,position = position_nudge(x=0.2),fill="white",outlier.alpha = 0)+
  scale_x_discrete(labels = c("NPA lncRNA","PAA lncRNA","PAS lncRNA","Protein-coding"))+
  coord_flip()+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-associated lncRNA"))+
  theme_half_open()+
  labs(x =NULL, y ="Number of tissues",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  ggsignif::geom_signif(annotations=c("***","***","***"),y_position=c(33,35,33),tip_length = 0,
              xmin = c(1.2,1.2,3.2),xmax = c(1.8,2.8,3.8))
ggsave("/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.lnc_tissuenum.pdf", 
       hs,width = 7, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#对比三类基因表达的广泛程度(小鼠)
mm_gene_type_tissueexp <- read.delim("~/PG/RNAseq/Mus_musculus/Mus_musculus.gene_type_tissueexp.txt")%>%
  group_by(geneid,type)%>%
  summarise(tnum=sum(exp))
#plot
mm_gene_type_tissueexp <- data.frame(mm_gene_type_tissueexp)%>%
  filter(tnum!=0)
mm_gene_type_tissueexp$type <- factor(mm_gene_type_tissueexp$type,
                                      levels = rev(c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                                                     "Non-pseudogene-associated lncRNA")))
wp1 <- signif(wilcox.test(filter(mm_gene_type_tissueexp,type=="Pseudogene-associated sense lncRNA")$tnum,
                          filter(mm_gene_type_tissueexp,type=="Non-pseudogene-associated lncRNA")$tnum)[["p.value"]],
              2)
wp12 <- signif(wilcox.test(filter(mm_gene_type_tissueexp,type=="Pseudogene-associated antisense lncRNA")$tnum,
                           filter(mm_gene_type_tissueexp,type=="Non-pseudogene-associated lncRNA")$tnum)[["p.value"]],
               2)
wp2 <- signif(wilcox.test(filter(mm_gene_type_tissueexp,type=="Pseudogene-associated sense lncRNA")$tnum,
                          filter(mm_gene_type_tissueexp,type=="Protein-coding")$tnum)[["p.value"]],
              2)
mm <- ggplot(data = mm_gene_type_tissueexp,aes(x=type,y=tnum))+
  geom_half_violin(width=1,side = "r",position = position_nudge(x=0.2),aes(color=type,fill=type))+
  geom_point(data = filter(mm_gene_type_tissueexp,type=="Protein-coding"),stroke=0,
             alpha=0.8/(18264/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(mm_gene_type_tissueexp,type=="Non-pseudogene-associated lncRNA"),stroke=0,
             alpha=0.8/(25582/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(mm_gene_type_tissueexp,type=="Pseudogene-associated sense lncRNA"),stroke=0,
             alpha=0.8,size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(mm_gene_type_tissueexp,type=="Pseudogene-associated antisense lncRNA"),stroke=0,
             alpha=0.8,size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_boxplot(fatten=3,width=0.05,position = position_nudge(x=0.2),fill="white",outlier.alpha = 0)+
  scale_y_continuous(breaks = c(0,2,4,6))+
  scale_x_discrete(labels = c("NPA lncRNA","PAA lncRNA","PAS lncRNA","Protein-coding"))+
  coord_flip()+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-associated lncRNA"))+
  theme_half_open()+
  labs(x =NULL, y ="Number of tissues",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  ggsignif::geom_signif(annotations=c("N.S.","N.S.","***"),y_position=c(6.5,7,6.5),tip_length = 0,
                        xmin = c(1.2,1.2,3.2),xmax = c(1.8,2.8,3.8))
ggsave("/home/zhzhang/PG/RNAseq/plot/Mus_musculus.lnc_tissuenum.pdf", 
       mm,width = 7, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```

