### 五.表观基因组特征分析
##### 1.基因TSS附近表观修饰
```r
ATAC bwwenjian（hg19）
##################################
cat "/home/zhzhang/PG/epigen/epimap/epimap_ATAC.txt"|while read i
do
file=${i}
wget https://epigenome.wustl.edu/epimap/data/averagetracks_pergroup/${file}
echo "${file}" >> /home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.all_ATAC.log
done
##################################
bigWigMerge /home/zhzhang/PG/epigen/epimap/ATAC/*.bigWig /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.bedGraph
liftOver /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.bedGraph /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.hg38.bedGraph /home/zhzhang/PG/epigen/epimap/ATAC/unmap.txt
cat /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.hg38.bedGraph|sed 's/^chr//g'|bedtools sort|grep -w -Ff "/home/zhzhang/PG/epigen/Homo_sapiens.GRCh38.chr.txt"|awk '$1!="M" {print $0}'|grep -v "alt"|grep -v "random"|grep -v "Un" > /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.grch38.bedGraph
bedtools merge -i /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.grch38.bedGraph -d -1 -c 4 -o mean > /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.grch38.merge.bedGraph
bedGraphToBigWig /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.grch38.merge.bedGraph "/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt" /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.grch38.merge.bw
#deeptools computeMatrix 
computeMatrix reference-point --referencePoint TSS -R /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pgdlncgeneTSS.bed /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.npgdlncgeneTSS.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed -S /home/zhzhang/PG/epigen/epimap/ATAC/all_ATACsum.grch38.merge.bw -b 5000 -a 5000 --binSize 10 -o /home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.ATAC_3geneTSS.matrix.gz -p 64 --skipZeros
#deeptools 可视化
zcat /home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.ATAC_3geneTSS.matrix.gz|sed 's/nan/0/g' > /home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.ATAC_3geneTSS.matrix.ok.txt
gzip /home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.ATAC_3geneTSS.matrix.ok.txt
plotProfile -m "/home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.ATAC_3geneTSS.matrix.ok.txt.gz" -o /home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.ATAC_3geneTSS.png --dpi 750 --legendLocation upper-right






bigWigToBedGraph "/home/zhzhang/PG/epigen/epimap/ATAC/average_ATAC-seq_imputed_Endocrine.bigWig" "/home/zhzhang/PG/epigen/epimap/ATAC/average_ATAC-seq_imputed_Endocrine.bedGraph"
liftOver /home/zhzhang/PG/epigen/epimap/ATAC/average_ATAC-seq_imputed_Endocrine.bedGraph /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/epigen/epimap/ATAC/average_ATAC-seq_imputed_Endocrine.hg38.bedGraph /home/zhzhang/PG/epigen/epimap/ATAC/unmap.txt
cat /home/zhzhang/PG/epigen/epimap/ATAC/average_ATAC-seq_imputed_Endocrine.hg38.bedGraph|sed 's/^chr//g'|bedtools sort|grep -w -Ff "/home/zhzhang/PG/epigen/Homo_sapiens.GRCh38.chr.txt"|awk '$1!="M" {print $0}'|grep -v "alt"|grep -v "random"|grep -v "Un"|bedtools merge -d -1 -c 4 -o mean > /home/zhzhang/PG/epigen/epimap/ATAC/average_ATAC-seq_imputed_Endocrine.hg38.ok.bedGraph
bedGraphToBigWig /home/zhzhang/PG/epigen/epimap/ATAC/average_ATAC-seq_imputed_Endocrine.hg38.ok.bedGraph "/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt" /home/zhzhang/PG/epigen/epimap/ATAC/average_ATAC-seq_imputed_Endocrine.hg38.ok.bw
#deeptools computeMatrix 
computeMatrix reference-point --referencePoint TSS -R /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pgdlncgeneTSS.bed /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.npgdlncgeneTSS.bed /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed -S /home/zhzhang/PG/epigen/epimap/ATAC/average_ATAC-seq_imputed_Endocrine.hg38.ok.bw -b 5000 -a 5000 --binSize 10 -o /home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.EndocrineATAC_3geneTSS.matrix.gz -p 64 --skipZeros
#deeptools 可视化
plotProfile -m /home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.EndocrineATAC_3geneTSS.matrix.gz -o /home/zhzhang/PG/epigen/epimap/ATAC/Homo_sapiens.EndocrineATAC_3geneTSS.png --dpi 750 --legendLocation upper-right



```
##### TSS bed
```r
#准备基因TSS bed文件
cp "/home/zhzhang/PG/Evolution/conserve/Homo_sapiens.allgene.bed" /home/zhzhang/PG/epigen/annotationbed/
#各类型基因id文件
#paslnc基因ID
grep -w "Pseudogene-associated sense" "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Homo_sapiens.geneid_class.txt"|awk '{print $1}' > /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgene.id.txt
#paalnc基因ID
grep -w "Pseudogene-associated antisense" "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Homo_sapiens.geneid_class.txt"|awk '{print $1}' > /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgene.id.txt
#npalnc基因ID
grep -w "Non-pseudogene-associated" "/share/home/zhzhang24/PG/RNAseqdata/newGTF/Homo_sapiens.geneid_class.txt"|awk '{print $1}' > /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgene.id.txt
#蛋白编码基因ID
grep -w Protein-coding "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.geneid_class.txt"|awk '{print $1}' > /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pcgene.id.txt
#根据各类基因IDgrep出各类型基因TSS bed
#paslnc基因bed
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.allgene.bed"|grep -w -Ff /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgene.id.txt|awk '$6=="+" {print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} $6=="-" {print $1"\t"$3"\t"$3+1"\t"$4"\t"$5"\t"$6}'|bedtools sort > /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.allgene.bed"|grep -w -Ff /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgene.id.txt|bedtools sort > /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgene.bed
#paalnc基因bed
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.allgene.bed"|grep -w -Ff /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgene.id.txt|awk '$6=="+" {print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} $6=="-" {print $1"\t"$3"\t"$3+1"\t"$4"\t"$5"\t"$6}'|bedtools sort > /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.allgene.bed"|grep -w -Ff /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgene.id.txt|bedtools sort > /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgene.bed
#npalnc基因bed
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.allgene.bed"|grep -w -Ff /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgene.id.txt|awk '$6=="+" {print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} $6=="-" {print $1"\t"$3"\t"$3+1"\t"$4"\t"$5"\t"$6}'|bedtools sort > /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.allgene.bed"|grep -w -Ff /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgene.id.txt|bedtools sort > /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgene.bed
#蛋白编码基因bed
cat "/home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.allgene.bed"|grep -w -Ff /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pcgene.id.txt|awk '$6=="+" {print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} $6=="-" {print $1"\t"$3"\t"$3+1"\t"$4"\t"$5"\t"$6}'|bedtools sort > /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed
cat "/home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.allgene.bed"|grep -w -Ff /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pcgene.id.txt|bedtools sort > /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pcgene.bed



#生成作为对照的基因间区
#生成随机的20000个6kb区域bed文件
bedtools random -n 20000 -l 6000 -seed 1024 -g "/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt" > /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.20000random_6kb_region.bed
bedtools shuffle -i /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.20000random_6kb_region.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -incl /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.intergenic.bed -noOverlapping|bedtools sort|awk '{print $1"\t"$2"\t"$3"\tintergenic_"$4"\tintergenic\t"$6}' > /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.20000random_6kb_intergenic.bed
#将基因间区的中点作为TSS阴性对照
cat /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.20000random_6kb_intergenic.bed|awk '{print $1"\t"$2+3000"\t"$2+3001"\t"$4"\t"$5"\t"$6}'|bedtools sort > /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed
#将基因间区的中心1000bp作为genebody阴性对照
cat /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.20000random_6kb_intergenic.bed|awk '{print $1"\t"$2+2500"\t"$2+3500"\t"$4"\t"$5"\t"$6}'|bedtools sort > /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.intergenicbody.bed



```
##### roadmap数据库下载数据
```r
#roadmap数据库下载数据(127 consolidated epigenomes：来源于183个样本，类型相同的样本数据被合并产生111个合并表观基因组，加上ENCODE项目的16个表观基因组)
#127 epigenomes染色质15种状态信息bed下载（人类）hg38
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/all_hg38lift.mnemonics.bedFiles.tgz
tar -xzvf "/home/zhzhang/PG/epigen/all_hg38lift.mnemonics.bedFiles.tgz" -C /home/zhzhang/PG/epigen/state/
gzip -d /home/zhzhang/PG/epigen/state/*
mv /home/zhzhang/PG/epigen/state/*.bed /home/zhzhang/PG/epigen/state/allstate/



#全部表观测序 peak bed下载（人类）hg19【先下载全部call 窄峰的组蛋白peak（Dnase仅下载MACS2的peak不下载hotspot的），再分别处理】
cd /home/zhzhang/PG/epigen/all_narrowPeakhistone_hg19/
##################################
grep -v hotspot /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i}
wget https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/${file}
echo "${file}" >> /home/zhzhang/PG/epigen/Homo_sapiens.all_narrowPeakhistone_hg19.log
done
##################################
gzip -d /home/zhzhang/PG/epigen/all_narrowPeakhistone_hg19/*.gz



#全部甲基化分数(37 WGBS和51 RRBS)bw下载（人类）hg19
cd /home/zhzhang/PG/epigen/dnamethylation/WGBS/
##################################
cat "/home/zhzhang/PG/epigen/dnamethylation/WGBSfile.txt"|while read i
do
file=${i}
wget https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/${file}
echo "${file}" >> /home/zhzhang/PG/epigen/dnamethylation/Homo_sapiens.all_WGBS.log
done
##################################

cd /home/zhzhang/PG/epigen/dnamethylation/RRBS/
##################################
cat "/home/zhzhang/PG/epigen/dnamethylation/RRBSfile.txt"|while read i
do
file=${i}
wget https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/RRBS/FractionalMethylation_bigwig/${file}
echo "${file}" >> /home/zhzhang/PG/epigen/dnamethylation/Homo_sapiens.all_RRBS.log
done
##################################



```


##### DNA甲基化
```r
#bwmerge生成全部表观基因组每个bin甲基化比率（Fractional methylation）的加和（bdg文件），第四列除以样本数，变成甲基化比率在所有表观基因组里的均值。最后lift over转化到hg38，bedtools merge用-c mean合并lift后重叠的区域，最后转为bw文件
bigWigMerge /home/zhzhang/PG/epigen/dnamethylation/WGBS/*.bigwig /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.bedGraph
bigWigMerge /home/zhzhang/PG/epigen/dnamethylation/RRBS/*.bigwig /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/allRRBS.bedGraph
cat /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.bedGraph|awk '{print $1"\t"$2"\t"$3"\t"$4/37}' > /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.mean_Fractional_methylation.bedGraph
cat /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/allRRBS.bedGraph|awk '{print $1"\t"$2"\t"$3"\t"$4/51}' > /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/allRRBS.mean_Fractional_methylation.bedGraph
liftOver /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.mean_Fractional_methylation.bedGraph /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.mean_Fractional_methylation.hg38.bedGraph /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/unmap.txt
liftOver /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/allRRBS.mean_Fractional_methylation.bedGraph /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/allRRBS.mean_Fractional_methylation.hg38.bedGraph /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/unmap.txt
cat /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.mean_Fractional_methylation.hg38.bedGraph|sed 's/^chr//g'|bedtools sort|awk '$1!="M" {print $0}'|grep -v "alt"|grep -v "random"|grep -v "Un"|bedtools merge -d -1 -c 4 -o mean > /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.mean_Fractional_methylation.merge.hg38.bedGraph
cat /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/allRRBS.mean_Fractional_methylation.hg38.bedGraph|sed 's/^chr//g'|bedtools sort|awk '$1!="M" {print $0}'|grep -v "alt"|grep -v "random"|grep -v "Un"|bedtools merge -d -1 -c 4 -o mean > /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/allRRBS.mean_Fractional_methylation.merge.hg38.bedGraph
bedGraphToBigWig /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.mean_Fractional_methylation.merge.hg38.bedGraph "/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt" /home/zhzhang/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.mean_Fractional_methylation.merge.hg38.bw
bedGraphToBigWig /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/allRRBS.mean_Fractional_methylation.merge.hg38.bedGraph "/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt" /home/zhzhang/PG/epigen/dnamethylation/RRBS/mergebw/allRRBS.mean_Fractional_methylation.merge.hg38.bw


```
```r
#生成三类基因TSS上下游各2kb区域每200bp一个bin的bed文件[五列：chr start end geneid bin相对于TSS上下游各2kb区域的start和end的中值]
rm /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
#蛋白编码基因
a=100
start=2000
end=1800
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2-start"\t"$2-end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start-200]
end=$[$end-200]
done
a=2100
start=0
end=200
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2+start"\t"$2+end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start+200]
end=$[$end+200]
done
#paslnc基因(照搬上述，改下TSS文件路径)
a=100
start=2000
end=1800
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2-start"\t"$2-end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start-200]
end=$[$end-200]
done
a=2100
start=0
end=200
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2+start"\t"$2+end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start+200]
end=$[$end+200]
done
#paalnc基因(照搬上述，改下TSS文件路径)
a=100
start=2000
end=1800
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2-start"\t"$2-end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start-200]
end=$[$end-200]
done
a=2100
start=0
end=200
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2+start"\t"$2+end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start+200]
end=$[$end+200]
done
#npalnc基因(照搬上述，改下TSS文件路径)
a=100
start=2000
end=1800
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2-start"\t"$2-end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start-200]
end=$[$end-200]
done
a=2100
start=0
end=200
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2+start"\t"$2+end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start+200]
end=$[$end+200]
done
#基因间区(照搬上述，改下TSS文件路径)
a=100
start=2000
end=1800
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2-start"\t"$2-end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start-200]
end=$[$end-200]
done
a=2100
start=0
end=200
for i in $(seq 1 10)
do
cat "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed"|awk -v a=$a -v start=$start -v end=$end '{print $1"\t"$2+start"\t"$2+end"\t"$4"\t"a}' >> /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed
a=$[$a+200]
start=$[$start+200]
end=$[$end+200]
done
#第四列修改为不重复（防止bigWigAverageOverBed报错）
cat /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.bed|awk '{print $1"\t"$2"\t"$3"\t"$4"__"$5}' > /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.4nodup.bed

#bigWigAverageOverBed确定每个基因TSS+—2kb区域每个bin的平均甲基化分数
micromamba run -n SEQ bigWigAverageOverBed /share/home/zhzhang24/PG/epigen/dnamethylation/WGBS/mergebw/allWGBS.mean_Fractional_methylation.merge.hg38.bw /share/home/zhzhang24/PG/epigen/dnamethylation/Homo_sapiens.3geneTSSud2kb_200bpbins.4nodup.bed /share/home/zhzhang24/PG/epigen/dnamethylation/WGBS/mergebw/Homo_sapiens.3geneTSSud2kb_200bpbins.WGBSmethylation.tab

#结果处理简化（仅保留第一列：每个bin的信息属于哪个基因以及bin相对位置；和第六列，bin中被检测到甲基化分数的碱基(cpg岛碱基)的平均甲基化分数。不用全部碱基的均值是因为去除cpg岛分布的影响，此处想表达bin甲基化分数低不是因为甲基化岛少，而是对甲基化岛的甲基化水平调控导致的。因为DNA甲基化的对象是cpg岛）
cat /share/home/zhzhang24/PG/epigen/dnamethylation/WGBS/mergebw/Homo_sapiens.3geneTSSud2kb_200bpbins.WGBSmethylation.tab |awk '{print $1"\t"$6}'|sed 's/__/\t/g' > /share/home/zhzhang24/PG/epigen/dnamethylation/WGBS/mergebw/Homo_sapiens.3geneTSSud2kb_200bpbins.WGBSmethylation.sim.tab


#传至实验室服务器
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/epigen/dnamethylation/WGBS/mergebw/Homo_sapiens.3geneTSSud2kb_200bpbins.WGBSmethylation.sim.tab zhzhang@122.205.95.67:/home/zhzhang/PG/epigen/dnamethy/



```
```r
#(a每个基因TSS+—2kb区域每个bin的平均甲基化分数file，b输出的图，
#c输出每类基因TSS+—2kb区域每类bin的甲基化分数均值，d输出两类lnc每个bin甲基化分数wilcox显著性检验p值)
#WGBS数据
a <- "~/PG/epigen/dnamethy/Homo_sapiens.3geneTSSud2kb_200bpbins.WGBSmethylation.sim.tab"
b <- "/home/zhzhang/PG/epigen/dnamethy/Homo_sapiens.3geneTSSud2kb_200bpbins.WGBSmethylation.pdf"
c <- "/home/zhzhang/PG/epigen/dnamethy/Homo_sapiens.3geneTSSud2kb_200bpbins.WGBSmethylation.tj.txt"
d <- "/home/zhzhang/PG/epigen/dnamethy/Homo_sapiens.2lncTSSud2kb_200bpbins.WGBSmethylation.wilcoxpvalue.txt"
#导入每个基因TSS+—2kb区域每个bin的平均甲基化分数
methylation <- read.delim(a, header=FALSE)
colnames(methylation) <- c("geneid","binpos","methy")
#导入基因分类文件
geneid_class <- read.delim("~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt")
#合并信息
he <- left_join(methylation,geneid_class,by="geneid")
he$type[is.na(he$type)==T] <- "Random intergenic"
he$type <- factor(he$type,levels = c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                                     "Non-pseudogene-associated lncRNA","Random intergenic"))
he$binpos <- factor(he$binpos,levels=seq(100,3900,by=200))
#去除没有被WGBS/RRBS测序检测到甲基化分数的bin（覆盖碱基数为0的bin,不含有cpg岛永远不会被甲基化，去除是为了防止甲基化岛的分布对甲基化分数的影响）
he_rm <- filter(he,methy!=0)
#统计每类基因TSS+—2kb区域每类bin的甲基化分数均值
meavg <- group_by(he_rm,type,binpos)%>%
  summarise(mean=mean(methy),median=median(methy))
data.table::fwrite(meavg,file =c,sep = '\t',row.names = F,quote = F,col.names = T)
#plot(bin的中点坐标作为x轴)
plot <- ggplot(data=he_rm,aes(x=binpos,y=methy))+
  geom_boxplot(aes(fill=type),position="dodge",
               color="black",size=0.1,outlier.alpha = 0,fatten=10,alpha=0.8)+
  geom_line(data=meavg,aes(color=type,y=mean,group=type),linetype="dashed",size=0.6)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                             "Non-pseudogene-associated lncRNA","Random intergenic"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-associated sense lncRNA","Pseudogene-associated antisense lncRNA",
                              "Non-pseudogene-associated lncRNA","Random intergenic"))+
  theme_half_open()+
  scale_x_discrete(breaks=c(100,700,1300,1900,2100,2700,3300,3900),
                   labels = c("[-2k,-1.8k]","[-1.4k,-1.2k]","[-0.8k,-0.6k]",
                              "[-0.2k,0]","[0,0.2k]",
                              "[0.6k,0.8k]","[1.2k,1.4k]","[1.8k,2k]"))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),
                     labels = c("0","25","50","75","100"))+
  labs(y ="Fractional methylation (%)",x ="Distance from TSS",fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 19),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),legend.position = "none")+
  theme(axis.text.x = element_text(angle =30)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave(b,plot,width = 10, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#显著性检验(两类lnc基因的TSS+—2kb区域每个bin的甲基化分数之间检验中位数差异是否显著)
wildata <- data.frame(binpos=seq(100,3900,by=200))
for (i in c(1:20)) {
  wildata[i,2] <- wilcox.test(filter(he_rm,type=="Pseudogene-associated sense lncRNA" & binpos==wildata[i,1])$methy,
                              filter(he_rm,type=="Non-pseudogene-associated lncRNA" & binpos==wildata[i,1])$methy)[["p.value"]]
  wildata[i,3] <- wilcox.test(filter(he_rm,type=="Pseudogene-associated antisense lncRNA" & binpos==wildata[i,1])$methy,
                              filter(he_rm,type=="Non-pseudogene-associated lncRNA" & binpos==wildata[i,1])$methy)[["p.value"]]
  wildata[i,4] <- wilcox.test(filter(he_rm,type=="Pseudogene-associated antisense lncRNA" & binpos==wildata[i,1])$methy,
                              filter(he_rm,type=="Pseudogene-associated sense lncRNA" & binpos==wildata[i,1])$methy)[["p.value"]]
}
colnames(wildata)[2] <- "PASvsNPApvalue"
colnames(wildata)[3] <- "PAAvsNPApvalue"
colnames(wildata)[4] <- "PASvsPAApvalue"
data.table::fwrite(wildata,file =d,sep = '\t',row.names = F,quote = F,col.names = T)


```
##### ATAC Cpeak
```r
#cpeak(biorivix s6table)/ATACdb v1.03数据库下载
#https://www.licpathway.net/ATACdb/index.php
#1493 samples染色质可及区域bed下载（人类）hg38
cd /home/zhzhang/PG/epigen/atac/
wget https://www.licpathway.net/ATACdb/download/packages/Accessible_chromatin_region_all.bed --no-check-certificate

```
```r
#cpeak数据（hg38）处理，提取出实验检测到的，分为持家的和特异性的
grep -w "observed" "/home/zhzhang/PG/epigen/atac/TableS6.cPeaks.gff3"|grep 'housekeeping=TRUE'|awk '{print $1"\t"$4"\t"$5}'|sed 's/^chr//g' > /home/zhzhang/PG/epigen/atac/cPeaks.observed.housekeeping.bed
grep -w "observed" "/home/zhzhang/PG/epigen/atac/TableS6.cPeaks.gff3"|grep 'housekeeping=FALSE'|awk '{print $1"\t"$4"\t"$5}'|sed 's/^chr//g' > /home/zhzhang/PG/epigen/atac/cPeaks.observed.specify.bed


#转换为bedgraph文件
bedtools genomecov -i /home/zhzhang/PG/epigen/atac/cPeaks.observed.housekeeping.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /home/zhzhang/PG/epigen/atac/cPeaks.observed.housekeeping.bedgraph
bedtools genomecov -i /home/zhzhang/PG/epigen/atac/cPeaks.observed.specify.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /home/zhzhang/PG/epigen/atac/cPeaks.observed.specify.bedgraph
#转换为bw文件
bedGraphToBigWig /home/zhzhang/PG/epigen/atac/cPeaks.observed.housekeeping.bedgraph /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt /home/zhzhang/PG/epigen/atac/cPeaks.observed.housekeeping.bw
bedGraphToBigWig /home/zhzhang/PG/epigen/atac/cPeaks.observed.specify.bedgraph /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt /home/zhzhang/PG/epigen/atac/cPeaks.observed.specify.bw
#deeptools computeMatrix 
computeMatrix reference-point --referencePoint TSS -R /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed" /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed -S /share/home/zhzhang24/PG/epigen/atac/cPeaks.observed.housekeeping.bw -b 3000 -a 3000 --binSize 10 -o /share/home/zhzhang24/PG/epigen/atac/Homo_sapiens.atachousekeeping_3geneTSS.matrix.gz -p 28
#deeptools 可视化
plotProfile -m /share/home/zhzhang24/PG/epigen/atac/Homo_sapiens.atachousekeeping_3geneTSS.matrix.gz -o /share/home/zhzhang24/PG/epigen/atac/Homo_sapiens.atachousekeeping_3geneTSS.png --dpi 1200 --legendLocation "upper-right" --colors "#B04150" "#00A087" "#628255" "#4DBBD5" "#8491B4" --samplesLabel "" --refPointLabel TSS -y "ATAC cPeak count frequency" --plotHeight 7 --plotWidth 8 --regionsLabel "Protein-coding" "PAS lncRNA" "PAA lncRNA" "NPA lncRNA" "Random intergenic" --outFileNameData /share/home/zhzhang24/PG/epigen/atac/Homo_sapiens.atachousekeeping_3geneTSS.profile.txt
#Homo_sapiens._3geneTSS.profile.txt删掉第一行导入R

```


##### Dnase（可及性）
```r
#DNase (53)
##################################
grep "DNase.macs2.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
cat /home/zhzhang/PG/epigen/all_narrowPeakhistone_hg19/${file} |awk '{print $1"\t"$2"\t"$3}' > /home/zhzhang/PG/epigen/dnase/hg19/${file}.bed
liftOver /home/zhzhang/PG/epigen/dnase/hg19/${file}.bed /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/epigen/dnase/hg38/${file}.bed /home/zhzhang/PG/epigen/dnase/hg19/unmap.txt
sed -i 's/^chr//g' /home/zhzhang/PG/epigen/dnase/hg38/${file}.bed
echo "${file}" >> /home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.log
done
##################################
grep "DNase.macs2.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
bedtools sort -i /home/zhzhang/PG/epigen/dnase/hg38/${file}.bed|bedtools merge >> /home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.bed
echo "${file}" >> /home/zhzhang/PG/epigen/dnase/Homo_sapiens.all.dnase.log
done
#去除转换后非染色体peak
bedtools sort -i /home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.bed |awk '$1!="M" {print $0}'|grep -v "alt"|grep -v "random"|grep -v "Un" > /home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.merge.bed
#转换为bedgraph文件,并将第四列的peak数量转换为所有样本中的平均值
bedtools genomecov -i /home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.merge.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.merge.bedgraph
cat "/home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.merge.bedgraph"|awk '{print $1"\t"$2"\t"$3"\t"$4/53}' > /home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.meanpeakcount.bedgraph
#转换为bw文件
bedGraphToBigWig /home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.meanpeakcount.bedgraph /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt /home/zhzhang/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.meanpeakcount.bw
#deeptools computeMatrix 
computeMatrix reference-point --referencePoint TSS -R /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed" /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed -S /share/home/zhzhang24/PG/epigen/dnase/Homo_sapiens.dnase.allpeak.meanpeakcount.bw -b 3000 -a 3000 --binSize 10 -o /share/home/zhzhang24/PG/epigen/dnase/Homo_sapiens.dnase_3geneTSS.matrix.gz -p 28
#deeptools 可视化
plotProfile -m /share/home/zhzhang24/PG/epigen/dnase/Homo_sapiens.dnase_3geneTSS.matrix.gz -o /share/home/zhzhang24/PG/epigen/dnase/Homo_sapiens.dnase_3geneTSS.png --dpi 1200 --legendLocation "upper-right" --colors "#B04150" "#00A087" "#628255" "#4DBBD5" "#8491B4" --samplesLabel "" --refPointLabel TSS -y "Dnase peak count frequency" --plotHeight 7 --plotWidth 8 --regionsLabel "Protein-coding" "PAS lncRNA" "PAA lncRNA" "NPA lncRNA" "Random intergenic" --outFileNameData /share/home/zhzhang24/PG/epigen/dnase/Homo_sapiens.dnase_3geneTSS.profile.txt
#Homo_sapiens.dnase_3geneTSS.profile.txt删掉第一行导入R


```


##### H3K27ac（转录激活）
```r
#H3K27AC (98)
##################################汇总全部98 epigenomes h3k27ac peak
grep "H3K27ac.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
cat /home/zhzhang/PG/epigen/all_narrowPeakhistone_hg19/${file} |awk '{print $1"\t"$2"\t"$3}' > /home/zhzhang/PG/epigen/h3k27ac/hg19/${file}.bed
liftOver /home/zhzhang/PG/epigen/h3k27ac/hg19/${file}.bed /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/epigen/h3k27ac/hg38/${file}.bed /home/zhzhang/PG/epigen/h3k27ac/hg19/unmap.txt
sed -i 's/^chr//g' /home/zhzhang/PG/epigen/h3k27ac/hg38/${file}.bed
echo "${file}" >> /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.log
done
##################################
grep "H3K27ac.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
bedtools sort -i /home/zhzhang/PG/epigen/h3k27ac/hg38/${file}.bed|bedtools merge >> /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.allpeak.bed
echo "${file}" >> /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.all.h3k27ac.log
done
#去除转换后非染色体peak
bedtools sort -i /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.allpeak.bed|awk '$1!="M" {print $0}'|grep -v "alt"|grep -v "random"|grep -v "Un" > /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.allpeak.merge.bed
#转换为bedgraph文件,并将第四列的peak数量转换为所有样本中的平均值
bedtools genomecov -i /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.allpeak.merge.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.allpeak.merge.bedgraph
cat "/home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.allpeak.merge.bedgraph"|awk '{print $1"\t"$2"\t"$3"\t"$4/98}' > /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.allpeak.meanpeakcount.bedgraph
#转换为bw文件
bedGraphToBigWig /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.allpeak.meanpeakcount.bedgraph /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt /home/zhzhang/PG/epigen/h3k27ac/Homo_sapiens.h3k27ac.allpeak.meanpeakcount.bw

epi="h3k27ac"
#deeptools computeMatrix 
computeMatrix reference-point --referencePoint TSS -R /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed" /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed -S /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}.allpeak.meanpeakcount.bw -b 3000 -a 3000 --binSize 10 -o /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.matrix.gz -p 28
#deeptools 可视化
plotProfile -m /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.matrix.gz -o /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.png --dpi 1200 --legendLocation "upper-right" --colors "#B04150" "#00A087" "#628255" "#4DBBD5" "#8491B4" --samplesLabel "" --refPointLabel TSS -y "${epi} peak count frequency" --plotHeight 7 --plotWidth 8 --regionsLabel "Protein-coding" "PAS lncRNA" "PAA lncRNA" "NPA lncRNA" "Random intergenic" --outFileNameData /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.profile.txt
#Homo_sapiens.h3k27ac_3geneTSS.profile.txt删掉第一行导入R

```


##### H3K9ac（转录激活）
```r
#H3K9ac (62)
##################################
grep "H3K9ac.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
cat /home/zhzhang/PG/epigen/all_narrowPeakhistone_hg19/${file} |awk '{print $1"\t"$2"\t"$3}' > /home/zhzhang/PG/epigen/h3k9ac/hg19/${file}.bed
liftOver /home/zhzhang/PG/epigen/h3k9ac/hg19/${file}.bed /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/epigen/h3k9ac/hg38/${file}.bed /home/zhzhang/PG/epigen/h3k9ac/hg19/unmap.txt
sed -i 's/^chr//g' /home/zhzhang/PG/epigen/h3k9ac/hg38/${file}.bed
echo "${file}" >> /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.log
done
##################################
grep "H3K9ac.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
bedtools sort -i /home/zhzhang/PG/epigen/h3k9ac/hg38/${file}.bed|bedtools merge >> /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.allpeak.bed
echo "${file}" >> /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.all.h3k9ac.log
done
#去除转换后非染色体peak
bedtools sort -i /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.allpeak.bed|awk '$1!="M" {print $0}'|grep -v "alt"|grep -v "random"|grep -v "Un" > /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.allpeak.merge.bed
#转换为bedgraph文件,并将第四列的peak数量转换为所有样本中的平均值
bedtools genomecov -i /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.allpeak.merge.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.allpeak.merge.bedgraph
cat "/home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.allpeak.merge.bedgraph"|awk '{print $1"\t"$2"\t"$3"\t"$4/62}' > /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.allpeak.meanpeakcount.bedgraph
#转换为bw文件
bedGraphToBigWig /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.allpeak.meanpeakcount.bedgraph /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt /home/zhzhang/PG/epigen/h3k9ac/Homo_sapiens.h3k9ac.allpeak.meanpeakcount.bw

epi="h3k9ac"
#deeptools computeMatrix 
computeMatrix reference-point --referencePoint TSS -R /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed" /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed -S /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}.allpeak.meanpeakcount.bw -b 3000 -a 3000 --binSize 10 -o /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.matrix.gz -p 28
#deeptools 可视化
plotProfile -m /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.matrix.gz -o /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.png --dpi 1200 --legendLocation "upper-right" --colors "#B04150" "#00A087" "#628255" "#4DBBD5" "#8491B4" --samplesLabel "" --refPointLabel TSS -y "${epi} peak count frequency" --plotHeight 7 --plotWidth 8 --regionsLabel "Protein-coding" "PAS lncRNA" "PAA lncRNA" "NPA lncRNA" "Random intergenic" --outFileNameData /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.profile.txt
#Homo_sapiens.h3k9ac_3geneTSS.profile.txt删掉第一行导入R



```


##### H3K4me3 (通过转录暂停-释放过程促进基因表达)
```r
#h3k4me3 (127)
##################################
grep "H3K4me3.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
cat /home/zhzhang/PG/epigen/all_narrowPeakhistone_hg19/${file} |awk '{print $1"\t"$2"\t"$3}' > /home/zhzhang/PG/epigen/h3k4me3/hg19/${file}.bed
liftOver /home/zhzhang/PG/epigen/h3k4me3/hg19/${file}.bed /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/epigen/h3k4me3/hg38/${file}.bed /home/zhzhang/PG/epigen/h3k4me3/hg19/unmap.txt
sed -i 's/^chr//g' /home/zhzhang/PG/epigen/h3k4me3/hg38/${file}.bed
echo "${file}" >> /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.log
done
##################################
grep "H3K4me3.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
bedtools sort -i /home/zhzhang/PG/epigen/h3k4me3/hg38/${file}.bed|bedtools merge >> /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.allpeak.bed
echo "${file}" >> /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.all.h3k4me3.log
done
#去除转换后非染色体peak
bedtools sort -i /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.allpeak.bed|awk '$1!="M" {print $0}'|grep -v "alt"|grep -v "random"|grep -v "Un" > /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.allpeak.merge.bed
#转换为bedgraph文件,并将第四列的peak数量转换为所有样本中的平均值
bedtools genomecov -i /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.allpeak.merge.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.allpeak.merge.bedgraph
cat "/home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.allpeak.merge.bedgraph"|awk '{print $1"\t"$2"\t"$3"\t"$4/127}' > /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.allpeak.meanpeakcount.bedgraph
#转换为bw文件
bedGraphToBigWig /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.allpeak.meanpeakcount.bedgraph /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt /home/zhzhang/PG/epigen/h3k4me3/Homo_sapiens.h3k4me3.allpeak.meanpeakcount.bw

epi="h3k4me3"
#deeptools computeMatrix 
computeMatrix reference-point --referencePoint TSS -R /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed" /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed -S /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}.allpeak.meanpeakcount.bw -b 3000 -a 3000 --binSize 10 -o /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.matrix.gz -p 28
#deeptools 可视化
plotProfile -m /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.matrix.gz -o /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.png --dpi 1200 --legendLocation "upper-right" --colors "#B04150" "#00A087" "#628255" "#4DBBD5" "#8491B4" --samplesLabel "" --refPointLabel TSS -y "${epi} peak count frequency" --plotHeight 7 --plotWidth 8 --regionsLabel "Protein-coding" "PAS lncRNA" "PAA lncRNA" "NPA lncRNA" "Random intergenic" --outFileNameData /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.profile.txt
#Homo_sapiens.h3k4me3_3geneTSS.profile.txt删掉第一行导入R


```


##### H3K4me2（转录激活）
```r
#h3k4me2 (24)
##################################
grep "H3K4me2.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
cat /home/zhzhang/PG/epigen/all_narrowPeakhistone_hg19/${file} |awk '{print $1"\t"$2"\t"$3}' > /home/zhzhang/PG/epigen/h3k4me2/hg19/${file}.bed
liftOver /home/zhzhang/PG/epigen/h3k4me2/hg19/${file}.bed /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/epigen/h3k4me2/hg38/${file}.bed /home/zhzhang/PG/epigen/h3k4me2/hg19/unmap.txt
sed -i 's/^chr//g' /home/zhzhang/PG/epigen/h3k4me2/hg38/${file}.bed
echo "${file}" >> /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.log
done
##################################
grep "H3K4me2.narrowPeak.gz" /home/zhzhang/PG/epigen/narrowPeak.filename.txt|while read i
do
file=${i%.gz}
bedtools sort -i /home/zhzhang/PG/epigen/h3k4me2/hg38/${file}.bed|bedtools merge >> /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.allpeak.bed
echo "${file}" >> /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.all.h3k4me2.log
done
#去除转换后非染色体peak
bedtools sort -i /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.allpeak.bed|awk '$1!="M" {print $0}'|grep -v "alt"|grep -v "random"|grep -v "Un" > /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.allpeak.merge.bed
#转换为bedgraph文件,并将第四列的peak数量转换为所有样本中的平均值
bedtools genomecov -i /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.allpeak.merge.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.allpeak.merge.bedgraph
cat "/home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.allpeak.merge.bedgraph"|awk '{print $1"\t"$2"\t"$3"\t"$4/24}' > /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.allpeak.meanpeakcount.bedgraph
#转换为bw文件
bedGraphToBigWig /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.allpeak.meanpeakcount.bedgraph /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt /home/zhzhang/PG/epigen/h3k4me2/Homo_sapiens.h3k4me2.allpeak.meanpeakcount.bw

epi="h3k4me2"
#deeptools computeMatrix 
computeMatrix reference-point --referencePoint TSS -R /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed" /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed -S /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}.allpeak.meanpeakcount.bw -b 3000 -a 3000 --binSize 10 -o /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.matrix.gz -p 28
#deeptools 可视化
plotProfile -m /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.matrix.gz -o /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.png --dpi 1200 --legendLocation "upper-right" --colors "#B04150" "#00A087" "#628255" "#4DBBD5" "#8491B4" --samplesLabel "" --refPointLabel TSS -y "${epi} peak count frequency" --plotHeight 7 --plotWidth 8 --regionsLabel "Protein-coding" "PAS lncRNA" "PAA lncRNA" "NPA lncRNA" "Random intergenic" --outFileNameData /share/home/zhzhang24/PG/epigen/${epi}/Homo_sapiens.${epi}_3geneTSS.profile.txt
#Homo_sapiens.h3k4me2_3geneTSS.profile.txt删掉第一行导入R



```


##### profile图
```r
#a输入profile数据矩阵路径，b输出profile图路径,c输入图y轴名称
#ATAC
a <- "~/PG/epigen/Homo_sapiens.atachousekeeping_3geneTSS.profile.txt"
b <- "/home/zhzhang/PG/epigen/Homo_sapiens.atachousekeeping_3geneTSS.pdf"
c <- "ATAC-seq housekeeping\ncPeak count frequency"
#DNASE
a <- "~/PG/epigen/Homo_sapiens.dnase_3geneTSS.profile.txt"
b <- "/home/zhzhang/PG/epigen/Homo_sapiens.dnase_3geneTSS.pdf"
c <- "DHS\ncount frequency"
#H3K27AC
a <- "~/PG/epigen/Homo_sapiens.h3k27ac_3geneTSS.profile.txt"
b <- "/home/zhzhang/PG/epigen/Homo_sapiens.h3k27ac_3geneTSS.pdf"
c <- "H3K27ac\npeak count frequency"
#H3K9ac
a <- "~/PG/epigen/Homo_sapiens.h3k9ac_3geneTSS.profile.txt"
b <- "/home/zhzhang/PG/epigen/Homo_sapiens.h3k9ac_3geneTSS.pdf"
c <- "H3K9ac\npeak count frequency"
#H3K4me3 
a <- "~/PG/epigen/Homo_sapiens.h3k4me3_3geneTSS.profile.txt"
b <- "/home/zhzhang/PG/epigen/Homo_sapiens.h3k4me3_3geneTSS.pdf"
c <- "H3K4me3\npeak count frequency"
#H3K4me2
a <- "~/PG/epigen/Homo_sapiens.h3k4me2_3geneTSS.profile.txt"
b <- "/home/zhzhang/PG/epigen/Homo_sapiens.h3k4me2_3geneTSS.pdf"
c <- "H3K4me2\npeak count frequency"
#导入profile数据矩阵（一行是一类样本，一列是一个bin，将每个bin列名改为其在6kb区间的起始位点）
profile <- read.delim(a)%>%
  select(-1)%>%
  column_to_rownames("X")
colnames(profile) <- seq(0,length.out=600,by=10)
#矩阵转换为长数据
profile <- t(profile)%>%
  data.frame()
pc <- mutate(select(profile,1),type="Protein-coding")%>%
  rownames_to_column("bin")
colnames(pc)[2] <- "fre"
pglnc <- mutate(select(profile,2),type="PAS lncRNA")%>%
  rownames_to_column("bin")
colnames(pglnc)[2] <- "fre"
pAAlnc <- mutate(select(profile,3),type="PAA lncRNA")%>%
  rownames_to_column("bin")
colnames(pAAlnc)[2] <- "fre"
npglnc <- mutate(select(profile,4),type="NPA lncRNA")%>%
  rownames_to_column("bin")
colnames(npglnc)[2] <- "fre"
inter <- mutate(select(profile,5),type="Random intergenic")%>%
  rownames_to_column("bin")
colnames(inter)[2] <- "fre"
profile_plot <- rbind(pc,pglnc,pAAlnc,npglnc,inter)
profile_plot$type <- factor(profile_plot$type,levels = c("Random intergenic","NPA lncRNA","PAA lncRNA",
                                                         "PAS lncRNA","Protein-coding"))
profile_plot$bin <- as.numeric(profile_plot$bin)
#plot(bin的起始坐标作为x轴)
profile_plot=mutate(profile_plot,bin=bin-1000)%>%
  filter(bin<=4000&bin>=0)
phs <- ggplot(data = profile_plot,aes(x=bin,y=fre))+
  geom_vline(xintercept=2000,size=1,linetype=2,color="gray")+
  geom_line(aes(color=type,group=type),size=1)+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA","Random intergenic"))+
  theme_half_open()+
  scale_x_continuous(expand=c(0,0),limits=c(0,4000),
                     breaks=c(0,1000,2000,3000,4000),
                     labels = c("-2","-1","0","1","2"))+
  labs(y =c,x ="Distance from TSS (kb)",fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 17),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),legend.position = "none",
        plot.margin = margin(t=15,r=15,b=0,l=5, unit = "pt"))
ggsave(b,phs,width = 3.5, height = 3.5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#ATAC/DNASE
ggsave(b,phs,width = 4.5, height = 3,dpi=1200, units = "in", device='pdf',bg = "transparent")


```




##### 2.基因TSS附近染色质修饰因子结合分布
##### TFBS下载
```r
#下载每个转录因子非冗余结合位点数据bigbed文件（GTRD数据库）(HG38),转换成bed,最后合并为一个bed（需要下载的文件名汇总于/home/zhzhang/PG/TFBS/物种/bigbed/物种TFBS.info.txt）[chr start end TFname]（所有TF所有潜在的结合位点）
#人
################"/home/zhzhang/PG/TFBS/Homo_sapiens/bigbed/downhsTFBS.sh"内容
#!/bin/bash
cd /home/zhzhang/PG/TFBS/Homo_sapiens/bigbed/
cat /home/zhzhang/PG/TFBS/Homo_sapiens/bigbed/hsTFBS.info.txt|while read i
do
file=${i%.bb}
wget http://gtrd.biouml.org:8888/egrid/bigBeds/hg38/ChIP-seq/Meta-clusters_by_TF/${file}.bb
bigBedToBed /home/zhzhang/PG/TFBS/Homo_sapiens/bigbed/${file}.bb /home/zhzhang/PG/TFBS/Homo_sapiens/bed/${file}.bed
TFname=${file%_ChIP-seq_Meta-clusters}
cat /home/zhzhang/PG/TFBS/Homo_sapiens/bed/${file}.bed|awk -v TFname="${TFname}" '{print $1"\t"$2"\t"$3"\t"TFname"_"$4}' >> /home/zhzhang/PG/TFBS/Homo_sapiens/allbed/Homo_sapiens.allTFBS.bed
echo "${TFname}" >> /home/zhzhang/PG/TFBS/Homo_sapiens/allbed/Homo_sapiens.allTFBS.log
done
###############################################################
nohup "/home/zhzhang/PG/TFBS/Homo_sapiens/bigbed/downhsTFBS.sh" &


```


##### 对比基因TSS附近染色质修饰因子结合位点count frequency
```r
cat "/home/zhzhang/PG/TFBS/Homo_sapiens/epiaTF/h3k27ac/TFname.txt"|while read i
cat "/home/zhzhang/PG/TFBS/Homo_sapiens/epiaTF/h3k4me3/TFname.txt"|while read i
cat "/share/home/zhzhang24/PG/TFBS/Homo_sapiens/epiaTF/NFKB/TFname.txt"|while read i
do
arr=($i)
qian=${arr[0]}
hou=${arr[1]}
#TFBS bed去掉chr
cat /share/home/zhzhang24/PG/TFBS/Homo_sapiens/bed/${qian}_ChIP-seq_Meta-clusters.bed|sed 's/^chr//g' > /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.bed
#TFBS bed转换为bedgraph文件
micromamba run -n SEQ bedtools genomecov -i /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.bed -g /share/home/zhzhang24/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.bedgraph
#转换为bw文件
micromamba run -n SEQ bedGraphToBigWig /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.bedgraph /share/home/zhzhang24/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.bw
done



module load arm/deeptools/3.5.5
cat "/share/home/zhzhang24/PG/TFBS/Homo_sapiens/epiaTF/h3k27ac/TFname.txt"|while read i
cat "/share/home/zhzhang24/PG/TFBS/Homo_sapiens/epiaTF/h3k4me3/TFname.txt"|while read i
do
arr=($i)
qian=${arr[0]}
hou=${arr[1]}
#deeptools computeMatrix 
computeMatrix reference-point --referencePoint TSS -R /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paslncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.paalncgeneTSS.bed" "/share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.npalncgeneTSS.bed" /share/home/zhzhang24/PG/epigen/annotationbed/Homo_sapiens.intergenicTSS.bed -S /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.bw -b 3000 -a 3000 --binSize 10 -o /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.matrix.gz -p 28
#deeptools 可视化
plotProfile -m /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.matrix.gz -o /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.png --dpi 1200 --legendLocation "upper-right" --colors "#B04150" "#00A087" "#628255" "#4DBBD5" "#8491B4" --samplesLabel "" --refPointLabel TSS -y "TFBS count frequency" --plotHeight 7 --plotWidth 8 --regionsLabel "Protein-coding" "PAS lncRNA" "PAA lncRNA" "NPA lncRNA" "Random intergenic" --outFileNameData /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.3geneTSS.profile.txt
sed "1d" -i /share/home/zhzhang24/PG/TFBS/Homo_sapiens/nochrbed/${hou}.3geneTSS.profile.txt
done


```
```r
#profile图
#
he <- data.frame()
TF <- c("RBBP5","SETD1A","WDR5","CBP","P300","ASH2L","CXXC1","HCFC1","GCN5","PCAF")
for (i in 1:length(TF)) {
  TFname <- TF[i]
  #导入profile数据矩阵（一行是一类样本，一列是一个bin，将每个bin列名改为其在6kb区间的起始位点）
  profile <- read.delim(paste("/home/zhzhang/PG/TFBS/Homo_sapiens/nochrbed/",
                              TFname,".3geneTSS.profile.txt",
                              sep = ""))%>%
    select(-1)%>%
    column_to_rownames("X")
  colnames(profile) <- seq(0,length.out=600,by=10)
  #矩阵转换为长数据
  profile <- t(profile)%>%
    data.frame()
  pc <- mutate(select(profile,1),type="Protein-coding")%>%
    rownames_to_column("bin")
  colnames(pc)[2] <- "fre"
  pglnc <- mutate(select(profile,2),type="PAS lncRNA")%>%
    rownames_to_column("bin")
  colnames(pglnc)[2] <- "fre"
  pAAlnc <- mutate(select(profile,3),type="PAA lncRNA")%>%
    rownames_to_column("bin")
  colnames(pAAlnc)[2] <- "fre"
  npglnc <- mutate(select(profile,4),type="NPA lncRNA")%>%
    rownames_to_column("bin")
  colnames(npglnc)[2] <- "fre"
  inter <- mutate(select(profile,5),type="Random intergenic")%>%
    rownames_to_column("bin")
  colnames(inter)[2] <- "fre"
  profile_plot <- rbind(pc,pglnc,pAAlnc,npglnc,inter)
  profile_plot$type <- factor(profile_plot$type,levels = c("Random intergenic","NPA lncRNA","PAA lncRNA",
                                                           "PAS lncRNA","Protein-coding"))
  profile_plot$bin <- as.numeric(profile_plot$bin)
  profile_plot <- mutate(profile_plot,TF=TFname)
  he <- rbind(he,profile_plot)
}
he$TF <- factor(he$TF,levels = c("RBBP5","SETD1A","WDR5","CBP","P300","ASH2L","CXXC1","HCFC1","GCN5","PCAF"))
#plot(bin的起始坐标作为x轴)
he=mutate(he,bin=bin-2000)%>%
  filter(bin<=2000&bin>=0)
phs <- ggplot(data = he,aes(x=bin,y=fre))+
  geom_vline(xintercept=1000,size=1,linetype=2,color="gray")+
  geom_line(aes(color=type,group=type),size=0.5,alpha=0.8)+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                     limits=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA","Random intergenic"))+
  theme_half_open()+
  scale_x_continuous(expand=c(0,0),limits=c(0,2000),
                     breaks=c(0,500,1000,1500,2000),
                     labels = c("-1","-0.5","0","0.5","1"))+
  labs(y ="TFBS count frequency",x ="Distance from TSS (kb)",fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 17),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),legend.position = "none",
        plot.margin = margin(t=15,r=15,b=0,l=5, unit = "pt"))+
  facet_wrap(~TF,nrow=2,strip.position="top",as.table=F,scales="free_y")+
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13))
ggsave("/home/zhzhang/PG/TFBS/Homo_sapiens/nochrbed/10epiaTF.3geneTSS.profile.pdf",
       phs,width = 15, height = 7,dpi=1200, units = "in", device='pdf',bg = "transparent")



```


##### 3.近端(Proximal) Promoter bed产生
```r
#近端启动子(TSS上下游500bp，来源于FANTOM5标准：https://www.sciencedirect.com/science/article/pii/S0167779917300562?via%3Dihub#bib0050)
#人类
cat /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pcgeneTSS.bed|awk '{print $1"\t"$2-500"\t"$2+500"\t"$4"\tProximal-promoter\t"$6}' > /home/zhzhang/PG/epigen/promoterbed/Homo_sapiens/Homo_sapiens.all.Proximalpromoter.bed
cat /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pgdlncgeneTSS.bed|awk '{print $1"\t"$2-500"\t"$2+500"\t"$4"\tProximal-promoter\t"$6}' >> /home/zhzhang/PG/epigen/promoterbed/Homo_sapiens/Homo_sapiens.all.Proximalpromoter.bed
cat /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.npgdlncgeneTSS.bed|awk '{print $1"\t"$2-500"\t"$2+500"\t"$4"\tProximal-promoter\t"$6}' >> /home/zhzhang/PG/epigen/promoterbed/Homo_sapiens/Homo_sapiens.all.Proximalpromoter.bed
#小鼠
cp "/home/zhzhang/PG/Evolution/conserve/Mus_musculus.allgene.bed" /home/zhzhang/PG/epigen/promoterbed/Mus_musculus/
cat "/home/zhzhang/PG/epigen/promoterbed/Mus_musculus/Mus_musculus.allgene.bed"|awk '$6=="+" {print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6} $6=="-" {print $1"\t"$3"\t"$3+1"\t"$4"\t"$5"\t"$6}'|bedtools sort > /home/zhzhang/PG/epigen/promoterbed/Mus_musculus/Mus_musculus.allgeneTSS.bed
cat /home/zhzhang/PG/epigen/promoterbed/Mus_musculus/Mus_musculus.allgeneTSS.bed|awk '{print $1"\t"$2-500"\t"$2+500"\t"$4"\tProximal-promoter\t"$6}' > /home/zhzhang/PG/epigen/promoterbed/Mus_musculus/Mus_musculus.all.Proximalpromoter.bed




```
##### 4.统计三类基因两类启动子中激活调控chromHMM状态(1,2)的平均占有率(平均指在127 epis)
```r
ls /home/zhzhang/PG/epigen/state/allstate/ > /home/zhzhang/PG/epigen/state/allstate.log
sed -i 's/^chr//g' /home/zhzhang/PG/epigen/state/allstate/*.bed
#127 epis激活调控chromHMM状态(1,2)bed信息合并
##################################
cat "/home/zhzhang/PG/epigen/state/allstate.log"|while read i
do
file=${i%_15_coreMarks_hg38lift_mnemonics.bed}
grep -E "1_TssA|2_TssAFlnk" /home/zhzhang/PG/epigen/state/allstate/${file}_15_coreMarks_hg38lift_mnemonics.bed |awk -v a=${file} '{print $1"\t"$2"\t"$3"\t"$4"\t"a}' >> /home/zhzhang/PG/epigen/state/Homo_sapiens.state1a2.all.bed
echo "${file}" >> /home/zhzhang/PG/epigen/state/Homo_sapiens.state1a2.all.log
done
##################################
#127 epis抑制调控chromHMM状态(10,11)bed信息合并
##################################
cat "/home/zhzhang/PG/epigen/state/allstate.log"|while read i
do
file=${i%_15_coreMarks_hg38lift_mnemonics.bed}
grep -E "10_TssBiv|11_BivFlnk" /home/zhzhang/PG/epigen/state/allstate/${file}_15_coreMarks_hg38lift_mnemonics.bed |awk -v a=${file} '{print $1"\t"$2"\t"$3"\t"$4"\t"a}' >> /home/zhzhang/PG/epigen/state/Homo_sapiens.state10a11.all.bed
echo "${file}" >> /home/zhzhang/PG/epigen/state/Homo_sapiens.state10a11.all.log
done
##################################


#基因近端启动子与激活调控chromHMM状态(1,2)bed交集
bedtools intersect -a /home/zhzhang/PG/epigen/promoterbed/Homo_sapiens/Homo_sapiens.all.Proximalpromoter.bed -b /home/zhzhang/PG/epigen/state/Homo_sapiens.state1a2.all.bed -wo > /home/zhzhang/PG/epigen/state/Homo_sapiens.Proximalpromoter_intersect_state1a2.bed
#基因间区与激活调控chromHMM状态(1,2)bed交集
bedtools intersect -a /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed -b /home/zhzhang/PG/epigen/state/Homo_sapiens.state1a2.all.bed -wo > /home/zhzhang/PG/epigen/state/Homo_sapiens.intergenic_intersect_state1a2.bed
#基因远端启动子与激活调控chromHMM状态(1,2)bed交集
bedtools intersect -a /home/zhzhang/PG/epigen/promoterbed/Homo_sapiens/Homo_sapiens.all.Distalpromoter.bed -b /home/zhzhang/PG/epigen/state/Homo_sapiens.state1a2.all.bed -wo > /home/zhzhang/PG/epigen/state/Homo_sapiens.Distalpromoter_intersect_state1a2.bed
#基因近端启动子与抑制调控chromHMM状态(10,11)bed交集
bedtools intersect -a /home/zhzhang/PG/epigen/promoterbed/Homo_sapiens/Homo_sapiens.all.Proximalpromoter.bed -b /home/zhzhang/PG/epigen/state/Homo_sapiens.state10a11.all.bed -wo > /home/zhzhang/PG/epigen/state/Homo_sapiens.Proximalpromoter_intersect_state10a11.bed
#基因间区与抑制调控chromHMM状态(10,11)bed交集
bedtools intersect -a /home/zhzhang/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.bed -b /home/zhzhang/PG/epigen/state/Homo_sapiens.state10a11.all.bed -wo > /home/zhzhang/PG/epigen/state/Homo_sapiens.intergenic_intersect_state10a11.bed
#基因远端启动子与抑制调控chromHMM状态(10,11)bed交集
bedtools intersect -a /home/zhzhang/PG/epigen/promoterbed/Homo_sapiens/Homo_sapiens.all.Distalpromoter.bed -b /home/zhzhang/PG/epigen/state/Homo_sapiens.state10a11.all.bed -wo > /home/zhzhang/PG/epigen/state/Homo_sapiens.Distalpromoter_intersect_state10a11.bed


scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/epigen/state/Homo_sapiens.Distalpromoter_intersect_state1a2.bed /home/zhzhang/PG/epigen/state/
scp -P 22 zhzhang@211.69.141.147:/home/zhzhang/PG/epigen/state/Homo_sapiens.Distalpromoter_intersect_state1a2.bed /home/zhzhang/PG/epigen/state/
scp -P 22 zhzhang@211.69.141.147:"/home/zhzhang/PG/epigen/state/Homo_sapiens.intergenic_intersect_state10a11.bed" /home/zhzhang/PG/epigen/state/
scp -P 22 zhzhang@211.69.141.147:"/home/zhzhang/PG/epigen/state/Homo_sapiens.intergenic_intersect_state1a2.bed" /home/zhzhang/PG/epigen/state/
scp -P 22 zhzhang@211.69.141.147:"/home/zhzhang/PG/epigen/state/Homo_sapiens.Proximalpromoter_intersect_state10a11.bed" /home/zhzhang/PG/epigen/state/
scp -P 2022 zhzhang@211.69.141.147:"/home/zhzhang/PG/epigen/state/Homo_sapiens.Proximalpromoter_intersect_state1a2.bed" /home/zhzhang/PG/epigen/state/


```
```r
#函数根据基因两类启动子和基因间区与指定state交集，计算每个注释中state的平均占有率
#a输入近端启动子与state交集，b输入远端启动子与state交集，e输入基因间区与state交集
#c输入基因分类文件，d输入基因间区分类文件
getstatemeanratio <- function(a,b,e,c,d){
  #导入基因分类文件
  geneid_class <- read.delim(c)%>%
    filter(type!="Interference lncRNA")
  #导入基因间区分类
  intergenic <- read.delim(d, header=FALSE)%>%
    select(1)%>%
    mutate(type="Random intergenic")
  colnames(intergenic)[1] <- "geneid"
  #计算每个基因近端启动子中state的平均占有率
  pro_intersect_state <- data.table::fread(a, header=FALSE)%>%
    data.frame()%>%
    group_by(V4)%>%
    summarise(len=sum(V12))%>%
    mutate(meanratio=len/1000/127)%>%
    select(1,3)%>%
    mutate(type2="Proximal promoter")
  colnames(pro_intersect_state)[1] <- "geneid"
  pro_intersect_state <- left_join(geneid_class,pro_intersect_state,by="geneid")
  pro_intersect_state$meanratio[is.na(pro_intersect_state$meanratio)==T] <- 0
  pro_intersect_state$type2[is.na(pro_intersect_state$type2)==T] <- "Proximal promoter"
  #计算每个基因远端启动子中state的平均占有率
  dis_intersect_state <- data.table::fread(b, header=FALSE)%>%
    data.frame()%>%
    group_by(V4)%>%
    summarise(len=sum(V12))%>%
    mutate(meanratio=len/3000/127)%>%
    select(1,3)%>%
    mutate(type2="Distal promoter")
  colnames(dis_intersect_state)[1] <- "geneid"
  dis_intersect_state <- left_join(geneid_class,dis_intersect_state,by="geneid")
  dis_intersect_state$meanratio[is.na(dis_intersect_state$meanratio)==T] <- 0
  dis_intersect_state$type2[is.na(dis_intersect_state$type2)==T] <- "Distal promoter"
  #计算每个基因间区中state的平均占有率
  intergenic_intersect_state <- read.delim(e, header=FALSE)%>%
    group_by(V4)%>%
    summarise(len=sum(V12))%>%
    mutate(meanratio=len/3000/127)%>%
    select(1,3)%>%
    mutate(type2="Random intergenic")
  colnames(intergenic_intersect_state)[1] <- "geneid"
  intergenic_intersect_state <- left_join(intergenic,intergenic_intersect_state,by="geneid")
  intergenic_intersect_state$meanratio[is.na(intergenic_intersect_state$meanratio)==T] <- 0
  intergenic_intersect_state$type2[is.na(intergenic_intersect_state$type2)==T] <- "Random intergenic"
  #合并
  allhe <- rbind(pro_intersect_state,dis_intersect_state,intergenic_intersect_state)
  return(allhe)
}

#人state1a2
hs_state1a2_meanratio <- getstatemeanratio(a="/home/zhzhang/PG/epigen/state/Homo_sapiens.Proximalpromoter_intersect_state1a2.bed",
                                           b="/home/zhzhang/PG/epigen/state/Homo_sapiens.Distalpromoter_intersect_state1a2.bed",
                                           e="~/PG/epigen/state/Homo_sapiens.intergenic_intersect_state1a2.bed",
                                           c="~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                                           d="~/PG/Evolution/conserve/Homo_sapiens.20000random_3kb_intergenic.phastCons.txt")%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA",
                        type=="Random intergenic" ~ "Random intergenic"))
data.table::fwrite(hs_state1a2_meanratio,
                   file ="/home/zhzhang/PG/epigen/state/Homo_sapiens.3gene2promoter.state1a2.meanratio.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#统计被state1a2注释到的注释的比例（与state1a2有交集的比例）
hs_state1a2_meanratio_fortj1 <- mutate(hs_state1a2_meanratio,num=1,state=1)
hs_state1a2_meanratio_fortj1$state[hs_state1a2_meanratio_fortj1$meanratio==0] <- 0
tj1 <- group_by(hs_state1a2_meanratio_fortj1,type,type2)%>%
  summarise(allnum=sum(num),statenum=sum(state))%>%
  mutate(ratio=statenum/allnum)
data.table::fwrite(tj1,file ="/home/zhzhang/PG/epigen/state/SPhs.3gene2promoter.own_state1a2_ratio.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#检验两类lncRNA近端启动子被state1a2注释到的比例差异是否显著，p-value = 0.2059/p-value = 0.5044
tj1_forfisher <- mutate(tj1,nonum=allnum-statenum)
fisher.test(rbind(tj1_forfisher[2,c(4,6)],tj1_forfisher[6,c(4,6)]))
fisher.test(rbind(tj1_forfisher[4,c(4,6)],tj1_forfisher[6,c(4,6)]))
#统计被state1a2注释到的注释中，state1a2平均占有率的均值/中位数
tj2 <- group_by(filter(hs_state1a2_meanratio,meanratio!=0),type,type2)%>%
  summarise(median=median(meanratio),mean=mean(meanratio))
data.table::fwrite(tj2,file ="/home/zhzhang/PG/epigen/state/SPhs.3gene2promoter.state1a2.meanratio.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot近端启动子被state1a2注释到的比例(100个里有几个与state1a2重叠)
tj1$type <- factor(tj1$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA",
                                       "NPA lncRNA","Random intergenic"))
phspro_over <- ggplot(data=filter(tj1,type2=="Proximal promoter" | type2=="Random intergenic"),
                      aes(x=type,y=ratio*100))+
  geom_col(width=0.5,aes(fill=type))+
  geom_signif(annotations=c("N.S.","N.S."),y_position=c(60,50),tip_length = 0.01,
              xmin = 2:3,xmax = c(4,4))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                             "NPA lncRNA","Random intergenic"))+
  theme_half_open()+
  scale_y_continuous(limits = c(0,100),breaks = c(0,25,50,75,100))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA","Random\nintergenic"))+
  labs(y = "Proportion of proximal promoter\noverlapping active promoter state (%)",
       x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 16),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =40)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave("/home/zhzhang/PG/epigen/state/Homo_sapiens.3Proximalpromoter.Proportion_overlap_state1a2.pdf",
       phspro_over,width = 5, height = 5.8,dpi=1200, units = "in", device='pdf',bg = "transparent")
#plot近端启动子state1a2平均占有率（近端启动子中与state1a2重叠的区域长度/近端启动子长度）
hs_state1a2_meanratio$type <- factor(hs_state1a2_meanratio$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA",
                                                                           "NPA lncRNA","Random intergenic"))
phspro <- ggplot(data=filter(filter(hs_state1a2_meanratio,meanratio!=0),type2=="Proximal promoter" | type2=="Random intergenic"),
                 aes(x=type,y=meanratio))+
  geom_boxplot(fatten = 3,outlier.alpha = 0,width=0.5,notch=T,aes(fill=type))+
  geom_signif(map_signif_level=T,y_position=c(1,1.1,1.2,1,1.1),tip_length = 0.01,
              comparisons = list(c("Protein-coding","PAS lncRNA"),
                                 c("PAS lncRNA","PAA lncRNA"),
                                 c("PAS lncRNA","NPA lncRNA"),
                                 c("PAA lncRNA","NPA lncRNA"),
                                 c("NPA lncRNA","Random intergenic")
              ))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465","#8491B4"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA",
                             "NPA lncRNA","Random intergenic"))+
  theme_half_open()+
  scale_y_continuous(limits = c(0,1.3),breaks = c(0,0.25,0.5,0.75,1))+
  scale_x_discrete(labels = c("Protein-coding","PAS lncRNA","PAA lncRNA",
                              "NPA lncRNA","Random\nintergenic"))+
  labs(y = "Mean ratio of regions\nwith active promoter state",x =NULL,fill = NULL,color = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none") +
  theme(axis.text.x = element_text(angle =40)) +
  theme(axis.text.x = element_text(hjust=1,vjust = 1))
ggsave("/home/zhzhang/PG/epigen/state/Homo_sapiens.3Proximalpromoter.state1a2_meanratio.pdf",
       phspro,width = 5, height = 5.5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```


```r
#127 epis异染色质chromHMM状态(9)bed信息合并
##################################
cat "/home/zhzhang/PG/epigen/state/allstate.log"|while read i
do
file=${i%_15_coreMarks_hg38lift_mnemonics.bed}
grep -w "9_Het" /home/zhzhang/PG/epigen/state/allstate/${file}_15_coreMarks_hg38lift_mnemonics.bed |awk -v a=${file} '{print $1"\t"$2"\t"$3"\t"$4"\t"a}' >> /home/zhzhang/PG/epigen/state/Homo_sapiens.het.all.bed
echo "${file}" >> /home/zhzhang/PG/epigen/state/Homo_sapiens.het.all.log
done
##################################
#去除线粒体染色体
grep -v "M" /home/zhzhang/PG/epigen/state/Homo_sapiens.het.all.bed|bedtools sort > /home/zhzhang/PG/epigen/state/Homo_sapiens.het.all.sort.bed

#转换为bedgraph文件,并将第四列的peak数量转换为所有样本中的平均值
bedtools genomecov -i /home/zhzhang/PG/epigen/state/Homo_sapiens.het.all.sort.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /home/zhzhang/PG/epigen/state/Homo_sapiens.het.all.sort.bedgraph
cat /home/zhzhang/PG/epigen/state/Homo_sapiens.het.all.sort.bedgraph|awk '{print $1"\t"$2"\t"$3"\t"$4/127}' > /home/zhzhang/PG/epigen/state/Homo_sapiens.het.all.sort.meancount.bedgraph
#提取平均覆盖度>0.1的区域作为组成性constitutive异染色质区域，转换为bedgraph文件
cat /home/zhzhang/PG/epigen/state/Homo_sapiens.het.all.sort.meancount.bedgraph|awk '$4>=0.1{print $1"\t"$2"\t"$3}'|bedtools sort|bedtools merge > /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehet.bed
bedtools genomecov -i /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehet.bed -g /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt -bga > /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehet.bedgraph
#转换为bw文件
bedGraphToBigWig /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehet.bedgraph /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chrsize.txt /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehet.bw
#deeptools computeMatrix 
computeMatrix scale-regions -R /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pcgene.bed /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.pgdlncgene.bed /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.npgdlncgene.bed /home/zhzhang/PG/epigen/annotationbed/Homo_sapiens.intergenicbody.bed -S /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehet.bw -b 10000 -a 10000 --binSize 10 -o /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehethet_3geneud.matrix.gz -p 64
#deeptools 可视化
plotProfile -m /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehethet_3geneud.matrix.gz -o /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehethet_3geneud.png --dpi 1200 --legendLocation "upper-right" --colors "#B04150" "#00A087" "#4DBBD5" "#8491B4" --samplesLabel "" --startLabel TSS --endLabel TES -y "Constitutivehet heterochromatin region frequency" --plotHeight 7 --plotWidth 14 --regionsLabel "Protein-coding" "Pseudogene-derived lncRNA" "Non-pseudogene-derived lncRNA" "Random intergenic" --outFileNameData /home/zhzhang/PG/epigen/state/Homo_sapiens.constitutivehethet_3geneud.profile.txt
#Homo_sapiens._3geneTSS.profile.txt删掉第一行导入R




```



