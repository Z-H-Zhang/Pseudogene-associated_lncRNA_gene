### 二.lncRNA注释的补充
##### 1.subreads生成CCS序列
```r
#小鼠
cd /home/zhzhang/PG/ISOseqdata/Mus_musculus/subreadsbam/
for i in PMouse_*.subreads.bam
do
sample=${i%.subreads.bam}
ccs -j 64 --min-passes 1 --min-length 300 --min-rq 0.9 --report-file /home/zhzhang/PG/ISOseqdata/Mus_musculus/ccsbam/${sample}.txt /home/zhzhang/PG/ISOseqdata/Mus_musculus/subreadsbam/${sample}.subreads.bam /home/zhzhang/PG/ISOseqdata/Mus_musculus/ccsbam/${sample}.ccs.bam
done

#鸡
cd /home/zhzhang/PG/ISOseqdata/Gallus_gallus/subreadsbam/
for i in PC*.subreads.bam
do
sample=${i%.subreads.bam}
ccs -j 64 --min-passes 1 --min-length 300 --min-rq 0.9 --report-file /home/zhzhang/PG/ISOseqdata/Gallus_gallus/ccsbam/${sample}.txt /home/zhzhang/PG/ISOseqdata/Gallus_gallus/subreadsbam/${sample}.subreads.bam /home/zhzhang/PG/ISOseqdata/Gallus_gallus/ccsbam/${sample}.ccs.bam
done




#斑马鱼
cd /home/zhzhang/PG/ISOseqdata/Danio_rerio/subreadsbam/
for i in PZ*.subreads.bam
do
sample=${i%.subreads.bam}
/home/zhzhang/software/smrtlink_7/smrtcmds/bin/ccs --numThreads 64 --minPasses 1 --minLength 300 --noPolish --reportFile /home/zhzhang/PG/ISOseqdata/Danio_rerio/ccsbam/${sample}.txt /home/zhzhang/PG/ISOseqdata/Danio_rerio/subreadsbam/${sample}.subreads.bam /home/zhzhang/PG/ISOseqdata/Danio_rerio/ccsbam/${sample}.ccs.bam
done


```
##### 2.去除测序引物序列，生成full-lenght reads 
```r
#小鼠
cd /home/zhzhang/PG/ISOseqdata/Mus_musculus/ccsbam/
for i in PMouse_*.ccs.bam
do
sample=${i%.ccs.bam}
lima /home/zhzhang/PG/ISOseqdata/Mus_musculus/ccsbam/${sample}.ccs.bam /home/zhzhang/PG/ISOseqdata/IsoSeqPrimers.fasta /home/zhzhang/PG/ISOseqdata/Mus_musculus/flbam/${sample}.fl.bam -j 64 --isoseq
done

#鸡
cd /home/zhzhang/PG/ISOseqdata/Gallus_gallus/ccsbam/
for i in PC*.ccs.bam
do
sample=${i%.ccs.bam}
lima /home/zhzhang/PG/ISOseqdata/Gallus_gallus/ccsbam/${sample}.ccs.bam /home/zhzhang/PG/ISOseqdata/IsoSeqPrimers.fasta /home/zhzhang/PG/ISOseqdata/Gallus_gallus/flbam/${sample}.fl.bam -j 64 --isoseq
done


#斑马鱼
cd /home/zhzhang/PG/ISOseqdata/Danio_rerio/ccsbam/
for i in PZ*.ccs.bam
do
sample=${i%.ccs.bam}
/home/zhzhang/software/smrtlink_7/smrtcmds/bin/lima /home/zhzhang/PG/ISOseqdata/Danio_rerio/ccsbam/${sample}.ccs.bam /home/zhzhang/PG/ISOseqdata/IsoSeqPrimers.fasta /home/zhzhang/PG/ISOseqdata/Danio_rerio/flbam/${sample}.fl.bam -j 28 --isoseq --no-pbi
done



```
##### 3.去除嵌合体以及polyA tail，生成full-lenght non-chemeric 序列
```r
#小鼠
cd /home/zhzhang/PG/ISOseqdata/Mus_musculus/flbam/
for i in PMouse_*.fl.primer_5p--primer_3p.bam
do
sample=${i%.fl.primer_5p--primer_3p.bam}
isoseq3 refine -j 64 --require-polya /home/zhzhang/PG/ISOseqdata/Mus_musculus/flbam/${sample}.fl.primer_5p--primer_3p.bam /home/zhzhang/PG/ISOseqdata/IsoSeqPrimers.fasta /home/zhzhang/PG/ISOseqdata/Mus_musculus/flncbam/${sample}.flnc.bam
done

#鸡
cd /home/zhzhang/PG/ISOseqdata/Gallus_gallus/flbam/
for i in PC*.fl.primer_5p--primer_3p.bam
do
sample=${i%.fl.primer_5p--primer_3p.bam}
isoseq3 refine -j 64 --require-polya /home/zhzhang/PG/ISOseqdata/Gallus_gallus/flbam/${sample}.fl.primer_5p--primer_3p.bam /home/zhzhang/PG/ISOseqdata/IsoSeqPrimers.fasta /home/zhzhang/PG/ISOseqdata/Gallus_gallus/flncbam/${sample}.flnc.bam
done


#斑马鱼
cd /home/zhzhang/PG/ISOseqdata/Danio_rerio/flbam/
for i in PZ*.fl.primer_5p--primer_3p.bam
do
sample=${i%.fl.primer_5p--primer_3p.bam}
/home/zhzhang/software/smrtlink_7/smrtcmds/bin/isoseq3 refine --require-polya /home/zhzhang/PG/ISOseqdata/Danio_rerio/flbam/${sample}.fl.primer_5p--primer_3p.bam /home/zhzhang/PG/ISOseqdata/IsoSeqPrimers.fasta /home/zhzhang/PG/ISOseqdata/Danio_rerio/flncbam/${sample}.flnc.bam
done





```


##### 4.聚类生成full-lenght转录本序列
```r
#小鼠
ls /home/zhzhang/PG/ISOseqdata/Mus_musculus/flncbam/*.flnc.bam > /home/zhzhang/PG/ISOseqdata/Mus_musculus/flncbam/mmflnc.fofn
isoseq3 cluster -j 64 --verbose /home/zhzhang/PG/ISOseqdata/Mus_musculus/flncbam/mmflnc.fofn /home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.unpolished.transcripts.bam
#合并全部转录本序列
gzip -d "/home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.unpolished.transcripts.hq.fasta.gz"
gzip -d "/home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.unpolished.transcripts.lq.fasta.gz"
cp /home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.unpolished.transcripts.hq.fasta /home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.unpolished.transcripts.fasta
cat /home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.unpolished.transcripts.lq.fasta >> /home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.unpolished.transcripts.fasta
#筛选长度>200的序列
seqkit seq /home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.unpolished.transcripts.fasta -m 200 > /home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.lt200.transcripts.fasta

#鸡
ls /home/zhzhang/PG/ISOseqdata/Gallus_gallus/flncbam/*.flnc.bam > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/flncbam/ggflnc.fofn
isoseq3 cluster -j 64 --verbose /home/zhzhang/PG/ISOseqdata/Gallus_gallus/flncbam/ggflnc.fofn /home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.unpolished.transcripts.bam
#合并全部转录本序列
gzip -d "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.unpolished.transcripts.hq.fasta.gz"
gzip -d "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.unpolished.transcripts.lq.fasta.gz"
cp /home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.unpolished.transcripts.hq.fasta /home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.unpolished.transcripts.fasta
cat /home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.unpolished.transcripts.lq.fasta >> /home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.unpolished.transcripts.fasta
#筛选长度>=200的序列
seqkit seq /home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.unpolished.transcripts.fasta -m 200 > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.lt200.transcripts.fasta



#斑马鱼
ls /home/zhzhang/PG/ISOseqdata/Danio_rerio/flncbam/*.flnc.bam > /home/zhzhang/PG/ISOseqdata/Danio_rerio/flncbam/drflnc.fofn
/home/zhzhang/software/smrtlink_7/smrtcmds/bin/isoseq3 cluster -j 60 --verbose /home/zhzhang/PG/ISOseqdata/Danio_rerio/flncbam/drflnc.fofn /home/zhzhang/PG/ISOseqdata/Danio_rerio/fl_transcripts/drALL.unpolished.transcripts.bam
#筛选长度>=200的序列
seqkit seq "/home/zhzhang/PG/ISOseqdata/Danio_rerio/fl_transcripts/drALL.unpolished.transcripts.fasta" -m 200 > /home/zhzhang/PG/ISOseqdata/Danio_rerio/fl_transcripts/drALL.lt200.transcripts.fasta





```
##### 5.全长转录本序列利用RNA-seq数据纠错
```r
#小鼠
source /home/zhzhang/miniconda3/bin/activate SEQ
lordec-correct -2 /home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_brain_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_brain_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_brain_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_brain_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_cerebellum_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_cerebellum_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_cerebellum_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_cerebellum_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_gonad_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_gonad_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_gonad_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_gonad_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_gut_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_gut_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_gut_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_gut_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_heart_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_heart_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_heart_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Mus_musculus/qc_fq/Mouse_heart_male_1_R2_val_2.fq.gz -k 19 -s 3 -T 28 -i /home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.lt200.transcripts.fasta -o /home/zhzhang/PG/ISOseqdata/Mus_musculus/corrected_transcripts/mmALL.corrected.transcripts.fasta
#筛选长度>=200的序列
seqkit seq /home/zhzhang/PG/ISOseqdata/Mus_musculus/corrected_transcripts/mmALL.corrected.transcripts.fasta -m 200 > /home/zhzhang/PG/ISOseqdata/Mus_musculus/corrected_transcripts/mmALL.lt200.corrected.transcripts.fasta

#鸡
source /home/zhzhang/miniconda3/bin/activate SEQ
lordec-correct -2 /home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_brain_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_brain_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_brain_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_brain_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_cerebellum_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_cerebellum_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_cerebellum_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_cerebellum_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_gonad_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_gonad_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_gonad_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_gonad_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_gut_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_gut_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_gut_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_gut_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_heart_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_heart_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_heart_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Gallus_gallus/qc_fq/Chicken_heart_male_1_R2_val_2.fq.gz -k 19 -s 3 -T 28 -i /home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.lt200.transcripts.fasta -o /home/zhzhang/PG/ISOseqdata/Gallus_gallus/corrected_transcripts/ggALL.corrected.transcripts.fasta
#筛选长度>=200的序列
seqkit seq /home/zhzhang/PG/ISOseqdata/Gallus_gallus/corrected_transcripts/ggALL.corrected.transcripts.fasta -m 200 > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/corrected_transcripts/ggALL.lt200.corrected.transcripts.fasta


#斑马鱼
source /home/zhzhang/miniconda3/bin/activate SEQ
lordec-correct -2 /home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_brain_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_brain_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_brain_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_brain_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_gonad_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_gonad_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_gonad_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_gonad_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_gut_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_gut_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_gut_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_gut_male_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_heart_female_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_heart_female_1_R2_val_2.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_heart_male_1_R1_val_1.fq.gz,/home/zhzhang/PG/RNAseqdata/Danio_rerio/qc_fq/Zebrafish_heart_male_1_R2_val_2.fq.gz -k 19 -s 3 -T 28 -i /home/zhzhang/PG/ISOseqdata/Danio_rerio/fl_transcripts/drALL.lt200.transcripts.fasta -o /home/zhzhang/PG/ISOseqdata/Danio_rerio/corrected_transcripts/drALL.corrected.transcripts.fasta
#筛选长度>=200的序列
seqkit seq /home/zhzhang/PG/ISOseqdata/Danio_rerio/corrected_transcripts/drALL.corrected.transcripts.fasta -m 200 > /home/zhzhang/PG/ISOseqdata/Danio_rerio/corrected_transcripts/drALL.lt200.corrected.transcripts.fasta




```


##### 6.全长转录本序列比对
```r
#小鼠
source /home/zhzhang/miniconda3/bin/activate SEQ
minimap2 -t 64 -ax splice --secondary=no -C5 -O6,24 -B4 -uf /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa /home/zhzhang/PG/ISOseqdata/Mus_musculus/corrected_transcripts/mmALL.lt200.corrected.transcripts.fasta |samtools sort -o /home/zhzhang/PG/ISOseqdata/Mus_musculus/ct_to_ref/mmALL_ct_to_ref.bam

#鸡
source /home/zhzhang/miniconda3/bin/activate SEQ
minimap2 -t 64 -ax splice --secondary=no -C5 -O6,24 -B4 -uf /home/zhzhang/PG/refgenome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.chr.fa /home/zhzhang/PG/ISOseqdata/Gallus_gallus/corrected_transcripts/ggALL.lt200.corrected.transcripts.fasta |samtools sort -o /home/zhzhang/PG/ISOseqdata/Gallus_gallus/ct_to_ref/ggALL_ct_to_ref.bam

#斑马鱼
source /home/zhzhang/miniconda3/bin/activate SEQ
minimap2 -t 64 -ax splice --secondary=no -C5 -O6,24 -B4 -uf /home/zhzhang/PG/refgenome/Danio_rerio.GRCz11.dna.chr.fa /home/zhzhang/PG/ISOseqdata/Danio_rerio/corrected_transcripts/drALL.lt200.corrected.transcripts.fasta |samtools sort -o /home/zhzhang/PG/ISOseqdata/Danio_rerio/ct_to_ref/drALL_ct_to_ref.bam




```


##### 7.collapse全长转录本注释
```r
#小鼠
source /home/zhzhang/miniconda3/bin/activate anaCogent
collapse_isoforms_by_sam.py --input /home/zhzhang/PG/ISOseqdata/Mus_musculus/corrected_transcripts/mmALL.lt200.corrected.transcripts.fasta -b /home/zhzhang/PG/ISOseqdata/Mus_musculus/ct_to_ref/mmALL_ct_to_ref.bam -o /home/zhzhang/PG/ISOseqdata/Mus_musculus/gtf/mmALL --dun-merge-5-shorter -c 0.85 -i 0.95 --cpus 28
#删除5端降解转录本
get_abundance_post_collapse.py /home/zhzhang/PG/ISOseqdata/Mus_musculus/gtf/mmALL.collapsed /home/zhzhang/PG/ISOseqdata/Mus_musculus/fl_transcripts/mmALL.unpolished.transcripts.cluster_report.csv
filter_away_subset.py /home/zhzhang/PG/ISOseqdata/Mus_musculus/gtf/mmALL.collapsed

#鸡
source /home/zhzhang/miniconda3/bin/activate anaCogent
collapse_isoforms_by_sam.py --input /home/zhzhang/PG/ISOseqdata/Gallus_gallus/corrected_transcripts/ggALL.lt200.corrected.transcripts.fasta -b /home/zhzhang/PG/ISOseqdata/Gallus_gallus/ct_to_ref/ggALL_ct_to_ref.bam -o /home/zhzhang/PG/ISOseqdata/Gallus_gallus/gtf/ggALL --dun-merge-5-shorter -c 0.85 -i 0.95 --cpus 28
#删除5端降解转录本
get_abundance_post_collapse.py /home/zhzhang/PG/ISOseqdata/Gallus_gallus/gtf/ggALL.collapsed /home/zhzhang/PG/ISOseqdata/Gallus_gallus/fl_transcripts/ggALL.unpolished.transcripts.cluster_report.csv
filter_away_subset.py /home/zhzhang/PG/ISOseqdata/Gallus_gallus/gtf/ggALL.collapsed


#斑马鱼
source /home/zhzhang/miniconda3/bin/activate anaCogent
collapse_isoforms_by_sam.py --input /home/zhzhang/PG/ISOseqdata/Danio_rerio/corrected_transcripts/drALL.lt200.corrected.transcripts.fasta -b /home/zhzhang/PG/ISOseqdata/Danio_rerio/ct_to_ref/drALL_ct_to_ref.bam -o /home/zhzhang/PG/ISOseqdata/Danio_rerio/gtf/drALL --dun-merge-5-shorter -c 0.85 -i 0.95 --cpus 28




```
##### 
##### 8.novel genes筛选
```r
#排除掉参考注释中的假基因(原注释中假基因没有明确分类lnc还是pc，没有用，又会导致测到的假基因来源的转录本会被归纳到kown范畴，而不是novel。最终导致假基因转录本的损失)
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.chr.gtf"|grep -v pseudogene|awk '$1!= "MT" {print $0}' > /home/zhzhang/PG/HPG/GTF_PG/Gallus_gallus.GRCg7b.108.chr.rmpg.gtf
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/Homo_sapiens.GRCh38.108.chr.gtf"|grep -v pseudogene|awk '$1!= "MT" {print $0}' > /home/zhzhang/PG/HPG/GTF_PG/Homo_sapiens.GRCh38.108.chr.rmpg.gtf
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.chr.gtf"|grep -v pseudogene|awk '$1!= "MT" {print $0}' > /home/zhzhang/PG/HPG/GTF_PG/Mus_musculus.GRCm39.108.chr.rmpg.gtf
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/Danio_rerio.GRCz11.108.chr.gtf"|grep -v pseudogene|awk '$1!= "MT" {print $0}' > "/home/zhzhang/PG/HPG/GTF_PG/Danio_rerio.GRCz11.108.chr.rmpg.gtf"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/Macaca_mulatta.Mmul_10.108.chr.gtf"|grep -v pseudogene|awk '$1!= "MT" {print $0}' > /home/zhzhang/PG/HPG/GTF_PG/Macaca_mulatta.Mmul_10.108.chr.rmpg.gtf
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/Rattus_norvegicus.mRatBN7.2.108.chr.gtf"|grep -v pseudogene|awk '$1!= "MT" {print $0}' > /home/zhzhang/PG/HPG/GTF_PG/Rattus_norvegicus.mRatBN7.2.108.chr.rmpg.gtf
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/Oryctolagus_cuniculus.OryCun2.0.108.chr.gtf"|grep -v pseudogene|awk '$1!= "MT" {print $0}' > /home/zhzhang/PG/HPG/GTF_PG/Oryctolagus_cuniculus.OryCun2.0.108.chr.rmpg.gtf
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/Monodelphis_domestica.ASM229v1.108.chr.gtf"|grep -v pseudogene|awk '$1!= "MT" {print $0}' > /home/zhzhang/PG/HPG/GTF_PG/Monodelphis_domestica.ASM229v1.108.chr.rmpg.gtf






#从参考基因组注释GTF文件，提取已注释的蛋白编码基因/lncRNA基因/其他ncRNA基因（除PG外所有注释的基因）bed文件【chr start end geneid gene_biotype strand】
#鸡
cat /home/zhzhang/PG/HPG/GTF_PG/Gallus_gallus.GRCg7b.108.chr.rmpg.gtf|awk '$3=="gene" && $17=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$18"\t"$7} $3=="gene" && $15=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$16"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.ref_gene.bed
#人
cat /home/zhzhang/PG/HPG/GTF_PG/Homo_sapiens.GRCh38.108.chr.rmpg.gtf|awk '$3=="gene" && $17=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$18"\t"$7} $3=="gene" && $15=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$16"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.ref_gene.bed
#小鼠
cat /home/zhzhang/PG/HPG/GTF_PG/Mus_musculus.GRCm39.108.chr.rmpg.gtf|awk '$3=="gene" && $17=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$18"\t"$7} $3=="gene" && $15=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$16"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.ref_gene.bed
#斑马鱼
cat /home/zhzhang/PG/HPG/GTF_PG/Danio_rerio.GRCz11.108.chr.rmpg.gtf|awk '$3=="gene" && $17=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$18"\t"$7} $3=="gene" && $15=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$16"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.ref_gene.bed
#猕猴
cat "/home/zhzhang/PG/HPG/GTF_PG/Macaca_mulatta.Mmul_10.108.chr.rmpg.gtf"|awk '$3=="gene" && $17=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$18"\t"$7} $3=="gene" && $15=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$16"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Macaca_mulatta/novel_lnc/Macaca_mulatta.ref_gene.bed
#大鼠
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/Rattus_norvegicus.mRatBN7.2.108.chr.gtf"|awk '$3=="gene" && $17=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$18"\t"$7} $3=="gene" && $15=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$16"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Rattus_norvegicus/novel_lnc/Rattus_norvegicus.ref_gene.bed
#兔
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/Oryctolagus_cuniculus.OryCun2.0.108.chr.gtf"|awk '$3=="gene" && $17=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$18"\t"$7} $3=="gene" && $15=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$16"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Oryctolagus_cuniculus/novel_lnc/Oryctolagus_cuniculus.ref_gene.bed
#负鼠
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/Monodelphis_domestica.ASM229v1.108.chr.gtf"|awk '$3=="gene" && $17=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$18"\t"$7} $3=="gene" && $15=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$16"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Monodelphis_domestica/novel_lnc/Monodelphis_domestica.ref_gene.bed



#从三代转录本gtf中提取三代转录本bed【chr start end geneid transcriptid strand】
#斑马鱼
cat "/home/zhzhang/PG/ISOseqdata/Danio_rerio/gtf/drALL.collapsed.gff"|awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_transcripts.bed
#鸡
cat "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/gtf/ggALL.collapsed.filtered.gff"|awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_transcripts.bed
#小鼠
cat "/home/zhzhang/PG/ISOseqdata/Mus_musculus/gtf/mmALL.collapsed.filtered.gff"|awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_transcripts.bed
#人
#数据来源于https://www.nature.com/articles/s41586-022-05035-y#data-availability（https://gtexportal.org/home/datasets）https://storage.googleapis.com/gtex_analysis_v9/long_read_data/flair_filter_transcripts.gtf.gz
#去掉chrm转录本，染色体只保留数字，与参考注释保持格式一致
cat "/home/zhzhang/PG/ISOseqdata/Homo_sapiens/gtf/flair_filter_transcripts.gtf"|awk '$1!="chrM" {print $0}'|sed "s/^chr//g" > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/gtf/flair_filter_transcripts.chrmatchref.gtf
#转录本id中的"_chr数字:数字"替换成""
sed -i "s/_chr.:[0-9]*//g;s/_chr..:[0-9]*//g" /home/zhzhang/PG/ISOseqdata/Homo_sapiens/gtf/flair_filter_transcripts.chrmatchref.gtf
#转录本id中的"_ENSG00000127527.13"替换成"-ENSG"
sed -i "s/_ENSG[0-9]*\.[0-9]\{1,2\}/-ENSG/g" /home/zhzhang/PG/ISOseqdata/Homo_sapiens/gtf/flair_filter_transcripts.chrmatchref.gtf
#提取
cat /home/zhzhang/PG/ISOseqdata/Homo_sapiens/gtf/flair_filter_transcripts.chrmatchref.gtf|awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_transcripts.bed




#筛选出与ref基因没有交集的三代新基因
#斑马鱼
#转录本与ref基因进行bedtools交集，获得与ref基因交集的geneid
bedtools intersect -a /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_transcripts.bed -b /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.ref_gene.bed -wa -s |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_transcripts_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取novelgene的转录本及其外显子的GTF注释
cat "/home/zhzhang/PG/ISOseqdata/Danio_rerio/gtf/drALL.collapsed.gff"|grep -v -w -Ff /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_transcripts_overlap_refgene.id.txt > /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.gtf


#鸡
#转录本与ref基因进行bedtools交集，获得与ref基因交集的geneid
bedtools intersect -a /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_transcripts.bed -b /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.ref_gene.bed -wa -s |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_transcripts_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取novelgene的转录本及其外显子的GTF注释
cat "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/gtf/ggALL.collapsed.gff"|grep -v -w -Ff /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_transcripts_overlap_refgene.id.txt > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.gtf

#小鼠
#转录本与ref基因进行bedtools交集，获得与ref基因交集的geneid
bedtools intersect -a /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_transcripts.bed -b /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.ref_gene.bed -wa -s |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_transcripts_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取novelgene的转录本及其外显子的GTF注释
cat "/home/zhzhang/PG/ISOseqdata/Mus_musculus/gtf/mmALL.collapsed.gff"|grep -v -w -Ff /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_transcripts_overlap_refgene.id.txt > /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.gtf

#人
#转录本与ref基因进行bedtools交集，获得与ref基因交集的geneid
bedtools intersect -a /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_transcripts.bed -b /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.ref_gene.bed -wa -s |awk '{print $4}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_transcripts_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取novelgene的转录本及其外显子的GTF注释
cat "/home/zhzhang/PG/ISOseqdata/Homo_sapiens/gtf/flair_filter_transcripts.chrmatchref.gtf"|grep -v -w -Ff /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_transcripts_overlap_refgene.id.txt > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.gtf


```


##### 9.novel lncRNA筛选-不具有蛋白编码能力
```r
#获取novel  gene的转录本序列
#鸡
gffread /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.gtf -g "/home/zhzhang/PG/refgenome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.chr.fa" -w /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.fa
#小鼠
gffread /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.gtf -g "/home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa" -w /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.fa
#人
gffread /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.gtf -g "/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa" -w /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.fa
#斑马鱼
gffread /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.gtf -g "/home/zhzhang/PG/refgenome/Danio_rerio.GRCz11.dna.chr.fa" -w /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.fa



#对每条novel  gene的转录本序列进行编码潜力检测
1.PLEK
python /home/zhzhang/software/PLEK.1.2/PLEK.py -fasta transcript.fa -out plek_output.txt -thread 28
#鸡
python /home/zhzhang/software/PLEK.1.2/PLEK.py -fasta /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.fa -out /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.plek.txt -thread 28
#小鼠
python /home/zhzhang/software/PLEK.1.2/PLEK.py -fasta /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.fa -out /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.plek.txt -thread 28
#人
python /home/zhzhang/software/PLEK.1.2/PLEK.py -fasta /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.fa -out /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.plek.txt -thread 28
#斑马鱼
python /home/zhzhang/software/PLEK.1.2/PLEK.py -fasta /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.fa -out /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.plek.txt -thread 28


2.LGC
python /home/zhzhang/software/LGC-1.0.py input.fasta LGC_output.txt -p 线程数 > LGC.log 2>&1
#鸡
python /home/zhzhang/software/LGC-1.0.py /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.fa /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.lgc.txt -p 28
#小鼠
python /home/zhzhang/software/LGC-1.0.py /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.fa /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.lgc.txt -p 28
#人
python /home/zhzhang/software/LGC-1.0.py /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.fa /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.lgc.txt -p 28
#斑马鱼
python /home/zhzhang/software/LGC-1.0.py /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.fa /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.lgc.txt -p 28


3.CPC2
python /home/zhzhang/software/CPC2-master/bin/CPC2.py -i transcript.fa -o CPC2_result.txt > cpc2.log 2>&1
#鸡
python /home/zhzhang/software/CPC2-master/bin/CPC2.py -i /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.fa -o /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.cpc2.txt
#小鼠
python /home/zhzhang/software/CPC2-master/bin/CPC2.py -i /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.fa -o /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.cpc2.txt
#人
python /home/zhzhang/software/CPC2-master/bin/CPC2.py -i /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.fa -o /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.cpc2.txt
#斑马鱼
python /home/zhzhang/software/CPC2-master/bin/CPC2.py -i /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.fa -o /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.cpc2.txt


4.CNCI
python /home/zhzhang/software/CNCI-master/CNCI.py -f transcript.fasta -o CNCI_result_dir -m ve -p 线程数 > cnci.log 2>&1
#-f指定转录本文件，可以是fasta格式，也可以是gtf格式，如果是gtf格式，需要同时指定-g和-d参数,-d参数需要指定2bit格式的参考基因组序列索引前缀路径;-p参数指定并行的CPU个数;-m指定使用的模型,ve代表脊椎动物，pl代表植物；-o指定输出结果的目录
#鸡
python /home/zhzhang/software/CNCI-master/CNCI.py -f /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.fa -o /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.CNCI -m ve -p 28
mv "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.CNCI/CNCI.index" /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.CNCI.txt
#小鼠
python /home/zhzhang/software/CNCI-master/CNCI.py -f /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.fa -o /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.CNCI -m ve -p 28
mv "/home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.CNCI/CNCI.index" /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.CNCI.txt
#人
python /home/zhzhang/software/CNCI-master/CNCI.py -f /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.fa -o /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.CNCI -m ve -p 28
mv "/home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.CNCI/CNCI.index" /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.CNCI.txt
#斑马鱼
python /home/zhzhang/software/CNCI-master/CNCI.py -f /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.fa -o /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.CNCI -m ve -p 28
mv "/home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.CNCI/CNCI.index" /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.CNCI.txt


#筛选出全部转录本都被四个软件预测non-coding（都是lncRNA）的gene，即潜在的新lncRNA基因
##获得转录本ID和基因ID对照信息表
#鸡
cat "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.gtf"|awk '$3=="transcript" {print $10"\t"$12}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.transcript_gene_id.txt
#小鼠
cat "/home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.gtf"|awk '$3=="transcript" {print $10"\t"$12}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.transcript_gene_id.txt
#人
cat "/home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.gtf"|awk '$3=="transcript" {print $12"\t"$10}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.transcript_gene_id.txt
#斑马鱼
cat "/home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.gtf"|awk '$3=="transcript" {print $10"\t"$12}'|sed "s/\"//g;s/;//g" > /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.transcript_gene_id.txt

##获得全部转录本isoform都被四个软件预测non-coding的geneID
#鸡
a <- "/home/zhzhang/PG/isoseq/Gallus_gallus/Gallus_gallus.3G_novel_transcripts.transcript_gene_id.txt"
b <- "/home/zhzhang/PG/isoseq/Gallus_gallus/Gallus_gallus.3G_novel_transcripts.CNCI.txt"
c <- "/home/zhzhang/PG/isoseq/Gallus_gallus/Gallus_gallus.3G_novel_transcripts.cpc2.txt"
d <- "/home/zhzhang/PG/isoseq/Gallus_gallus/Gallus_gallus.3G_novel_transcripts.plek.txt"
e <- "/home/zhzhang/PG/isoseq/Gallus_gallus/Gallus_gallus.3G_novel_transcripts.lgc.txt"
out <- "/home/zhzhang/PG/isoseq/Gallus_gallus/Gallus_gallus.3G_novel_plncgene.id.txt"
#小鼠
a <- "/home/zhzhang/PG/isoseq/Mus_musculus/Mus_musculus.3G_novel_transcripts.transcript_gene_id.txt"
b <- "/home/zhzhang/PG/isoseq/Mus_musculus/Mus_musculus.3G_novel_transcripts.CNCI.txt"
c <- "/home/zhzhang/PG/isoseq/Mus_musculus/Mus_musculus.3G_novel_transcripts.cpc2.txt"
d <- "/home/zhzhang/PG/isoseq/Mus_musculus/Mus_musculus.3G_novel_transcripts.plek.txt"
e <- "/home/zhzhang/PG/isoseq/Mus_musculus/Mus_musculus.3G_novel_transcripts.lgc.txt"
out <- "/home/zhzhang/PG/isoseq/Mus_musculus/Mus_musculus.3G_novel_plncgene.id.txt"
#人
a <- "/home/zhzhang/PG/isoseq/Homo_sapiens/Homo_sapiens.3G_novel_transcripts.transcript_gene_id.txt"
b <- "/home/zhzhang/PG/isoseq/Homo_sapiens/Homo_sapiens.3G_novel_transcripts.CNCI.txt"
c <- "/home/zhzhang/PG/isoseq/Homo_sapiens/Homo_sapiens.3G_novel_transcripts.cpc2.txt"
d <- "/home/zhzhang/PG/isoseq/Homo_sapiens/Homo_sapiens.3G_novel_transcripts.plek.txt"
e <- "/home/zhzhang/PG/isoseq/Homo_sapiens/Homo_sapiens.3G_novel_transcripts.lgc.txt"
out <- "/home/zhzhang/PG/isoseq/Homo_sapiens/Homo_sapiens.3G_novel_plncgene.id.txt"
#斑马鱼
a <- "/home/zhzhang/PG/isoseq/Danio_rerio/Danio_rerio.3G_novel_transcripts.transcript_gene_id.txt"
b <- "/home/zhzhang/PG/isoseq/Danio_rerio/Danio_rerio.3G_novel_transcripts.CNCI.txt"
c <- "/home/zhzhang/PG/isoseq/Danio_rerio/Danio_rerio.3G_novel_transcripts.cpc2.txt"
d <- "/home/zhzhang/PG/isoseq/Danio_rerio/Danio_rerio.3G_novel_transcripts.plek.txt"
e <- "/home/zhzhang/PG/isoseq/Danio_rerio/Danio_rerio.3G_novel_transcripts.lgc.txt"
out <- "/home/zhzhang/PG/isoseq/Danio_rerio/Danio_rerio.3G_novel_plncgene.id.txt"

#导入转录本id-基因id对照表
transcript_gene_id <- read.delim(a, header=FALSE)
colnames(transcript_gene_id) <- c("tid","gid")
#导入CNCI结果
CNCI <- read.delim(b)%>%
  select(1,2)
colnames(CNCI) <- c("tid","CNCI")
#导入cpc2结果
CPC2 <- read.delim(c)%>%
  select(1,8)
colnames(CPC2) <- c("tid","CPC2")
#导入plek结果
PLEK <- read.delim(d, header=FALSE)%>%
  separate(V3,c("no","tid"),sep = '>')%>%
  select(4,1)
colnames(PLEK) <- c("tid","PLEK")
#导入cpc2结果
LGC <- read.delim(e, header=FALSE, comment.char="#")%>%
  select(1,5)
colnames(LGC) <- c("tid","LGC")
#合并
he <- left_join(transcript_gene_id,CNCI,by="tid")%>%
  left_join(CPC2,by="tid")%>%
  left_join(PLEK,by="tid")%>%
  left_join(LGC,by="tid")%>%
  mutate(lnc=0,tnum=1)
#被四个软件判断为noncoding的转录本，lnc属性设定为1（是lncRNA转录本）
he$lnc[he$CNCI=="noncoding" & he$CPC2=="noncoding" & he$PLEK=="Non-coding" & he$LGC=="Non-coding"] <- 1
#统计基因的转录本数以及lnc转录本数，数目相同的为全部转录本都是lnc的基因
lncgene <- group_by(he,gid)%>%
  summarise(lnctnum=sum(lnc),alltnum=sum(tnum))%>%
  filter(lnctnum==alltnum)%>%
  select(1)
#保存潜在的lncRNA基因ID
data.table::fwrite(lncgene,file =out,sep = '\t',row.names = F,quote = F,col.names = F)
##根据geneID提取潜在的3Gnovel lncrna gene对应的转录本及其exon注释GTF
cat "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.gtf"|grep -w -Ff "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.id.txt" > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.gtf
cat "/home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.gtf"|grep -w -Ff /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.id.txt > /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.gtf
cat "/home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.gtf"|grep -w -Ff /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.id.txt > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.gtf
cat "/home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.gtf"|grep -w -Ff "/home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.id.txt" > /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.gtf



```


##### 10.novel lncRNA筛选-不具有pfam蛋白结构域
```r
###筛选出全部转录本都不具有pfam蛋白结构域的基因
#获取转录本序列
#斑马鱼
gffread /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.gtf -g "/home/zhzhang/PG/refgenome/Danio_rerio.GRCz11.dna.chr.fa" -w /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.fa
#鸡
gffread /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.gtf -g "/home/zhzhang/PG/refgenome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.chr.fa" -w /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.fa
#小鼠
gffread /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.gtf -g "/home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa" -w /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.fa
#人
gffread /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.gtf -g "/home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa" -w /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.fa


#翻译成蛋白序列
#斑马鱼
transeq /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.fa /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.pep.fa -frame=F -clean Y -trim Y
#鸡
transeq /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.fa /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.pep.fa -frame=F -clean Y -trim Y
#小鼠
transeq /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.fa /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.pep.fa -frame=F -clean Y -trim Y
#人
transeq /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.fa /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.pep.fa -frame=F -clean Y -trim Y



#获得具有pfam蛋白结构域的转录本相关信息
#斑马鱼
interproscan.sh -cpu 26 -dp -appl Pfam -dra -f TSV -i /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.pep.fa -o /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.ownpfam.tsv
#鸡
interproscan.sh -cpu 26 -dp -appl Pfam -dra -f TSV -i /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.pep.fa -o /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.ownpfam.tsv
#小鼠
interproscan.sh -cpu 26 -dp -appl Pfam -dra -f TSV -i /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.pep.fa -o /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.ownpfam.tsv
#人
interproscan.sh -cpu 26 -dp -appl Pfam -dra -f TSV -i /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.pep.fa -o /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.ownpfam.tsv

#根据interproscan 输出结果，提取具有pfam蛋白结构域（E-value < 1e-5）的转录本对应的geneID
#鸡
cat /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.ownpfam.tsv|awk -F "\t" '$9<1e-5 {print $1}'|awk -F "_" '{print $1}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.ownpfam.tid.txt
cat "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_transcripts.transcript_gene_id.txt"|grep -w -Ff /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.ownpfam.tid.txt|awk '{print $2}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.ownpfam.gid.txt
#小鼠
cat /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.ownpfam.tsv|awk -F "\t" '$9<1e-5 {print $1}'|awk -F "_" '{print $1}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.ownpfam.tid.txt
cat "/home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_transcripts.transcript_gene_id.txt"|grep -w -Ff /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.ownpfam.tid.txt|awk '{print $2}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.ownpfam.gid.txt
#人
cat /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.ownpfam.tsv|awk -F "\t" '$9<1e-5 {print $1}'|awk -F "_" '{print $1}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.ownpfam.tid.txt
cat "/home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_transcripts.transcript_gene_id.txt"|grep -w -Ff /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.ownpfam.tid.txt|awk '{print $2}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.ownpfam.gid.txt
#斑马鱼
cat /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.ownpfam.tsv|awk -F "\t" '$9<1e-5 {print $1}'|awk -F "_" '{print $1}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.ownpfam.tid.txt
cat "/home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_transcripts.transcript_gene_id.txt"|grep -w -Ff /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.ownpfam.tid.txt|awk '{print $2}'|sort|uniq > /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.ownpfam.gid.txt

#根据具有pfam蛋白结构域的转录本对应的基因ID，反向匹配出全部转录本都不具有pfam蛋白结构域的基因GTF注释，即3G novel lncrna gene(三种翻译方式均不产生pfam结构域的转录本)
#鸡
cat "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.gtf"|grep -w -v -Ff "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_plncgene.ownpfam.gid.txt" > /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_lncRNA.gtf
#小鼠
cat "/home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.gtf"|grep -w -v -Ff /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_plncgene.ownpfam.gid.txt > /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_lncRNA.gtf
#人
cat "/home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.gtf"|grep -w -v -Ff /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_plncgene.ownpfam.gid.txt > /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_lncRNA.gtf
#斑马鱼
cat "/home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.gtf"|grep -w -v -Ff "/home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_plncgene.ownpfam.gid.txt" > /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_lncRNA.gtf

```


##### 11.下载处理基于RNA-seq的lncRNA注释, 筛选出与Ref-gene以及3代novel-lncRNA没有交集的二代lncRNA GENE
```r
#二代RNA-seq数据组装的lncRNA基因集（包括已知lncRNA和novel lncRNA），外显子水平gtf注释来源于（https://www.nature.com/articles/s41586-019-1341-x#Sec9）
#二代组装的lnc基因注释liftover到新版基因组
##lnc exon水平gtf转bed【最后一列是外显子顺序】，染色体号转为符合UCSC chain文件的形式
#鸡
cat "/home/zhzhang/PG/2G_lnc/chicken.lncRNA.gtf"|awk '$1~/^[0-9]/ {print "chr"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14} $1~/^[WZ]/ {print "chr"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14} $1~/^LGE.*/ {print "chr"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14} $1~/^AADN.*/ {print "chrUn_"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14} $1~/^JH.*/ {print "chrUn_"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14}'|sed "s/\;//g" > /home/zhzhang/PG/2G_lnc/chicken.lncRNA.bed
#人
cat "/home/zhzhang/PG/2G_lnc/human.lncRNA.gtf"|awk '{print "chr"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14}'|sed "s/\;//g" > /home/zhzhang/PG/2G_lnc/human.lncRNA.bed
#小鼠
cat "/home/zhzhang/PG/2G_lnc/mouse.lncRNA.gtf"|awk '$1~/^[0-9]/ {print "chr"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14} $1~/^[XY]/ {print "chr"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14}'|sed "s/\;//g" > "/home/zhzhang/PG/2G_lnc/mouse.lncRNA.bed"
##外显子bed生成转录本bed
#a是二代lnc exon bed路径，b是输出的lnc 转录本bed路径
lncexon2t <- function(a,b){
  #导入二代lnc exon bed
  lncRNA <- read.delim(a, header=FALSE)
  #lnc 转录本bed
  lncRNAt <- group_by(lncRNA,V1,V4,V6)%>%
    summarise(start=min(V2),end=max(V3))%>%
    mutate(a="transcript")%>%
    select(1,4,5,2,6,3)
  #储存
  data.table::fwrite(lncRNAt,file =b,sep = '\t',row.names = F,quote = F,col.names = F)
}
#鸡
lncexon2t(a="~/PG/2G_lnc/chicken.lncRNA.bed",
          b="~/PG/2G_lnc/chicken.lncRNA.transcript.bed")
#人
lncexon2t(a="~/PG/2G_lnc/human.lncRNA.bed",
          b="~/PG/2G_lnc/human.lncRNA.transcript.bed")
#小鼠
lncexon2t(a="~/PG/2G_lnc/mouse.lncRNA.bed",
          b="~/PG/2G_lnc/mouse.lncRNA.transcript.bed")
##liftover转换lnc转录本到新版基因组。获得无法正确转换的转录本id(为了避免出现以下情况：lift over后转录本结构出现错误（例如：liftover后一个转录本的外显子出现在不同的染色体上，外显子顺序改变，或外显子之间间距过大）)
#鸡
liftOver /home/zhzhang/PG/2G_lnc/chicken.lncRNA.transcript.bed /home/zhzhang/software/liftover/galGal4ToGalGal6.over.chain.gz /home/zhzhang/PG/2G_lnc/chicken.lncRNA_GalGal6.transcript.bed /home/zhzhang/PG/2G_lnc/lncTgalGal4ToGalGal6.unmap.txt
liftOver /home/zhzhang/PG/2G_lnc/chicken.lncRNA_GalGal6.transcript.bed "/home/zhzhang/software/liftover/galGal6ToGRCg7b.over.chain.gz" /home/zhzhang/PG/2G_lnc/chicken.lncRNA_GRCg7b.transcript.bed /home/zhzhang/PG/2G_lnc/lncTgalGal6ToGRCg7b.unmap.txt
rm /home/zhzhang/PG/2G_lnc/chicken.lncRNA_GalGal6.transcript.bed
rm /home/zhzhang/PG/2G_lnc/chicken.lncRNA_GRCg7b.transcript.bed
cat /home/zhzhang/PG/2G_lnc/lncTgalGal4ToGalGal6.unmap.txt|awk '$5=="transcript" {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Gallus_gallus.unmaplnctid.txt
cat /home/zhzhang/PG/2G_lnc/lncTgalGal6ToGRCg7b.unmap.txt|awk '$5=="transcript" {print $4}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Gallus_gallus.unmaplnctid.txt
#人
liftOver /home/zhzhang/PG/2G_lnc/human.lncRNA.transcript.bed /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/2G_lnc/human.lncRNA_hg38.transcript.bed /home/zhzhang/PG/2G_lnc/lncThg19ToHg38.unmap.txt
rm /home/zhzhang/PG/2G_lnc/human.lncRNA_hg38.transcript.bed
cat /home/zhzhang/PG/2G_lnc/lncThg19ToHg38.unmap.txt|awk '$5=="transcript" {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Homo_sapiens.unmaplnctid.txt
#小鼠
liftOver /home/zhzhang/PG/2G_lnc/mouse.lncRNA.transcript.bed /home/zhzhang/software/liftover/mm10ToMm39.over.chain.gz /home/zhzhang/PG/2G_lnc/mouse.lncRNA_mm39.transcript.bed /home/zhzhang/PG/2G_lnc/lncTmm10ToMm39.unmap.txt
rm /home/zhzhang/PG/2G_lnc/mouse.lncRNA_mm39.transcript.bed
cat /home/zhzhang/PG/2G_lnc/lncTmm10ToMm39.unmap.txt|awk '$5=="transcript" {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Mus_musculus.unmaplnctid.txt
##筛选出可以转换至新版基因组的lnc转录本的外显子 bed文件，并转换至新版基因组
#鸡
cat "/home/zhzhang/PG/2G_lnc/chicken.lncRNA.bed"|grep -w -v -Ff "/home/zhzhang/PG/2G_lnc/Gallus_gallus.unmaplnctid.txt" > /home/zhzhang/PG/2G_lnc/chicken.lncRNA.ok.bed
liftOver /home/zhzhang/PG/2G_lnc/chicken.lncRNA.ok.bed /home/zhzhang/software/liftover/galGal4ToGalGal6.over.chain.gz /home/zhzhang/PG/2G_lnc/chicken.lncRNA_GalGal6.bed /home/zhzhang/PG/2G_lnc/galGal4ToGalGal6.unmap.txt
liftOver /home/zhzhang/PG/2G_lnc/chicken.lncRNA_GalGal6.bed "/home/zhzhang/software/liftover/galGal6ToGRCg7b.over.chain.gz" /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.bed /home/zhzhang/PG/2G_lnc/galGal6ToGRCg7b.unmap.txt
#NCBI染色体号转正常数字（ENSEMBLE格式）
cat /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.bed|sed "s/NC_052532.1/1/g;s/NC_052533.1/2/g;s/NC_052534.1/3/g;s/NC_052535.1/4/g;s/NC_052536.1/5/g;s/NC_052537.1/6/g;s/NC_052538.1/7/g;s/NC_052539.1/8/g;s/NC_052540.1/9/g;s/NC_052541.1/10/g;s/NC_052542.1/11/g;s/NC_052543.1/12/g;s/NC_052544.1/13/g;s/NC_052545.1/14/g;s/NC_052546.1/15/g;s/NC_052547.1/16/g;s/NC_052548.1/17/g;s/NC_052549.1/18/g;s/NC_052550.1/19/g;s/NC_052551.1/20/g;s/NC_052552.1/21/g;s/NC_052553.1/22/g;s/NC_052554.1/23/g;s/NC_052555.1/24/g;s/NC_052556.1/25/g;s/NC_052557.1/26/g;s/NC_052558.1/27/g;s/NC_052559.1/28/g;s/NC_052560.1/29/g;s/NC_052561.1/30/g;s/NC_052562.1/31/g;s/NC_052563.1/32/g;s/NC_052564.1/33/g;s/NC_052565.1/34/g;s/NC_052566.1/35/g;s/NC_052567.1/36/g;s/NC_052568.1/37/g;s/NC_052569.1/38/g;s/NC_052570.1/39/g;s/NC_052571.1/W/g;s/NC_052572.1/Z/g" > /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.chrENS.bed
#人
cat "/home/zhzhang/PG/2G_lnc/human.lncRNA.bed"|grep -w -v -Ff "/home/zhzhang/PG/2G_lnc/Homo_sapiens.unmaplnctid.txt" > /home/zhzhang/PG/2G_lnc/human.lncRNA.ok.bed
liftOver /home/zhzhang/PG/2G_lnc/human.lncRNA.ok.bed /home/zhzhang/software/liftover/hg19ToHg38.over.chain.gz /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.bed /home/zhzhang/PG/2G_lnc/hg19ToHg38.unmap.txt
#UCSC染色体号转正常数字（ENSEMBLE格式）
cat /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.bed |sed "s/chr//g" > /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.chrENS.bed
#小鼠
cat "/home/zhzhang/PG/2G_lnc/mouse.lncRNA.bed"|grep -w -v -Ff "/home/zhzhang/PG/2G_lnc/Mus_musculus.unmaplnctid.txt" > /home/zhzhang/PG/2G_lnc/mouse.lncRNA.ok.bed
liftOver /home/zhzhang/PG/2G_lnc/mouse.lncRNA.ok.bed /home/zhzhang/software/liftover/mm10ToMm39.over.chain.gz /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.bed /home/zhzhang/PG/2G_lnc/mm10ToMm39.unmap.txt
#UCSC染色体号转正常数字（ENSEMBLE格式）
cat /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.bed |sed "s/chr//g" > /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.chrENS.bed
##删除至少一个外显子匹配不到新版基因组核染色体的lnc转录本
#鸡
#获取至少一个外显子匹配不到新版基因组核染色体的lnc转录本的ID（lift over中损失的exon，以及因为转换后匹配到非核染色体上损失的exon）
cat "/home/zhzhang/PG/2G_lnc/galGal4ToGalGal6.unmap.txt"|awk '$4~/^Chi/ {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Gallus_gallus.unmaplncid.txt
cat "/home/zhzhang/PG/2G_lnc/galGal6ToGRCg7b.unmap.txt"|awk '$4~/^Chi/ {print $4}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Gallus_gallus.unmaplncid.txt
#删除上述转录本ID对应的全部exon
cat /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.chrENS.bed |grep -w -v -Ff /home/zhzhang/PG/2G_lnc/Gallus_gallus.unmaplncid.txt > /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.1.bed
#人
#获取至少一个外显子匹配不到新版基因组核染色体的lnc转录本的ID（lift over中损失的exon，以及因为转换后匹配到非核染色体上损失的exon）
cat "/home/zhzhang/PG/2G_lnc/hg19ToHg38.unmap.txt"|awk '$4~/^Hum/ {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Homo_sapiens.unmaplncid.txt
cat "/home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.chrENS.bed"|awk '$1~/_/ {print $4}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Homo_sapiens.unmaplncid.txt
#删除上述转录本ID对应的全部exon
cat /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.chrENS.bed |grep -w -v -Ff /home/zhzhang/PG/2G_lnc/Homo_sapiens.unmaplncid.txt > /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.1.bed
#小鼠
#获取至少一个外显子匹配不到新版基因组核染色体的lnc转录本的ID（lift over中损失的exon，以及因为转换后匹配到非核染色体上损失的exon）
cat "/home/zhzhang/PG/2G_lnc/mm10ToMm39.unmap.txt"|awk '$4~/^Mou/ {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Mus_musculus.unmaplncid.txt
#删除上述转录本ID对应的全部exon
cat /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.chrENS.bed |grep -w -v -Ff /home/zhzhang/PG/2G_lnc/Mus_musculus.unmaplncid.txt > /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.1.bed


#斑马鱼二代RNA-seq数据组装的lncRNA基因集（包括已知lncRNA和novel lncRNA）来源于fishget数据库，基因组版本一致不用liftover；
#提取参考基因组中没有的lncRNA基因ID（长度>=200，所以为lncRNA GENEid）
grep "non-coding" "/home/zhzhang/PG/2G_lnc/zebrafish_gene_summary.txt"|grep -v ENSDAR|awk '$4-$3 >= 200 {print $1}' > /home/zhzhang/PG/2G_lnc/zebrafish.lncRNA.ID.txt
#根据上述ID提取染色体上的lncRNA gtf
grep -w -Ff /home/zhzhang/PG/2G_lnc/zebrafish.lncRNA.ID.txt "/home/zhzhang/PG/2G_lnc/zebrafish.gtf"|grep "^NC"|sed 's/NC_007112.7/1/g;s/NC_007113.7/2/g;s/NC_007114.7/3/g;s/NC_007115.7/4/g;s/NC_007116.7/5/g;s/NC_007117.7/6/g;s/NC_007118.7/7/g;s/NC_007119.7/8/g;s/NC_007120.7/9/g;s/NC_007121.7/10/g;s/NC_007122.7/11/g;s/NC_007123.7/12/g;s/NC_007124.7/13/g;s/NC_007125.7/14/g;s/NC_007126.7/15/g;s/NC_007127.7/16/g;s/NC_007128.7/17/g;s/NC_007129.7/18/g;s/NC_007130.7/19/g;s/NC_007131.7/20/g;s/NC_007132.7/21/g;s/NC_007133.7/22/g;s/NC_007134.7/23/g;s/NC_007135.7/24/g;s/NC_007136.7/25/g' > "/home/zhzhang/PG/2G_lnc/zebrafish.lncRNA.gtf"
#转为bed
cat "/home/zhzhang/PG/2G_lnc/zebrafish.lncRNA.gtf"|awk '{print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7"\t"$14}'|sed 's/\"//g;s/\;//g' > /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA.1.bed


###筛选出与Ref-gene以及3代novel-lncRNA没有交集的二代lnc GENE
#斑马鱼
#二代lnc外显子bed与ref基因进行bedtools交集，获得与ref基因交集的二代lncRNA geneid
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA.1.bed -b /home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.ref_gene.bed -wa -s |awk '{print $5}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA_overlap_refgene.id.txt
#二代lnc外显子bed与三代novel lnc基因进行bedtools交集，获得与三代novel lnc基因交集的二代lncRNA geneid
cat "/home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_lncRNA.gtf"|awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/;//g;s/\"//g" > /home/zhzhang/PG/2G_lnc/Danio_rerio.3G_novel_lncRNAt.bed
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA.1.bed -b /home/zhzhang/PG/2G_lnc/Danio_rerio.3G_novel_lncRNAt.bed -wa -s |awk '{print $5}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取与Ref-gene以及3代novel-lncRNA没有交集的二代lncRNA GENE的外显子的bed注释
cat /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA.1.bed|grep -v -w -Ff /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA_overlap_refgene.id.txt > /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA.2.bed
#二代lncRNA GENE exon bed注释转换为gtf格式
cat /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA.2.bed|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$5"\"; transcript_id \""$4"\";"}' > /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA.gtf

#鸡
#二代lnc外显子bed与ref基因进行bedtools交集，获得与ref基因交集的二代lncRNA geneid
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.1.bed -b /home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.ref_gene.bed -wa -s |awk '{print $5}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA_overlap_refgene.id.txt
#二代lnc外显子bed与三代novel lnc基因进行bedtools交集，获得与三代novel lnc基因交集的二代lncRNA geneid
cat "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_lncRNA.gtf"|awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/;//g;s/\"//g" > /home/zhzhang/PG/2G_lnc/Gallus_gallus.3G_novel_lncRNAt.bed
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.1.bed -b /home/zhzhang/PG/2G_lnc/Gallus_gallus.3G_novel_lncRNAt.bed -wa -s |awk '{print $5}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取与Ref-gene以及3代novel-lncRNA没有交集的二代lncRNA GENE的外显子的bed注释
cat /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.1.bed|grep -v -w -Ff /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA_overlap_refgene.id.txt > /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.2.bed
#二代lncRNA GENE exon bed注释转换为gtf格式
cat /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.2.bed|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$5"\"; transcript_id \""$4"\";"}' > /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.gtf

#小鼠
#二代lnc外显子bed与ref基因进行bedtools交集，获得与ref基因交集的二代lncRNA geneid
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.1.bed -b /home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.ref_gene.bed -wa -s |awk '{print $5}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA_overlap_refgene.id.txt
#二代lnc外显子bed与三代novel lnc基因进行bedtools交集，获得与三代novel lnc基因交集的二代lncRNA geneid
cat "/home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_lncRNA.gtf"|awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/;//g;s/\"//g" > /home/zhzhang/PG/2G_lnc/Mus_musculus.3G_novel_lncRNAt.bed
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.1.bed -b /home/zhzhang/PG/2G_lnc/Mus_musculus.3G_novel_lncRNAt.bed -wa -s |awk '{print $5}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取与Ref-gene以及3代novel-lncRNA没有交集的二代lncRNA GENE的外显子的bed注释
cat /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.1.bed|grep -v -w -Ff /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA_overlap_refgene.id.txt > /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.2.bed
#二代lncRNA GENE exon bed注释转换为gtf格式
cat /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.2.bed|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$5"\"; transcript_id \""$4"\";"}' > /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.gtf

#人
#二代lnc外显子bed与ref基因进行bedtools交集，获得与ref基因交集的二代lncRNA geneid
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.1.bed -b /home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.ref_gene.bed -wa -s |awk '{print $5}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA_overlap_refgene.id.txt
#二代lnc外显子bed与三代novel lnc基因进行bedtools交集，获得与三代novel lnc基因交集的二代lncRNA geneid
cat "/home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_lncRNA.gtf"|awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/;//g;s/\"//g" > /home/zhzhang/PG/2G_lnc/Homo_sapiens.3G_novel_lncRNAt.bed
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.1.bed -b /home/zhzhang/PG/2G_lnc/Homo_sapiens.3G_novel_lncRNAt.bed -wa -s |awk '{print $5}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取与Ref-gene以及3代novel-lncRNA没有交集的二代lncRNA GENE的外显子的bed注释
cat /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.1.bed|grep -v -w -Ff /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA_overlap_refgene.id.txt > /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.2.bed
#二代lncRNA GENE exon bed注释转换为gtf格式
cat /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.2.bed|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$5"\"; transcript_id \""$4"\";"}' > /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.gtf


```


```r
#兔和负鼠的二代RNA-seq数据组装的lncRNA基因集，基因组版本一致不用liftover；
#gtf转为bed
#兔
cat "/home/zhzhang/PG/2G_lnc/rabbit.lncRNA.gtf"|awk '$1~/^[0-9]/ {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14} $1~/^[XY]/ {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14}'|sed 's/\"//g;s/\;//g' > /home/zhzhang/PG/2G_lnc/Oryctolagus_cuniculus.2G_lncRNA.1.bed
#负鼠
cat "/home/zhzhang/PG/2G_lnc/opossum.lncRNA.gtf"|awk '$1~/^[0-9]/ {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14} $1~/^[XY]/ {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14}'|sed 's/\"//g;s/\;//g' > /home/zhzhang/PG/2G_lnc/Monodelphis_domestica.2G_lncRNA.1.bed



#猕猴和大鼠的二代RNA-seq数据组装的lncRNA基因集需要lift over至新版基因组
##lnc exon水平gtf转bed【最后一列是外显子顺序】，染色体号转为符合UCSC chain文件的形式
#猕猴
cat "/home/zhzhang/PG/2G_lnc/macaque.lncRNA.gtf"|awk '$1!~/^1099/ {print "chr"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14}'|sed "s/\;//g" > /home/zhzhang/PG/2G_lnc/macaque.lncRNA.bed
#大鼠
cat "/home/zhzhang/PG/2G_lnc/rat.lncRNA.gtf"|awk '$1!~/^AABR/ {print "chr"$1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7"\t"$14}'|sed "s/\;//g" > /home/zhzhang/PG/2G_lnc/rat.lncRNA.bed
##外显子bed生成转录本bed
#a是二代lnc exon bed路径，b是输出的lnc 转录本bed路径
lncexon2t <- function(a,b){
  #导入二代lnc exon bed
  lncRNA <- read.delim(a, header=FALSE)
  #lnc 转录本bed
  lncRNAt <- group_by(lncRNA,V1,V4,V6)%>%
    summarise(start=min(V2),end=max(V3))%>%
    mutate(a="transcript")%>%
    select(1,4,5,2,6,3)
  #储存
  data.table::fwrite(lncRNAt,file =b,sep = '\t',row.names = F,quote = F,col.names = F)
}
#猕猴
lncexon2t(a="~/PG/2G_lnc/macaque.lncRNA.bed",
          b="~/PG/2G_lnc/macaque.lncRNA.transcript.bed")
#大鼠
lncexon2t(a="~/PG/2G_lnc/rat.lncRNA.bed",
          b="~/PG/2G_lnc/rat.lncRNA.transcript.bed")
##liftover转换lnc转录本到新版基因组。获得无法正确转换的转录本id(为了避免出现以下情况：lift over后转录本结构出现错误（例如：liftover后一个转录本的外显子出现在不同的染色体上，外显子顺序改变，或外显子之间间距过大）)
#猕猴
liftOver /home/zhzhang/PG/2G_lnc/macaque.lncRNA.transcript.bed /home/zhzhang/software/liftover/rheMac3ToRheMac8.over.chain.gz /home/zhzhang/PG/2G_lnc/macaque.lncRNA_RheMac8.transcript.bed /home/zhzhang/PG/2G_lnc/lncTrheMac3ToRheMac8.unmap.txt
liftOver /home/zhzhang/PG/2G_lnc/macaque.lncRNA_RheMac8.transcript.bed "/home/zhzhang/software/liftover/rheMac8ToRheMac10.over.chain.gz" /home/zhzhang/PG/2G_lnc/macaque.lncRNA_RheMac10.transcript.bed /home/zhzhang/PG/2G_lnc/lncTrheMac8ToRheMac10.unmap.txt
rm /home/zhzhang/PG/2G_lnc/macaque.lncRNA_RheMac8.transcript.bed
rm /home/zhzhang/PG/2G_lnc/macaque.lncRNA_RheMac10.transcript.bed
cat /home/zhzhang/PG/2G_lnc/lncTrheMac3ToRheMac8.unmap.txt|awk '$5=="transcript" {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Macaca_mulatta.unmaplnctid.txt
cat /home/zhzhang/PG/2G_lnc/lncTrheMac8ToRheMac10.unmap.txt|awk '$5=="transcript" {print $4}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Macaca_mulatta.unmaplnctid.txt
#大鼠
liftOver /home/zhzhang/PG/2G_lnc/rat.lncRNA.transcript.bed "/home/zhzhang/software/liftover/rn5ToRn7.over.chain.gz" /home/zhzhang/PG/2G_lnc/mouse.lncRNA_Rn7.transcript.bed /home/zhzhang/PG/2G_lnc/lncTrn5ToRn7.unmap.txt
rm /home/zhzhang/PG/2G_lnc/mouse.lncRNA_Rn7.transcript.bed
cat /home/zhzhang/PG/2G_lnc/lncTrn5ToRn7.unmap.txt|awk '$5=="transcript" {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.unmaplnctid.txt
##筛选出可以转换至新版基因组的lnc转录本的外显子 bed文件，并转换至新版基因组
#猕猴
cat "/home/zhzhang/PG/2G_lnc/macaque.lncRNA.bed"|grep -w -v -Ff /home/zhzhang/PG/2G_lnc/Macaca_mulatta.unmaplnctid.txt > /home/zhzhang/PG/2G_lnc/macaque.lncRNA.ok.bed
liftOver /home/zhzhang/PG/2G_lnc/macaque.lncRNA.ok.bed "/home/zhzhang/software/liftover/rheMac3ToRheMac8.over.chain.gz" /home/zhzhang/PG/2G_lnc/macaque.lncRNA_RheMac8.bed /home/zhzhang/PG/2G_lnc/rheMac3ToRheMac8.unmap.txt
liftOver /home/zhzhang/PG/2G_lnc/macaque.lncRNA_RheMac8.bed "/home/zhzhang/software/liftover/rheMac8ToRheMac10.over.chain.gz" /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.bed /home/zhzhang/PG/2G_lnc/rheMac8ToRheMac10.unmap.txt
#NCBI染色体号转正常数字（ENSEMBLE格式）
cat /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.bed|sed "s/chr//g" > /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.chrENS.bed
#大鼠
cat "/home/zhzhang/PG/2G_lnc/rat.lncRNA.bed"|grep -w -v -Ff /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.unmaplnctid.txt > /home/zhzhang/PG/2G_lnc/rat.lncRNA.ok.bed
liftOver /home/zhzhang/PG/2G_lnc/rat.lncRNA.ok.bed /home/zhzhang/software/liftover/rn5ToRn7.over.chain.gz /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.bed /home/zhzhang/PG/2G_lnc/rn5ToRn7.unmap.txt
#UCSC染色体号转正常数字（ENSEMBLE格式）
cat /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.bed |sed "s/chr//g" > /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.chrENS.bed
##删除至少一个外显子匹配不到新版基因组核染色体的lnc转录本
#猕猴
#获取至少一个外显子匹配不到新版基因组核染色体的lnc转录本的ID（lift over中损失的exon，以及因为转换后匹配到非核染色体上损失的exon）
cat /home/zhzhang/PG/2G_lnc/rheMac3ToRheMac8.unmap.txt|awk '$4~/^Mac/ {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Macaca_mulatta.unmaplncid.txt
cat /home/zhzhang/PG/2G_lnc/rheMac8ToRheMac10.unmap.txt|awk '$4~/^Mac/ {print $4}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Macaca_mulatta.unmaplncid.txt
cat "/home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.chrENS.bed"|awk '$1~/_/ {print $4}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Macaca_mulatta.unmaplncid.txt
#删除上述转录本ID对应的全部exon
cat /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.chrENS.bed |grep -w -v -Ff /home/zhzhang/PG/2G_lnc/Macaca_mulatta.unmaplncid.txt > /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.1.1.bed
#大鼠
#获取至少一个外显子匹配不到新版基因组核染色体的lnc转录本的ID（lift over中损失的exon，以及因为转换后匹配到非核染色体上损失的exon）
cat /home/zhzhang/PG/2G_lnc/rn5ToRn7.unmap.txt|awk '$4~/^Rat/ {print $4}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.unmaplncid.txt
cat /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.chrENS.bed|awk '$1~/_/ {print $4}'|sort|uniq >> /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.unmaplncid.txt
#删除上述转录本ID对应的全部exon
cat /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.chrENS.bed |grep -w -v -Ff /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.unmaplncid.txt > /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.1.1.bed
##删除liftover后一个转录本中不同exon位于不同chr的转录本
#remove transcripts (after liftover ,exons in different chr)
library(dplyr)
#Macaca_mulatta
a <- "~/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.1.1.bed"
b <- "~/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.1.1.rm.txt"
#Rattus_norvegicus
a <- "/home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.1.1.bed"
b <- "/home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.1.1.rm.txt"
#input
twoG_lncRNA <- read.delim(a, header=FALSE)%>%
  select(1,4)%>%
  distinct(V1,V4)%>%
  group_by(V4)%>%
  summarise(num=n())%>%
  filter(num>1)%>%
  select(1)
#OUT
data.table::fwrite(twoG_lncRNA,file = b,col.names = F,row.names = F,quote = F)
#猕猴
cat /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.1.1.bed |grep -w -v -Ff /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.1.1.rm.txt > /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.1.bed
#大鼠
cat /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.1.1.bed |grep -w -v -Ff /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.1.1.rm.txt > /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.1.bed




###筛选出与Ref-gene没有交集的二代lnc GENE
#猕猴
#二代lnc外显子bed与ref基因进行bedtools交集，获得与ref基因交集的二代lncRNA geneid
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.1.bed -b /home/zhzhang/PG/ISOseqdata/Macaca_mulatta/novel_lnc/Macaca_mulatta.ref_gene.bed -wa -s |awk '{print $5}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取与Ref-gene以及3代novel-lncRNA没有交集的二代lncRNA GENE的外显子的bed注释
cat /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.1.bed|grep -v -w -Ff /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA_overlap_refgene.id.txt > /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.2.bed
#二代lncRNA GENE exon bed注释转换为gtf格式
cat /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.2.bed|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$5"\"; transcript_id \""$4"\";"}' > /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.gtf

#大鼠
#二代lnc外显子bed与ref基因进行bedtools交集，获得与ref基因交集的二代lncRNA geneid
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.1.bed -b /home/zhzhang/PG/ISOseqdata/Rattus_norvegicus/novel_lnc/Rattus_norvegicus.ref_gene.bed -wa -s |awk '{print $5}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取与Ref-gene以及3代novel-lncRNA没有交集的二代lncRNA GENE的外显子的bed注释
cat /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.1.bed|grep -v -w -Ff /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA_overlap_refgene.id.txt > /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.2.bed
#二代lncRNA GENE exon bed注释转换为gtf格式
cat /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.2.bed|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$5"\"; transcript_id \""$4"\";"}' > /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.gtf

#兔
#二代lnc外显子bed与ref基因进行bedtools交集，获得与ref基因交集的二代lncRNA geneid
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Oryctolagus_cuniculus.2G_lncRNA.1.bed -b /home/zhzhang/PG/ISOseqdata/Oryctolagus_cuniculus/novel_lnc/Oryctolagus_cuniculus.ref_gene.bed -wa -s |awk '{print $5}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Oryctolagus_cuniculus.2G_lncRNA_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取与Ref-gene以及3代novel-lncRNA没有交集的二代lncRNA GENE的外显子的bed注释
cat /home/zhzhang/PG/2G_lnc/Oryctolagus_cuniculus.2G_lncRNA.1.bed|grep -v -w -Ff /home/zhzhang/PG/2G_lnc/Oryctolagus_cuniculus.2G_lncRNA_overlap_refgene.id.txt > /home/zhzhang/PG/2G_lnc/Oryctolagus_cuniculus.2G_lncRNA.2.bed
#二代lncRNA GENE exon bed注释转换为gtf格式
cat /home/zhzhang/PG/2G_lnc/Oryctolagus_cuniculus.2G_lncRNA.2.bed|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$5"\"; transcript_id \""$4"\";"}' > /home/zhzhang/PG/2G_lnc/Oryctolagus_cuniculus.2G_lncRNA.gtf

#负鼠
#二代lnc外显子bed与ref基因进行bedtools交集，获得与ref基因交集的二代lncRNA geneid
bedtools intersect -a /home/zhzhang/PG/2G_lnc/Monodelphis_domestica.2G_lncRNA.1.bed -b /home/zhzhang/PG/ISOseqdata/Monodelphis_domestica/novel_lnc/Monodelphis_domestica.ref_gene.bed -wa -s |awk '{print $5}'|sort|uniq > /home/zhzhang/PG/2G_lnc/Monodelphis_domestica.2G_lncRNA_overlap_refgene.id.txt
#根据筛选出的geneID，反向提取与Ref-gene以及3代novel-lncRNA没有交集的二代lncRNA GENE的外显子的bed注释
cat /home/zhzhang/PG/2G_lnc/Monodelphis_domestica.2G_lncRNA.1.bed|grep -v -w -Ff /home/zhzhang/PG/2G_lnc/Monodelphis_domestica.2G_lncRNA_overlap_refgene.id.txt > /home/zhzhang/PG/2G_lnc/Monodelphis_domestica.2G_lncRNA.2.bed
#二代lncRNA GENE exon bed注释转换为gtf格式
cat /home/zhzhang/PG/2G_lnc/Monodelphis_domestica.2G_lncRNA.2.bed|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$5"\"; transcript_id \""$4"\";"}' > /home/zhzhang/PG/2G_lnc/Monodelphis_domestica.2G_lncRNA.gtf






```


