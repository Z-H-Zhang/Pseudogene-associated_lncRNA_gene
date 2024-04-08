### FIG3
### in expression
##### 1.public RNA-seq数据下载
```r
#人类数据（https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-1733/sdrf）
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
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna.chromosome.${i}.fa >> /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa
done
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa" >> /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa" >> /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa


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
STAR --runThreadN 28 --runMode genomeGenerate --genomeDir /home/zhzhang/PG/refgenome/STAR_Homo_sapiens_GRCh38/ --genomeFastaFiles /home/zhzhang/PG/refgenome/Homo_sapiens.GRCh38.dna.chr.fa --sjdbGTFfile /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/Homo_sapiens.GRCh38.108.chr.gtf --sjdbOverhang 149
STAR --runThreadN 28 --runMode genomeGenerate --genomeDir /home/zhzhang/PG/refgenome/STAR_Mus_musculus_GRCm39/ --genomeFastaFiles /home/zhzhang/PG/refgenome/Mus_musculus.GRCm39.dna.chr.fa --sjdbGTFfile /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.chr.gtf --sjdbOverhang 149

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
cat "/home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.3Gnovel_lncRNA_exon.bed"|awk '{print $1"\t3G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Gallus_gallus.GRCg7b.108.chr.rmpg.novellncRNA.gtf
#人
cp "/home/zhzhang/PG/HPG/GTF_PG/Homo_sapiens.GRCh38.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.3Gnovel_lncRNA_exon.bed"|awk '{print $1"\t3G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf
#小鼠
cp "/home/zhzhang/PG/HPG/GTF_PG/Mus_musculus.GRCm39.108.chr.rmpg.gtf" /home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.2Gnovel_lncRNA_exon.bed"|awk '{print $1"\t2G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf
cat "/home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.3Gnovel_lncRNA_exon.bed"|awk '{print $1"\t3G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf

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
cat "/home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.3Gnovel_lncRNA_exon.bed"|awk '{print $1"\t3G_lnc\texon\t"$2"\t"$3"\t.\t"$6"\t.\t""gene_id \""$4"\"; transcript_id \""$5"\";"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Danio_rerio.GRCz11.108.chr.rmpg.novellncRNA.gtf


```
##### 7.表达定量+定量结果整理为表达矩阵+标准化TPM
```r
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


```
##### 8.基因分类文件(蛋白编码/lnc的geneid-分类对照list文件)
```r
#lncrna基因部分
cp "/home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.lncRNA_class.txt" /home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.geneid_class.txt
cp "/home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.lncRNA_class.txt" /home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.geneid_class.txt

#蛋白编码基因部分
cat "/home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.GRCh38.108.chr.rmpg.novellncRNA.gtf"|awk '$3== "gene" && $17=="gene_biotype" {print $10"\t"$18} $3== "gene" && $15=="gene_biotype" {print $10"\t"$16}'|sed "s/\"//g;s/\;//g"|grep "protein_coding"|awk '{print $1"\t""Protein-coding"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Homo_sapiens.geneid_class.txt
cat "/home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.GRCm39.108.chr.rmpg.novellncRNA.gtf"|awk '$3== "gene" && $17=="gene_biotype" {print $10"\t"$18} $3== "gene" && $15=="gene_biotype" {print $10"\t"$16}'|sed "s/\"//g;s/\;//g"|grep "protein_coding"|awk '{print $1"\t""Protein-coding"}' >> /home/zhzhang/PG/RNAseqdata/newGTF/Mus_musculus.geneid_class.txt


```
##### \[9.\]两类lncRNA与蛋白编码基因的表达比例(基因在某样本中TPM>=1视为在该样本中表达)
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
      filter(type!="Interference lncRNA")%>%
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
mm_spexp <- three_gene_exppercent(a = "~/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
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
allsp_spexp$type <- factor(allsp_spexp$type,levels = c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))
#储存
data.table::fwrite(allsp_spexp,file ="/home/zhzhang/PG/RNAseq/ALLSP_3gene_spexppercent.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#画图
#人
ph <- ggplot(data = filter(allsp_spexp,sp=="Human"),aes(x=threshold,y=spexpratio*100))+
  geom_point(aes(color=type),size=2)+
  geom_line(aes(group=type,color=type))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  scale_y_continuous(limits = c(70,105))+
  labs(x = "Expression cutoff (TPM)", y ="Expressed genes (%)",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "top", legend.direction = "vertical")+
  theme(legend.position = c(0.1, 0.9))
ggsave("/home/zhzhang/PG/RNAseq/plot/Human.pglnc_sp_exp_percent.png", 
       ph,width = 4.5, height = 4.5,dpi=1200, units = "in", device='png',bg = "transparent")
#小鼠
pm <- ggplot(data = filter(allsp_spexp,sp=="Mouse"),aes(x=threshold,y=spexpratio*100))+
  geom_point(aes(color=type),size=2)+
  geom_line(aes(group=type,color=type))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  scale_y_continuous(limits = c(20,95))+
  labs(x = "Expression cutoff (TPM)", y ="Expressed genes (%)",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "top", legend.direction = "vertical")+
  theme(legend.position = c(0.1, 0.9))
ggsave("/home/zhzhang/PG/RNAseq/plot/Mouse.pglnc_sp_exp_percent.png", 
       pm,width = 4.5, height = 4.5,dpi=1200, units = "in", device='png',bg = "transparent")
#2.函数输出每个样本中三类基因的表达比例。
#（a输入基因表达TPM矩阵，b输入基因分类文件，c输入物种名）
three_gene_sampleexppercent <- function(a,b,c){
  #导入基因表达TPM矩阵
  allsample_TPM <- read.delim(a, row.names=1)
  #导入基因ID分类
  geneid_class <- read.delim(b)%>%
    filter(type!="Interference lncRNA")
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
  #计算每个样本中pglnc基因的表达比例
  pglnc_gene <- filter(gene_exp,type=="Pseudogene-derived lncRNA")%>%
    select(-1,-2)
  pglnc_gene_sampleexpratio <- data.frame(expratio=apply(pglnc_gene,2,sum)/nrow(pglnc_gene),
                                          type="Pseudogene-derived lncRNA")%>%
    rownames_to_column("sample")
  #计算每个样本中npglnc基因的表达比例
  npglnc_gene <- filter(gene_exp,type=="Non-pseudogene-derived lncRNA")%>%
    select(-1,-2)
  npglnc_gene_sampleexpratio <- data.frame(expratio=apply(npglnc_gene,2,sum)/nrow(npglnc_gene),
                                           type="Non-pseudogene-derived lncRNA")%>%
    rownames_to_column("sample")
  #合并
  gene_sampleexpratio <- rbind(coding_gene_sampleexpratio,pglnc_gene_sampleexpratio,npglnc_gene_sampleexpratio)
  #添加物种信息
  gene_sampleexpratio <- mutate(gene_sampleexpratio,sp=c)
  return(gene_sampleexpratio)
}
#计算各物种三类基因表达比例
#小鼠
mm_sampleexp <- three_gene_sampleexppercent(a = "~/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                                            b = "~/PG/RNAseq/Mus_musculus/Mus_musculus.geneid_class.txt",
                                            c = "Mouse")
#人
hs_sampleexp <- three_gene_sampleexppercent(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                                            b = "~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt",
                                            c = "Human")
#合并
allsp_sampleexp <- rbind(mm_sampleexp,hs_sampleexp)
allsp_sampleexp$type <- factor(allsp_sampleexp$type,levels = c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))
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
allsp_tissueexp <- select(allsp_tissueexp,-sex)%>%
  group_by(sp,tissue,type)%>%
  summarise(mean=mean(expratio),median=median(expratio),sd=sd(expratio))
data.table::fwrite(allsp_tissueexp,file ="/home/zhzhang/PG/RNAseq/ALLSP_3gene_tissueexppercent.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot展示每种组织三类基因的平均表达比例
pmouse <- ggplot(data = filter(allsp_tissueexp,sp=="Mouse"),aes(x=tissue,y=mean*100))+
  geom_col(width=0.6,position="dodge",aes(fill=type))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Mouse",type=="Protein-coding"),
                position=position_nudge(x=-0.2),aes(ymax=(mean+sd)*100,ymin=(mean-sd)*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Mouse",type=="Pseudogene-derived lncRNA"),
                aes(ymax=(mean+sd)*100,ymin=(mean-sd)*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Mouse",type=="Non-pseudogene-derived lncRNA"),
                position=position_nudge(x=0.2),aes(ymax=(mean+sd)*100,ymin=(mean-sd)*100))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,20,40,60),labels = c("0","20","40","60"))+
  scale_x_discrete(labels = c("Brain","Cerebellum","Colon","Heart","Ovary","Testis"))+
  labs(x = NULL, y ="Expression proportion (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")
ggsave("/home/zhzhang/PG/RNAseq/plot/Mus_musculus.tissuelnc_exp_percent.png", 
       pmouse,width = 8, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")
phuman <- ggplot(data = filter(allsp_tissueexp,sp=="Human"),aes(x=tissue,y=mean*100))+
  geom_col(width=0.6,position="dodge",aes(fill=type))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Human",type=="Protein-coding"),
                position=position_nudge(x=-0.2),aes(ymax=(mean+sd)*100,ymin=(mean-sd)*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Human",type=="Pseudogene-derived lncRNA"),
                aes(ymax=(mean+sd)*100,ymin=(mean-sd)*100))+
  geom_errorbar(width=0.1,data = filter(allsp_tissueexp,sp=="Human",type=="Non-pseudogene-derived lncRNA"),
                position=position_nudge(x=0.2),aes(ymax=(mean+sd)*100,ymin=(mean-sd)*100))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,25,50,75,100),labels = c("0","25","50","75","100"))+
  scale_x_discrete(labels = str_to_title(distinct(filter(data.frame(allsp_tissueexp),sp=="Human"),tissue)$tissue))+
  labs(x = NULL, y ="Expression proportion (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  theme(axis.text.x = element_text(angle =30))
ggsave("/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.tissuelnc_exp_percent.pdf", 
       phuman,width = 18, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#ttest验证睾丸中两类lnc表达比例均值差异性是否显著
#小鼠，p-value = 0.005173
t.test(allsp_sampleexp$expratio[45:47],allsp_sampleexp$expratio[74:76])
#人，p-value = 4.791e-06
t.test(allsp_sampleexp$expratio[462:469],allsp_sampleexp$expratio[662:669])


```
##### \[10.\]两类lncRNA与蛋白编码基因（表达的）在物种/各组织中的最大表达量
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
mm_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Mus_musculus/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Mus_musculus/Mus_musculus.spexp_geneid_class.txt",
                      c = "Mouse")
#人类
hs_TPM <- tgenemaxTPM(a = "~/PG/RNAseq/Homo_sapiens/allsample_TPM.txt",
                      b = "/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.spexp_geneid_class.txt",
                      c = "Human")
#不同物种三类基因表达量数据合并
ALLsp_TPM <- rbind(mm_TPM,hs_TPM)
ALLsp_TPM$type <- factor(ALLsp_TPM$type,levels = c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))
#统计储存
tj <- group_by(ALLsp_TPM,sp,type)%>%
  summarise(mean=mean(log2(maxTPM)),median=median(log2(maxTPM)))
data.table::fwrite(tj,file ="/home/zhzhang/PG/RNAseq/ALLSP_3gene_maxexplog2TPM.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
pmm <- ggplot(data = filter(ALLsp_TPM,sp=="Mouse"),aes(x=type,y=log2(maxTPM)))+
  geom_boxplot(fatten = 3,notch = T,width=0.5,outlier.alpha = 0,aes(fill=type))+
  geom_signif(map_signif_level=T,y_position = c(10,6.5),tip_length = 0.01,
              comparisons=list(c("Protein-coding","Pseudogene-derived lncRNA"),
                               c("Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  scale_y_continuous(breaks = c(0,2,4,6,8,10))+
  coord_cartesian(ylim = c(0, 11))+
  scale_x_discrete(labels = c("Protein-coding","Pseudogene-derived\nlncRNA","Non-pseudogene-\nderived lncRNA"))+
  labs(x = NULL, y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =30))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/RNAseq/plot/Mus_musculus.lnc_maxexp.png", 
       pmm,width = 4.5, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")
phs <- ggplot(data = filter(ALLsp_TPM,sp=="Human"),aes(x=type,y=log2(maxTPM)))+
  geom_boxplot(fatten = 3,notch = T,width=0.5,outlier.alpha = 0,aes(fill=type))+
  geom_signif(map_signif_level=T,y_position = c(11.5,8),tip_length = 0.01,
              comparisons=list(c("Protein-coding","Pseudogene-derived lncRNA"),
                               c("Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA")))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 13))+
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10,12.5),labels = c("0","2.5","5","7.5","10","12.5"))+
  scale_x_discrete(labels = c("Protein-coding","Pseudogene-derived\nlncRNA","Non-pseudogene-\nderived lncRNA"))+
  labs(x = NULL, y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(angle =30))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave("/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.lnc_maxexp.png", 
       phs,width = 4.5, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")


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
allsp_sampleexp <- rbind(hs_tissumaxexp,mm_tissumaxexp)%>%
  filter(type!="Interference lncRNA")
allsp_sampleexp$type <- factor(allsp_sampleexp$type,levels = c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))
#统计储存
allsp_sampleexp_tj <- group_by(allsp_sampleexp,sp,tissue,type)%>%
  summarise(median=median(exp))
data.table::fwrite(allsp_sampleexp_tj,file ="/home/zhzhang/PG/RNAseq/ALLSP_alltissue_3gene_maxTPM.tj.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#wilcox test p值计算(pgdlnc vs npgdlnc)
pdata <- distinct(allsp_sampleexp_tj,sp,tissue)
for (i in 1:nrow(pdata)) {
  forsp <- pdata$sp[i]
  fortissue <- pdata$tissue[i]
  pva <- wilcox.test(filter(allsp_sampleexp,sp==forsp & tissue==fortissue & type=="Pseudogene-derived lncRNA")$exp,
                     filter(allsp_sampleexp,sp==forsp & tissue==fortissue & type=="Non-pseudogene-derived lncRNA")$exp)[["p.value"]]
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
}
colnames(pdata)[c(3,4)] <- c("pvalue","anno")
data.table::fwrite(pdata,file ="/home/zhzhang/PG/RNAseq/ALLSP_alltissue_3gene_maxTPM.diffpvalue.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot每种组织，三类基因最大表达量
pmouse <- ggplot(data = filter(allsp_sampleexp,sp=="Mouse"),aes(x=tissue,y=log2(exp)))+
  geom_boxplot(fatten = 3,notch = T,width=0.5,outlier.alpha = 0,aes(fill=type))+
  geom_signif(annotations = filter(pdata,sp=="Mouse")$anno,y_position = c(11),
              xmin = c(1:6),xmax = c(1.166:6.166),tip_length = 0)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
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
  geom_boxplot(fatten = 3,notch = T,width=0.5,outlier.alpha = 0,aes(fill=type))+
  geom_signif(annotations = filter(pdata,sp=="Human")$anno,y_position = c(11),
              xmin = c(1:32),xmax = c(1.166:32.166),tip_length = 0)+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  coord_cartesian(ylim = c(0, 11))+
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10,12.5),labels = c("0","2.5","5","7.5","10","12.5"))+
  scale_x_discrete(labels = str_to_title(distinct(filter(data.frame(allsp_sampleexp),sp=="Human"),tissue)$tissue))+
  labs(x = NULL, y =expression("M"*"a"*"x"*"i"*"m"*"a"*"l"~"e"*"x"*"p"*"r"*"e"*"s"*"s"*"i"*"o"*"n"~"("*"l"*"o"*"g"[2]*"T"*"P"*"M"*")"),
       fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  theme(axis.text.x = element_text(hjust=1,vjust = 1))+
  theme(axis.text.x = element_text(angle =30))
ggsave("/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.tissue_3gene_maxexp.pdf", 
       phuman,width = 18, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```



##### \[12.\]lnc组织特异性和组织分布情况
```r
#函数输出每个基因在每种组织中是否表达(在某组织至少一个样本中表达即为在该组织中表达)
#（a输入基因表达TPM矩阵，b输入基因分类文件，c输出路径）
three_gene_sampleexppercent <- function(a,b,c){
  #导入基因表达TPM矩阵
  allsample_TPM <- read.delim(a, row.names=1)
  #导入基因ID分类
  geneid_class <- read.delim(b)%>%
    filter(type!="Interference lncRNA")
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
hs_gene_type_tissueexp$type <- factor(hs_gene_type_tissueexp$type,
                                      levels = c("Non-pseudogene-derived lncRNA","Pseudogene-derived lncRNA","Protein-coding"))
#储存
  data.table::fwrite(hs_gene_type_tissueexp,file ="/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.gene_tissuenum.txt",sep = '\t',row.names = F,quote = F,col.names = T)
#plot
library(gghalves)
wp1 <- signif(wilcox.test(filter(hs_gene_type_tissueexp,type=="Pseudogene-derived lncRNA")$tnum,
                          filter(hs_gene_type_tissueexp,type=="Non-pseudogene-derived lncRNA")$tnum)[["p.value"]],
              2)
wp2 <- signif(wilcox.test(filter(hs_gene_type_tissueexp,type=="Pseudogene-derived lncRNA")$tnum,
                          filter(hs_gene_type_tissueexp,type=="Protein-coding")$tnum)[["p.value"]],
              2)
hs <- ggplot(data = hs_gene_type_tissueexp,aes(x=type,y=tnum))+
  geom_half_violin(width=1,side = "r",position = position_nudge(x=0.2),aes(color=type,fill=type))+
  geom_point(data = filter(hs_gene_type_tissueexp,type=="Protein-coding"),stroke=0,
             alpha=0.8/(18264/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(hs_gene_type_tissueexp,type=="Non-pseudogene-derived lncRNA"),stroke=0,
             alpha=0.8/(25582/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(hs_gene_type_tissueexp,type=="Pseudogene-derived lncRNA"),stroke=0,
             alpha=0.8,size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_boxplot(fatten=3,width=0.05,position = position_nudge(x=0.2),fill="white",outlier.alpha = 0)+
  scale_x_discrete(labels = c("Non-pseudogene-\nderived lncRNA","Pseudogene-\nderived lncRNA","Protein-coding"))+
  coord_flip()+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  labs(x =NULL, y ="Number of tissues",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  geom_signif(annotations=c(wp1,wp2),y_position=33,tip_length = 0,
              xmin = c(1.2,2.2),xmax = c(1.8,2.8))
ggsave("/home/zhzhang/PG/RNAseq/plot/Homo_sapiens.lnc_tissuenum.png", 
       hs,width = 7, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")
#对比三类基因表达的广泛程度(小鼠)
mm_gene_type_tissueexp <- read.delim("~/PG/RNAseq/Mus_musculus/Mus_musculus.gene_type_tissueexp.txt")%>%
  group_by(geneid,type)%>%
  summarise(tnum=sum(exp))
mm_gene_type_tissueexp$type <- factor(mm_gene_type_tissueexp$type,
                                      levels = c("Non-pseudogene-derived lncRNA","Pseudogene-derived lncRNA","Protein-coding"))
#plot
library(gghalves)
wp1 <- signif(wilcox.test(filter(mm_gene_type_tissueexp,type=="Pseudogene-derived lncRNA")$tnum,
                          filter(mm_gene_type_tissueexp,type=="Non-pseudogene-derived lncRNA")$tnum)[["p.value"]],
              2)
wp2 <- signif(wilcox.test(filter(mm_gene_type_tissueexp,type=="Pseudogene-derived lncRNA")$tnum,
                          filter(mm_gene_type_tissueexp,type=="Protein-coding")$tnum)[["p.value"]],
              2)
mm <- ggplot(data = mm_gene_type_tissueexp,aes(x=type,y=tnum))+
  geom_half_violin(width=1,side = "r",position = position_nudge(x=0.2),aes(color=type,fill=type))+
  geom_point(data = filter(mm_gene_type_tissueexp,type=="Protein-coding"),stroke=0,
             alpha=0.8/(18264/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(mm_gene_type_tissueexp,type=="Non-pseudogene-derived lncRNA"),stroke=0,
             alpha=0.8/(25582/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(mm_gene_type_tissueexp,type=="Pseudogene-derived lncRNA"),stroke=0,
             alpha=0.8,size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_boxplot(fatten=3,width=0.05,position = position_nudge(x=0.2),fill="white",outlier.alpha = 0)+
  scale_x_discrete(labels = c("Non-pseudogene-\nderived lncRNA","Pseudogene-\nderived lncRNA","Protein-coding"))+
  coord_flip()+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                    limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#FAA465","#8491B4"),
                     limits=c("Protein-coding","Pseudogene-derived lncRNA","Non-pseudogene-derived lncRNA"))+
  theme_half_open()+
  labs(x =NULL, y ="Number of tissues",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  geom_signif(annotations=c(wp1,wp2),y_position=7,tip_length = 0,
              xmin = c(1.2,2.2),xmax = c(1.8,2.8))
ggsave("/home/zhzhang/PG/RNAseq/plot/Mus_musculus.lnc_tissuenum.png", 
       mm,width = 7, height = 5,dpi=1200, units = "in", device='png',bg = "transparent")



```



