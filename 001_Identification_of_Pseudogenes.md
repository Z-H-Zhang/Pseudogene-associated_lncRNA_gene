# Pg-lnc
# PG
### 一.鉴定假基因
##### 1.rpeatmasked(rm)基因组序列下载和unmasked染色体序列下载
```r
#斑马鱼
#下载rmdna
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/dna
wget https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna_rm.primary_assembly.fa.gz
for i in $(seq 1 25)
do
wget https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome.${i}.fa.gz
done
#解压
gzip -d /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/dna/*
#rmdna去除非必要部分
cat -n "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/dna/Danio_rerio.GRCz11.dna_rm.primary_assembly.fa"|grep ">"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/dna/Danio_rerio.GRCz11.dna_rm.primary_assembly.fa"|sed -n "1,22418402p" > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/dna/Danio_rerio.GRCz11.dna_rm.mainchr.fa"





#鸡
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/
wget https://ftp.ensembl.org/pub/release-108/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_rm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.W.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.Z.fa.gz
for i in $(seq 1 39)
do
wget https://ftp.ensembl.org/pub/release-108/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.${i}.fa.gz
done
#解压
gzip -d /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/*
#rmdna去除非必要部分
cat -n "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_rm.toplevel.fa"|grep ">"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_rm.toplevel.fa"|sed -n "1,17352103p" > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_rm.mainchr.fa"





#小鼠
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/
wget https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.X.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.Y.fa.gz
for i in $(seq 1 19)
do
wget https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.${i}.fa.gz
done
#解压
gzip -d /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/*
#rmdna去除非必要部分
cat -n "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa"|grep ">"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa"|sed -n "1,41041397p" > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.mainchr.fa"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa"|sed -n "41041671,45390549p" >> "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.mainchr.fa"



#人
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
for i in $(seq 1 22)
do
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz
done
#解压
gzip -d /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/*
#rmdna去除非必要部分
cat -n "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"|grep ">"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"|sed -n "1,47916726p" > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna_rm.mainchr.fa"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"|sed -n "47917005,51471479p" >> "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna_rm.mainchr.fa"





#猕猴
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/
wget https://ftp.ensembl.org/pub/release-108/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna_rm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.X.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.Y.fa.gz
for i in $(seq 1 20)
do
wget https://ftp.ensembl.org/pub/release-108/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.${i}.fa.gz
done
#解压
gzip -d /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/*
#rmdna去除非必要部分
cat -n "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna_rm.toplevel.fa"|grep ">"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna_rm.toplevel.fa"|sed -n "1,47566307p" > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna_rm.mainchr.fa"



#大鼠
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/
wget https://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_rm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.X.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.Y.fa.gz
for i in $(seq 1 20)
do
wget https://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.${i}.fa.gz
done
#解压
gzip -d /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/*
#rmdna去除非必要部分
cat -n "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_rm.toplevel.fa"|grep ">"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_rm.toplevel.fa"|sed -n "1,43891257p" > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_rm.mainchr.fa"



#兔
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/dna/
wget https://ftp.ensembl.org/pub/release-108/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna_rm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.X.fa.gz
for i in $(seq 1 21)
do
wget https://ftp.ensembl.org/pub/release-108/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.${i}.fa.gz
done
#解压
gzip -d /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/dna/*
#rmdna去除非必要部分
cat -n "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna_rm.toplevel.fa"|grep ">"
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna_rm.toplevel.fa|sed -n "1,37462566p" > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna_rm.mainchr.fa"



#负鼠
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/dna/
wget https://ftp.ensembl.org/pub/release-108/fasta/monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna_rm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna.primary_assembly.X.fa.gz
for i in $(seq 1 8)
do
wget https://ftp.ensembl.org/pub/release-108/fasta/monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna.primary_assembly.${i}.fa.gz
done
#解压
gzip -d /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/dna/*
#rmdna去除非必要部分
cat -n "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna_rm.toplevel.fa"|grep ">"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna_rm.toplevel.fa"|sed -n "1,58372899p" > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna_rm.mainchr.fa"



```
##### 2.蛋白质序列下载
```r
#斑马鱼
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/pep/
wget https://ftp.ensembl.org/pub/release-108/fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz
gzip -d Danio_rerio.GRCz11.pep.all.fa.gz
#查看蛋白所属的基因类型（能编码蛋白的基因的类型）
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa|awk '{print $6}'|sort|uniq
###（除polymorphic_pseudogene外，其余五类基因的exon坐标需要加入外显子坐标文件）
gene_biotype:IG_C_gene
gene_biotype:polymorphic_pseudogene
gene_biotype:protein_coding
gene_biotype:TR_D_gene
gene_biotype:TR_J_gene
gene_biotype:TR_V_gene




#鸡
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/pep/
wget https://ftp.ensembl.org/pub/release-108/fasta/gallus_gallus/pep/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa.gz
gzip -d Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa.gz
#查看蛋白所属的基因类型（能编码蛋白的基因的类型）
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/pep/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa"|awk '{print $6}'|sort|uniq
###（protein_coding基因的exon坐标需要加入外显子坐标文件）
gene_biotype:protein_coding



#小鼠
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/pep/
wget https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz
gzip -d Mus_musculus.GRCm39.pep.all.fa.gz
#查看蛋白所属的基因类型（能编码蛋白的基因的类型）
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa"|awk '{print $6}'|sort|uniq
###（十类基因的exon坐标需要加入外显子坐标文件）
gene_biotype:IG_C_gene
gene_biotype:IG_D_gene
gene_biotype:IG_J_gene
gene_biotype:IG_LV_gene
gene_biotype:IG_V_gene
gene_biotype:protein_coding
gene_biotype:TR_C_gene
gene_biotype:TR_D_gene
gene_biotype:TR_J_gene
gene_biotype:TR_V_gene




#人
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/pep/
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gzip -d Homo_sapiens.GRCh38.pep.all.fa.gz
#查看蛋白所属的基因类型（能编码蛋白的基因的类型）
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/pep/Homo_sapiens.GRCh38.pep.all.fa"|awk '{print $6}'|sort|uniq
###（九类基因的exon坐标需要加入外显子坐标文件）
gene_biotype:IG_C_gene
gene_biotype:IG_D_gene
gene_biotype:IG_J_gene
gene_biotype:IG_V_gene
gene_biotype:protein_coding
gene_biotype:TR_C_gene
gene_biotype:TR_D_gene
gene_biotype:TR_J_gene
gene_biotype:TR_V_gene




#猕猴
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/pep/
wget https://ftp.ensembl.org/pub/release-108/fasta/macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa.gz
gzip -d Macaca_mulatta.Mmul_10.pep.all.fa.gz
#查看蛋白所属的基因类型（能编码蛋白的基因的类型）
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa"|awk '{print $6}'|sort|uniq
###（除polymorphic_pseudogene外，其余五类基因的exon坐标需要加入外显子坐标文件）
gene_biotype:IG_C_gene
gene_biotype:IG_D_gene
gene_biotype:IG_V_gene
gene_biotype:protein_coding
gene_biotype:TR_C_gene
gene_biotype:TR_J_gene
gene_biotype:TR_V_gene


#大鼠
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/pep/
wget https://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/pep/Rattus_norvegicus.mRatBN7.2.pep.all.fa.gz
gzip -d Rattus_norvegicus.mRatBN7.2.pep.all.fa.gz
#查看蛋白所属的基因类型（能编码蛋白的基因的类型）
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/pep/Rattus_norvegicus.mRatBN7.2.pep.all.fa"|awk '{print $6}'|sort|uniq
###（除polymorphic_pseudogene外，其余五类基因的exon坐标需要加入外显子坐标文件）
gene_biotype:IG_V_gene
gene_biotype:protein_coding
gene_biotype:TR_C_gene
gene_biotype:TR_J_gene
gene_biotype:TR_V_gene


#兔
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/pep/
wget https://ftp.ensembl.org/pub/release-108/fasta/oryctolagus_cuniculus/pep/Oryctolagus_cuniculus.OryCun2.0.pep.all.fa.gz
gzip -d Oryctolagus_cuniculus.OryCun2.0.pep.all.fa.gz
#查看蛋白所属的基因类型（能编码蛋白的基因的类型）
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/pep/Oryctolagus_cuniculus.OryCun2.0.pep.all.fa"|awk '{print $6}'|sort|uniq
###（除polymorphic_pseudogene外，其余五类基因的exon坐标需要加入外显子坐标文件）
gene_biotype:IG_C_gene
gene_biotype:IG_V_gene
gene_biotype:protein_coding
gene_biotype:TR_C_gene
gene_biotype:TR_J_gene
gene_biotype:TR_V_gene

#负鼠
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/pep/
wget https://ftp.ensembl.org/pub/release-108/fasta/monodelphis_domestica/pep/Monodelphis_domestica.ASM229v1.pep.all.fa.gz
gzip -d Monodelphis_domestica.ASM229v1.pep.all.fa.gz
#查看蛋白所属的基因类型（能编码蛋白的基因的类型）
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/pep/Monodelphis_domestica.ASM229v1.pep.all.fa"|awk '{print $6}'|sort|uniq
###（除polymorphic_pseudogene外，其余五类基因的exon坐标需要加入外显子坐标文件）
gene_biotype:IG_V_gene
gene_biotype:protein_coding
gene_biotype:TR_V_gene

```
##### 3.基因组注释下载及蛋白编码基因exon坐标提取
```r
#斑马鱼【GRCz11
#下载
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/
wget https://ftp.ensembl.org/pub/release-108/gtf/danio_rerio/Danio_rerio.GRCz11.108.chr.gtf.gz
gzip -d "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/Danio_rerio.GRCz11.108.chr.gtf.gz"
#提取全部染色体上的蛋白编码基因外显子坐标
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/Danio_rerio.GRCz11.108.chr.gtf"|grep -E "protein_coding|IG_C_gene|TR_D_gene|TR_J_gene|TR_V_gene"|sort -n -k 1 -k 4 -k 5|awk '$3=="exon" {print NR"\t"$1"\t"$4"\t"$5}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/Danio_rerio.GRCz11.108.exLocs.txt
#从全部染色体外显子坐标文件中提取出每个单独的染色体外显子坐标
for i in $(seq 1 25)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/Danio_rerio.GRCz11.108.exLocs.txt|awk -v a=${i} '$2==a {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/chr${i}_exLocs
done




#鸡【GRCg7b
#下载
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/
wget https://ftp.ensembl.org/pub/release-108/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.chr.gtf.gz
gzip -d "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.chr.gtf.gz"
#提取全部染色体上的蛋白编码基因外显子坐标
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.chr.gtf"|grep protein_coding|sort -n -k 1 -k 4 -k 5|awk '$3=="exon" {print NR"\t"$1"\t"$4"\t"$5}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.exLocs.txt
#从全部染色体外显子坐标文件中提取出每个单独的染色体外显子坐标
for i in $(seq 1 39)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.exLocs.txt|awk -v a=${i} '$2==a {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/chr${i}_exLocs
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.exLocs.txt|awk '$2=="W" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/chrW_exLocs
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.exLocs.txt|awk '$2=="Z" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/chrZ_exLocs



#小鼠【GRCm39
#下载
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/
wget https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.chr.gtf.gz
gzip -d "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.chr.gtf.gz"
#提取全部染色体上的蛋白编码基因外显子坐标
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.chr.gtf"|grep -E "protein_coding|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene"|sort -n -k 1 -k 4 -k 5|awk '$3=="exon" {print NR"\t"$1"\t"$4"\t"$5}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.exLocs.txt
#从全部染色体外显子坐标文件中提取出每个单独的染色体外显子坐标
for i in $(seq 1 19)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.exLocs.txt|awk -v a=${i} '$2==a {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/chr${i}_exLocs
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.exLocs.txt|awk '$2=="X" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/chrX_exLocs
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/Mus_musculus.GRCm39.108.exLocs.txt|awk '$2=="Y" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/chrY_exLocs





#人【GRCh38
#下载
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz
gzip -d "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/Homo_sapiens.GRCh38.108.chr.gtf.gz"
#提取全部染色体上的蛋白编码基因外显子坐标
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/Homo_sapiens.GRCh38.108.chr.gtf"|grep -E "protein_coding|IG_C_gene|IG_D_gene|IG_J_gene|IG_V_gene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene"|sort -n -k 1 -k 4 -k 5|awk '$3=="exon" {print NR"\t"$1"\t"$4"\t"$5}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/Homo_sapiens.GRCh38.108.exLocs.txt
#从全部染色体外显子坐标文件中提取出每个单独的染色体外显子坐标
for i in $(seq 1 22)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/Homo_sapiens.GRCh38.108.exLocs.txt|awk -v a=${i} '$2==a {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/chr${i}_exLocs
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/Homo_sapiens.GRCh38.108.exLocs.txt|awk '$2=="X" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/chrX_exLocs
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/Homo_sapiens.GRCh38.108.exLocs.txt|awk '$2=="Y" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/chrY_exLocs





#猕猴
#下载
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/
wget https://ftp.ensembl.org/pub/release-108/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.108.chr.gtf.gz
gzip -d "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/Macaca_mulatta.Mmul_10.108.chr.gtf.gz"
#提取全部染色体上的蛋白编码基因外显子坐标
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/Macaca_mulatta.Mmul_10.108.chr.gtf"|grep -E "protein_coding|IG_C_gene|IG_D_gene|IG_V_gene|TR_C_gene|TR_J_gene|TR_V_gene"|sort -n -k 1 -k 4 -k 5|awk '$3=="exon" {print NR"\t"$1"\t"$4"\t"$5}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/Macaca_mulatta.Mmul_10.108.exLocs.txt
#从全部染色体外显子坐标文件中提取出每个单独的染色体外显子坐标
for i in $(seq 1 20)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/Macaca_mulatta.Mmul_10.108.exLocs.txt|awk -v a=${i} '$2==a {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/chr${i}_exLocs
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/Macaca_mulatta.Mmul_10.108.exLocs.txt|awk '$2=="X" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/chrX_exLocs
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/Macaca_mulatta.Mmul_10.108.exLocs.txt|awk '$2=="Y" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/chrY_exLocs



#大鼠
#下载
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/
wget https://ftp.ensembl.org/pub/release-108/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.108.chr.gtf.gz
gzip -d "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/Rattus_norvegicus.mRatBN7.2.108.chr.gtf.gz"
#提取全部染色体上的蛋白编码基因外显子坐标
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/Rattus_norvegicus.mRatBN7.2.108.chr.gtf"|grep -E "protein_coding|IG_V_gene|TR_C_gene|TR_J_gene|TR_V_gene"|sort -n -k 1 -k 4 -k 5|awk '$3=="exon" {print NR"\t"$1"\t"$4"\t"$5}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/Rattus_norvegicus.mRatBN7.2.108.exLocs.txt
#从全部染色体外显子坐标文件中提取出每个单独的染色体外显子坐标
for i in $(seq 1 20)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/Rattus_norvegicus.mRatBN7.2.108.exLocs.txt|awk -v a=${i} '$2==a {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/chr${i}_exLocs
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/Rattus_norvegicus.mRatBN7.2.108.exLocs.txt|awk '$2=="X" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/chrX_exLocs
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/Rattus_norvegicus.mRatBN7.2.108.exLocs.txt|awk '$2=="Y" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/chrY_exLocs



#兔
#下载
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/
wget https://ftp.ensembl.org/pub/release-108/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.108.chr.gtf.gz
gzip -d "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/Oryctolagus_cuniculus.OryCun2.0.108.chr.gtf.gz"
#提取全部染色体上的蛋白编码基因外显子坐标
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/Oryctolagus_cuniculus.OryCun2.0.108.chr.gtf"|grep -E "protein_coding|IG_V_gene|IG_C_gene|TR_C_gene|TR_J_gene|TR_V_gene"|sort -n -k 1 -k 4 -k 5|awk '$3=="exon" {print NR"\t"$1"\t"$4"\t"$5}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/Oryctolagus_cuniculus.OryCun2.0.108.exLocs.txt
#从全部染色体外显子坐标文件中提取出每个单独的染色体外显子坐标
for i in $(seq 1 21)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/Oryctolagus_cuniculus.OryCun2.0.108.exLocs.txt|awk -v a=${i} '$2==a {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/chr${i}_exLocs
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/Oryctolagus_cuniculus.OryCun2.0.108.exLocs.txt|awk '$2=="X" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/chrX_exLocs



#负鼠
#下载
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/
wget https://ftp.ensembl.org/pub/release-108/gtf/monodelphis_domestica/Monodelphis_domestica.ASM229v1.108.chr.gtf.gz
gzip -d Monodelphis_domestica.ASM229v1.108.chr.gtf.gz
#提取全部染色体上的蛋白编码基因外显子坐标
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/Monodelphis_domestica.ASM229v1.108.chr.gtf"|grep -E "protein_coding|IG_V_gene|TR_V_gene"|sort -n -k 1 -k 4 -k 5|awk '$3=="exon" {print NR"\t"$1"\t"$4"\t"$5}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/Monodelphis_domestica.ASM229v1.108.exLocs.txt
#从全部染色体外显子坐标文件中提取出每个单独的染色体外显子坐标
for i in $(seq 1 8)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/Monodelphis_domestica.ASM229v1.108.exLocs.txt|awk -v a=${i} '$2==a {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/chr${i}_exLocs
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/Monodelphis_domestica.ASM229v1.108.exLocs.txt|awk '$2=="X" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/chrX_exLocs


```


##### 4.假基因候选鉴定
```r
#斑马鱼
source /home/zhzhang/miniconda3/bin/activate pseudopipe
cd /home/zhzhang/software/pgenes/pseudopipe/bin
./pseudopipe_pro.sh /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/dna/Danio_rerio.GRCz11.dna_rm.mainchr.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome.%s.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Danio_rerio/mysql/chr%s_exLocs 0



#鸡
source /home/zhzhang/miniconda3/bin/activate pseudopipe
cd /home/zhzhang/software/pgenes/pseudopipe/bin
./pseudopipe_pro.sh /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_rm.mainchr.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.primary_assembly.%s.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/pep/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.pep.all.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Gallus_gallus/mysql/chr%s_exLocs 0


#小鼠
source /home/zhzhang/miniconda3/bin/activate pseudopipe
cd /home/zhzhang/software/pgenes/pseudopipe/bin
./pseudopipe_pro.sh /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.mainchr.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.%s.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Mus_musculus/mysql/chr%s_exLocs 0


#人
source /home/zhzhang/miniconda3/bin/activate pseudopipe
cd /home/zhzhang/software/pgenes/pseudopipe/bin
./pseudopipe_pro.sh /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna_rm.mainchr.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/dna/Homo_sapiens.GRCh38.dna.chromosome.%s.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/pep/Homo_sapiens.GRCh38.pep.all.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens_n/mysql/chr%s_exLocs 0



#猕猴
source /home/zhzhang/miniconda3/bin/activate pseudopipe
cd /home/zhzhang/software/pgenes/pseudopipe/bin
./pseudopipe_pro.sh /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna_rm.mainchr.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.primary_assembly.%s.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Macaca_mulatta/mysql/chr%s_exLocs 0
#大鼠
source /home/zhzhang/miniconda3/bin/activate pseudopipe
cd /home/zhzhang/software/pgenes/pseudopipe/bin
./pseudopipe_pro.sh /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_rm.mainchr.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.%s.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/pep/Rattus_norvegicus.mRatBN7.2.pep.all.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Rattus_norvegicus/mysql/chr%s_exLocs 0
#兔
source /home/zhzhang/miniconda3/bin/activate pseudopipe
cd /home/zhzhang/software/pgenes/pseudopipe/bin
./pseudopipe_pro.sh /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna_rm.mainchr.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.chromosome.%s.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/pep/Oryctolagus_cuniculus.OryCun2.0.pep.all.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Oryctolagus_cuniculus/mysql/chr%s_exLocs 0
#负鼠
source /home/zhzhang/miniconda3/bin/activate pseudopipe
cd /home/zhzhang/software/pgenes/pseudopipe/bin
./pseudopipe_pro.sh /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna_rm.mainchr.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/dna/Monodelphis_domestica.ASM229v1.dna.primary_assembly.%s.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/pep/Monodelphis_domestica.ASM229v1.pep.all.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Monodelphis_domestica/mysql/chr%s_exLocs 0


```
##### 5.假基因外显子信息合并
```r
#斑马鱼
#负链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/minus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/minus/pexons/ALL.pgexon.gff"
for i in $(seq 2 25)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/minus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/minus/pexons/ALL.pgexon.gff
done
#正链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/plus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/plus/pexons/ALL.pgexon.gff"
for i in $(seq 2 25)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/plus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/plus/pexons/ALL.pgexon.gff
done
#正负链外显子合并
cp /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/minus/pexons/ALL.pgexon.gff /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_pgexon.txt
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/plus/pexons/ALL.pgexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_pgexon.txt



#青鳉
#负链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/minus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/minus/pexons/ALL.pgexon.gff"
for i in $(seq 2 24)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/minus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/minus/pexons/ALL.pgexon.gff
done
#正链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/plus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/plus/pexons/ALL.pgexon.gff"
for i in $(seq 2 24)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/plus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/plus/pexons/ALL.pgexon.gff
done
#正负链外显子合并
cp /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/minus/pexons/ALL.pgexon.gff /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/Oryzias_latipes_pgexon.txt
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/plus/pexons/ALL.pgexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/Oryzias_latipes_pgexon.txt








#鸡
#负链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/minus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/minus/pexons/ALL.pgexon.gff"
for i in $(seq 2 39)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/minus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/minus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/minus/pexons/W.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/minus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/minus/pexons/Z.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/minus/pexons/ALL.pgexon.gff
#正链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/plus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/plus/pexons/ALL.pgexon.gff"
for i in $(seq 2 39)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/plus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/plus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/plus/pexons/W.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/plus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/plus/pexons/Z.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/plus/pexons/ALL.pgexon.gff
#正负链外显子合并
cp /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/minus/pexons/ALL.pgexon.gff /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_pgexon.txt
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/plus/pexons/ALL.pgexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_pgexon.txt



#小鼠
#负链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/minus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/minus/pexons/ALL.pgexon.gff"
for i in $(seq 2 19)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/minus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/minus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/minus/pexons/X.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/minus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/minus/pexons/Y.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/minus/pexons/ALL.pgexon.gff
#正链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/plus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/plus/pexons/ALL.pgexon.gff"
for i in $(seq 2 19)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/plus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/plus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/plus/pexons/X.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/plus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/plus/pexons/Y.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/plus/pexons/ALL.pgexon.gff
#正负链外显子合并
cp /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/minus/pexons/ALL.pgexon.gff /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_pgexon.txt
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/plus/pexons/ALL.pgexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_pgexon.txt



#人
#负链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/minus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/minus/pexons/ALL.pgexon.gff"
for i in $(seq 2 22)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/minus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/minus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/minus/pexons/X.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/minus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/minus/pexons/Y.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/minus/pexons/ALL.pgexon.gff
#正链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/plus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/plus/pexons/ALL.pgexon.gff"
for i in $(seq 2 22)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/plus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/plus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/plus/pexons/X.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/plus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/plus/pexons/Y.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/plus/pexons/ALL.pgexon.gff
#正负链外显子合并
cp /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/minus/pexons/ALL.pgexon.gff /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_pgexon.txt
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/plus/pexons/ALL.pgexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_pgexon.txt


```
##### 6.筛选出高质量假基因，并生成HPG\_ID
```r
#假基因文件传至实验室服务器/home/zhzhang/PG/PP_wd
library(tidyverse)
#pg输入假基因文件，hpg输出高质量假基因文件,sp物种名,pid2gid输入蛋白ID和基因ID对照关系文件
filterpg <- function(pg,hpg,sp,pid2gid){
  #导入假基因
  pgenes <- read.delim(pg)
  #筛选(类型，比对覆盖率和比对长度aa，e值，序列一致性)
  h_pg <- mutate(pgenes,len=(end-start)/3)%>%
    filter(frac>=0.1 & len>=30 & expect<=0.00001 & ident>=0.2)%>%
    arrange(chr,start,end)
  #增加假基因ID
  hpgnum <- nrow(h_pg)
  h_pg <- mutate(h_pg,geneid=paste(sp,"PG",c(1:hpgnum),sep = "_"))
  #根据假基因母本蛋白ID，以及蛋白ID-基因ID对照关系，增加假基因的母本基因ID
  PID2GID <- read.delim(pid2gid,header=FALSE)
  colnames(PID2GID) <- c("query","parentgene_id")
  h_pg <- left_join(h_pg,PID2GID,by="query")
  #储存高质量假基因
  data.table::fwrite(h_pg,file = hpg,sep = '\t',row.names = F,quote = F,col.names = T)
  return(h_pg)
}


#斑马鱼
filterpg("/home/zhzhang/PG/PP_wd/Danio_rerio_pgenes.txt",
         "/home/zhzhang/PG/pg_message/Danio_rerio_hpg.txt",
         "DAR",
         "/home/zhzhang/PG/PP_wd/Danio_rerio.GRCz11.PID2GID.txt")

#青鳉
filterpg("/home/zhzhang/PG/PP_wd/Oryzias_latipes_pgenes.txt",
         "/home/zhzhang/PG/pg_message/Oryzias_latipes_hpg.txt",
         "ORL",
         "/home/zhzhang/PG/PP_wd/Oryzias_latipes.ASM223467v1.PID2GID.txt")

#鸡
filterpg("/home/zhzhang/PG/PP_wd/Gallus_gallus_pgenes.txt",
         "/home/zhzhang/PG/pg_message/Gallus_gallus_hpg.txt",
         "GAL",
         "/home/zhzhang/PG/PP_wd/Gallus_gallus.GRCg7b.PID2GID.txt")

#小鼠
filterpg("/home/zhzhang/PG/PP_wd/Mus_musculus_pgenes.txt",
         "/home/zhzhang/PG/pg_message/Mus_musculus_hpg.txt",
         "MUS",
         "/home/zhzhang/PG/PP_wd/Mus_musculus.GRCm39.PID2GID.txt")

#人
filterpg("/home/zhzhang/PG/PP_wd/Homo_sapiens_pgenes.txt",
         "/home/zhzhang/PG/pg_message/Homo_sapiens_hpg.txt",
         "HOS",
         "/home/zhzhang/PG/PP_wd/Homo_sapiens.GRCh38.PID2GID.txt")

#猕猴
filterpg("/home/zhzhang/PG/PP_wd/Macaca_mulatta_pgenes.txt",
         "/home/zhzhang/PG/pg_message/Macaca_mulatta_hpg.txt",
         "MMU",
         "/home/zhzhang/PG/PP_wd/Macaca_mulatta.Mmul_10.PID2GID.txt")

#大鼠
filterpg("/home/zhzhang/PG/PP_wd/Rattus_norvegicus_pgenes.txt",
         "/home/zhzhang/PG/pg_message/Rattus_norvegicus_hpg.txt",
         "RNO",
         "/home/zhzhang/PG/PP_wd/Rattus_norvegicus.mRatBN7.2.PID2GID.txt")

#兔子
filterpg("/home/zhzhang/PG/PP_wd/Oryctolagus_cuniculus_pgenes.txt",
         "/home/zhzhang/PG/pg_message/Oryctolagus_cuniculus_hpg.txt",
         "OCU",
         "/home/zhzhang/PG/PP_wd/Oryctolagus_cuniculus.OryCun2.0.PID2GID.txt")

#负鼠
filterpg("/home/zhzhang/PG/PP_wd/Monodelphis_domestica_pgenes.txt",
         "/home/zhzhang/PG/pg_message/Monodelphis_domestica_hpg.txt",
         "MOD",
         "/home/zhzhang/PG/PP_wd/Monodelphis_domestica.ASM229v1.PID2GID.txt")



#DUP和PSSD替换为全称
#斑马鱼
sed -i "s/DUP/Duplicated/g;s/PSSD/Processed/g;s/FRAG/Fragment/g" /home/zhzhang/PG/pg_message/Danio_rerio_hpg.txt
#青鳉
sed -i "s/DUP/Duplicated/g;s/PSSD/Processed/g;s/FRAG/Fragment/g" /home/zhzhang/PG/pg_message/Oryzias_latipes_hpg.txt
#鸡
sed -i "s/DUP/Duplicated/g;s/PSSD/Processed/g;s/FRAG/Fragment/g" /home/zhzhang/PG/pg_message/Gallus_gallus_hpg.new.txt
#小鼠
sed -i "s/DUP/Duplicated/g;s/PSSD/Processed/g;s/FRAG/Fragment/g" /home/zhzhang/PG/pg_message/Mus_musculus_hpg.new.txt
#人
sed -i "s/DUP/Duplicated/g;s/PSSD/Processed/g;s/FRAG/Fragment/g" /home/zhzhang/PG/pg_message/Homo_sapiens_hpg.new.txt
#猕猴
sed -i "s/DUP/Duplicated/g;s/PSSD/Processed/g;s/FRAG/Fragment/g" "/home/zhzhang/PG/pg_message/Macaca_mulatta_hpg.txt"
#大鼠
sed -i "s/DUP/Duplicated/g;s/PSSD/Processed/g;s/FRAG/Fragment/g" "/home/zhzhang/PG/pg_message/Rattus_norvegicus_hpg.txt"
#兔子
sed -i "s/DUP/Duplicated/g;s/PSSD/Processed/g;s/FRAG/Fragment/g" "/home/zhzhang/PG/pg_message/Oryctolagus_cuniculus_hpg.txt"
#负鼠
sed -i "s/DUP/Duplicated/g;s/PSSD/Processed/g;s/FRAG/Fragment/g" "/home/zhzhang/PG/pg_message/Monodelphis_domestica_hpg.txt"




#生成的高质量假基因文件，传至集群服务器：/home/zhzhang/PG/HPG/和/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/物种/pgenes/各一份
```
##### 7.高质量假基因生成bed格式
```r
#斑马鱼
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg.bed"


#青鳉
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/Oryzias_latipes_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryzias_latipes/pgenes/Oryzias_latipes_hpg.bed"


#鸡
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg.bed"


#小鼠
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg.bed"


#人
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg.bed"


#猕猴
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg.bed"

#大鼠
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg.bed"

#兔子
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg.bed"

#负鼠
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg.bed"




```

