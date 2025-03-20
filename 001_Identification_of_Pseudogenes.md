### 一.PseudoPipe鉴定八个物种假基因
##### 1.rpeatmasked(rm)基因组序列下载和unmasked染色体序列下载（还需解压和rm内去掉MT等非染色体）
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
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz
for i in $(seq 1 22)
do
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz
done
#解压
gzip -d /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/*
#rmdna去除非必要部分
cat -n "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"|grep ">"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"|sed -n "1,47916726p" > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.mainchr.fa"
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.primary_assembly.fa"|sed -n "47917005,51471479p" >> "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.mainchr.fa"











#sup
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
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/pep/
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
gzip -d Homo_sapiens.GRCh38.pep.all.fa.gz
#查看蛋白所属的基因类型（能编码蛋白的基因的类型）
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa"|awk '{print $6}'|sort|uniq
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




#sup
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
cd /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/
wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz
gzip -d "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/Homo_sapiens.GRCh38.108.chr.gtf.gz"
#提取全部染色体上的蛋白编码基因外显子坐标
cat "/home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/Homo_sapiens.GRCh38.108.chr.gtf"|grep -E "protein_coding|IG_C_gene|IG_D_gene|IG_J_gene|IG_V_gene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene"|sort -n -k 1 -k 4 -k 5|awk '$3=="exon" {print NR"\t"$1"\t"$4"\t"$5}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/Homo_sapiens.GRCh38.108.exLocs.txt
#从全部染色体外显子坐标文件中提取出每个单独的染色体外显子坐标
for i in $(seq 1 22)
do
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/Homo_sapiens.GRCh38.108.exLocs.txt|awk -v a=${i} '$2==a {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/chr${i}_exLocs
done
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/Homo_sapiens.GRCh38.108.exLocs.txt|awk '$2=="X" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/chrX_exLocs
cat /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/Homo_sapiens.GRCh38.108.exLocs.txt|awk '$2=="Y" {print $1"\t"$2"\t"$3"\t"$4}' > /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/chrY_exLocs








#sup
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


gene_biotype:IG_V_gene
gene_biotype:protein_coding
gene_biotype:TR_V_gene
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
./pseudopipe_pro.sh /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.mainchr.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.%s.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/Homo_sapiens/mysql/chr%s_exLocs 0




#sup
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


#猕猴
#负链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/minus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/minus/pexons/ALL.pgexon.gff"
for i in $(seq 2 20)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/minus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/minus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/minus/pexons/X.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/minus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/minus/pexons/Y.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/minus/pexons/ALL.pgexon.gff
#正链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/plus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/plus/pexons/ALL.pgexon.gff"
for i in $(seq 2 20)
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/plus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/plus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/plus/pexons/X.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/plus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/plus/pexons/Y.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/plus/pexons/ALL.pgexon.gff
#正负链外显子合并
cp /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/minus/pexons/ALL.pgexon.gff /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_pgexon.txt
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/plus/pexons/ALL.pgexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_pgexon.txt


#负鼠兔子大鼠
#注释掉Y
sp="Monodelphis_domestica"
num="8"
#注释掉Y
sp="Oryctolagus_cuniculus"
num="21"
#不注释
sp="Rattus_norvegicus"
num="20"
#负链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/minus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/minus/pexons/ALL.pgexon.gff"
for i in $(seq 2 ${num})
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/minus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/minus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/minus/pexons/X.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/minus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/minus/pexons/Y.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/minus/pexons/ALL.pgexon.gff
#正链合并
cp "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/plus/pexons/1.pexon.gff" "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/plus/pexons/ALL.pgexon.gff"
for i in $(seq 2 ${num})
do
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/plus/pexons/${i}.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/plus/pexons/ALL.pgexon.gff
done
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/plus/pexons/X.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/plus/pexons/ALL.pgexon.gff
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/plus/pexons/Y.pexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/plus/pexons/ALL.pgexon.gff
#正负链外显子合并
cp /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/minus/pexons/ALL.pgexon.gff /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/${sp}_pgexon.txt
tail -n +2 /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/plus/pexons/ALL.pgexon.gff >> /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/${sp}_pgexon.txt



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
##### 7.高质量假基因和全部假基因外显子生成bed格式
```r
#斑马鱼
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg.bed"
#假基因外显子生成bed格式【chr start end 母蛋白ID exon strand】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_pgexon.txt" | awk '{print $1"\t"$5"\t"$6"\t"$8"\t""exon""\t"$7}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_pgexon.bed"


#鸡
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg.bed"
#假基因外显子生成bed格式【chr start end 母蛋白ID exon strand】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_pgexon.txt" | awk '{print $1"\t"$5"\t"$6"\t"$8"\t""exon""\t"$7}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_pgexon.bed"


#小鼠
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg.bed"
#假基因外显子生成bed格式【chr start end 母蛋白ID exon strand】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_pgexon.txt" | awk '{print $1"\t"$5"\t"$6"\t"$8"\t""exon""\t"$7}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_pgexon.bed"


#人
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg.bed"
#假基因外显子生成bed格式【chr start end 母蛋白ID exon strand】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_pgexon.txt" | awk '{print $1"\t"$5"\t"$6"\t"$8"\t""exon""\t"$7}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_pgexon.bed"


#猕猴
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg.bed"
#假基因外显子生成bed格式【chr start end 母蛋白ID exon strand】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_pgexon.txt" | awk '{print $1"\t"$5"\t"$6"\t"$8"\t""exon""\t"$7}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_pgexon.bed"

#大鼠
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg.bed"
#假基因外显子生成bed格式【chr start end 母蛋白ID exon strand】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_pgexon.txt" | awk '{print $1"\t"$5"\t"$6"\t"$8"\t""exon""\t"$7}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_pgexon.bed"

#兔子
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg.bed"
#假基因外显子生成bed格式【chr start end 母蛋白ID exon strand】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_pgexon.txt" | awk '{print $1"\t"$5"\t"$6"\t"$8"\t""exon""\t"$7}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_pgexon.bed"

#负鼠
#高质量假基因生成bed格式【chr start end 母蛋白ID gene strand 假基因ID】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg.txt" | awk '{print $1"\t"$2"\t"$3"\t"$5"\t""gene""\t"$4"\t"$16}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg.bed"
#假基因外显子生成bed格式【chr start end 母蛋白ID exon strand】
tail -n +2 "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_pgexon.txt" | awk '{print $1"\t"$5"\t"$6"\t"$8"\t""exon""\t"$7}' > "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_pgexon.bed"


```
##### 8.Bedtools筛选出高质量假基因的外显子
```r
#斑马鱼
#与某个高质量假基因位置有交集且母蛋白与该假基因一致的外显子，就是该高质量假基因的外显子。根据这条准则筛选出高质量假基因的外显子，以及外显子所对应的高质量假基因ID【输出bed：chr start end 所属假基因ID exon strand】
bedtools intersect -a /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_pgexon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg.bed -wo -s |awk '$4==$10 {print $1"\t"$2"\t"$3"\t"$13"\t"$5"\t"$6}'> /home/zhzhang/PG/HPG/Danio_rerio_hpgexon.bed


#鸡
#与某个高质量假基因位置有交集且母蛋白与该假基因一致的外显子，就是该高质量假基因的外显子。根据这条准则筛选出高质量假基因的外显子，以及外显子所对应的高质量假基因ID【输出bed：chr start end 所属假基因ID exon strand】
bedtools intersect -a /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_pgexon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg.bed -wo -s |awk '$4==$10 {print $1"\t"$2"\t"$3"\t"$13"\t"$5"\t"$6}'> /home/zhzhang/PG/HPG/Gallus_gallus_hpgexon.bed



#小鼠
#与某个高质量假基因位置有交集且母蛋白与该假基因一致的外显子，就是该高质量假基因的外显子。根据这条准则筛选出高质量假基因的外显子，以及外显子所对应的高质量假基因ID【输出bed：chr start end 所属假基因ID exon strand】
bedtools intersect -a /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_pgexon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg.bed -wo -s |awk '$4==$10 {print $1"\t"$2"\t"$3"\t"$13"\t"$5"\t"$6}'> /home/zhzhang/PG/HPG/Mus_musculus_hpgexon.bed


#人
#与某个高质量假基因位置有交集且母蛋白与该假基因一致的外显子，就是该高质量假基因的外显子。根据这条准则筛选出高质量假基因的外显子，以及外显子所对应的高质量假基因ID【输出bed：chr start end 所属假基因ID exon strand】
bedtools intersect -a /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_pgexon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg.bed -wo -s |awk '$4==$10 {print $1"\t"$2"\t"$3"\t"$13"\t"$5"\t"$6}'> /home/zhzhang/PG/HPG/Homo_sapiens_hpgexon.bed

#other 4
#与某个高质量假基因位置有交集且母蛋白与该假基因一致的外显子，就是该高质量假基因的外显子。根据这条准则筛选出高质量假基因的外显子，以及外显子所对应的高质量假基因ID【输出bed：chr start end 所属假基因ID exon strand】
#猕猴负鼠兔子大鼠
sp="Macaca_mulatta"
sp="Monodelphis_domestica"
sp="Oryctolagus_cuniculus"
sp="Rattus_norvegicus"
bedtools intersect -a /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/${sp}_pgexon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${sp}/pgenes/${sp}_hpg.bed -wo -s |awk '$4==$10 {print $1"\t"$2"\t"$3"\t"$13"\t"$5"\t"$6}'> /home/zhzhang/PG/HPG/${sp}_hpgexon.bed


```


##### 9.对比鉴定的假基因和ENSEMBLE注释的假基因
```r
#ensemble pseudogene region bed file
ls /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/ |while read i
do
gtf=`ls /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/${i}/mysql/*.gtf`
grep pseudogene ${gtf}|awk '$3=="gene" && $17=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$18"\t"$7} $3=="gene" && $15=="gene_biotype" {print $1"\t"$4"\t"$5"\t"$10"\t"$16"\t"$7}'|sed "s/\"//g;s/;//g"|bedtools sort > /home/zhzhang/PG/PseudoPipe_wd/compare/${i}.ensemble.pseudogene.bed
done
grep -v rRNA_pseudogene "/home/zhzhang/PG/PseudoPipe_wd/compare/Homo_sapiens.ensemble.pseudogene.bed" > /home/zhzhang/PG/PseudoPipe_wd/compare/Homo_sapiens.ensemble.pseudogene.1bed
rm /home/zhzhang/PG/PseudoPipe_wd/compare/Homo_sapiens.ensemble.pseudogene.bed
mv /home/zhzhang/PG/PseudoPipe_wd/compare/Homo_sapiens.ensemble.pseudogene.1bed /home/zhzhang/PG/PseudoPipe_wd/compare/Homo_sapiens.ensemble.pseudogene.bed


#ensemble假基因总数统计
wc -l /home/zhzhang/PG/PseudoPipe_wd/compare/*ensemble.pseudogene.bed
#被覆盖的数量
ls /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/ |grep -v Oryzias_latipes|while read i
do
bedtools intersect -a "/home/zhzhang/PG/PseudoPipe_wd/compare/${i}.ensemble.pseudogene.bed" -b "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg.bed" -wa -s -f 0.5|sort -u|wc -l
done
#新的数量
ls /home/zhzhang/PG/PseudoPipe_wd/ppipe_input/ |grep -v Oryzias_latipes|while read i
do
bedtools intersect -a "/home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/${i}/pgenes/${i}_hpg.bed" -b "/home/zhzhang/PG/PseudoPipe_wd/compare/${i}.ensemble.pseudogene.bed" -wa -s -F 0.1 -v|sort -u|wc -l
done
#(sp ensemble_num cover_num)
Zebrafish 324 251 
Chicken 48 29
Human 14728 12398
Macaque 749 271
Opossum 828 352
Mouse 13652 11142
Rabbit 477 349
Rat 913 683











```

