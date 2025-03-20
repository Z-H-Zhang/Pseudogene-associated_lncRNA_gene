### 七.人类多组织单细胞转录组数据进行衰老相关分析
#### 1.数据下载
```r
#创建参考基因组
cd /home/zhzhang/PG/scRNAseq/crref/
cellranger mkref --genome=grch38_108 --fasta="/home/zhzhang/PG/scRNAseq/crref/Homo_sapiens.GRCh38.dna.chrall.fa" --genes="/home/zhzhang/PG/scRNAseq/crref/Homo_sapiens.GRCh38.108.chrall.rmpg.novellncRNA.gtf" --nthreads=64 --memgb=254


#ovary
#下sra数据
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255690
#SRA生成fastq，压缩，以及标准命名
cat "/home/zhzhang/PG/aging/ovary/ovarysrr.txt"|while read i
do
srr=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
fasterq-dump -e 64 -m 2048GB --split-3 -O /home/zhzhang/PG/aging/ovary/fastq/${sample}/ -o ${sample} /home/zhzhang/PG/aging/ovary/srr/${srr}
pigz -p 64 /home/zhzhang/PG/aging/ovary/fastq/${sample}/${sample}_1.fastq
pigz -p 64 /home/zhzhang/PG/aging/ovary/fastq/${sample}/${sample}_2.fastq
mv /home/zhzhang/PG/aging/ovary/fastq/${sample}/${sample}_1.fastq.gz /home/zhzhang/PG/aging/ovary/fastq/${sample}/${sample}_S1_L001_R1_001.fastq.gz
mv /home/zhzhang/PG/aging/ovary/fastq/${sample}/${sample}_2.fastq.gz /home/zhzhang/PG/aging/ovary/fastq/${sample}/${sample}_S1_L001_R2_001.fastq.gz
done



#testis
#下sra数据
https://www.ebi.ac.uk/ena/browser/view/PRJNA757777
#SRA生成fastq，压缩，以及标准命名
ls /home/zhzhang/PG/scRNAseq/sra/ > /home/zhzhang/PG/scRNAseq/samplename.txt
cat /home/zhzhang/PG/scRNAseq/samplename.txt|while read i
do
fasterq-dump -e 64 -m 2048GB --split-3 -O /home/zhzhang/PG/scRNAseq/fastq/${i}/ -o ${i} /home/zhzhang/PG/scRNAseq/sra/${i}
pigz -p 64 /home/zhzhang/PG/scRNAseq/fastq/${i}/${i}_1.fastq
pigz -p 64 /home/zhzhang/PG/scRNAseq/fastq/${i}/${i}_2.fastq
mv /home/zhzhang/PG/scRNAseq/fastq/${i}/${i}_1.fastq /home/zhzhang/PG/scRNAseq/fastq/${i}/${i}_S1_L001_R1_001.fastq.gz
mv /home/zhzhang/PG/scRNAseq/fastq/${i}/${i}_2.fastq /home/zhzhang/PG/scRNAseq/fastq/${i}/${i}_S1_L001_R2_001.fastq.gz
done




#skin
#fastq下载
#R1（10,12 21,24）
cat "/home/zhzhang/PG/aging/skin/skinfastq.txt"|sed -n '21,24p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}_f1.fastq.gz" -O /home/zhzhang/PG/aging/skin/fastq/${sample}/${sample}_S1_L00${lane}_R1_001.fastq.gz
done
#R2（12 21）
cat "/home/zhzhang/PG/aging/skin/skinfastq.txt"|sed -n '12p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}_r2.fastq.gz" -O /home/zhzhang/PG/aging/skin/fastq/${sample}/${sample}_S1_L00${lane}_R2_001.fastq.gz
done




#muscle
#SRR
cat "/home/zhzhang/PG/aging/muscle/musclefastq.txt"|sed -n '15p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
wget -t 0 -c "${fastq}" -P /home/zhzhang/PG/aging/muscle/fastq/${sample}/
done
#SRA生成fastq，以及标准命名
source /home/zhzhang/miniconda3/bin/activate parafq
cat "/home/zhzhang/PG/aging/muscle/musclefastq.txt"|while read i
do
srr=`echo "${i}"|awk '{print $3}'`
sample=`echo "${i}"|awk '{print $2}'`
parallel-fastq-dump -t 28 -O /home/zhzhang/PG/aging/muscle/fastq/${sample}/ --split-3 --gzip -s /home/zhzhang/PG/aging/muscle/fastq/${sample}/${srr}
mv /home/zhzhang/PG/aging/muscle/fastq/${sample}/${srr}_1.fastq.gz /home/zhzhang/PG/aging/muscle/fastq/${sample}/${sample}_S1_L001_R1_001.fastq.gz
mv /home/zhzhang/PG/aging/muscle/fastq/${sample}/${srr}_2.fastq.gz /home/zhzhang/PG/aging/muscle/fastq/${sample}/${sample}_S1_L001_R2_001.fastq.gz
done




#bonemarrow
#SRR
cat "/home/zhzhang/PG/aging/bonemarrow/bonemarrowfastq.txt"|sed -n '1,8p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
wget -t 0 -c "${fastq}" -P /home/zhzhang/PG/aging/bonemarrow/fastq/${sample}/
done
#SRA生成fastq，以及标准命名
source /home/zhzhang/miniconda3/bin/activate parafq
cat "/home/zhzhang/PG/aging/bonemarrow/bonemarrowfastq.txt"|while read i
do
srr=`echo "${i}"|awk '{print $3}'`
sample=`echo "${i}"|awk '{print $2}'`
parallel-fastq-dump -t 32 -O /home/zhzhang/PG/aging/bonemarrow/fastq/${sample}/ --split-3 --gzip -s /home/zhzhang/PG/aging/bonemarrow/fastq/${sample}/${srr}
mv /home/zhzhang/PG/aging/bonemarrow/fastq/${sample}/${srr}_1.fastq.gz /home/zhzhang/PG/aging/bonemarrow/fastq/${sample}/${sample}_S1_L001_R1_001.fastq.gz
mv /home/zhzhang/PG/aging/bonemarrow/fastq/${sample}/${srr}_2.fastq.gz /home/zhzhang/PG/aging/bonemarrow/fastq/${sample}/${sample}_S1_L001_R2_001.fastq.gz
done



#heart
#fastq下载
#R1 8
zt="2\n3\n10\n21"
for ii in $(seq 1 4)
do
num=`echo -e ${zt}|sed -n "${ii}p"`
#
cat "/home/zhzhang/PG/aging/heart/heartfastq.txt"|sed -n "9p"|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}_1.fastq.gz" -O /home/zhzhang/PG/aging/heart/fastq/${sample}/${sample}_S1_L00${lane}_R1_001.fastq.gz
done
#
done
#R2
zt="2\n9\n14\n21"
for ii in $(seq 1 4)
do
num=`echo -e ${zt}|sed -n "${ii}p"`
#
cat "/home/zhzhang/PG/aging/heart/heartfastq.txt"|sed -n "${num}p"|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}_2.fastq.gz" -O /home/zhzhang/PG/aging/heart/fastq/${sample}/${sample}_S1_L00${lane}_R2_001.fastq.gz
done
#
done


#lung
#fastq下载
#R1 18
cat "/home/zhzhang/PG/aging/lung/lungfastq.txt"|sed -n '18p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}_1.fastq.gz" -O /home/zhzhang/PG/aging/lung/fastq/${sample}/${sample}_S1_L00${lane}_R1_001.fastq.gz
done
#R2 18
cat "/home/zhzhang/PG/aging/lung/lungfastq.txt"|sed -n '18p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}_2.fastq.gz" -O /home/zhzhang/PG/aging/lung/fastq/${sample}/${sample}_S1_L00${lane}_R2_001.fastq.gz
done


#brain
#fastq下载
#R1 （7-8）
cat "/home/zhzhang/PG/aging/brain/brainfastq.txt"|sed -n '7,8p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}_1.fastq.gz" -O /home/zhzhang/PG/aging/brain/fastq/${sample}/${sample}_S1_L00${lane}_R1_001.fastq.gz
done
#R2 （4 7 9）
cat "/home/zhzhang/PG/aging/brain/brainfastq.txt"|sed -n '9p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}_2.fastq.gz" -O /home/zhzhang/PG/aging/brain/fastq/${sample}/${sample}_S1_L00${lane}_R2_001.fastq.gz
done
#
cat "/home/zhzhang/PG/aging/brain/brainfastq.txt"|sed -n '1,9p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
seqkit stats -j 60 -T /home/zhzhang/PG/aging/brain/fastq/${sample}/${sample}_S1_L00${lane}_R1_001.fastq.gz
seqkit stats -j 60 -T /home/zhzhang/PG/aging/brain/fastq/${sample}/${sample}_S1_L00${lane}_R2_001.fastq.gz
done




#liver
#bam/srr下载(仅Y5是srr)
#BAM
cat "/home/zhzhang/PG/aging/liver/liverfastq.txt"|sed -n '4p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}" -O /home/zhzhang/PG/aging/liver/fastq/${sample}/${sample}.bam
done
#srr
cat "/home/zhzhang/PG/aging/liver/liverfastq.txt"|sed -n '8,13p'|while read i
do
fastq=`echo "${i}"|awk '{print $1}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
wget -t 0 -c "${fastq}" -P /home/zhzhang/PG/aging/liver/fastq/${sample}/
done
#SRA生成fastq，以及标准命名
source /home/zhzhang/miniconda3/bin/activate parafq
cat "/home/zhzhang/PG/aging/liver/liverfastq.txt"|sed -n '8,13p'|while read i
do
srr=`echo "${i}"|awk '{print $1}'|awk -F '/' '{print $6}'`
sample=`echo "${i}"|awk '{print $2}'`
lane=`echo "${i}"|awk '{print $3}'`
parallel-fastq-dump -t 32 -O /home/zhzhang/PG/aging/liver/fastq/${sample}/ --split-3 --gzip -s /home/zhzhang/PG/aging/liver/fastq/${sample}/${srr}
mv /home/zhzhang/PG/aging/liver/fastq/${sample}/${srr}_1.fastq.gz /home/zhzhang/PG/aging/liver/fastq/${sample}/${sample}_S1_L00${lane}_R1_001.fastq.gz
mv /home/zhzhang/PG/aging/liver/fastq/${sample}/${srr}_2.fastq.gz /home/zhzhang/PG/aging/liver/fastq/${sample}/${sample}_S1_L00${lane}_R2_001.fastq.gz
mv /home/zhzhang/PG/aging/liver/fastq/${sample}/${srr}_3.fastq.gz /home/zhzhang/PG/aging/liver/fastq/${sample}/${sample}_S1_L00${lane}_I1_001.fastq.gz
done
#bam生成fastq
cat "/home/zhzhang/PG/aging/liver/liverfastq.txt"|sed -n '1,7p'|while read i
do
sample=`echo "${i}"|awk '{print $2}'`
bamtofastq_linux --nthreads=64 --reads-per-fastq=50000000000000000 /home/zhzhang/PG/aging/liver/fastq/${sample}/${sample}.bam /home/zhzhang/PG/aging/liver/fastq/${sample}/tmp/
done
#bam生成fastq的改名
cat "/home/zhzhang/PG/aging/liver/liverfastq.txt"|sed -n '1,7p'|while read i
do
sample=`echo "${i}"|awk '{print $2}'`
#fastq文件夹
random=`ls /home/zhzhang/PG/aging/liver/fastq/${sample}/tmp/`
#fastq文件夹数量
randomnum=`echo "${random}"|wc -l`
#
if [ "$randomnum" == "1" ]
then
random1="/home/zhzhang/PG/aging/liver/fastq/${sample}/tmp/${random}"
ls ${random1} > /home/zhzhang/PG/aging/liver/fastq/jiu.txt
ls ${random1}|sed "s/bamtofastq/${sample}/g" > /home/zhzhang/PG/aging/liver/fastq/xin.txt
filenum=`wc -l /home/zhzhang/PG/aging/liver/fastq/jiu.txt|awk '{print $1}'`
for ii in $(seq 1 ${filenum})
do
jiufile=`cat /home/zhzhang/PG/aging/liver/fastq/jiu.txt|sed -n "${ii}p"`
xinfile=`cat /home/zhzhang/PG/aging/liver/fastq/xin.txt|sed -n "${ii}p"`
mv "${random1}/${jiufile}" "/home/zhzhang/PG/aging/liver/fastq/${sample}/${xinfile}"
done
fi
#
if [ "$randomnum" == "2" ]
then
pre1=`echo "${random}"|sed -n '1p'`
pre2=`echo "${random}"|sed -n '2p'`
random1="/home/zhzhang/PG/aging/liver/fastq/${sample}/tmp/${pre1}"
random2="/home/zhzhang/PG/aging/liver/fastq/${sample}/tmp/${pre2}"
ls ${random1} > /home/zhzhang/PG/aging/liver/fastq/jiu.txt
ls ${random1}|sed "s/bamtofastq/${sample}/g" > /home/zhzhang/PG/aging/liver/fastq/xin.txt
filenum=`wc -l /home/zhzhang/PG/aging/liver/fastq/jiu.txt|awk '{print $1}'`
for iii in $(seq 1 ${filenum})
do
jiufile=`cat /home/zhzhang/PG/aging/liver/fastq/jiu.txt|sed -n "${iii}p"`
xinfile=`cat /home/zhzhang/PG/aging/liver/fastq/xin.txt|sed -n "${iii}p"`
mv "${random1}/${jiufile}" "/home/zhzhang/PG/aging/liver/fastq/${sample}/${xinfile}"
done
ls ${random2} > /home/zhzhang/PG/aging/liver/fastq/jiu.txt
ls ${random2}|sed "s/bamtofastq_S1/${sample}_S2/g" > /home/zhzhang/PG/aging/liver/fastq/xin.txt
filenum=`wc -l /home/zhzhang/PG/aging/liver/fastq/jiu.txt|awk '{print $1}'`
for iiii in $(seq 1 ${filenum})
do
jiufile=`cat /home/zhzhang/PG/aging/liver/fastq/jiu.txt|sed -n "${iiii}p"`
xinfile=`cat /home/zhzhang/PG/aging/liver/fastq/xin.txt|sed -n "${iiii}p"`
mv "${random2}/${jiufile}" "/home/zhzhang/PG/aging/liver/fastq/${sample}/${xinfile}"
done
fi
done
rm -r /home/zhzhang/PG/aging/liver/fastq/*/tmp/
#liver O3的s2改为s1 L9 L10 .Y2的s2改为s1 L9 L10



```
#### 2.cellranger count
```r
#cellranger产生矩阵
#组织变量
tissue="ovary"
tissue="testis" #file:cat /home/zhzhang/PG/scRNAseq/samplename.txt|while read i
tissue="skin"
tissue="heart"
tissue="liver"
tissue="lung"
tissue="brain" #--chemistry='ARC-v1'
tissue="bonemarrow"
tissue="muscle"
#cellranger count
cd /home/zhzhang/PG/aging/${tissue}/count/
awk '{print $2}' "/home/zhzhang/PG/aging/${tissue}/${tissue}fastq.txt"|sort|uniq|while read i
do
sample="${i}"
cellranger count --id=${sample} --transcriptome=/home/zhzhang/PG/scRNAseq/crref/grch38_108/ --fastqs=/home/zhzhang/PG/aging/${tissue}/fastq/${sample}/ --sample=${sample} --nosecondary --localcores=64 --localmem=2048
done


```


#### 3.细胞聚类
##### ovary
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)
plan("multiprocess", workers = 28)
options(future.globals.maxSize= 2500000000000)
# Load the dataset
O1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/ovary/count/Ovary_O1/outs/filtered_feature_bc_matrix/")
O2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/ovary/count/Ovary_O2/outs/filtered_feature_bc_matrix/")
O3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/ovary/count/Ovary_O3/outs/filtered_feature_bc_matrix/")
M1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/ovary/count/Ovary_M1/outs/filtered_feature_bc_matrix/")
M2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/ovary/count/Ovary_M2/outs/filtered_feature_bc_matrix/")
M3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/ovary/count/Ovary_M3/outs/filtered_feature_bc_matrix/")
Y1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/ovary/count/Ovary_Y1/outs/filtered_feature_bc_matrix/")
Y2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/ovary/count/Ovary_Y2/outs/filtered_feature_bc_matrix/")
Y3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/ovary/count/Ovary_Y3/outs/filtered_feature_bc_matrix/")
#创建对象
O1 <- CreateSeuratObject(counts = O1, project = "O1", min.cells = 3, min.features = 200)
O2 <- CreateSeuratObject(counts = O2, project = "O2", min.cells = 3, min.features = 200)
O3 <- CreateSeuratObject(counts = O3, project = "O3", min.cells = 3, min.features = 200)
M1 <- CreateSeuratObject(counts = M1, project = "M1", min.cells = 3, min.features = 200)
M2 <- CreateSeuratObject(counts = M2, project = "M2", min.cells = 3, min.features = 200)
M3 <- CreateSeuratObject(counts = M3, project = "M3", min.cells = 3, min.features = 200)
Y1 <- CreateSeuratObject(counts = Y1, project = "Y1", min.cells = 3, min.features = 200)
Y2 <- CreateSeuratObject(counts = Y2, project = "Y2", min.cells = 3, min.features = 200)
Y3 <- CreateSeuratObject(counts = Y3, project = "Y3", min.cells = 3, min.features = 200)

#merge
ALLsample <- merge(O1,y = c(O2,O3,M1,M2,M3,Y1,Y2,Y3), add.cell.ids = c("O1","O2","O3","M1","M2","M3","Y1","Y2","Y3"), project = "ALLsample")

#细胞质控
#VlnPlot(ALLsample, features = c("nFeature_RNA", "nCount_RNA"), ncol =2, group.by = "orig.ident", pt.size = 0)
ALLsample <- PercentageFeatureSet(ALLsample,pattern = "^MT-",col.name = "percent_mt")
ALLsample <- subset(ALLsample, subset = nFeature_RNA > 500 & percent_mt < 50)
#标准化归一化
ALLsample <- NormalizeData(ALLsample, assay ="RNA",normalization.method = "LogNormalize", scale.factor = 10000)
ALLsample <- FindVariableFeatures(ALLsample, selection.method = "vst",nfeatures=5000)
ALLsample <- ScaleData(ALLsample)
#降维
ALLsample <- RunPCA(ALLsample,npcs = 50)
#DimPlot(object = ALLsample, dims=c(1,2),reduction = "pca")
#ElbowPlot(ALLsample, ndims = 50, reduction = "pca")
ALLsample <- RunHarmony(ALLsample,reduction.use = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#聚类
ALLsample <- FindNeighbors(ALLsample, dims = 1:20,reduction = 'harmony')
ALLsample <- FindClusters(ALLsample, resolution = 0.8)
#clustree(ALLsample)
#DimPlot(object = ALLsample, reduction = "pca",dim=c(1,2),label = TRUE,group.by="orig.ident")
#可视化降维
ALLsample <- RunUMAP(ALLsample, dims = 1:20,reduction = 'harmony')
#DimPlot(object = ALLsample, reduction = "umap",dim=c(1,2),group.by="seurat_clusters",label = TRUE)
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/ovary/seurat/RNAnormal_ALLsample.rds")
#marker
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,return.thresh = 0.01)
data.table::fwrite(ALLsample.markers,file ="/home/zhzhang/PG/aging/ovary/seurat/ALLsample.allcluster.markergene.byRNA.txt",sep = '\t',row.names = F,quote = F,col.names = T)


```
##### testis
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)
plan("multiprocess", workers = 28)
options(future.globals.maxSize= 2500000000000)

# Load the dataset
O1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/O1/outs/filtered_feature_bc_matrix/")
O2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/O2/outs/filtered_feature_bc_matrix/")
O3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/O3/outs/filtered_feature_bc_matrix/")
O4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/O4/outs/filtered_feature_bc_matrix/")
O5 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/O5/outs/filtered_feature_bc_matrix/")
O6 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/O6/outs/filtered_feature_bc_matrix/")
O7 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/O7/outs/filtered_feature_bc_matrix/")
O8 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/O8/outs/filtered_feature_bc_matrix/")
Y1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/Y1/outs/filtered_feature_bc_matrix/")
Y2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/Y2/outs/filtered_feature_bc_matrix/")
Y3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/Y3/outs/filtered_feature_bc_matrix/")
Y4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/testis/count/Y4/outs/filtered_feature_bc_matrix/")
#创建对象
O1 <- CreateSeuratObject(counts = O1, project = "O1", min.cells = 3, min.features = 200)
O2 <- CreateSeuratObject(counts = O2, project = "O2", min.cells = 3, min.features = 200)
O3 <- CreateSeuratObject(counts = O3, project = "O3", min.cells = 3, min.features = 200)
O4 <- CreateSeuratObject(counts = O4, project = "O4", min.cells = 3, min.features = 200)
O5 <- CreateSeuratObject(counts = O5, project = "O5", min.cells = 3, min.features = 200)
O6 <- CreateSeuratObject(counts = O6, project = "O6", min.cells = 3, min.features = 200)
O7 <- CreateSeuratObject(counts = O7, project = "O7", min.cells = 3, min.features = 200)
O8 <- CreateSeuratObject(counts = O8, project = "O8", min.cells = 3, min.features = 200)
Y1 <- CreateSeuratObject(counts = Y1, project = "Y1", min.cells = 3, min.features = 200)
Y2 <- CreateSeuratObject(counts = Y2, project = "Y2", min.cells = 3, min.features = 200)
Y3 <- CreateSeuratObject(counts = Y3, project = "Y3", min.cells = 3, min.features = 200)
Y4 <- CreateSeuratObject(counts = Y4, project = "Y4", min.cells = 3, min.features = 200)

#merge
ALLsample <- merge(O1,y = c(O2,O3,O4,O5,O6,O7,O8,Y1,Y2,Y3,Y4), add.cell.ids = c("O1","O2","O3","O4","O5","O6","O7","O8","Y1","Y2","Y3","Y4"), project = "ALLsample")

#细胞质控
#VlnPlot(ALLsample, features = c("nFeature_RNA", "nCount_RNA"), ncol =2, group.by = "orig.ident", pt.size = 0)
ALLsample <- PercentageFeatureSet(ALLsample,pattern = "^MT-",col.name = "percent_mt")
ALLsample <- subset(ALLsample, subset = nFeature_RNA > 500 & percent_mt < 50)
#标准化归一化
ALLsample <- NormalizeData(ALLsample, assay ="RNA",normalization.method = "LogNormalize", scale.factor = 10000)
ALLsample <- FindVariableFeatures(ALLsample, selection.method = "vst",nfeatures=5000)
ALLsample <- ScaleData(ALLsample)
#降维
ALLsample <- RunPCA(ALLsample,npcs = 50)
#DimPlot(object = ALLsample, dims=c(1,2),reduction = "pca")
#ElbowPlot(ALLsample, ndims = 50, reduction = "pca")
ALLsample <- RunHarmony(ALLsample,reduction.use = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#聚类
ALLsample <- FindNeighbors(ALLsample, dims = 1:20,reduction = 'harmony')
ALLsample <- FindClusters(ALLsample)
#clustree(ALLsample)
#DimPlot(object = ALLsample, reduction = "pca",dim=c(1,2),label = TRUE,group.by="orig.ident")
#可视化降维
ALLsample <- RunUMAP(ALLsample, dims = 1:20,reduction = 'harmony')
#DimPlot(object = ALLsample, reduction = "umap",dim=c(1,2),group.by="seurat_clusters",label = TRUE)
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/testis/seurat/RNAnormal_ALLsample.rds")
#marker
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,return.thresh = 0.01)
data.table::fwrite(ALLsample.markers,file ="/home/zhzhang/PG/aging/testis/seurat/ALLsample.allcluster.markergene.byRNA.txt",sep = '\t',row.names = F,quote = F,col.names = T)


```
##### skin
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)
plan("multiprocess", workers = 28)
options(future.globals.maxSize= 2500000000000)
# Load the dataset
O1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/skin/count/O1/outs/filtered_feature_bc_matrix/")
O2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/skin/count/O2/outs/filtered_feature_bc_matrix/")
O3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/skin/count/O3/outs/filtered_feature_bc_matrix/")
Y1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/skin/count/Y1/outs/filtered_feature_bc_matrix/")
Y2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/skin/count/Y2/outs/filtered_feature_bc_matrix/")
Y3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/skin/count/Y3/outs/filtered_feature_bc_matrix/")
#创建对象
O1 <- CreateSeuratObject(counts = O1, project = "O1", min.cells = 3, min.features = 200)
O2 <- CreateSeuratObject(counts = O2, project = "O2", min.cells = 3, min.features = 200)
O3 <- CreateSeuratObject(counts = O3, project = "O3", min.cells = 3, min.features = 200)
Y1 <- CreateSeuratObject(counts = Y1, project = "Y1", min.cells = 3, min.features = 200)
Y2 <- CreateSeuratObject(counts = Y2, project = "Y2", min.cells = 3, min.features = 200)
Y3 <- CreateSeuratObject(counts = Y3, project = "Y3", min.cells = 3, min.features = 200)

#merge
ALLsample <- merge(O1,y = c(O2,O3,Y1,Y2,Y3), add.cell.ids = c("O1","O2","O3","Y1","Y2","Y3"), project = "ALLsample")

#细胞质控
#VlnPlot(ALLsample, features = c("nFeature_RNA", "nCount_RNA"), ncol =2, group.by = "orig.ident", pt.size = 0)
ALLsample <- PercentageFeatureSet(ALLsample,pattern = "^MT-",col.name = "percent_mt")
ALLsample <- subset(ALLsample, subset = nFeature_RNA > 500 & percent_mt < 50)
#标准化归一化
ALLsample <- NormalizeData(ALLsample, assay ="RNA",normalization.method = "LogNormalize", scale.factor = 10000)
ALLsample <- FindVariableFeatures(ALLsample, selection.method = "vst",nfeatures=5000)
ALLsample <- ScaleData(ALLsample)
#降维
ALLsample <- RunPCA(ALLsample,npcs = 50)
#DimPlot(object = ALLsample, dims=c(1,2),reduction = "pca")
#ElbowPlot(ALLsample, ndims = 50, reduction = "pca")
ALLsample <- RunHarmony(ALLsample,reduction.use = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#聚类
ALLsample <- FindNeighbors(ALLsample, dims = 1:20,reduction = 'harmony')
ALLsample <- FindClusters(ALLsample)
#clustree(ALLsample)
#DimPlot(object = ALLsample, reduction = "pca",dim=c(1,2),label = TRUE,group.by="orig.ident")
#可视化降维
ALLsample <- RunUMAP(ALLsample, dims = 1:20,reduction = 'harmony')
#DimPlot(object = ALLsample, reduction = "umap",dim=c(1,2),group.by="seurat_clusters",label = TRUE)
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/skin/seurat/RNAnormal_ALLsample.rds")
#marker
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,return.thresh = 0.01)
data.table::fwrite(ALLsample.markers,file ="/home/zhzhang/PG/aging/skin/seurat/ALLsample.allcluster.markergene.byRNA.txt",sep = '\t',row.names = F,quote = F,col.names = T)



```
##### liver
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)
plan("multiprocess", workers = 28)
options(future.globals.maxSize= 2500000000000)
# Load the dataset
O1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/liver/count/O1/outs/filtered_feature_bc_matrix/")
O2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/liver/count/O2/outs/filtered_feature_bc_matrix/")
O3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/liver/count/O3/outs/filtered_feature_bc_matrix/")
Y1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/liver/count/Y1/outs/filtered_feature_bc_matrix/")
Y2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/liver/count/Y2/outs/filtered_feature_bc_matrix/")
Y3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/liver/count/Y3/outs/filtered_feature_bc_matrix/")
Y4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/liver/count/Y4/outs/filtered_feature_bc_matrix/")
Y5 <- Read10X(data.dir = "/home/zhzhang/PG/aging/liver/count/Y5/outs/filtered_feature_bc_matrix/")
#创建对象
O1 <- CreateSeuratObject(counts = O1, project = "O1", min.cells = 3, min.features = 200)
O2 <- CreateSeuratObject(counts = O2, project = "O2", min.cells = 3, min.features = 200)
O3 <- CreateSeuratObject(counts = O3, project = "O3", min.cells = 3, min.features = 200)
Y1 <- CreateSeuratObject(counts = Y1, project = "Y1", min.cells = 3, min.features = 200)
Y2 <- CreateSeuratObject(counts = Y2, project = "Y2", min.cells = 3, min.features = 200)
Y3 <- CreateSeuratObject(counts = Y3, project = "Y3", min.cells = 3, min.features = 200)
Y4 <- CreateSeuratObject(counts = Y4, project = "Y4", min.cells = 3, min.features = 200)
Y5 <- CreateSeuratObject(counts = Y5, project = "Y5", min.cells = 3, min.features = 200)
#merge
ALLsample <- merge(O1,y = c(O2,O3,Y1,Y2,Y3,Y4,Y5), add.cell.ids = c("O1","O2","O3","Y1","Y2","Y3","Y4","Y5"), project = "ALLsample")

#细胞质控
#VlnPlot(ALLsample, features = c("nFeature_RNA", "nCount_RNA"), ncol =2, group.by = "orig.ident", pt.size = 0)
ALLsample <- PercentageFeatureSet(ALLsample,pattern = "^MT-",col.name = "percent_mt")
ALLsample <- subset(ALLsample, subset = nFeature_RNA > 500 & percent_mt < 50)
#标准化归一化
ALLsample <- NormalizeData(ALLsample, assay ="RNA",normalization.method = "LogNormalize", scale.factor = 10000)
ALLsample <- FindVariableFeatures(ALLsample, selection.method = "vst",nfeatures=5000)
ALLsample <- ScaleData(ALLsample)
#降维
ALLsample <- RunPCA(ALLsample,npcs = 50)
#DimPlot(object = ALLsample, dims=c(1,2),reduction = "pca")
#ElbowPlot(ALLsample, ndims = 50, reduction = "pca")
ALLsample <- RunHarmony(ALLsample,reduction.use = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#聚类
ALLsample <- FindNeighbors(ALLsample, dims = 1:20,reduction = 'harmony')
ALLsample <- FindClusters(ALLsample)
#clustree(ALLsample)
#DimPlot(object = ALLsample, reduction = "pca",dim=c(1,2),label = TRUE,group.by="orig.ident")
#可视化降维
ALLsample <- RunUMAP(ALLsample, dims = 1:20,reduction = 'harmony')
#DimPlot(object = ALLsample, reduction = "umap",dim=c(1,2),group.by="seurat_clusters",label = TRUE)
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/liver/seurat/RNAnormal_ALLsample.rds")
#marker
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,return.thresh = 0.01)
data.table::fwrite(ALLsample.markers,file ="/home/zhzhang/PG/aging/liver/seurat/ALLsample.allcluster.markergene.byRNA.txt",sep = '\t',row.names = F,quote = F,col.names = T)



```
##### heart
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)
plan("multiprocess", workers = 28)
options(future.globals.maxSize= 2500000000000)

# Load the dataset
O1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/heart/count/O1/outs/filtered_feature_bc_matrix/")
O2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/heart/count/O2/outs/filtered_feature_bc_matrix/")
O3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/heart/count/O3/outs/filtered_feature_bc_matrix/")
O4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/heart/count/O4/outs/filtered_feature_bc_matrix/")
O5 <- Read10X(data.dir = "/home/zhzhang/PG/aging/heart/count/O5/outs/filtered_feature_bc_matrix/")
O6 <- Read10X(data.dir = "/home/zhzhang/PG/aging/heart/count/O6/outs/filtered_feature_bc_matrix/")
Y1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/heart/count/Y1/outs/filtered_feature_bc_matrix/")
Y2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/heart/count/Y2/outs/filtered_feature_bc_matrix/")
Y3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/heart/count/Y3/outs/filtered_feature_bc_matrix/")
#创建对象
O1 <- CreateSeuratObject(counts = O1, project = "O1", min.cells = 3, min.features = 200)
O2 <- CreateSeuratObject(counts = O2, project = "O2", min.cells = 3, min.features = 200)
O3 <- CreateSeuratObject(counts = O3, project = "O3", min.cells = 3, min.features = 200)
O4 <- CreateSeuratObject(counts = O4, project = "O4", min.cells = 3, min.features = 200)
O5 <- CreateSeuratObject(counts = O5, project = "O5", min.cells = 3, min.features = 200)
O6 <- CreateSeuratObject(counts = O6, project = "O6", min.cells = 3, min.features = 200)
Y1 <- CreateSeuratObject(counts = Y1, project = "Y1", min.cells = 3, min.features = 200)
Y2 <- CreateSeuratObject(counts = Y2, project = "Y2", min.cells = 3, min.features = 200)
Y3 <- CreateSeuratObject(counts = Y3, project = "Y3", min.cells = 3, min.features = 200)

#merge
ALLsample <- merge(O1,y = c(O2,O3,O4,O5,O6,Y1,Y2,Y3), add.cell.ids = c("O1","O2","O3","O4","O5","O6","Y1","Y2","Y3"), project = "ALLsample")

#细胞质控
#VlnPlot(ALLsample, features = c("nFeature_RNA", "nCount_RNA"), ncol =2, group.by = "orig.ident", pt.size = 0)
ALLsample <- PercentageFeatureSet(ALLsample,pattern = "^MT-",col.name = "percent_mt")
ALLsample <- subset(ALLsample, subset = nFeature_RNA > 500 & percent_mt < 50)
#标准化归一化
ALLsample <- NormalizeData(ALLsample, assay ="RNA",normalization.method = "LogNormalize", scale.factor = 10000)
ALLsample <- FindVariableFeatures(ALLsample, selection.method = "vst",nfeatures=5000)
ALLsample <- ScaleData(ALLsample)
#降维
ALLsample <- RunPCA(ALLsample,npcs = 50)
#DimPlot(object = ALLsample, dims=c(1,2),reduction = "pca")
#ElbowPlot(ALLsample, ndims = 50, reduction = "pca")
ALLsample <- RunHarmony(ALLsample,reduction.use = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#聚类
ALLsample <- FindNeighbors(ALLsample, dims = 1:20,reduction = 'harmony')
ALLsample <- FindClusters(ALLsample)
#clustree(ALLsample)
#DimPlot(object = ALLsample, reduction = "pca",dim=c(1,2),label = TRUE,group.by="orig.ident")
#可视化降维
ALLsample <- RunUMAP(ALLsample, dims = 1:20,reduction = 'harmony')
#DimPlot(object = ALLsample, reduction = "umap",dim=c(1,2),group.by="seurat_clusters",label = TRUE)
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/heart/seurat/RNAnormal_ALLsample.rds")
#marker
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,return.thresh = 0.01)
data.table::fwrite(ALLsample.markers,file ="/home/zhzhang/PG/aging/heart/seurat/ALLsample.allcluster.markergene.byRNA.txt",sep = '\t',row.names = F,quote = F,col.names = T)



```
##### lung
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)
plan("multiprocess", workers = 28)
options(future.globals.maxSize= 2500000000000)

# Load the dataset
O1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/O1/outs/filtered_feature_bc_matrix/")
O2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/O2/outs/filtered_feature_bc_matrix/")
O3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/O3/outs/filtered_feature_bc_matrix/")
O4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/O4/outs/filtered_feature_bc_matrix/")
O5 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/O5/outs/filtered_feature_bc_matrix/")
O6 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/O6/outs/filtered_feature_bc_matrix/")
Y1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/Y1/outs/filtered_feature_bc_matrix/")
Y2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/Y2/outs/filtered_feature_bc_matrix/")
Y3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/Y3/outs/filtered_feature_bc_matrix/")
Y4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/Y4/outs/filtered_feature_bc_matrix/")
Y5 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/Y5/outs/filtered_feature_bc_matrix/")
Y6 <- Read10X(data.dir = "/home/zhzhang/PG/aging/lung/count/Y6/outs/filtered_feature_bc_matrix/")
#创建对象
O1 <- CreateSeuratObject(counts = O1, project = "O1", min.cells = 3, min.features = 200)
O2 <- CreateSeuratObject(counts = O2, project = "O2", min.cells = 3, min.features = 200)
O3 <- CreateSeuratObject(counts = O3, project = "O3", min.cells = 3, min.features = 200)
O4 <- CreateSeuratObject(counts = O4, project = "O4", min.cells = 3, min.features = 200)
O5 <- CreateSeuratObject(counts = O5, project = "O5", min.cells = 3, min.features = 200)
O6 <- CreateSeuratObject(counts = O6, project = "O6", min.cells = 3, min.features = 200)
Y1 <- CreateSeuratObject(counts = Y1, project = "Y1", min.cells = 3, min.features = 200)
Y2 <- CreateSeuratObject(counts = Y2, project = "Y2", min.cells = 3, min.features = 200)
Y3 <- CreateSeuratObject(counts = Y3, project = "Y3", min.cells = 3, min.features = 200)
Y4 <- CreateSeuratObject(counts = Y4, project = "Y4", min.cells = 3, min.features = 200)
Y5 <- CreateSeuratObject(counts = Y5, project = "Y5", min.cells = 3, min.features = 200)
Y6 <- CreateSeuratObject(counts = Y6, project = "Y6", min.cells = 3, min.features = 200)
#merge
ALLsample <- merge(O1,y = c(O2,O3,O4,O5,O6,Y1,Y2,Y3,Y4,Y5,Y6), add.cell.ids = c("O1","O2","O3","O4","O5","O6","Y1","Y2","Y3","Y4","Y5","Y6"), project = "ALLsample")

#细胞质控
#VlnPlot(ALLsample, features = c("nFeature_RNA", "nCount_RNA"), ncol =2, group.by = "orig.ident", pt.size = 0)
ALLsample <- PercentageFeatureSet(ALLsample,pattern = "^MT-",col.name = "percent_mt")
ALLsample <- subset(ALLsample, subset = nFeature_RNA > 500 & percent_mt < 50)
#标准化归一化
ALLsample <- NormalizeData(ALLsample, assay ="RNA",normalization.method = "LogNormalize", scale.factor = 10000)
ALLsample <- FindVariableFeatures(ALLsample, selection.method = "vst",nfeatures=5000)
ALLsample <- ScaleData(ALLsample)
#降维
ALLsample <- RunPCA(ALLsample,npcs = 50)
#DimPlot(object = ALLsample, dims=c(1,2),reduction = "pca")
#ElbowPlot(ALLsample, ndims = 50, reduction = "pca")
ALLsample <- RunHarmony(ALLsample,reduction.use = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#聚类
ALLsample <- FindNeighbors(ALLsample, dims = 1:20,reduction = 'harmony')
ALLsample <- FindClusters(ALLsample)
#clustree(ALLsample)
#DimPlot(object = ALLsample, reduction = "pca",dim=c(1,2),label = TRUE,group.by="orig.ident")
#可视化降维
ALLsample <- RunUMAP(ALLsample, dims = 1:20,reduction = 'harmony')
#DimPlot(object = ALLsample, reduction = "umap",dim=c(1,2),group.by="seurat_clusters",label = TRUE)
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/lung/seurat/RNAnormal_ALLsample.rds")
#marker
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,return.thresh = 0.01)
data.table::fwrite(ALLsample.markers,file ="/home/zhzhang/PG/aging/lung/seurat/ALLsample.allcluster.markergene.byRNA.txt",sep = '\t',row.names = F,quote = F,col.names = T)



```
##### brain
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)
plan("multiprocess", workers = 28)
options(future.globals.maxSize= 2500000000000)
# Load the dataset
O1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/brain/count/O1/outs/filtered_feature_bc_matrix/")
O2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/brain/count/O2/outs/filtered_feature_bc_matrix/")
O3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/brain/count/O3/outs/filtered_feature_bc_matrix/")
O4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/brain/count/O4/outs/filtered_feature_bc_matrix/")
Y1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/brain/count/Y1/outs/filtered_feature_bc_matrix/")
Y2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/brain/count/Y2/outs/filtered_feature_bc_matrix/")
Y3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/brain/count/Y3/outs/filtered_feature_bc_matrix/")
#创建对象
O1 <- CreateSeuratObject(counts = O1, project = "O1", min.cells = 3, min.features = 200)
O2 <- CreateSeuratObject(counts = O2, project = "O2", min.cells = 3, min.features = 200)
O3 <- CreateSeuratObject(counts = O3, project = "O3", min.cells = 3, min.features = 200)
O4 <- CreateSeuratObject(counts = O4, project = "O4", min.cells = 3, min.features = 200)
Y1 <- CreateSeuratObject(counts = Y1, project = "Y1", min.cells = 3, min.features = 200)
Y2 <- CreateSeuratObject(counts = Y2, project = "Y2", min.cells = 3, min.features = 200)
Y3 <- CreateSeuratObject(counts = Y3, project = "Y3", min.cells = 3, min.features = 200)

#merge
ALLsample <- merge(O1,y = c(O2,O3,O4,Y1,Y2,Y3), add.cell.ids = c("O1","O2","O3","O4","Y1","Y2","Y3"), project = "ALLsample")

#细胞质控
#VlnPlot(ALLsample, features = c("nFeature_RNA", "nCount_RNA"), ncol =2, group.by = "orig.ident", pt.size = 0)
ALLsample <- PercentageFeatureSet(ALLsample,pattern = "^MT-",col.name = "percent_mt")
ALLsample <- subset(ALLsample, subset = nFeature_RNA > 500 & percent_mt < 50)
#标准化归一化
ALLsample <- NormalizeData(ALLsample, assay ="RNA",normalization.method = "LogNormalize", scale.factor = 10000)
ALLsample <- FindVariableFeatures(ALLsample, selection.method = "vst",nfeatures=5000)
ALLsample <- ScaleData(ALLsample)
#降维
ALLsample <- RunPCA(ALLsample,npcs = 50)
#DimPlot(object = ALLsample, dims=c(1,2),reduction = "pca")
#ElbowPlot(ALLsample, ndims = 50, reduction = "pca")
ALLsample <- RunHarmony(ALLsample,reduction.use = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#聚类
ALLsample <- FindNeighbors(ALLsample, dims = 1:20,reduction = 'harmony')
ALLsample <- FindClusters(ALLsample,resolution=c(0.8,0.2))
#clustree(ALLsample)
#DimPlot(object = ALLsample, reduction = "pca",dim=c(1,2),label = TRUE,group.by="orig.ident")
#可视化降维
ALLsample <- RunUMAP(ALLsample, dims = 1:20,reduction = 'harmony')
#DimPlot(object = ALLsample, reduction = "umap",dim=c(1,2),group.by="seurat_clusters",label = TRUE)
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/brain/seurat/RNAnormal_ALLsample.rds")
#marker
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5,return.thresh = 0.01)
data.table::fwrite(ALLsample.markers,file ="/home/zhzhang/PG/aging/brain/seurat/ALLsample.allcluster.markergene.byRNA.txt",sep = '\t',row.names = F,quote = F,col.names = T)






```
##### BM
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)
plan("multiprocess", workers = 28)
options(future.globals.maxSize= 2500000000000)
# Load the dataset
O1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/bonemarrow/count/O1/outs/filtered_feature_bc_matrix/")
O2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/bonemarrow/count/O2/outs/filtered_feature_bc_matrix/")
O3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/bonemarrow/count/O3/outs/filtered_feature_bc_matrix/")
Y1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/bonemarrow/count/Y1/outs/filtered_feature_bc_matrix/")
Y2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/bonemarrow/count/Y2/outs/filtered_feature_bc_matrix/")
Y3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/bonemarrow/count/Y3/outs/filtered_feature_bc_matrix/")
Y4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/bonemarrow/count/Y4/outs/filtered_feature_bc_matrix/")
Y5 <- Read10X(data.dir = "/home/zhzhang/PG/aging/bonemarrow/count/Y5/outs/filtered_feature_bc_matrix/")
#创建对象
O1 <- CreateSeuratObject(counts = O1, project = "O1", min.cells = 3, min.features = 200)
O2 <- CreateSeuratObject(counts = O2, project = "O2", min.cells = 3, min.features = 200)
O3 <- CreateSeuratObject(counts = O3, project = "O3", min.cells = 3, min.features = 200)
Y1 <- CreateSeuratObject(counts = Y1, project = "Y1", min.cells = 3, min.features = 200)
Y2 <- CreateSeuratObject(counts = Y2, project = "Y2", min.cells = 3, min.features = 200)
Y3 <- CreateSeuratObject(counts = Y3, project = "Y3", min.cells = 3, min.features = 200)
Y4 <- CreateSeuratObject(counts = Y4, project = "Y4", min.cells = 3, min.features = 200)
Y5 <- CreateSeuratObject(counts = Y5, project = "Y5", min.cells = 3, min.features = 200)
#merge
ALLsample <- merge(O1,y = c(O2,O3,Y1,Y2,Y3,Y4,Y5), add.cell.ids = c("O1","O2","O3","Y1","Y2","Y3","Y4","Y5"), project = "ALLsample")

#细胞质控
#VlnPlot(ALLsample, features = c("nFeature_RNA", "nCount_RNA"), ncol =2, group.by = "orig.ident", pt.size = 0)
ALLsample <- PercentageFeatureSet(ALLsample,pattern = "^MT-",col.name = "percent_mt")
ALLsample <- subset(ALLsample, subset = nFeature_RNA > 500 & percent_mt < 50)
#标准化归一化
ALLsample <- NormalizeData(ALLsample, assay ="RNA",normalization.method = "LogNormalize", scale.factor = 10000)
ALLsample <- FindVariableFeatures(ALLsample, selection.method = "vst",nfeatures=5000)
ALLsample <- ScaleData(ALLsample)
#降维
ALLsample <- RunPCA(ALLsample,npcs = 50)
#DimPlot(object = ALLsample, dims=c(1,2),reduction = "pca")
#ElbowPlot(ALLsample, ndims = 50, reduction = "pca")
ALLsample <- RunHarmony(ALLsample,reduction.use = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#聚类
ALLsample <- FindNeighbors(ALLsample, dims = 1:20,reduction = 'harmony')
ALLsample <- FindClusters(ALLsample)
#clustree(ALLsample)
#DimPlot(object = ALLsample, reduction = "pca",dim=c(1,2),label = TRUE,group.by="orig.ident")
#可视化降维
ALLsample <- RunUMAP(ALLsample, dims = 1:20,reduction = 'harmony')
#DimPlot(object = ALLsample, reduction = "umap",dim=c(1,2),group.by="seurat_clusters",label = TRUE)
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/bonemarrow/seurat/RNAnormal_ALLsample.rds")
#marker
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,return.thresh = 0.01)
data.table::fwrite(ALLsample.markers,file ="/home/zhzhang/PG/aging/bonemarrow/seurat/ALLsample.allcluster.markergene.byRNA.txt",sep = '\t',row.names = F,quote = F,col.names = T)



```
##### muscle
```r
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(future)
plan("multiprocess", workers = 28)
options(future.globals.maxSize= 2500000000000)

# Load the dataset
O1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O1/outs/filtered_feature_bc_matrix/")
O2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O2/outs/filtered_feature_bc_matrix/")
O3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O3/outs/filtered_feature_bc_matrix/")
O4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O4/outs/filtered_feature_bc_matrix/")
O5 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O5/outs/filtered_feature_bc_matrix/")
O6 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O6/outs/filtered_feature_bc_matrix/")
O7 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O7/outs/filtered_feature_bc_matrix/")
O8 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O8/outs/filtered_feature_bc_matrix/")
O9 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O9/outs/filtered_feature_bc_matrix/")
O10 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O10/outs/filtered_feature_bc_matrix/")
O11 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/O11/outs/filtered_feature_bc_matrix/")
Y1 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/Y1/outs/filtered_feature_bc_matrix/")
Y2 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/Y2/outs/filtered_feature_bc_matrix/")
Y3 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/Y3/outs/filtered_feature_bc_matrix/")
Y4 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/Y4/outs/filtered_feature_bc_matrix/")
Y5 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/Y5/outs/filtered_feature_bc_matrix/")
Y6 <- Read10X(data.dir = "/home/zhzhang/PG/aging/muscle/count/Y6/outs/filtered_feature_bc_matrix/")
#创建对象
O1 <- CreateSeuratObject(counts = O1, project = "O1", min.cells = 3, min.features = 200)
O2 <- CreateSeuratObject(counts = O2, project = "O2", min.cells = 3, min.features = 200)
O3 <- CreateSeuratObject(counts = O3, project = "O3", min.cells = 3, min.features = 200)
O4 <- CreateSeuratObject(counts = O4, project = "O4", min.cells = 3, min.features = 200)
O5 <- CreateSeuratObject(counts = O5, project = "O5", min.cells = 3, min.features = 200)
O6 <- CreateSeuratObject(counts = O6, project = "O6", min.cells = 3, min.features = 200)
O7 <- CreateSeuratObject(counts = O7, project = "O7", min.cells = 3, min.features = 200)
O8 <- CreateSeuratObject(counts = O8, project = "O8", min.cells = 3, min.features = 200)
O9 <- CreateSeuratObject(counts = O9, project = "O9", min.cells = 3, min.features = 200)
O10 <- CreateSeuratObject(counts = O10, project = "O10", min.cells = 3, min.features = 200)
O11 <- CreateSeuratObject(counts = O11, project = "O11", min.cells = 3, min.features = 200)
Y1 <- CreateSeuratObject(counts = Y1, project = "Y1", min.cells = 3, min.features = 200)
Y2 <- CreateSeuratObject(counts = Y2, project = "Y2", min.cells = 3, min.features = 200)
Y3 <- CreateSeuratObject(counts = Y3, project = "Y3", min.cells = 3, min.features = 200)
Y4 <- CreateSeuratObject(counts = Y4, project = "Y4", min.cells = 3, min.features = 200)
Y5 <- CreateSeuratObject(counts = Y5, project = "Y5", min.cells = 3, min.features = 200)
Y6 <- CreateSeuratObject(counts = Y6, project = "Y6", min.cells = 3, min.features = 200)
#merge
ALLsample <- merge(O1,y = c(O2,O3,O4,O5,O6,O7,O8,O9,O10,O11,Y1,Y2,Y3,Y4,Y5,Y6), add.cell.ids = c("O1","O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","Y1","Y2","Y3","Y4","Y5","Y6"), project = "ALLsample")

#细胞质控
#VlnPlot(ALLsample, features = c("nFeature_RNA", "nCount_RNA"), ncol =2, group.by = "orig.ident", pt.size = 0)
ALLsample <- PercentageFeatureSet(ALLsample,pattern = "^MT-",col.name = "percent_mt")
ALLsample <- subset(ALLsample, subset = nFeature_RNA > 500 & percent_mt < 50)
#标准化归一化
ALLsample <- NormalizeData(ALLsample, assay ="RNA",normalization.method = "LogNormalize", scale.factor = 10000)
ALLsample <- FindVariableFeatures(ALLsample, selection.method = "vst",nfeatures=5000)
ALLsample <- ScaleData(ALLsample)
#降维
ALLsample <- RunPCA(ALLsample,npcs = 50)
#DimPlot(object = ALLsample, dims=c(1,2),reduction = "pca")
#ElbowPlot(ALLsample, ndims = 50, reduction = "pca")
ALLsample <- RunHarmony(ALLsample,reduction.use = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
#聚类
ALLsample <- FindNeighbors(ALLsample, dims = 1:20,reduction = 'harmony')
ALLsample <- FindClusters(ALLsample)
#clustree(ALLsample)
#DimPlot(object = ALLsample, reduction = "pca",dim=c(1,2),label = TRUE,group.by="orig.ident")
#可视化降维
ALLsample <- RunUMAP(ALLsample, dims = 1:20,reduction = 'harmony')
#DimPlot(object = ALLsample, reduction = "umap",dim=c(1,2),group.by="seurat_clusters",label = TRUE)
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/muscle/seurat/RNAnormal_ALLsample.rds")
#marker
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1,return.thresh = 0.01)
data.table::fwrite(ALLsample.markers,file ="/home/zhzhang/PG/aging/muscle/seurat/ALLsample.allcluster.markergene.byRNA.txt",sep = '\t',row.names = F,quote = F,col.names = T)






```


#### 4.细胞注释
##### ovary
```r
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/aging/ovary/seurat/RNAnormal_ALLsample.rds /home/zhzhang/PG/aging/ovary/seurat/
scp -P 2022 "/home/zhz/GreenHub.Setup.2.2.0.exe" zhzhang@211.69.141.147:/home/zhzhang/


library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
#导入R
ALLsample <- readRDS("~/PG/aging/ovary/seurat/RNAnormal_ALLsample.rds")
DimPlot(object = ALLsample, reduction = "umap",
        group.by="seurat_clusters",label = TRUE,raster=FALSE)
#手工注释细胞类型
cluster2celltype <- c("0"="T&S",
                      "1"="SMC", 
                      "2"="T&S",
                      "3"= "T&S", 
                      "4"= "T&S", 
                      "5"= "T&S",
                      "6"= "EC", 
                      "7"= "EC", 
                      "8"= "SMC",
                      "9"= "SMC",
                      "10"= "EC",
                      "11"= "NK",
                      "12"= "SMC",
                      "13"= "MONO",
                      "14"= "T-cell",
                      "15"= "GC",
                      "16"= "NK",
                      "17"= "GC",
                      "18"= "T&S",
                      "19"= "GC",
                      "20"= "GC",
                      "21"= "EC",
                      "22"= "OO"
)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
#plot
ALLsample@meta.data[['cell_type']] <- factor(ALLsample@meta.data[['cell_type']],
                                             levels =c("GC","OO","T&S","SMC",
                                                       "EC","MONO","NK","T-cell"))
umapcell <- DimPlot(object = ALLsample, reduction = "umap",
                    group.by="cell_type",label = TRUE,raster=FALSE)+
  scale_color_manual(values = c("#ECE4B2","#654B85","#D37823",
                                "#BB4833","#2296B1","#1B608E",
                                "#ECA790","#4D8341"))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15))+
  labs(title = NULL,x="UMAP-1",y="UMAP-2")+
  theme(legend.text = element_text(size = 20),legend.key.size = unit(25,"pt"))+
  guides(color=guide_legend(override.aes = list(size=5)))
ggsave("/home/zhzhang/PG/aging/ovary/seurat/Human.ovary.celltype_umap.pdf", 
       umapcell,width = 8, height = 6.6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#储存
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/ovary/seurat/celltype_ALLsample.rds")






#颗粒细胞granulosa cells（GC: GSTA1, AMH , HSD17B1）
#卵母细胞oocytes(OO: TUBB8, ZP3 and FIGLA)
#卵泡膜和基质细胞theca and stroma cells (T&S, DCN and STAR OGN COL1A2)
#平滑肌细胞smooth muscle cells (SMC, ACTA2 and MUSTN1 MYH11)
#内皮细胞endothelial cells (EC, TM4SF1 and VWF KRT18)
#单核细胞 monocytes (MONO, TYROBP and IFI30)
#NK细胞 natural killer cells (NK, CCL5 and NKG7) 
#T细胞 (T-cell，IL7R and KLRB1) 


```
##### testis
```r
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/aging/testis/seurat/RNAnormal_ALLsample.rds /home/zhzhang/PG/aging/testis/seurat/

library(dplyr)
library(Seurat)
library(patchwork)
#导入R
ALLsample <- readRDS("/home/zhzhang/PG/aging/testis/seurat/RNAnormal_ALLsample.rds")
DimPlot(object = ALLsample, reduction = "umap",group.by="seurat_clusters",label = TRUE)
#手工注释细胞类型
cluster2celltype <- c("0"="PMC",
                      "1"="SMC", 
                      "2"="ES",
                      "3"= "PMC", 
                      "4"= "EC", 
                      "5"= "ES",
                      "6"= "Late SPC", 
                      "7"= "Leydig", 
                      "8"= "RS",
                      "9"= "ES",
                      "10"= "RS",
                      "11"= "EC",
                      "12"= "Diff.SPG",
                      "13"= "Early SPC",
                      "14"= "EC",
                      "15"= "ES",
                      "16"= "SSC",
                      "17"= "RS",
                      "18"= "Sertoli",
                      "19"= "Late SPC",
                      "20"= "Diff.SPG",
                      "21"= "Sertoli",
                      "22"= "Macro",
                      "23"= "EC",
                      "24"= "Macro",
                      "25"= "T-cell",
                      "26"= "Sertoli",
                      "27"= "EC"
)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
#plot
ALLsample@meta.data[['cell_type']] <- factor(ALLsample@meta.data[['cell_type']],
                                             levels =c("SSC","Diff.SPG","Early SPC","Late SPC",
                                                       "RS","ES",
                                                       "Sertoli","Leydig","PMC","SMC","EC",
                                                       "Macro","T-cell"))
umapcell <- DimPlot(object = ALLsample, reduction = "umap",
                    group.by="cell_type",label = TRUE,raster=FALSE)+
  scale_color_manual(values = c("#2A6434","#50A15C","#7DBB78",
                                "#ACD18E","#A6C9DD","#337CA5",
                                "#E0DC95","#984E27",
                                "#B86577","#E2A0A4","#B58F7A",
                                "#BFB5D7","#D3A665"))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15))+
  labs(title = NULL,x="UMAP-1",y="UMAP-2")+
  theme(legend.text = element_text(size = 20),legend.key.size = unit(25,"pt"))+
  guides(color=guide_legend(override.aes = list(size=5)))
ggsave("/home/zhzhang/PG/aging/testis/seurat/Human.testis.celltype_umap.pdf", 
       umapcell,width = 8, height = 6.6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#储存
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/testis/seurat/celltype_ALLsample.rds")



#支持细胞Sertoli cells（Sertoli ,CITED1，FATE1，SOX9，PRND，CLDN11，WFDC2，BEX1，BEX2）
#巨噬细胞macrophages （Macro, CD14，CD163，LYZ，C1QA,C1QC,C1QB）
#内皮细胞endothelial cells （EC, NOTCH4,CD34,VWF，TIE1，PECAM1）
#间质细胞Leydig cells（Leydig,CFD，LUM）
#管周肌样细胞peritubular cells （PMC, AEBP1，PTCH1，WFDC1）
#平滑肌细胞smooth muscle cells (SMC, NOTCH3，FABP4)
#T细胞（T-cell, CD52 ID2）

#精原干细胞spermatogonial stem cells (SSC ,PIWIL4 GFRA1 UTF1)
#分化的精原细胞differentiating spermatogonia （Diff.SPG, KIT MKI67 DCAF4L1）
#早期精母细胞early primary spermatocytes (Early SPC, DMC1，SPO11)
#晚期精母细胞late primary spermatocytes （Late SPC, PIWIL1，SPDYA，ADAM2）
#球形精子细胞round spermatids （RS, ACRV1，CATSPER3，SPACA1）
#长型精子细胞elongating/elongated spermatids (ES, PRM3,PRM2,ETNK2,TSSK6,LELP1,DCUN1D1）

```
##### skin
```r
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/aging/skin/seurat/RNAnormal_ALLsample.rds /home/zhzhang/PG/aging/skin/seurat/

library(dplyr)
library(Seurat)
library(patchwork)
#导入R
ALLsample <- readRDS("/home/zhzhang/PG/aging/skin/seurat/RNAnormal_ALLsample.rds")
DimPlot(object = ALLsample, reduction = "umap",group.by="seurat_clusters",label = TRUE)
#手工注释细胞类型
cluster2celltype <- c("0"="BC",
                      "1"="SC", 
                      "2"="VHF",
                      "3"= "SC", 
                      "4"= "SC", 
                      "5"= "SC",
                      "6"= "MC", 
                      "7"= "BC", 
                      "8"= "BC",
                      "9"= "GC",
                      "10"= "BC",
                      "11"= "FB",
                      "12"= "SC",
                      "13"= "ME",
                      "14"= "VHF",
                      "15"= "MC",
                      "16"= "IC",
                      "17"= "IC",
                      "18"= "VHF",
                      "19"= "IC",
                      "20"= "GC",
                      "21"= "PC",
                      "22"= "VHF",
                      "23"= "EC",
                      "24"= "PC"
)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
#plot
ALLsample@meta.data[['cell_type']] <- factor(ALLsample@meta.data[['cell_type']],
                                             levels =c("BC","MC","VHF","ME",
                                                       "SC","GC",
                                                       "EC","IC","FB","PC"))
umapcell <- DimPlot(object = ALLsample, reduction = "umap",
                    group.by="cell_type",label = TRUE,raster=FALSE)+
  scale_color_manual(values = c("#D657A4","#E7B22F","#CE3E1E",
                                "#73D81E","#720354","#9B1AF7",
                                "#528ACE","#304C3B",
                                "#4B7DF8","#83FDFC"))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15))+
  labs(title = NULL,x="UMAP-1",y="UMAP-2")+
  theme(legend.text = element_text(size = 20),legend.key.size = unit(25,"pt"))+
  guides(color=guide_legend(override.aes = list(size=5)))
ggsave("/home/zhzhang/PG/aging/skin/seurat/Human.skin.celltype_umap.pdf", 
       umapcell,width = 8, height = 6.6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#储存
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/skin/seurat/celltype_ALLsample.rds")





#基底细胞basal cells（BC, KRT14, KRT5, KRT15,COL17A1）
#有丝分裂细胞mitotic cell (MC, MKI67,TK1)
#绒毛毛囊细胞vellus hair follicle cell (VHF, SOX9, KRT6B, and SFRP1)
#黑色素细胞melanocyte (ME, TYR and DCT)
#棘细胞spinous cell (SC, KRT10 KRT1 KRTDAP OVOL1)
#颗粒细胞granular cell  (GC，FLG)
#内皮细胞endothelial cells (EC, CCL21 CLDN5)
#免疫细胞 immune cell (IC, PTPRC CD74)
#成纤维细胞 fibroblast (FB, PDGFRA DCN)
#周细胞pericyte (PC, RGS5 ACTA2 PDGFRB)


```
##### heart
```r
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/aging/heart/seurat/RNAnormal_ALLsample.rds /home/zhzhang/PG/aging/heart/seurat/


library(dplyr)
library(Seurat)
library(patchwork)
#导入R
ALLsample <- readRDS("/home/zhzhang/PG/aging/heart/seurat/RNAnormal_ALLsample.rds")
DimPlot(object = ALLsample, reduction = "umap",group.by="seurat_clusters",label = TRUE)
#手工注释细胞类型
cluster2celltype <- c("0"="FB",
                      "1"="EC", 
                      "2"="PC",
                      "3"= "CM", 
                      "4"= "CM", 
                      "5"= "Myeloid",
                      "6"= "Endocardial", 
                      "7"= "CM", 
                      "8"= "Myeloid",
                      "9"= "EC",
                      "10"= "CM",
                      "11"= "Myeloid",
                      "12"= "T/NK-cell",
                      "13"= "SMC",
                      "14"= "FB",
                      "15"= "CM",
                      "16"= "Myeloid",
                      "17"= "Endocardial",
                      "18"= "Myeloid",
                      "19"= "FB",
                      "20"= "LEC",
                      "21"= "CM",
                      "22"= "Neuron",
                      "23"= "CM",
                      "24"= "CM",
                      "25"="FB",
                      "26"="CM", 
                      "27"="PC",
                      "28"= "Mast", 
                      "29"= "CM", 
                      "30"= "EC",
                      "31"= "FB", 
                      "32"= "CM", 
                      "33"= "Myeloid",
                      "34"= "Epicardial",
                      "35"= "Adipo",
                      "36"= "T/NK-cell",
                      "37"= "Endocardial",
                      "38"= "Endocardial",
                      "39"= "CM",
                      "40"= "CM",
                      "41"= "Myeloid"
)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
#plot
ALLsample@meta.data[['cell_type']] <- factor(ALLsample@meta.data[['cell_type']],
                                             levels =c("CM","SMC","FB","EC",
                                                       "LEC","Neuron",
                                                       "Epicardial","Mast","Adipo","Myeloid","PC","Endocardial","T/NK-cell"))
umapcell <- DimPlot(object = ALLsample, reduction = "umap",
                    group.by="cell_type",label = TRUE,raster=FALSE)+
  scale_color_manual(values = c("#D29339","#7BC94E","#0095C7",
                                "#9C3440","#B42F70","#A0CF38",
                                "#F7E733","#EF6C33","#CF6EAC",
                                "#6EC558","#9A6DA6","#6DCFF5",
                                "#43549C"))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15))+
  labs(title = NULL,x="UMAP-1",y="UMAP-2")+
  theme(legend.text = element_text(size = 20),legend.key.size = unit(25,"pt"))+
  guides(color=guide_legend(override.aes = list(size=5)))
ggsave("/home/zhzhang/PG/aging/heart/seurat/Human.heart.celltype_umap.pdf", 
       umapcell,width = 8, height = 6.6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#储存
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/heart/seurat/celltype_ALLsample.rds")




#心肌细胞 cardiomyocytes （CM, RYR2 FGF12 MYBPC3 TNNT2）
#平滑肌细胞smooth muscle cells (SMC, MYL9 MYH11)
#成纤维细胞 fibroblast (FB, DCN LUM CCDC80 FN1 NEGR1 ABCA8 CDH19)
#内皮细胞endothelial cells (EC, VWF FLT1 EGFL7 PECAM1)
#淋巴内皮细胞lymphatic endothelial cells（LEC, CCL21 PROX1 FLT4 PDPN）
#神经元neurons (Neuron, NRXN1 NCAM1)
#心外膜细胞 epicardial cells (Epicardial, WWC1 PRG4 HAS1)
#肥大细胞mast cells (Mast, KIT CPA3 GATA2 FER)
#脂肪细胞 adipocytes (Adipo, PLIN1 DGAT2 GPD1)
#髓细胞 myeloid cells (Myeloid, CSF1R C1QC MERTK F13A1 MRC1) 
#周细胞pericyte (PC, RGS5 AGT KCNJ8 PDGFRB NOTCH3 DLC1 ABCC9)
#心内膜细胞 endocardial cells（Endocardial, NRG1 NRG3 PCDH15 CDH11 LEPR PKHD1L1 MYRIP）
#T cells and natural killer cells （T/NK-cell, T:CD3e CD2 ITK LTB;NK: KLRD1 ）



```
##### Brain
```r
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/aging/brain/seurat/RNAnormal_ALLsample.rds /home/zhzhang/PG/aging/brain/seurat/


library(dplyr)
library(Seurat)
library(patchwork)
#导入R
ALLsample <- readRDS("/home/zhzhang/PG/aging/brain/seurat/RNAnormal_ALLsample.rds")
DimPlot(object = ALLsample, reduction = "umap",group.by="seurat_clusters",label = TRUE)
#手工注释细胞类型
cluster2celltype <- c("0"="Ex",
                      "1"="Oli", 
                      "2"="Ast",
                      "3"= "EC", 
                      "4"= "In", 
                      "5"= "Ex",
                      "6"= "In", 
                      "7"= "OPC", 
                      "8"= "Mic",
                      "9"= "Ex",
                      "10"= "Ex",
                      "11"= "Ex",
                      "12"= "Ex",
                      "13"= "VLMC",
                      "14"= "Ast",
                      "15"= "Ast"
)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
#plot
ALLsample@meta.data[['cell_type']] <- factor(ALLsample@meta.data[['cell_type']],
                                             levels =c("Ex","In","Ast","Oli",
                                                       "OPC","Mic",
                                                       "EC","VLMC"))
umapcell <- DimPlot(object = ALLsample, reduction = "umap",
                    group.by="cell_type",label = TRUE,raster=FALSE)+
  scale_color_manual(values = c("#4EA666","#BC8A82","#72C081",
                                "#F29898","#9FD191","#C9DCF2",
                                "#CD97B8","#267EB6"))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15))+
  labs(title = NULL,x="UMAP-1",y="UMAP-2")+
  theme(legend.text = element_text(size = 20),legend.key.size = unit(25,"pt"))+
  guides(color=guide_legend(override.aes = list(size=5)))
ggsave("/home/zhzhang/PG/aging/brain/seurat/Human.brain.celltype_umap.pdf", 
       umapcell,width = 8, height = 6.6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#储存
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/brain/seurat/celltype_ALLsample.rds")



#Molecular and cellular evolution of the primate dorsolateral prefrontal cortex（https://www.science.org/doi/10.1126/science.abo7257#sec-11）
#Single-nucleus transcriptome analysis of human brain immune response in patients with severe COVID-19(https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00933-8#Sec14)

#兴奋性神经元excitatory neurons (Ex, SLC17A7 CAMK2A CNKSR2 SYT1 ) 
#抑制性神经元inhibitory neurons (In, GAD1 GAD2 SYT1 )
#星形胶质细胞 astrocytes (Ast, AQP4 GFAP ALDH1L1)
#少突胶质细胞 oligodendrocytes (Oli, MOG MOBP)
#少突胶质前体细胞oligodendrocyte precursor cells (OPC, PDGFRA VCAN)
#小胶质细胞microglia (Mic, P2RY12 APBB1IP)
#内皮细胞endothelial cells (EC, CLDN5 FLT1)
#血管软脑膜细胞vascular leptomeningeal cells (VLMC; COL1A2 CEMIP)






#Note that four of the non-neuronal subclasses are unique to the Ma-Sestan dataset: Immune, RB, PC, and SMC
#Immune cell, PTPRC:macrophages (F13A1+COLEC12+), myeloid cells (LSP1+LYZ+), B cells (EBF1+IGKC+) and T cells (SKAP1+CD247+)
#pericytes (PC,P2RY14 GRM8 ) 
#smooth muscle cells (SMC, ACTA2 )

```


##### lung
```r
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/aging/lung/seurat/RNAnormal_ALLsample.rds /home/zhzhang/PG/aging/lung/seurat/

library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
#导入R
ALLsample <- readRDS("/home/zhzhang/PG/aging/lung/seurat/RNAnormal_ALLsample.rds")
DimPlot(object = ALLsample, reduction = "umap",
        group.by="seurat_clusters",label = TRUE,raster=FALSE)
#手工注释细胞类型
cluster2celltype <- c("0"="Macro",
                      "1"="MONO", 
                      "2"="T-cell",
                      "3"= "AM", 
                      "4"= "AM", 
                      "5"= "AM",
                      "6"= "Macro", 
                      "7"= "AM", 
                      "8"= "MONO",
                      "9"= "NK",
                      "10"= "DC",
                      "11"= "B-cell",
                      "12"= "CEC",
                      "13"= "AT2",
                      "14"= "FB",
                      "15"= "Ciliated",
                      "16"= "Prolif",
                      "17"= "LEC",
                      "18"= "AT1",
                      "19"= "T-cell"
)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
#plot
ALLsample@meta.data[['cell_type']] <- factor(ALLsample@meta.data[['cell_type']],
                                             levels =c("AT1","AT2","Ciliated","FB","CEC","LEC",
"MONO","Macro","AM","DC","NK","T-cell","B-cell","Prolif"))
umapcell <- DimPlot(object = ALLsample, reduction = "umap",
                    group.by="cell_type",label = TRUE,raster=FALSE)+
  scale_color_manual(values = c("#D29339","#7BC94E","#0095C7",
                                "#9C3440","#B42F70","#A0CF38",
                                "#F7E733","#EF6C33","#CF6EAC",
                                "#6EC558","#9A6DA6","#6DCFF5",
                                "#43549C","#A2A2A6"))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15))+
  labs(title = NULL,x="UMAP-1",y="UMAP-2")+
  theme(legend.text = element_text(size = 20),legend.key.size = unit(25,"pt"))+
  guides(color=guide_legend(override.aes = list(size=5)))
ggsave("/home/zhzhang/PG/aging/lung/seurat/Human.lung.celltype_umap.pdf", 
       umapcell,width = 8, height = 6.6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#储存
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/lung/seurat/celltype_ALLsample.rds")




#https://www.science.org/doi/10.1126/sciadv.aba1983#acknowledgments
#A single-cell transcriptomic landscape of the lungs of patients with COVID-19 & A molecular single-cell lung atlas of lethal COVID-19 & Single-cell analysis reveals prognostic fibroblast subpopulations linked to molecular and immunological subtypes of lung cancer
上皮细胞（Epithelical cell）
#肺泡1型细胞alveolar type 1 cells （AT1:AGER EMP2 CAV1 CAV2 CLDN18 ;CLIC5 ）
#肺泡2型细胞alveolar type 2 cells （AT2:LRRK2 SFTPA1;SFTPB SFTPC ABCA3 ）
#纤毛细胞 ciliated cells (Ciliated:ANKRD66 DYDC2 C22orf15;FOXJ1 TP73 CCDC78)
基质细胞（Stromal cells）
#成纤维细胞 fibroblast (FB, PLA2G2A SCARA5 MFAP5 GAS1)
#毛细血管内皮细胞capillary endothelial cells (CEC, HPGD CA4 IL1RL1 FCN3 IL7R BTNL9 SLC6A4)
#淋巴内皮细胞lymphatic endothelial cells（LEC, NTS TFF3 SCN3B PROX1 TBX1）
免疫细胞（Immune cells）
#单核细胞 monocytes (MONO, S100A12 VCAN S100A8 S100A9 LILRB2 LILRA5 IL1R2 TLR2)
#巨噬细胞macrophages (Macro,RNASE1 PMP22 FPR3 SDS )
#肺泡巨噬细胞alveolar macrophages (AM,GPD1 PPIC AMIGO2 RBP4 RETN LSAMP)
#树突状细胞 dendritic cells (DC, SIPA1L3 CD74 HDAC9)
#NK细胞 natural killer cells (NK, S1PR5 FGFBP2 SH2D1B) 
#T细胞 (T-cell, CD3D THEMIS GPR171 ADAM19 IL7R)
#B细胞B-cells (B-cell, MS4A1 EBF1 BCL11A BLK JCHAIN)
#增殖细胞proliferative cell（Prolif，BIRC5 CDK1 HMGB2 MKI67 TOP2A）


#NO
#动脉内皮细胞 arterial endothelial cells (AEC, DKK2 HEY1 GJA5 IGFBP3 ARL15)
#静脉内皮细胞 venous endothelial cells (VEC, ACKR1 SULT1E1 HDAC9 ADGRG6 FAM155A)
#支气管周围血管内皮细胞 peribronchial vascular endothelial cells (PVEC, PLVAP SELE COL15A1 MPZL2 POSTN)
#肌成纤维细胞 Myofibroblast (MFB, BMP5 NKD2 LUM SCN7A MOXD1)
#平滑肌细胞smooth muscle cells (SMC, PLN DES TNNT2 ACTG2 ATP1A2)
#周细胞pericyte (PC, KCNK3 CDH6 COX4I2 NDUFA4L2 PDGFRB)
#棒状细胞 club cells （Club：SCGB3A2 CYP2B7P ITGA9 MGP ATL2; MUC4 CP CLIC6 SCGB3A1 SCGB1A1）
#杯状细胞No goblet cells （Goblet: BPIFB1 MUC5B LTF SCGB3A1 SCGB1A1 ; MUC5AC SPDEF）
#肺神经内分泌细胞 pulmonary neuroendocrine cells (PNEC: GRP CHGA SCG2 ASCL1 CALCA)
#肺离子细胞(MAY NOT HAVE) ionocytes (Ionocyte: ASCL3 FOXI1 ATP6V1G3 BSND  HEPACAM2)
#间皮细胞 mesothelial cells (Mesothelial :CALB2 VCAM1 MEDAG GAS1 HAS1)
#基底细胞NO basal cells （Basal: MIR205HG KRT5 EYA2 CYP24A1 KRT17 ;TP63 ）
```
##### liver
```r
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/aging/liver/seurat/RNAnormal_ALLsample.rds /home/zhzhang/PG/aging/liver/seurat/

library(dplyr)
library(Seurat)
library(patchwork)
#导入R
ALLsample <- readRDS("/home/zhzhang/PG/aging/liver/seurat/RNAnormal_ALLsample.rds")
DimPlot(object = ALLsample, reduction = "umap",group.by="seurat_clusters",label = TRUE)
#手工注释细胞类型
cluster2celltype <- c("0"="Hepatocyte",
                      "1"="Hepatocyte", 
                      "2"="Hepatocyte",
                      "3"= "Hepatocyte", 
                      "4"= "Hepatocyte", 
                      "5"= "EC",
                      "6"= "MONO", 
                      "7"= "Hepatocyte", 
                      "8"= "Stellate",
                      "9"= "T/NK-cell",
                      "10"= "T/NK-cell",
                      "11"= "EC", #cvLSEC
                      "12"= "Chol",
                      "13"= "Hepatocyte",
                      "14"= "KC",
                      "15"= "EC",
                      "16"= "EC",
                      "17"= "T/NK-cell",
                      "18"= "Hepatocyte",
                      "19"= "T/NK-cell",
                      "20"= "KC",
                      "21"= "Prolif",
                      "22"= "B-cell",
                      "23"= "MONO",
                      "24"= "Hepatocyte",
                      "25"= "Stellate",
                      "26"= "RBC",
                      "27"= "Chol",
                      "28"= "Hepatocyte",
                      "29"= "Hepatocyte"
)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
#plot
ALLsample@meta.data[['cell_type']] <- factor(ALLsample@meta.data[['cell_type']],
                                             levels =c("Hepatocyte","EC","Chol","Stellate",
"MONO","KC",
                                                       "T/NK-cell","B-cell",
                                                       "Prolif","RBC"))
umapcell <- DimPlot(object = ALLsample, reduction = "umap",
                    group.by="cell_type",label = TRUE,raster=FALSE)+
  scale_color_manual(values = c("#F9AE8F","#DE7B33","#D84902",
                                "#FFC450","#41B7C5","#2271B6",
                                "#41AE75","#FA8071",
                                "#010080","#890201"))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15))+
  labs(title = NULL,x="UMAP-1",y="UMAP-2")+
  theme(legend.text = element_text(size = 20),legend.key.size = unit(25,"pt"))+
  guides(color=guide_legend(override.aes = list(size=5)))
ggsave("/home/zhzhang/PG/aging/liver/seurat/Human.liver.celltype_umap.pdf", 
       umapcell,width = 8, height = 6.6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#储存
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/liver/seurat/celltype_ALLsample.rds")




#fig4a s1 （https://www.sciencedirect.com/science/article/pii/S0168827824000035?via=ihub）
#肝细胞hepatocyte (Hepatocyte, AGT ALB ALDOB APOA1 CYP1A2 CYP2A6 CYP2A7 ADH1A ADH1B ADH1C ADH4)
#单核细胞 monocytes (MONO, CD14 LYZ S100A4 S100A6 S100A8 S100A9)
#库普弗细胞kupffer cells（KC，CD68 VCAM1 SLC40A1 C1QA C1QB C1QC CD5L CTSB MARCO）
#内皮细胞endothelial cells （EC，C7 DNASE1L3 ENG GSN INMT LIFR PECAM1 RAMP3 TIMP3 VWF）
#胆管细胞cholangiocytes （Chol，ELF3 CFTR KRT18 KRT8 ANXA4）
#肝星形细胞stellate cells （Stellate，COL6A1 COLEC11 HGF PTH1R VIPR1）
#T cells and natural killer cells （T/NK-cell, CD247 NKG7 PARP8. T:CD3D CD3E TRAC TRBC2 IL7R LTB CD7 CTSW  NK:GNLY CD69 CMC1 GZMK KLRF1）
#B细胞（B-cell，CD79A CD79B IGHM IGHA1）
#血红细胞red blood cell (RBC, HBA1 HBA2 HBB HBD)
#增殖细胞proliferative cell（Prolif，BIRC5 CDK1 HMGB2 MKI67 TOP2A）


```
##### muscle
```r
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/aging/muscle/seurat/RNAnormal_ALLsample.rds /home/zhzhang/PG/aging/muscle/seurat/

library(dplyr)
library(Seurat)
library(patchwork)
#导入R
ALLsample <- readRDS("/home/zhzhang/PG/aging/muscle/seurat/RNAnormal_ALLsample.rds")
DimPlot(object = ALLsample, reduction = "umap",group.by="seurat_clusters",label = TRUE)
#手工注释细胞类型
cluster2celltype <- c("0"="MF-II",
                      "1"="MF-II", 
                      "2"="MF-I",
                      "3"= "MF-I", 
                      "4"= "MF-I", 
                      "5"= "FAP",
                      "6"= "MF-I", 
                      "7"= "FAP", 
                      "8"= "MF-I",
                      "9"= "MF-II",
                      "10"= "Myeloid",
                      "11"= "MF-II", 
                      "12"= "MuSC",
                      "13"= "EC",
                      "14"= "FAP",
                      "15"= "MF-I",
                      "16"= "SMC",
                      "17"= "NMJ",
                      "18"= "EC",
                      "19"= "T/NK-cell",
                      "20"= "EC",
                      "21"= "SMC",
                      "22"= "Mast",
                      "23"= "SMC",
                      "24"= "FAP",
                      "25"= "FAP"
)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
#plot
ALLsample@meta.data[['cell_type']] <- factor(ALLsample@meta.data[['cell_type']],
                                             levels =c("MF-I","MF-II","NMJ","MuSC",
"FAP","SMC","EC","Myeloid","Mast","T/NK-cell"))
umapcell <- DimPlot(object = ALLsample, reduction = "umap",
                    group.by="cell_type",label = TRUE,raster=FALSE)+
  scale_color_manual(values = c("#D29339","#7BC94E","#0095C7",
                                "#9C3440","#B42F70","#A0CF38",
                                "#F7E733","#EF6C33","#CF6EAC",
                                "#6EC558","#9A6DA6","#6DCFF5",
                                "#43549C"))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15))+
  labs(title = NULL,x="UMAP-1",y="UMAP-2")+
  theme(legend.text = element_text(size = 20),legend.key.size = unit(25,"pt"))+
  guides(color=guide_legend(override.aes = list(size=5)))
ggsave("/home/zhzhang/PG/aging/muscle/seurat/Human.muscle.celltype_umap.pdf", 
       umapcell,width = 8, height = 6.6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#储存
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/muscle/seurat/celltype_ALLsample.rds")






#https://www.nature.com/articles/s41586-024-07348-6#Sec2 (https://www.nature.com/articles/s41586-024-07348-6/figures/7)
#慢肌纤维type I myofibers (MF-I, TNNT1 MYH7 MYH7B TNNC1 TNNI1 ATP2A2) #15:TNFRSF12A DNAJA4 OTUD1 (https://www.nature.com/articles/s43587-024-00613-3/figures/10)
#快肌纤维type II myofibers (MF-II, TNNT3 MYH1 MYH2 TNNC2 TNNI2 ATP2A1)
#神经肌肉接头neuromuscular junction (NMJ, PKD2L2 GRIA2 BMPR1B EFNA5)
#肌肉干细胞muscle stem cells (MuSC, PAX7)
#纤维脂肪前体细胞fibro-adipogenic progenitors (FAP, PDGFRA PTGDS SCN7A MYOC)
#平滑肌细胞smooth muscle cells (SMC, PDGFRB ACTA2 MYH11 TAGLN MYL9 MYOCD)
#内皮细胞endothelial cells (EC, PECAM1 FLT1 CD34 CDH5 VWF JAM2)
#髓系细胞 myeloid cells (Myeloid, CSF1R MERTK F13A1 MS4A6A MRC1) 
#肥大细胞mast cells (Mast, TPSB2 MS4A2)
#T cell and natural killer cell (T/NK-cell, CD247 GNLY)



```
##### BM
```r
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/aging/bonemarrow/seurat/RNAnormal_ALLsample.rds /home/zhzhang/PG/aging/bonemarrow/seurat/

library(dplyr)
library(Seurat)
library(patchwork)
#导入R
ALLsample <- readRDS("/home/zhzhang/PG/aging/bonemarrow/seurat/RNAnormal_ALLsample.rds")
DimPlot(object = ALLsample, reduction = "umap",group.by="seurat_clusters",label = TRUE)
#手工注释细胞类型
cluster2celltype <- c("0"="HSC",
                      "1"="LMPP", 
                      "2"="MEP",
                      "3"= "pDC", 
                      "4"= "Ery", 
                      "5"= "GMP",
                      "6"= "MONO", 
                      "7"= "Ery", 
                      "8"= "Mk",
                      "9"= "Pro-B",
                      "10"= "HSC",
                      "11"= "GMP", 
                      "12"= "HSC",
                      "13"= "pDC",
                      "14"= "Basophil",
                      "15"= "HSC",
                      "16"= "CLP",
                      "17"= "pDC",
                      "18"= "CLP",
                      "19"= "Ery",
                      "20"= "Ery",
                      "21"= "T/NK-cell",
                      "22"= "HSC",
                      "23"= "CLP",
                      "24"= "T/NK-cell",
                      "25"= "EC",
                      "26"= "Ery",
                      "27"= "MONO"
)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
#plot
ALLsample@meta.data[['cell_type']] <- factor(ALLsample@meta.data[['cell_type']],
                                             levels =c("HSC","LMPP","GMP","MONO",
"pDC","CLP",
                                                       "Pro-B","T/NK-cell",
                                                       "MEP","Mk","Ery","Basophil","EC"))
umapcell <- DimPlot(object = ALLsample, reduction = "umap",
                    group.by="cell_type",label = TRUE,raster=FALSE)+
  scale_color_manual(values = c("#D29339","#7BC94E","#0095C7",
                                "#9C3440","#B42F70","#A0CF38",
                                "#F7E733","#EF6C33","#CF6EAC",
                                "#6EC558","#9A6DA6","#6DCFF5",
                                "#43549C","#A2A2A6"))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15))+
  labs(title = NULL,x="UMAP-1",y="UMAP-2")+
  theme(legend.text = element_text(size = 20),legend.key.size = unit(25,"pt"))+
  guides(color=guide_legend(override.aes = list(size=5)))
ggsave("/home/zhzhang/PG/aging/bonemarrow/seurat/Human.bonemarrow.celltype_umap.pdf", 
       umapcell,width = 8, height = 6.6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#储存
saveRDS(ALLsample, file = "/home/zhzhang/PG/aging/bonemarrow/seurat/celltype_ALLsample.rds")



#https://iiif.elifesciences.org/lax:79363%2Felife-79363-fig1-v2.tif/full/,1500/0/default.jpg
#造血干细胞hematopoietic stem cells (HSC, CRHBP AVP FOS CD164 MLLT3 MDK BST2)
#淋巴-髓系原始祖细胞lympho-myeloid primed progenitors (LMPP, SPINK2 SMIM24 CSF3R)
#粒细胞-单核细胞祖细胞 granulocyte-monocyte progenitors (GMP, CTSG PRTN3 MPO AZU1 ELANE CST7)
#单核细胞 monocytes (MONO, LYZ CSTA CST3 ANXA2)
#浆细胞样树突状细胞plasmacytoid dendritic cells (pDC, IRF7 IGKC IGLC2 SPIB STMN1 TUBB6 MYBL2 SNRNP25)
#常见淋巴祖细胞common lymphoid progenitors (CLP, IL7R DNTT)
#祖B细胞 pro-B cells (Pro-B, VPREB1 VPREB3 CD79B CD24 CD9 AKAP12 LAPTM5)
#T cells and natural killer cells (T/NK-cell, HOXA9 HOXA10 CXXC5 PRSS2 CD247 NKG7 CD3D CD3E TRAC)
#巨核细胞-红细胞祖细胞megakaryocyte-erythroid progenitors (MEP, FCER1A HPGDS CYTL1 TESPA1 GATA2)
#巨核细胞megakaryocytes (Mk, PLEK PBX1 GP1BA FLI1 VWF)
#红细胞erythroid (Ery, HBD KLF1 APOC1 BLVRB MYC HBB CA1 FAM178B PRDX2 REXO2) 
#嗜碱粒细胞basophil (Basophil, CLC TPSAB1 MS4A2 HPGD HDC)
#内皮细胞endothelial cells （EC，C7 DNASE1L3 ENG GSN RAMP3)


```










#### 5.各cell type三类基因表达比例
```r
#统一
#单细胞基因ID和基因名对照表
tail -n +2 /home/zhzhang/PG/scRNAseq/crref/grch38_108/star/geneInfo.tab|awk '{print $1"\t"$2}' > /home/zhzhang/PG/scRNAseq/crref/gene.ID.name.human.txt
scp -P 2022 zhzhang@211.69.141.147:/home/zhzhang/PG/scRNAseq/crref/gene.ID.name.human.txt /home/zhzhang/PG/scRNAseq/


#基因名genename-基因类型对照表
#R
#导入单细胞基因ID和基因名对照表
geneIDname <- read.delim("~/PG/scRNAseq/gene.ID.name.human.txt", header=FALSE)
colnames(geneIDname) <- c("geneid","genename")
#导入基因ID分类
geneid_class <- read.delim("~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt")
#合并产生基因名-基因类型对照表
geneid_he <- left_join(geneid_class,geneIDname,by="geneid")%>%
  select(3,2)%>%
  mutate(genename=stringr::str_replace_all(genename,"_","-"))%>%#因为单细胞seurat包会把基因name中的_替换成-，所以这里与seurat保持一致
  distinct(genename,.keep_all = T)
#储存
data.table::fwrite(geneid_he,file ="/home/zhzhang/PG/scRNAseq/genename.type.human.txt",sep = '\t',row.names = F,quote = F,col.names = T)


```
```r
#脚本部分
#生科院服务器，计算每个celltype细胞群中，每个基因的阳性细胞比例pct1
library(dplyr)
library(Seurat)
library(patchwork)
library(future)
#a输入的未注释rds，b保存的注释rds，c保存的每个celltype每个基因表达信息
#testis
a="/home/zhzhang/PG/aging/testis/seurat/RNAnormal_ALLsample.rds"
b="/home/zhzhang/PG/aging/testis/seurat/celltype_ALLsample.rds"
c="/home/zhzhang/PG/aging/testis/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
#ovary
a="/home/zhzhang/PG/aging/ovary/seurat/RNAnormal_ALLsample.rds"
b="/home/zhzhang/PG/aging/ovary/seurat/celltype_ALLsample.rds"
c="/home/zhzhang/PG/aging/ovary/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
#skin
a="/home/zhzhang/PG/aging/skin/seurat/RNAnormal_ALLsample.rds"
b="/home/zhzhang/PG/aging/skin/seurat/celltype_ALLsample.rds"
c="/home/zhzhang/PG/aging/skin/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
#liver
a="/home/zhzhang/PG/aging/liver/seurat/RNAnormal_ALLsample.rds"
b="/home/zhzhang/PG/aging/liver/seurat/celltype_ALLsample.rds"
c="/home/zhzhang/PG/aging/liver/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
#lung
a="/home/zhzhang/PG/aging/lung/seurat/RNAnormal_ALLsample.rds"
b="/home/zhzhang/PG/aging/lung/seurat/celltype_ALLsample.rds"
c="/home/zhzhang/PG/aging/lung/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
#heart
a="/home/zhzhang/PG/aging/heart/seurat/RNAnormal_ALLsample.rds"
b="/home/zhzhang/PG/aging/heart/seurat/celltype_ALLsample.rds"
c="/home/zhzhang/PG/aging/heart/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
#brain
a="/home/zhzhang/PG/aging/brain/seurat/RNAnormal_ALLsample.rds"
b="/home/zhzhang/PG/aging/brain/seurat/celltype_ALLsample.rds"
c="/home/zhzhang/PG/aging/brain/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
#bonemarrow
a="/home/zhzhang/PG/aging/bonemarrow/seurat/RNAnormal_ALLsample.rds"
b="/home/zhzhang/PG/aging/bonemarrow/seurat/celltype_ALLsample.rds"
c="/home/zhzhang/PG/aging/bonemarrow/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
#muscle
a="/home/zhzhang/PG/aging/muscle/seurat/RNAnormal_ALLsample.rds"
b="/home/zhzhang/PG/aging/muscle/seurat/celltype_ALLsample.rds"
c="/home/zhzhang/PG/aging/muscle/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
#
ALLsample <- readRDS(a)
#手工注释细胞类型(不同组织自行变化代码)
cluster2celltype <- c("0"="EC",
                      ......)
ALLsample@meta.data[['cell_type']] <- unname(cluster2celltype[ALLsample@meta.data$seurat_clusters])
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "cell_type"
saveRDS(ALLsample, file = b)
#
plan("multiprocess", workers = 64)
options(future.globals.maxSize= 2500000000000)
ALLsample.markers <- FindAllMarkers(object = ALLsample,assay="RNA",min.pct=0,logfc.threshold=0,only.pos = FALSE,return.thresh=1)
data.table::fwrite(ALLsample.markers,file =c,sep = '\t',row.names = F,quote = F,col.names = T)




```
```r
#对比每类细胞中三类基因的表达比例（在某细胞群中阳性细胞比例>=0.01即为在该细胞群中表达）
library(dplyr)
library(ggsignif)
library(patchwork)
library(cowplot)
#ovary
a <- "/home/zhzhang/PG/aging/ovary/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
level <- c("GC","OO","T&S","SMC","EC","MONO","NK","T-cell")
b <- "/home/zhzhang/PG/aging/ovary/seurat/genetype_celltype_expratio.tj.txt"
c <- "/home/zhzhang/PG/aging/ovary/seurat/lnctype_celltype_expratio.diffpvalue.tj.txt"
d <- "/home/zhzhang/PG/aging/ovary/seurat/genetype_celltype_expratio.pdf"
e <- "/home/zhzhang/PG/aging/ovary/seurat/celltype_expgenename_type.txt"
f <- "/home/zhzhang/PG/aging/ovary/seurat/genetype_expin_celltypenum.pdf"
#testis
a <- "/home/zhzhang/PG/aging/testis/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
level <- c("SSC","Diff.SPG","Early SPC","Late SPC",
                    "RS","ES",
                    "Sertoli","Leydig","PMC","SMC","EC",
                    "Macro","T-cell")
b <- "/home/zhzhang/PG/aging/testis/seurat/genetype_celltype_expratio.tj.txt"
c <- "/home/zhzhang/PG/aging/testis/seurat/lnctype_celltype_expratio.diffpvalue.tj.txt"
d <- "/home/zhzhang/PG/aging/testis/seurat/genetype_celltype_expratio.pdf"
e <- "/home/zhzhang/PG/aging/testis/seurat/celltype_expgenename_type.txt"
f <- "/home/zhzhang/PG/aging/testis/seurat/genetype_expin_celltypenum.pdf"
#skin
a <- "/home/zhzhang/PG/aging/skin/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
level <- c("BC","MC","VHF","ME","SC","GC","EC","IC","FB","PC")
b <- "/home/zhzhang/PG/aging/skin/seurat/genetype_celltype_expratio.tj.txt"
c <- "/home/zhzhang/PG/aging/skin/seurat/lnctype_celltype_expratio.diffpvalue.tj.txt"
d <- "/home/zhzhang/PG/aging/skin/seurat/genetype_celltype_expratio.pdf"
e <- "/home/zhzhang/PG/aging/skin/seurat/celltype_expgenename_type.txt"
f <- "/home/zhzhang/PG/aging/skin/seurat/genetype_expin_celltypenum.pdf"
#liver
a <- "/home/zhzhang/PG/aging/liver/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
level <- c("Hepatocyte","EC","Chol","Stellate","MONO","KC","T/NK-cell","B-cell","Prolif","RBC")
b <- "/home/zhzhang/PG/aging/liver/seurat/genetype_celltype_expratio.tj.txt"
c <- "/home/zhzhang/PG/aging/liver/seurat/lnctype_celltype_expratio.diffpvalue.tj.txt"
d <- "/home/zhzhang/PG/aging/liver/seurat/genetype_celltype_expratio.pdf"
e <- "/home/zhzhang/PG/aging/liver/seurat/celltype_expgenename_type.txt"
f <- "/home/zhzhang/PG/aging/liver/seurat/genetype_expin_celltypenum.pdf"
#lung
a <- "/home/zhzhang/PG/aging/lung/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
level <- c("AT1","AT2","Ciliated","FB","CEC","LEC","MONO","Macro","AM","DC","NK","T-cell","B-cell","Prolif")
b <- "/home/zhzhang/PG/aging/lung/seurat/genetype_celltype_expratio.tj.txt"
c <- "/home/zhzhang/PG/aging/lung/seurat/lnctype_celltype_expratio.diffpvalue.tj.txt"
d <- "/home/zhzhang/PG/aging/lung/seurat/genetype_celltype_expratio.pdf"
e <- "/home/zhzhang/PG/aging/lung/seurat/celltype_expgenename_type.txt"
f <- "/home/zhzhang/PG/aging/lung/seurat/genetype_expin_celltypenum.pdf"
#heart
a <- "/home/zhzhang/PG/aging/heart/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
level <- c("CM","SMC","FB","EC","LEC","Neuron","Epicardial","Mast","Adipo","Myeloid","PC","Endocardial","T/NK-cell")
b <- "/home/zhzhang/PG/aging/heart/seurat/genetype_celltype_expratio.tj.txt"
c <- "/home/zhzhang/PG/aging/heart/seurat/lnctype_celltype_expratio.diffpvalue.tj.txt"
d <- "/home/zhzhang/PG/aging/heart/seurat/genetype_celltype_expratio.pdf"
e <- "/home/zhzhang/PG/aging/heart/seurat/celltype_expgenename_type.txt"
f <- "/home/zhzhang/PG/aging/heart/seurat/genetype_expin_celltypenum.pdf"
#brain
a <- "/home/zhzhang/PG/aging/brain/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
level <- c("Ex","In","Ast","Oli","OPC","Mic","EC","VLMC")
b <- "/home/zhzhang/PG/aging/brain/seurat/genetype_celltype_expratio.tj.txt"
c <- "/home/zhzhang/PG/aging/brain/seurat/lnctype_celltype_expratio.diffpvalue.tj.txt"
d <- "/home/zhzhang/PG/aging/brain/seurat/genetype_celltype_expratio.pdf"
e <- "/home/zhzhang/PG/aging/brain/seurat/celltype_expgenename_type.txt"
f <- "/home/zhzhang/PG/aging/brain/seurat/genetype_expin_celltypenum.pdf"
#bonemarrow
a <- "/home/zhzhang/PG/aging/bonemarrow/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
level <- c("HSC","LMPP","GMP","MONO","pDC","CLP","Pro-B","T/NK-cell","MEP","Mk","Ery","Basophil","EC")
b <- "/home/zhzhang/PG/aging/bonemarrow/seurat/genetype_celltype_expratio.tj.txt"
c <- "/home/zhzhang/PG/aging/bonemarrow/seurat/lnctype_celltype_expratio.diffpvalue.tj.txt"
d <- "/home/zhzhang/PG/aging/bonemarrow/seurat/genetype_celltype_expratio.pdf"
e <- "/home/zhzhang/PG/aging/bonemarrow/seurat/celltype_expgenename_type.txt"
f <- "/home/zhzhang/PG/aging/bonemarrow/seurat/genetype_expin_celltypenum.pdf"
#muscle
a <- "/home/zhzhang/PG/aging/muscle/seurat/ALLsample.allcelltype.markergene.byRNA.txt"
level <- c("MF-I","MF-II","NMJ","MuSC","FAP","SMC","EC","Myeloid","Mast","T/NK-cell")
b <- "/home/zhzhang/PG/aging/muscle/seurat/genetype_celltype_expratio.tj.txt"
c <- "/home/zhzhang/PG/aging/muscle/seurat/lnctype_celltype_expratio.diffpvalue.tj.txt"
d <- "/home/zhzhang/PG/aging/muscle/seurat/genetype_celltype_expratio.pdf"
e <- "/home/zhzhang/PG/aging/muscle/seurat/celltype_expgenename_type.txt"
f <- "/home/zhzhang/PG/aging/muscle/seurat/genetype_expin_celltypenum.pdf"
#导入基因name分类
genename_class <- read.delim("/home/zhzhang/PG/scRNAseq/genename.type.human.txt")
#导入全部基因在各celltype的阳性细胞比例，并添加基因分类
allcelltype_genepct <- read.delim(a)%>%
  select(6,7,3)
colnames(allcelltype_genepct) <- c("celltype","genename","pct.1")
allcelltype_genepct <- left_join(allcelltype_genepct,genename_class,by="genename")
allcelltype_genepct <- allcelltype_genepct[is.na(allcelltype_genepct$type)==F,]
#基因在某细胞群中阳性细胞比例>=0.01即为在该细胞群中表达(后续储存表达的基因列表)
allcelltype_genepct <- mutate(allcelltype_genepct,exp=0)
allcelltype_genepct$exp[allcelltype_genepct$pct.1>=0.01] <- 1
#tj三类基因总数
genenum <- group_by(mutate(genename_class,num=1),type)%>%
  summarise(allnum=sum(num))
#统计三类基因表达比例
celltype_expratio <- group_by(allcelltype_genepct,celltype,type)%>%
  summarise(expnum=sum(exp))%>%
  left_join(genenum,by="type")%>%
  mutate(expratio=expnum/allnum)%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
celltype_expratio$type <- factor(celltype_expratio$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
celltype_expratio$celltype <- factor(celltype_expratio$celltype,
                                     levels = level)
data.table::fwrite(celltype_expratio,file =b,sep = '\t',row.names = F,quote = F,col.names = T)
#验证每个celltype中pglnc表达比例高于npglnc表达比例的显著性（fishertest）
celltype_expratio <- mutate(celltype_expratio,nonum=allnum-expnum)
celltype_pvalue <- data.frame("celltype"=level)
for (i in c(1:length(level))) {
  fishertest <- fisher.test(filter(celltype_expratio,celltype==celltype_pvalue$celltype[i])[c(1,4),c(3,6)])
  celltype_pvalue[i,2] <- fishertest[["p.value"]]
  fishertest <- fisher.test(filter(celltype_expratio,celltype==celltype_pvalue$celltype[i])[c(1,3),c(3,6)])
  celltype_pvalue[i,3] <- fishertest[["p.value"]]
}
colnames(celltype_pvalue) <- c("celltype","pvalue1","pvalue2")
celltype_pvalue=mutate(celltype_pvalue,anno1=case_when(pvalue1<0.001 ~ "***",
                                                      pvalue1<0.01 ~ "**",
                                                      pvalue1<0.05 ~ "*",
                                                      pvalue1>=0.05 ~ "N.S."),
                       anno2=case_when(pvalue2<0.001 ~ "***",
                                       pvalue2<0.01 ~ "**",
                                       pvalue2<0.05 ~ "*",
                                       pvalue2>=0.05 ~ "N.S."))
data.table::fwrite(celltype_pvalue,file =c,sep = '\t',row.names = F,quote = F,col.names = T)
#画图
ppcell <- ggplot(data = celltype_expratio,aes(x=celltype,y=expratio*100))+
  geom_col(aes(fill=type),width = 0.6,position = "dodge")+
  ggsignif::geom_signif(annotations=celltype_pvalue$anno1,y_position=100,tip_length = rep(100,8),
              xmin = c(0.925:length(level)),xmax = c(1.225:(length(level)+1))
              )+
  ggsignif::geom_signif(annotations=celltype_pvalue$anno2,y_position=90,tip_length = rep(100,8),
                        xmin = c(1.075:(length(level)+1)),xmax = c(1.225:(length(level)+1))
  )+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  theme_half_open()+
  scale_y_continuous(limits = c(0,100))+
  labs(x = NULL, y ="Expression proportion (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        legend.position = "none")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))
ggsave(d, 
       ppcell,width = 12, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
#输出每个celltype表达的genename及其type
ctexp_genenametype <- filter(allcelltype_genepct,exp==1)%>%
  select(1,2,4)%>%
  mutate(ctagn=paste(celltype,genename,sep = "----"))%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
data.table::fwrite(ctexp_genenametype,file =e,sep = '\t',row.names = F,quote = F,col.names = T)
#统计celltype-specify
tjcellsp <- group_by(ctexp_genenametype,genename,type)%>%
  summarise(celltypenum=n())
tjcellsp$type <- factor(tjcellsp$type,levels = rev(c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA")))
wilpvalue <- data.frame(pvalue=c(wilcox.test(filter(tjcellsp,type=="Protein-coding")$celltypenum,
                                             filter(tjcellsp,type=="PAS lncRNA")$celltypenum)[["p.value"]],
                                 wilcox.test(filter(tjcellsp,type=="PAS lncRNA")$celltypenum,
                                             filter(tjcellsp,type=="NPA lncRNA")$celltypenum)[["p.value"]],
                                 wilcox.test(filter(tjcellsp,type=="PAA lncRNA")$celltypenum,
                                             filter(tjcellsp,type=="NPA lncRNA")$celltypenum)[["p.value"]]))%>%
  mutate(str=case_when(pvalue<=0.001~"***",
                       pvalue<=0.01~"**",
                       pvalue<=0.05~"*",
                       T~"N.S."))
hs <- ggplot(data = tjcellsp,aes(x=type,y=celltypenum))+
  gghalves::geom_half_violin(width=1,side = "r",position = position_nudge(x=0.2),aes(color=type,fill=type))+
  geom_point(data = filter(tjcellsp,type=="Protein-coding"),stroke=0,
             alpha=0.8/(18264/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(tjcellsp,type=="NPA lncRNA"),stroke=0,
             alpha=0.8/(25582/952/5),size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(tjcellsp,type=="PAS lncRNA"),stroke=0,
             alpha=0.8,size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_point(data = filter(tjcellsp,type=="PAA lncRNA"),stroke=0,
             alpha=0.8,size=0.5,position = position_jitter(width=0.2),aes(color=type))+
  geom_boxplot(fatten=3,width=0.05,position = position_nudge(x=0.2),fill="white",outlier.alpha = 0)+
  scale_x_discrete(labels = c("NPA lncRNA","PAA lncRNA","PAS lncRNA","Protein-coding"))+
  scale_y_continuous(breaks = c(2,4,6,8,10,12,14,16,18,20))+
  coord_flip()+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  theme_half_open()+
  labs(x =NULL, y ="Number of cell types",fill = NULL,color=NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "none")+
  geom_signif(annotations=wilpvalue$str[c(1,3)],y_position=max(tjcellsp$celltypenum)+1,tip_length = 0,
              xmin = c(3.2,1.2),xmax = c(3.8,1.8))+
  geom_signif(annotations=wilpvalue$str[c(2)],y_position=max(tjcellsp$celltypenum)+2,tip_length = 0,
              xmin = c(1.2),xmax = c(2.8))
ggsave(f, 
       hs,width = 7, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
#### 6.组织水平比较表达
```r
#统计三类基因表达状况,tissue level
#导入基因name分类
genename_class <- read.delim("/home/zhzhang/PG/scRNAseq/genename.type.human.txt")
#tj三类基因总数
genenum <- group_by(mutate(genename_class,num=1),type)%>%
  summarise(allnum=sum(num))%>%
  data.frame()%>%
  mutate(type=case_when(type=="Protein-coding" ~ "Protein-coding",
                        type=="Non-pseudogene-associated lncRNA" ~ "NPA lncRNA",
                        type=="Pseudogene-associated sense lncRNA" ~ "PAS lncRNA",
                        type=="Pseudogene-associated antisense lncRNA" ~ "PAA lncRNA"))
#
tissuelist <- c("brain","skin","heart","muscle","bonemarrow","liver","lung","ovary","testis")
alltissue <- data.frame()
for (i in 1:length(tissuelist)) {
  #导入
  OYDEA <- read.delim(paste("/home/zhzhang/PG/aging/",tissuelist[i],"/seurat/celltype_expgenename_type.txt",sep=""))%>%
    select(2,3)%>%
    distinct(genename,type)%>%
    group_by(type)%>%
    summarise(num=n())%>%
    mutate(TISSUE=Hmisc::capitalize(tissuelist[i]))
  alltissue <- rbind(alltissue,OYDEA)
  gc()
}
#合并
he <- left_join(alltissue,genenum,by="type")
he$TISSUE[he$TISSUE=="Bonemarrow"]="Bone marrow"
he$TISSUE <- factor(he$TISSUE,levels = c("Brain","Skin", "Heart", "Muscle", "Bone marrow", "Liver"     
                                         , "Lung"  ,"Ovary"  ,"Testis"))
he$type <- factor(he$type,levels =c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
he <- arrange(he,TISSUE,type)%>%
  data.frame()%>%
  mutate(nonum=allnum-num,ratio=num*100/allnum)
#fishertest
newtl <- c("Brain","Skin", "Heart", "Muscle", "Bone marrow", "Liver", "Lung"  ,"Ovary"  ,"Testis")
tiztpv1 <- data.frame()
for (i in 1:length(newtl)) {
  fisherre <- fisher.test(filter(he,TISSUE==newtl[i] & type!="Protein-coding" & type!="PAA lncRNA")[,c("num","nonum")])
  pv <- fisherre[["p.value"]]
  fpv <- data.frame(tissue=newtl[i],pvalue=pv)
  tiztpv1 <- rbind(tiztpv1,fpv)
}
tiztpv1 <- mutate(tiztpv1,pstr=case_when(pvalue<0.001~"***",
                                         pvalue<0.01~"**",
                                         pvalue<0.05~"*",
                                         T~"N.S."))
#
tiztpv2 <- data.frame()
for (i in 1:length(newtl)) {
  fisherre <- fisher.test(filter(he,TISSUE==newtl[i] & type!="NPA lncRNA" & type!="PAA lncRNA")[,c("num","nonum")])
  pv <- fisherre[["p.value"]]
  fpv <- data.frame(tissue=newtl[i],pvalue=pv)
  tiztpv2 <- rbind(tiztpv2,fpv)
}
tiztpv2 <- mutate(tiztpv2,pstr=case_when(pvalue<0.001~"***",
                                         pvalue<0.01~"**",
                                         pvalue<0.05~"*",
                                         T~"N.S."))
#
tiztpv3 <- data.frame()
for (i in 1:length(newtl)) {
  fisherre <- fisher.test(filter(he,TISSUE==newtl[i] & type!="Protein-coding" & type!="PAS lncRNA")[,c("num","nonum")])
  pv <- fisherre[["p.value"]]
  fpv <- data.frame(tissue=newtl[i],pvalue=pv)
  tiztpv3 <- rbind(tiztpv3,fpv)
}
tiztpv3 <- mutate(tiztpv3,pstr=case_when(pvalue<0.001~"***",
                                         pvalue<0.01~"**",
                                         pvalue<0.05~"*",
                                         T~"N.S."))
#plot
ph3 <- ggplot(data =he,aes(x=TISSUE,y=ratio))+
  geom_col(width=0.6,position = "dodge",aes(color=type),fill="white")+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+  cowplot::theme_half_open()+
  scale_y_continuous(breaks = c(0,25,50,75,100))+
  ggsignif::geom_signif(annotations=tiztpv1$pstr,textsize=4,
                        y_position=c(45,40,50,40,42,37,47,47.5,83),tip_length = 0,size=0.4,
                        xmin = c(0.925:8.925),xmax = c(1.225:9.225))+
  ggsignif::geom_signif(annotations=tiztpv2$pstr,textsize=4,
                        y_position=c(80,80,80,80,80,80,80,83,95),tip_length = 0,size=0.4,
                        xmin = c(0.775:8.775),xmax = c(0.925:8.925))+
  ggsignif::geom_signif(annotations=tiztpv3$pstr,textsize=4,
                        y_position=c(37,32,42,32,34,29,39,39.5,75),tip_length = 0,size=0.4,
                        xmin = c(1.075:9.075),xmax = c(1.225:9.225))+
  labs(x =NULL, y ="Expression proportion (%)",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "none", legend.direction = "vertical")+
  theme(axis.text.x = element_text(hjust=0.5,vjust =0.5))+
  theme(axis.text.x = element_text(angle =90))
ggsave("/home/zhzhang/PG/aging/alltissue.3genetype_expratio.pdf", 
       ph3,width = 8, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")



```


#### 7\. celltype-level O/Y 中每类细胞的差异表达分析（两组比较时，在两组中celltype的num需均满足cellnum>500）
```r
#脚本部分
#生科院服务器
library(dplyr)
library(Seurat)
library(patchwork)
library(future)
#ovary
a <- "/home/zhzhang/PG/aging/ovary/seurat/celltype_ALLsample.rds"
level <- c("GC","T&S","SMC","EC","MONO","NK","T-cell")
b <- "/home/zhzhang/PG/aging/ovary/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
c <- "/home/zhzhang/PG/aging/ovary/seurat/ALLsample.group.celltype.cellnum.txt"
#testis
a <- "/home/zhzhang/PG/aging/testis/seurat/celltype_ALLsample.rds"
level <- c("SSC","Diff.SPG","Early SPC","Late SPC","RS","ES","Sertoli","Leydig","PMC","SMC","EC","Macro","T-cell")
b <- "/home/zhzhang/PG/aging/testis/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
c <- "/home/zhzhang/PG/aging/testis/seurat/ALLsample.group.celltype.cellnum.txt"
#skin
a <- "/home/zhzhang/PG/aging/skin/seurat/celltype_ALLsample.rds"
level <- c("BC","MC","VHF","ME","SC","GC","EC","IC","FB","PC")
b <- "/home/zhzhang/PG/aging/skin/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
c <- "/home/zhzhang/PG/aging/skin/seurat/ALLsample.group.celltype.cellnum.txt"
#liver
a <- "/home/zhzhang/PG/aging/liver/seurat/celltype_ALLsample.rds"
level <- c("Hepatocyte","EC","Chol","Stellate","MONO","KC","T/NK-cell","B-cell","Prolif","RBC")
b <- "/home/zhzhang/PG/aging/liver/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
c <- "/home/zhzhang/PG/aging/liver/seurat/ALLsample.group.celltype.cellnum.txt"
#lung
a <- "/home/zhzhang/PG/aging/lung/seurat/celltype_ALLsample.rds"
level <- c("AT1","AT2","Ciliated","FB","CEC","LEC","MONO","Macro","AM","DC","NK","T-cell","B-cell","Prolif")
b <- "/home/zhzhang/PG/aging/lung/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
c <- "/home/zhzhang/PG/aging/lung/seurat/ALLsample.group.celltype.cellnum.txt"
#heart
a <- "/home/zhzhang/PG/aging/heart/seurat/celltype_ALLsample.rds"
level <- c("CM","SMC","FB","EC","LEC","Neuron","Epicardial","Mast","Adipo","Myeloid","PC","Endocardial","T/NK-cell")
b <- "/home/zhzhang/PG/aging/heart/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
c <- "/home/zhzhang/PG/aging/heart/seurat/ALLsample.group.celltype.cellnum.txt"
#brain
a <- "/home/zhzhang/PG/aging/brain/seurat/celltype_ALLsample.rds"
level <- c("Ex","In","Ast","Oli","OPC","Mic","EC","VLMC")
b <- "/home/zhzhang/PG/aging/brain/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
c <- "/home/zhzhang/PG/aging/brain/seurat/ALLsample.group.celltype.cellnum.txt"
#bonemarrow
a <- "/home/zhzhang/PG/aging/bonemarrow/seurat/celltype_ALLsample.rds"
level <- c("HSC","LMPP","GMP","MONO","pDC","CLP","Pro-B","T/NK-cell","MEP","Mk","Ery","Basophil","EC")
b <- "/home/zhzhang/PG/aging/bonemarrow/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
c <- "/home/zhzhang/PG/aging/bonemarrow/seurat/ALLsample.group.celltype.cellnum.txt"
#muscle
a <- "/home/zhzhang/PG/aging/muscle/seurat/celltype_ALLsample.rds"
level <- c("MF-I","MF-II","NMJ","MuSC","FAP","SMC","EC","Myeloid","Mast","T/NK-cell")
b <- "/home/zhzhang/PG/aging/muscle/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
c <- "/home/zhzhang/PG/aging/muscle/seurat/ALLsample.group.celltype.cellnum.txt"
#input
ALLsample <- readRDS(a)
#包括分组的细胞类型注释
ALLsample@meta.data[['group_cell_type']] <- paste(do::Replace(data=ALLsample@meta.data$orig.ident,pattern=c("O.*:Old","Y.*:Young")),
                                                  ALLsample@meta.data$cell_type,
                                                  sep = "----")
#设置每个细胞分类标识的变量名
Idents(ALLsample) <- "group_cell_type"

#循环对每类细胞，进行O/Y之间的差异表达分析.
#O/Y的DEA结果储存在dataframe中
celltype <- level
OY <- data.frame()
plan("multiprocess", workers = 64)
options(future.globals.maxSize= 2500000000000)
for (i in c(1:length(level))) {
  #OY
  print(paste("now:OY",celltype[i],sep = "----"))
  difexpOY <- FindMarkers(ALLsample,assay="RNA",
                           ident.1=paste("Old",celltype[i],sep = "----"),
                           ident.2=paste("Young",celltype[i],sep = "----"),
                           only.pos = FALSE, min.pct = 0, logfc.threshold =0)%>%
    tibble::rownames_to_column("genename")%>%
    mutate(cell_type=celltype[i])
  OY <- rbind(OY,difexpOY)
}
#储存
data.table::fwrite(OY,
                   file =b,sep = '\t',row.names = F,quote = F,col.names = T)
#细胞频数
pinshu <- data.frame(table(ALLsample@meta.data$group_cell_type))
data.table::fwrite(pinshu,
                   file =c,sep = '\t',row.names = F,quote = F,col.names = T)


```


#### 8.每个celltype对比表达的几类基因，在ON/Y和OD/ON中差异表达的比例
```r
#a输入每种celltype表达的基因，b输入每个基因差异表达信息，c输出处理后的每个基因差异表达信息
#d输出差异表达比例，e输出差异表达比例差异p值
#Ovary
a <- "/home/zhzhang/PG/aging/ovary/seurat/celltype_expgenename_type.txt"
b <- "/home/zhzhang/PG/aging/ovary/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
level <- c("GC","T&S","SMC","EC","MONO","NK","T-cell")
c <- "/home/zhzhang/PG/aging/ovary/seurat/Celltype_expgene_OYDEA.txt"
d <- "/home/zhzhang/PG/aging/ovary/seurat/Celltype_expgene_OYDEA.3typegene_deratio.tj.txt"
e <- "/home/zhzhang/PG/aging/ovary/seurat/Celltype_expgene_OYDEA.2typelnc_deratio_diffpvalue.tj.txt"
Ti <- "Ovary"
p1 <- "/home/zhzhang/PG/aging/ovary/dea/Celltype_expgene.3typegene_downratio.pdf"
p2 <- "/home/zhzhang/PG/aging/ovary/dea/Celltype_expgene.3typegene_upratio.pdf"
#testis
a <- "/home/zhzhang/PG/aging/testis/seurat/celltype_expgenename_type.txt"
b <- "/home/zhzhang/PG/aging/testis/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
level <- c("SSC","Diff.SPG","Early SPC","Late SPC","RS","ES","Sertoli","Leydig","PMC","SMC","EC","Macro","T-cell")
c <- "/home/zhzhang/PG/aging/testis/seurat/Celltype_expgene_OYDEA.txt"
d <- "/home/zhzhang/PG/aging/testis/seurat/Celltype_expgene_OYDEA.3typegene_deratio.tj.txt"
e <- "/home/zhzhang/PG/aging/testis/seurat/Celltype_expgene_OYDEA.2typelnc_deratio_diffpvalue.tj.txt"
Ti <- "Testis"
p1 <- "/home/zhzhang/PG/aging/testis/dea/Celltype_expgene.3typegene_downratio.pdf"
p2 <- "/home/zhzhang/PG/aging/testis/dea/Celltype_expgene.3typegene_upratio.pdf"
#skin
a <- "/home/zhzhang/PG/aging/skin/seurat/celltype_expgenename_type.txt"
b <- "/home/zhzhang/PG/aging/skin/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
level <- c("BC","MC","VHF","ME","SC","GC","EC","IC","FB","PC")
c <- "/home/zhzhang/PG/aging/skin/seurat/Celltype_expgene_OYDEA.txt"
d <- "/home/zhzhang/PG/aging/skin/seurat/Celltype_expgene_OYDEA.3typegene_deratio.tj.txt"
e <- "/home/zhzhang/PG/aging/skin/seurat/Celltype_expgene_OYDEA.2typelnc_deratio_diffpvalue.tj.txt"
Ti <- "Skin"
p1 <- "/home/zhzhang/PG/aging/skin/dea/Celltype_expgene.3typegene_downratio.pdf"
p2 <- "/home/zhzhang/PG/aging/skin/dea/Celltype_expgene.3typegene_upratio.pdf"
#liver
a <- "/home/zhzhang/PG/aging/liver/seurat/celltype_expgenename_type.txt"
b <- "/home/zhzhang/PG/aging/liver/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
level <- c("Hepatocyte","EC","Chol","Stellate","MONO","KC","T/NK-cell","B-cell","Prolif","RBC")
c <- "/home/zhzhang/PG/aging/liver/seurat/Celltype_expgene_OYDEA.txt"
d <- "/home/zhzhang/PG/aging/liver/seurat/Celltype_expgene_OYDEA.3typegene_deratio.tj.txt"
e <- "/home/zhzhang/PG/aging/liver/seurat/Celltype_expgene_OYDEA.2typelnc_deratio_diffpvalue.tj.txt"
Ti <- "Liver"
p1 <- "/home/zhzhang/PG/aging/liver/dea/Celltype_expgene.3typegene_downratio.pdf"
p2 <- "/home/zhzhang/PG/aging/liver/dea/Celltype_expgene.3typegene_upratio.pdf"
#lung
a <- "/home/zhzhang/PG/aging/lung/seurat/celltype_expgenename_type.txt"
b <- "/home/zhzhang/PG/aging/lung/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
level <- c("AT1","AT2","Ciliated","FB","CEC","LEC","MONO","Macro","AM","DC","NK","T-cell","B-cell","Prolif")
c <- "/home/zhzhang/PG/aging/lung/seurat/Celltype_expgene_OYDEA.txt"
d <- "/home/zhzhang/PG/aging/lung/seurat/Celltype_expgene_OYDEA.3typegene_deratio.tj.txt"
e <- "/home/zhzhang/PG/aging/lung/seurat/Celltype_expgene_OYDEA.2typelnc_deratio_diffpvalue.tj.txt"
Ti <- "Lung"
p1 <- "/home/zhzhang/PG/aging/lung/dea/Celltype_expgene.3typegene_downratio.pdf"
p2 <- "/home/zhzhang/PG/aging/lung/dea/Celltype_expgene.3typegene_upratio.pdf"
#heart
a <- "/home/zhzhang/PG/aging/heart/seurat/celltype_expgenename_type.txt"
b <- "/home/zhzhang/PG/aging/heart/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
level <- c("CM","SMC","FB","EC","LEC","Neuron","Epicardial","Mast","Adipo","Myeloid","PC","Endocardial","T/NK-cell")
c <- "/home/zhzhang/PG/aging/heart/seurat/Celltype_expgene_OYDEA.txt"
d <- "/home/zhzhang/PG/aging/heart/seurat/Celltype_expgene_OYDEA.3typegene_deratio.tj.txt"
e <- "/home/zhzhang/PG/aging/heart/seurat/Celltype_expgene_OYDEA.2typelnc_deratio_diffpvalue.tj.txt"
Ti <- "Heart"
p1 <- "/home/zhzhang/PG/aging/heart/dea/Celltype_expgene.3typegene_downratio.pdf"
p2 <- "/home/zhzhang/PG/aging/heart/dea/Celltype_expgene.3typegene_upratio.pdf"
#brain
a <- "/home/zhzhang/PG/aging/brain/seurat/celltype_expgenename_type.txt"
b <- "/home/zhzhang/PG/aging/brain/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
level <- c("Ex","In","Ast","Oli","OPC","Mic","EC","VLMC")
c <- "/home/zhzhang/PG/aging/brain/seurat/Celltype_expgene_OYDEA.txt"
d <- "/home/zhzhang/PG/aging/brain/seurat/Celltype_expgene_OYDEA.3typegene_deratio.tj.txt"
e <- "/home/zhzhang/PG/aging/brain/seurat/Celltype_expgene_OYDEA.2typelnc_deratio_diffpvalue.tj.txt"
Ti <- "Brain"
p1 <- "/home/zhzhang/PG/aging/brain/dea/Celltype_expgene.3typegene_downratio.pdf"
p2 <- "/home/zhzhang/PG/aging/brain/dea/Celltype_expgene.3typegene_upratio.pdf"
#bonemarrow
a <- "/home/zhzhang/PG/aging/bonemarrow/seurat/celltype_expgenename_type.txt"
b <- "/home/zhzhang/PG/aging/bonemarrow/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
level <- c("HSC","LMPP","GMP","MONO","pDC","CLP","Pro-B","T/NK-cell","MEP","Mk","Ery","Basophil","EC")
c <- "/home/zhzhang/PG/aging/bonemarrow/seurat/Celltype_expgene_OYDEA.txt"
d <- "/home/zhzhang/PG/aging/bonemarrow/seurat/Celltype_expgene_OYDEA.3typegene_deratio.tj.txt"
e <- "/home/zhzhang/PG/aging/bonemarrow/seurat/Celltype_expgene_OYDEA.2typelnc_deratio_diffpvalue.tj.txt"
Ti <- "Bone marrow"
p1 <- "/home/zhzhang/PG/aging/bonemarrow/dea/Celltype_expgene.3typegene_downratio.pdf"
p2 <- "/home/zhzhang/PG/aging/bonemarrow/dea/Celltype_expgene.3typegene_upratio.pdf"
#muscle
a <- "/home/zhzhang/PG/aging/muscle/seurat/celltype_expgenename_type.txt"
b <- "/home/zhzhang/PG/aging/muscle/seurat/ALLsample.allcelltype.OY.DEA.byRNA.txt"
level <- c("MF-I","MF-II","NMJ","MuSC","FAP","SMC","EC","Myeloid","Mast","T/NK-cell")
c <- "/home/zhzhang/PG/aging/muscle/seurat/Celltype_expgene_OYDEA.txt"
d <- "/home/zhzhang/PG/aging/muscle/seurat/Celltype_expgene_OYDEA.3typegene_deratio.tj.txt"
e <- "/home/zhzhang/PG/aging/muscle/seurat/Celltype_expgene_OYDEA.2typelnc_deratio_diffpvalue.tj.txt"
Ti <- "Muscle"
p1 <- "/home/zhzhang/PG/aging/muscle/dea/Celltype_expgene.3typegene_downratio.pdf"
p2 <- "/home/zhzhang/PG/aging/muscle/dea/Celltype_expgene.3typegene_upratio.pdf"

#导入每个celltype表达的genname列表
cte_gnt <- read.delim(a)%>%
  select(4,3)
#导入O/Y 衰老过程中每类细胞的差异表达分析结果，并添加分类
allcelltype_ONY <- read.delim(b)%>%
  mutate(ctagn=paste(cell_type,genename,sep = "----"))%>%
  left_join(cte_gnt,by="ctagn")%>%
  select(-ctagn)
allcelltype_ONY <- allcelltype_ONY[is.na(allcelltype_ONY$type)==F,]
allcelltype_ONY$type <- factor(allcelltype_ONY$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
allcelltype_ONY$cell_type <- factor(allcelltype_ONY$cell_type,
                                    levels = level)
#筛选各celltype的O/Y差异表达基因(|log2FC|>0.25且padj<0.05 ，de赋为1)
allcelltype_ONYDEG <- mutate(allcelltype_ONY,pct=case_when(pct.1>=0.1|pct.2>=0.1~1,
                                                           T~0),
                             de=0,DEtype="NoDE")
allcelltype_ONYDEG$de[which(allcelltype_ONYDEG$pct==1&allcelltype_ONYDEG$p_val_adj<0.05 & abs(allcelltype_ONYDEG$avg_log2FC)>0.25)] <- 1
allcelltype_ONYDEG$DEtype[which(allcelltype_ONYDEG$de==1 & allcelltype_ONYDEG$avg_log2FC> 0.25)] <- "Up"
allcelltype_ONYDEG$DEtype[which(allcelltype_ONYDEG$de==1 & allcelltype_ONYDEG$avg_log2FC< -0.25)] <- "Down"
data.table::fwrite(allcelltype_ONYDEG,file =c,sep = '\t',row.names = F,quote = F,col.names = T)
#统计各celltype，表达的三类基因中，差异表达的比例
allcelltype_ONYDEG_ftj <- mutate(allcelltype_ONYDEG,
                                 up=case_when(DEtype=="Up"~1,
                                              T~0),
                                 down=case_when(DEtype=="Down"~1,
                                                T~0))
tj <- group_by(allcelltype_ONYDEG_ftj,cell_type,type)%>%
  summarise(denum=sum(de),num=n(),upnum=sum(up),downnum=sum(down))%>%
  mutate(node=num-denum,noup=num-upnum,nodown=num-downnum,
         deratio=denum*100/num,upratio=upnum*100/num,downratio=downnum*100/num)
data.table::fwrite(tj,file =d,sep = '\t',row.names = F,quote = F,col.names = T)
#验证每个celltype ON/Y中pglnc 差异表达比例高于npglnc的显著性（fisher）
ONY_celltype_pvalue <- data.frame("celltype"=level)
for (i in c(1:length(level))) {
  test <- fisher.test(filter(tj,cell_type==ONY_celltype_pvalue$celltype[i])[2:4,c(3,7)])
  ONY_celltype_pvalue[i,2] <- test[["p.value"]]
  test1 <- fisher.test(filter(tj,cell_type==ONY_celltype_pvalue$celltype[i])[2:4,c(5,8)])
  ONY_celltype_pvalue[i,3] <- test1[["p.value"]]
  test2 <- fisher.test(filter(tj,cell_type==ONY_celltype_pvalue$celltype[i])[2:4,c(6,9)])
  ONY_celltype_pvalue[i,4] <- test2[["p.value"]]
  test <- fisher.test(filter(tj,cell_type==ONY_celltype_pvalue$celltype[i])[3:4,c(3,7)])
  ONY_celltype_pvalue[i,5] <- test[["p.value"]]
  test1 <- fisher.test(filter(tj,cell_type==ONY_celltype_pvalue$celltype[i])[3:4,c(5,8)])
  ONY_celltype_pvalue[i,6] <- test1[["p.value"]]
  test2 <- fisher.test(filter(tj,cell_type==ONY_celltype_pvalue$celltype[i])[3:4,c(6,9)])
  ONY_celltype_pvalue[i,7] <- test2[["p.value"]]
}
colnames(ONY_celltype_pvalue) <- c("celltype","deratiopvalue","upratiopvalue","downratiopvalue",
                                   "deratiopvalue2","upratiopvalue2","downratiopvalue2")
ONY_celltype_pvalue <- mutate(ONY_celltype_pvalue,deratiopvalue=signif(deratiopvalue,1),
                              upratiopvalue=signif(upratiopvalue,1),
                              downratiopvalue=signif(downratiopvalue,1))%>%
  mutate(deratiostr=case_when(deratiopvalue<=0.001~"***",
                              deratiopvalue<=0.01~"**",
                              deratiopvalue<=0.05~"*",
                              T~"N.S."),
         upratiostr=case_when(upratiopvalue<=0.001~"***",
                              upratiopvalue<=0.01~"**",
                              upratiopvalue<=0.05~"*",
                              T~"N.S."),
         downratiostr=case_when(downratiopvalue<=0.001~"***",
                                downratiopvalue<=0.01~"**",
                                downratiopvalue<=0.05~"*",
                                T~"N.S."),
         deratiostr2=case_when(deratiopvalue2<=0.001~"***",
                              deratiopvalue2<=0.01~"**",
                              deratiopvalue2<=0.05~"*",
                              T~"N.S."),
         upratiostr2=case_when(upratiopvalue2<=0.001~"***",
                              upratiopvalue2<=0.01~"**",
                              upratiopvalue2<=0.05~"*",
                              T~"N.S."),
         downratiostr2=case_when(downratiopvalue2<=0.001~"***",
                                downratiopvalue2<=0.01~"**",
                                downratiopvalue2<=0.05~"*",
                                T~"N.S."))
data.table::fwrite(ONY_celltype_pvalue,file =e,sep = '\t',row.names = F,quote = F,col.names = T)

#组织内up/down ratio
#plot 差异表达比例对比
OYdown <- ggplot(data = tj,aes(x=cell_type,y=downratio))+
  geom_col(aes(fill=type),width = 0.6,position = "dodge")+
  geom_signif(annotations=ONY_celltype_pvalue$downratiostr,
              y_position=max(tj$downratio+4),tip_length = 0,
              xmin = c(0.925:length(level)),xmax = c(1.225:(length(level)+1)))+
  geom_signif(annotations=ONY_celltype_pvalue$downratiostr2,
              y_position=max(tj$downratio+2),tip_length = 0,
              xmin = c(1.075:(length(level)+1)),xmax = c(1.225:(length(level)+1)))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                     limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+  theme_half_open()+
  coord_cartesian(ylim = c(0,max(tj$downratio)+5))+
  labs(x = NULL, y ="Proportions (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        legend.position = "none")+
  theme(axis.text.x = element_text(angle =90))+
  theme(axis.text.x = element_text(vjust =0.5,hjust=0.5))+
  labs(title = paste(Ti," (Down)",sep=""))
ggsave(p1, 
       OYdown,width = 10, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")
OYup <- ggplot(data = tj,aes(x=cell_type,y=upratio))+
  geom_col(aes(fill=type),width = 0.6,position = "dodge")+
  geom_signif(annotations=ONY_celltype_pvalue$upratiostr,
              y_position=max(tj$upratio+4),tip_length = 0,
              xmin = c(0.925:length(level)),xmax = c(1.225:(length(level)+1)))+
  geom_signif(annotations=ONY_celltype_pvalue$upratiostr2,
              y_position=max(tj$upratio+2),tip_length = 0,
              xmin = c(1.075:(length(level)+1)),xmax = c(1.225:(length(level)+1)))+
  scale_fill_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+  theme_half_open()+
  theme_half_open()+
  coord_cartesian(ylim = c(0,max(tj$upratio)+5))+
  labs(x = NULL, y ="Proportions (%)",fill = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        legend.position = "none")+
  theme(axis.text.x = element_text(angle =90))+
  theme(axis.text.x = element_text(vjust =0.5,hjust=0.5))+
  labs(title = paste(Ti," (Up)",sep=""))
ggsave(p2, 
       OYup,width = 10, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
#### 9.组织水平（each tissue），比较三类基因的DE比例
```r
tissuelist <- c("brain","skin","heart","muscle","bonemarrow","liver","lung","ovary","testis")
tizt <- data.frame()
for (i in 1:length(tissuelist)) {
  #导入
  OYDEA <- read.delim(paste("/home/zhzhang/PG/aging/",tissuelist[i],"/seurat/Celltype_expgene_OYDEA.txt",sep=""))
  #统计基因在组织水平是否de
  tjzt <- group_by(OYDEA,genename,type)%>%
    summarise(denum=sum(de))%>%
    mutate(de=case_when(denum>=1~1,denum==0~0))%>%
    data.frame()%>%
    select(-denum)
  #统计组织水平上，3类基因中de基因和表达基因
  zongti <- group_by(tjzt,type)%>%
    summarise(allnum=n(),denum=sum(de))%>%
    mutate(nodenum=allnum-denum,deratio=denum*100/allnum,tissue=stringr::str_to_title(tissuelist[i]))
  tizt <- rbind(tizt,zongti)
  
}
tizt$type <- factor(tizt$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
tizt <- arrange(tizt,tissue,type)
data.table::fwrite(tizt,
                   file ="/home/zhzhang/PG/aging/alltissue.deratio.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#PLOT
tiztplot=tizt
tiztplot$tissue[tiztplot$tissue=="Bonemarrow"] <- "Bone marrow"
tiztplot$tissue <- factor(tiztplot$tissue ,
                      levels = rev( do::Replace(data=Hmisc::capitalize(tissuelist),pattern=c("Bonemarrow:Bone marrow"))))
ph <- ggplot(data =tiztplot,aes(x=tissue,y=deratio))+
  geom_point(aes(color=type),size=2)+
  geom_line(aes(group=type,color=type))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                     limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  scale_y_continuous(limits = c(0,50))+
  labs(x =NULL, y ="Proportions (%)",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 18),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")+
  coord_flip()
ggsave("/home/zhzhang/PG/aging/alltissue.deratio.pdf", 
       ph,width = 5, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#odd ratio plot
#fishertest【pas vs npa】
tiztpv <- data.frame()
for (i in 1:length(tissuelist)) {
  fisherre <- fisher.test(filter(tizt,tissue==stringr::str_to_title(tissuelist[i]) & type!="Protein-coding" & type!="PAA lncRNA")[,3:4])
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
                   file ="/home/zhzhang/PG/aging/alltissue.deratio.fisherOR.pasVSnpa.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#
tiztpvplot=tiztpv
tiztpvplot$tissue[tiztpvplot$tissue=="Bonemarrow"] <- "Bone marrow"
tiztpvplot$tissue <- factor(tiztpvplot$tissue ,
                        levels = rev( do::Replace(data=Hmisc::capitalize(tissuelist),pattern=c("Bonemarrow:Bone marrow"))))
ph1 <- ggplot(data =tiztpvplot,aes(x=tissue,y=OR))+
  geom_hline(yintercept=1,linetype=2,color="#FAA465")+
  geom_errorbar(aes(ymin=ORlow,ymax=ORhigh),width=0.1,position = position_nudge(x=-0.1))+
  geom_point(color="#7197AD",size=5,shape=18,position = position_nudge(x=-0.1))+
  geom_text(aes(label=pstr),position = position_nudge(x=0.1))+
  cowplot::theme_half_open()+
  scale_y_continuous(limits = c(0.15,max(tiztpvplot$ORhigh+0.2)))+
  labs(x =NULL, y ="Odds ratio",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 18),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")+
  coord_flip()
ggsave("/home/zhzhang/PG/aging/alltissue.deratio.fisherOR.pasVSnpa.pdf", 
       ph1,width = 4, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#fishertest【paa vs npa】
tiztpv <- data.frame()
for (i in 1:length(tissuelist)) {
  fisherre <- fisher.test(filter(tizt,tissue==stringr::str_to_title(tissuelist[i]) & type!="Protein-coding" & type!="PAS lncRNA")[,3:4])
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
                   file ="/home/zhzhang/PG/aging/alltissue.deratio.fisherOR.paaVSnpa.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#
tiztpvplot=tiztpv
tiztpvplot$tissue[tiztpvplot$tissue=="Bonemarrow"] <- "Bone marrow"
tiztpvplot$tissue <- factor(tiztpvplot$tissue ,
                            levels = rev( do::Replace(data=Hmisc::capitalize(tissuelist),pattern=c("Bonemarrow:Bone marrow"))))
ph1 <- ggplot(data =tiztpvplot,aes(x=tissue,y=OR))+
  geom_hline(yintercept=1,linetype=2,color="#FAA465")+
  geom_errorbar(aes(ymin=ORlow,ymax=ORhigh),width=0.1,position = position_nudge(x=-0.1))+
  geom_point(color="#628255",size=5,shape=18,position = position_nudge(x=-0.1))+
  geom_text(aes(label=pstr),position = position_nudge(x=0.1))+
  cowplot::theme_half_open()+
  scale_y_continuous(limits = c(0.15,max(tiztpvplot$ORhigh+0.2)))+
  labs(x =NULL, y ="Odds ratio",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 18),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")+
  coord_flip()
ggsave("/home/zhzhang/PG/aging/alltissue.deratio.fisherOR.paaVSnpa.pdf", 
       ph1,width = 4, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#fishertest【pas vs paa】
tiztpv <- data.frame()
for (i in 1:length(tissuelist)) {
  fisherre <- fisher.test(filter(tizt,tissue==stringr::str_to_title(tissuelist[i]) & type!="Protein-coding" & type!="NPA lncRNA")[,3:4])
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
                   file ="/home/zhzhang/PG/aging/alltissue.deratio.fisherOR.pasVSpaa.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#
tiztpvplot=tiztpv
tiztpvplot$tissue[tiztpvplot$tissue=="Bonemarrow"] <- "Bone marrow"
tiztpvplot$tissue <- factor(tiztpvplot$tissue ,
                            levels = rev( do::Replace(data=Hmisc::capitalize(tissuelist),pattern=c("Bonemarrow:Bone marrow"))))
ph1 <- ggplot(data =tiztpvplot,aes(x=tissue,y=OR))+
  geom_hline(yintercept=1,linetype=2,color="#628255")+
  geom_errorbar(aes(ymin=ORlow,ymax=ORhigh),width=0.1,position = position_nudge(x=-0.1))+
  geom_point(color="#7197AD",size=5,shape=18,position = position_nudge(x=-0.1))+
  geom_text(aes(label=pstr),position = position_nudge(x=0.1))+
  cowplot::theme_half_open()+
  scale_y_continuous(limits = c(0.1,max(tiztpvplot$ORhigh+0.2)))+
  labs(x =NULL, y ="Odds ratio",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 18),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")+
  coord_flip()
ggsave("/home/zhzhang/PG/aging/alltissue.deratio.fisherOR.pasVSpaa.pdf", 
       ph1,width = 4, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
#### 10\. DEG中，三类基因中组织共享性agingDEG的比例
```r
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
tjalltissue <- group_by(alltissue,genename,type)%>%
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
                   file ="/home/zhzhang/PG/aging/alltissue.tissueshared_degratio.inallDEG.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
fisher.test(tj1[c(1,3),c(3,5)])
fisher.test(tj1[c(1,2),c(3,5)])
fisher.test(tj1[c(3,4),c(3,5)])
#plot tissuesharedDE基因比例差异
cedif <- ggplot(data = tj1,aes(x=type,y=ratio))+
  geom_col(aes(fill=type),width = 0.5)+
  geom_signif(annotations=c("***","***","***"),y_position=c(80,52,46),tip_length = 0,
              xmin = c(1,2,3),xmax = c(2,4,4))+
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
ggsave("/home/zhzhang/PG/aging/alltissue.tissueshared_degratio.inallDEG.pdf", 
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
                   file ="/home/zhzhang/PG/aging/alltissue.tissuespecificgene_ratio.inallDEG.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
fisher.test(tj2[c(1,3),c(3,5)])
fisher.test(tj2[c(1,2),c(3,5)])
fisher.test(tj2[c(3,4),c(3,5)])
#plot exptissuespecific基因inDEG比例差异
cedif2 <- ggplot(data = tj2,aes(x=type,y=ratio))+
  geom_col(aes(fill=type),width = 0.5)+
  geom_signif(annotations=c("***","***","N.S."),y_position=c(12,23,21),tip_length = 0,
              xmin = c(1,2,3),xmax = c(2,4,4))+
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
ggsave("/home/zhzhang/PG/aging/alltissue.tissuespecificgene_ratio.inallDEG.pdf", 
       cedif2,width = 4, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")
#统计非组织特异性且DE的基因中，在几个组织中出现aging-related DE
tjzz1 <- group_by(fil_tjalltissue,expnum,type)%>%
  summarise(mediandenum=median(denum),meandenum=mean(denum))
tjzz1$type <- factor(tjzz1$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
tjzz1 <- arrange(tjzz1,expnum,type)
data.table::fwrite(tjzz1,
                   file ="/home/zhzhang/PG/aging/alltissue.detissuenum.allDEG.tj.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#wilcoxtest
tiztpv1 <- data.frame()
for (i in 2:9) {
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
for (i in 2:9) {
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
for (i in 2:9) {
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
                                 levels = c(1:9))
fil_tjalltissue$type <- factor(fil_tjalltissue$type,levels = c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))
tjzz1$expnum <- factor(tjzz1$expnum ,
                       levels = c(1:9))
ph3 <- ggplot(data =filter(fil_tjalltissue,expnum!=1),aes(x=expnum,y=denum))+
  geom_boxplot(width=0.6,outlier.alpha = 0,aes(color=type))+
  geom_point(data =filter(tjzz1,expnum!=1&type=="Protein-coding"),aes(x=expnum,y=meandenum),shape=23,color="#BB5A5D",position = position_nudge(x=-0.225))+
  geom_point(data =filter(tjzz1,expnum!=1&type=="PAS lncRNA"),aes(x=expnum,y=meandenum),shape=23,color="#7197AD",position = position_nudge(x=-0.075))+
  geom_point(data =filter(tjzz1,expnum!=1&type=="PAA lncRNA"),aes(x=expnum,y=meandenum),shape=23,color="#628255",position = position_nudge(x=0.075))+
  geom_point(data =filter(tjzz1,expnum!=1&type=="NPA lncRNA"),aes(x=expnum,y=meandenum),shape=23,color="#FAA465",position = position_nudge(x=0.225))+
  scale_color_manual(values=c("#BB5A5D","#7197AD","#628255","#FAA465"),
                    limits=c("Protein-coding","PAS lncRNA","PAA lncRNA","NPA lncRNA"))+
  cowplot::theme_half_open()+
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9))+
  ggsignif::geom_signif(annotations=tiztpv1$pstr,textsize=4,
                        y_position=c(2.2,2.2,2.3,2.3,3.5,4.3,6.2,9),tip_length = 0,size=0.4,
                        xmin = c(0.925:7.925),xmax = c(1.225:8.225))+
  ggsignif::geom_signif(annotations=tiztpv2$pstr,textsize=4,
                        y_position=c(3,3,3,3.5,4.5,5.3,7,9.5),tip_length = 0,size=0.4,
                        xmin = c(0.775:7.775),xmax = c(0.925:7.925))+
  ggsignif::geom_signif(annotations=tiztpv3$pstr,textsize=4,
                        y_position=c(1.4,1.4,1.5,1.5,2.2,3.5,5.2,8.2),tip_length = 0,size=0.4,
                        xmin = c(1.075:8.075),xmax = c(1.225:8.225))+
  labs(x ="Number of tissues with expression", y ="Number of tissues with\naging-related differential expression",colour = NULL)+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 17),
        legend.position = "none", legend.direction = "vertical")
ggsave("/home/zhzhang/PG/aging/alltissue.detissuenum.allDEG.pdf", 
       ph3,width = 12, height = 6,dpi=1200, units = "in", device='pdf',bg = "transparent")


```
#### 11\. JPX example
```r
#选基因
tissuelist <- c("brain","skin","heart","muscle","bonemarrow","liver","lung","ovary","testis")
alltissue <- data.frame()
for (i in 1:length(tissuelist)) {
  #导入
  OYDEA <- read.delim(paste("/home/zhzhang/PG/aging/",tissuelist[i],"/seurat/Celltype_expgene_OYDEA.txt",sep=""))
  #统计基因在特定的组织水平是否de
  tjzt <- group_by(OYDEA,genename,type)%>%
    summarise(dectnum=sum(de),expctnum=n())%>%
    mutate(tissuede=case_when(dectnum>=1~1,dectnum==0~0))%>%
    data.frame()
  alltissue <- rbind(alltissue,tjzt)
}
#统计在至少一个组织中表达的基因，表达的组织数以及差异表达的组织数
tjalltissue <- group_by(alltissue,genename,type)%>%
  summarise(exptsnum=n(),detsnum=sum(tissuede),dectnum=sum(dectnum),expctnum=sum(expctnum))%>%
  data.frame()
tjalltissue_filterpglnc <- filter(tjalltissue,type=="Pseudogene-derived lncRNA")





#三类基因，跨组织和celltype 表达变化的相关性
#计算相关性时只取两种celltype共同表达的基因
tissuelist <- c("brain","skin","heart","muscle","bonemarrow","liver","lung","ovary","testis")
alltissue <- data.frame()
for (i in 1:length(tissuelist)) {
  #导入
  OYDEA <- read.delim(paste("/home/zhzhang/PG/aging/",tissuelist[i],"/seurat/Celltype_expgene_OYDEA.txt",sep=""))%>%
    select(1,3,6,7,10,11)%>%
    mutate(TISSUE=Hmisc::capitalize(tissuelist[i]))
  #cell_type=paste(cell_type," (",Hmisc::capitalize(tissuelist[i]),")",sep="")
  alltissue <- rbind(alltissue,OYDEA)
  gc()
}
# JPX
JPX <- filter(alltissue,genename=="JPX")%>%
  filter(de==1)
JPX$TISSUE[JPX$TISSUE=="Bonemarrow"] <- "Bone marrow"
JPX$TISSUE <- factor(JPX$TISSUE,levels = c("Brain","Skin", "Heart", "Muscle", "Bone marrow", "Liver"     
                                           , "Lung"  ,"Ovary"  ,"Testis"))
JPX <- arrange(JPX,TISSUE)
JPX$cell_type <- factor(JPX$cell_type,levels = distinct(JPX,cell_type)$cell_type)

phuman <- ggplot(data = JPX,aes(x=cell_type,y=avg_log2FC))+
  geom_col(width=0.01,aes(fill=avg_log2FC))+
  geom_point(aes(size=-log10(p_val_adj),fill=avg_log2FC),shape=21,color="black")+
  geom_hline(yintercept=0,linetype=2,color="gray")+
  scale_fill_gradient2(low="#2E406E",high="#92272B",midpoint=0,mid="white",
                       breaks=c(-1,-0.5,0,0.5),limits=c(-1.1,0.5))+
  cowplot::theme_half_open()+
  scale_size(range = c(3,15),breaks = c(50,100,150,200))+
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5),
                     labels = c("-1","-0.5","0","0.5"))+
  coord_cartesian(ylim = c(-1.2, 0.5))+
  labs(x = NULL, y ="log2FC",fill = "log2FC",size="-log10(p.adj)")+
  theme(axis.title = element_text(size = 20),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.position = "right")+
  theme(axis.text.y = element_text(hjust=0.5,vjust =0.5))+
  theme(axis.text.y = element_text(angle =90))+
  theme(axis.text.x = element_text(hjust=0.5,vjust =0.5))+
  theme(axis.text.x = element_text(angle =90))
ggsave("/home/zhzhang/PG/aging/JPX.pdf", 
       phuman,width =10.5, height = 5,dpi=1200, units = "in", device='pdf',bg = "transparent")



```
#### 12\. 衰老相关lncRNA的guilt-by-association method功能富集
```r
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
humanNPA=filter(hstjalltissue,type=="NPA lncRNA" & denum>0)
#储存
data.table::fwrite(rbind(dplyr::select(humanPAS,1),
                         dplyr::select(humanPAA,1),
                         dplyr::select(humanNPA,1)),file ="/home/zhzhang/PG/aging/Homo_sapiens.agingPAlnc.id.txt",
                   sep = '\t',row.names = F,quote = F,col.names = F)

###pglnc富集的coexp基因功能富集
#导入基因分类
geneid_class <- read.delim("/home/zhzhang/PG/RNAseq/Homo_sapiens/Homo_sapiens.cbindsample.rmtestis.spexp_geneid_class.txt", header=FALSE)
colnames(geneid_class) <- c("geneid","type")
#lnc gene-coexp pcg
pglnc_pcgene <- read.delim("/home/zhzhang/PG/RNAseq/Homo_sapiens/addtarget_new/Homo_sapiens.spexp.lnc_pcgene.cordata.filter.txt", header=FALSE)%>%
  select(-3,-4)
colnames(pglnc_pcgene)=c("geneid","pcgid")
#提取GOID和description对照
library(GO.db)
goterms <- Term(GOTERM)
GOlist=as.data.frame(goterms)%>%
  dplyr::rename(Description=goterms)%>%
  rownames_to_column("GO")
#提取基因ID-GOID
k <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL")
ko= AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
         keys = k,
         columns = c("GO", "ONTOLOGY"),
         keytype="ENSEMBL")%>%
  filter(ONTOLOGY=="BP")%>%
  dplyr::select(-EVIDENCE,-ONTOLOGY)%>%
  dplyr::rename(pcgid=ENSEMBL)%>%
  left_join(GOlist,by="GO")
ko_rm=group_by(ko,GO)%>%
  summarise(num=n())%>%
  data.frame()%>%
  filter(num>=10 & num<=500)%>%
  dplyr::select(-num)%>%
  left_join(ko,by="GO")
#lncRNA gene 关联通路
lncgo=left_join(pglnc_pcgene,ko_rm,by="pcgid")%>%
  group_by(geneid,GO,Description)%>%
  summarise(num=n())%>%
  data.frame()
lncgo=lncgo[is.na(lncgo$GO)==F,]
lncgof=filter(lncgo,num>=2)
#富集注释
go2gene <- lncgof[, c(2, 1)]
go2name <- lncgof[, c(2, 3)]
#coexp富集结果
set.seed(10)
npglnc <- dplyr::select(humanNPA,1)%>%
  dplyr::sample_n(2000)
ego_npg <- clusterProfiler::enricher(npglnc$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name,pAdjustMethod = 'fdr', pvalueCutoff = 1, qvalueCutoff = 1,maxGSSize = 1500)@result%>%
  mutate(type="NPA lncRNA")
ego_pas <- clusterProfiler::enricher(humanPAS$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name,pAdjustMethod = 'fdr', pvalueCutoff = 1, qvalueCutoff = 1,maxGSSize = 1500)@result%>%
  mutate(type="PAS lncRNA")
ego_paa <- clusterProfiler::enricher(humanPAA$geneid, universe=filter(geneid_class,type!="Protein-coding")$geneid,TERM2GENE = go2gene, TERM2NAME = go2name,pAdjustMethod = 'fdr', pvalueCutoff = 1, qvalueCutoff = 1,maxGSSize = 1500)@result%>%
  mutate(type="PAA lncRNA")
#储存
data.table::fwrite(rbind(ego_pas,ego_paa,ego_npg)%>%dplyr::select(10,1:9),file ="/home/zhzhang/PG/aging/Homo_sapiens.agingPAlnc.gba_enrich.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
ego <- read.delim("/home/zhzhang/PG/aging/Homo_sapiens.agingPAlnc.gba_enrich.txt")

ego <- separate(ego,GeneRatio,c("GeneRatio1","GeneRatio2"))%>%
  separate(BgRatio,c("BgRatio1","BgRatio2"))%>%
  mutate(GeneRatio1=as.numeric(GeneRatio1),GeneRatio2=as.numeric(GeneRatio2),
         BgRatio1=as.numeric(BgRatio1),BgRatio2=as.numeric(BgRatio2))%>%
  mutate(ES=(GeneRatio1/GeneRatio2)/(BgRatio1/BgRatio2))
ego$type <- factor(ego$type,levels=rev(c("NPA lncRNA","PAA lncRNA",
                                     "PAS lncRNA")))
#
respect=c("telomere maintenance",
          "telomere capping",
          "cellular response to DNA damage stimulus",
          "double-strand break repair",
          "mismatch repair",
          "nucleotide-excision repair",
          "cellular senescence",
          "autophagy",
          "response to oxidative stress",
          "stem cell population maintenance",
          "Wnt signaling pathway",
          "TOR signaling",
          "I-kappaB kinase/NF-kappaB signaling",
          "regulation of inflammatory response",
          "immune system process")
#富集图
egoplot <- filter(ego,Description %in% respect)%>%
  mutate(Description=Hmisc::capitalize(Description))
egoplot$Description <- factor(egoplot$Description,levels = Hmisc::capitalize(respect))
pp1 <- ggplot(data = egoplot,aes(x=type,y=Description))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:14.5),color="gray")+
  geom_vline(xintercept=c(1.5:2.5),color="gray")+
  scale_fill_gradient2(low="#2E406E",high="#C12039",midpoint=-log10(0.05),mid="white")+
  scale_size(breaks = c(0.5,1,1.5,2,2.5,3),range = c(1,15))+
  cowplot::theme_half_open()+
  labs(x = NULL, y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave("/home/zhzhang/PG/aging/Homo_sapiens.agingPAlnc.gba_enrich.pdf", 
       pp1,width = 8, height = 9,dpi=1200, units = "in", device='pdf',bg = "transparent")




```
#### 13\. 三类lncRNA基因启动子TF富集
```r
#全部基因的近端/远端启动子，与TFBS交集.获取全部基因两类启动子具有的TFBS情况，传至实验室服务器/home/zhzhang/PG/TFBS/
#人
bedtools sort -i /home/zhzhang/PG/TFBS/Homo_sapiens/allbed/Homo_sapiens.allTFBS.bed|sed 's/^chr//g' > /home/zhzhang/PG/TFBS/Homo_sapiens/intersect/Homo_sapiens.allTFBS.sort.bed
bedtools sort -i /share/home/zhzhang24/PG/epigen/promoterbed/Homo_sapiens/Homo_sapiens.all.Proximalpromoter.bed > /share/home/zhzhang24/PG/TFBS/Homo_sapiens/proximal_intersect/Homo_sapiens.all.Proximalpromoter.sort.bed
#近端
bedtools intersect -a /share/home/zhzhang24/PG/TFBS/Homo_sapiens/proximal_intersect/Homo_sapiens.all.Proximalpromoter.sort.bed -b /share/home/zhzhang24/PG/TFBS/Homo_sapiens/intersect/Homo_sapiens.allTFBS.sort.bed -wo > /share/home/zhzhang24/PG/TFBS/Homo_sapiens/proximal_intersect/Homo_sapiens.proximalpromoter_intersect_TFBS.txt


#脚本跑
#函数统计每个基因启动子&基因间区具有的TFBS种类及对应结合比例的信息文件
#a输入每个基因启动子与TFBS交集文件，b输入每个基因间区与TFBS交集文件，
#c输出每个基因启动子以及基因间区拥有的TFBS种类及其结合长度比例统计结果
genepromoterTFBS <- function(a,c){
  #导入启动子与TFBS交集,获取每个启动子中每种TFBS的结合长度比例
  allgenepromoter_intersect_TFBS <- data.table::fread(a, header=FALSE)%>%
    data.frame()%>%
    select(4,10,11)%>%
    separate(V10,c("TF"),remove = T,sep = "_")%>%
    group_by(V4,TF)%>%
    summarise(len=sum(V11))%>%
    mutate(len=len/1000)
  colnames(allgenepromoter_intersect_TFBS) <- c("geneid","TFBS","num")
  #储存
  data.table::fwrite(allgenepromoter_intersect_TFBS,file =c,sep = '\t',row.names = F,quote = F,col.names = T)
}
#人近端启动子
genepromoterTFBS(a="/share/home/zhzhang24/PG/TFBS/Homo_sapiens/proximal_intersect/Homo_sapiens.proximalpromoter_intersect_TFBS.txt",
                 c="/share/home/zhzhang24/PG/TFBS/Homo_sapiens/proximal_intersect/Homo_sapiens.proximalpromoterAintergenic.TFBS.txt")


#传至实验室服务器/home/zhzhang/PG/TFBS/
rsync -P -u -r -e "ssh -p 5348" /share/home/zhzhang24/PG/TFBS/Homo_sapiens/proximal_intersect/Homo_sapiens.proximalpromoterAintergenic.TFBS.txt zhzhang@122.205.95.67:/home/zhzhang/PG/TFBS/Homo_sapiens/


```


```r
#导入每个基因id的分类
geneid_class <- read.delim("~/PG/RNAseq/Homo_sapiens/Homo_sapiens.geneid_class.txt")
#富集注释
allgeneRBP <- read.delim("~/PG/TFBS/Homo_sapiens/Homo_sapiens.proximalpromoterAintergenic.TFBS.txt")%>%
  select(-3)
anno <- mutate(allgeneRBP,GO_Description=TFBS)
colnames(anno)[2] <- "GO"
go2gene <- anno[, c(2, 1)]
go2name <- anno[, c(2, 3)]
#富集结果
set.seed(10)
npalnc <- filter(geneid_class,type=="Non-pseudogene-associated lncRNA")%>%
  select(1)%>%
  dplyr::sample_n(2000)
paslnc <- filter(geneid_class,type=="Pseudogene-associated sense lncRNA")%>%
  select(1)
paalnc <- filter(geneid_class,type=="Pseudogene-associated antisense lncRNA")%>%
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
  select(10,1:9)%>%
  select(-Description)
#导入ID和name对照表
idmapping <- read.delim("~/PG/TFBS/Homo_sapiens/idmapping.tsv")%>%
  filter(To!="SMN2")
colnames(idmapping)=c("ID","Description")
ego=left_join(ego,idmapping,by="ID")%>%
  select(1:2,10,3:9)
#储存
data.table::fwrite(ego,file ="/home/zhzhang/PG/TFBS/Homo_sapiens/Human.lncRNA.type.TFBS_enrich.txt",
                   sep = '\t',row.names = F,quote = F,col.names = T)
#plot
ego <- read.delim("/home/zhzhang/PG/TFBS/Homo_sapiens/Human.lncRNA.type.TFBS_enrich.txt")
#
ego <- separate(ego,GeneRatio,c("GeneRatio1","GeneRatio2"))%>%
  separate(BgRatio,c("BgRatio1","BgRatio2"))%>%
  mutate(GeneRatio1=as.numeric(GeneRatio1),GeneRatio2=as.numeric(GeneRatio2),
         BgRatio1=as.numeric(BgRatio1),BgRatio2=as.numeric(BgRatio2))%>%
  mutate(ES=(GeneRatio1/GeneRatio2)/(BgRatio1/BgRatio2))
ego$type <- factor(ego$type,levels=rev(c("NPA lncRNA","PAA lncRNA",
                                     "PAS lncRNA")))

#
respect=c("NRF1",#维持细胞氧化还原稳态，抵御炎性衰老【Nrf1 is an indispensable redox-determining factor for mitochondrial homeostasis by integrating multi-hierarchical regulatory networks，Regulation and Functions of the ER-Associated Nrf1 Transcription Factor】
          "E2F6",#介导细胞衰老相关的细胞周期阻滞【SenNet recommendations for detecting senescent cells in different tissues】
          "E2F1",
          "E2F4",
          "E2F3",
          "RB1",#介导细胞衰老相关的细胞周期阻滞【SenNet recommendations for detecting senescent cells in different tissues】
          "PML",#与细胞衰老相关,介导细胞衰老相关的细胞周期阻滞【Regulation of E2Fs and senescence by PML nuclear bodies
          "FOXO3",#Longevity Factor /geroprotective gene老年保护基因[Pharmaceutical and nutraceutical activation of FOXO3 for healthy longevity,Single-nucleus profiling unveils a geroprotective role of the FOXO3 in primate skeletal muscle aging
          "SIRT6",#Longevity Factor 【Emerging roles of SIRT6 in human diseases and its modulators
          "MTOR",#mTOR is a key modulator of ageing and age-related disease
          "SOX5",#rejuvenation factor缓解细胞衰老【Genome-wide CRISPR activation screening in senescent cells reveals SOX5 as a driver and therapeutic target of rejuvenation
          "SOX9",#与细胞衰老相关【Sox9 Accelerates Vascular Aging by Regulating Extracellular Matrix Composition and Stiffness
          "SP1",#Sp1 蛋白的水平随着年龄的增长而降低，与细胞衰老相关【DNA damage-induced degradation of Sp1 promotes cellular senescence
          "ZBTB7A",#putative aging modulators[The Gene-Regulatory Footprint of Aging Highlights Conserved Central Regulators] 
          "CXXC1","SREBF1","SREBF2"#【A single-nucleus transcriptomic atlas of primate liver aging uncovers the pro-senescence role of SREBP2 in hepatocytes
          )

#富集图
egoplot <- filter(ego,Description %in% respect)
egoplot$Description <- factor(egoplot$Description,levels = rev(respect))
pp1 <- ggplot(data = egoplot,aes(x=type,y=Description))+
  geom_point(aes(fill=-log10(p.adjust),size=ES),shape=21,color="black")+
  geom_hline(yintercept=c(1.5:16.5),color="gray")+
  geom_vline(xintercept=c(1.5:2.5),color="gray")+
  scale_fill_gradient(low="white",high="#C12039")+
  #scale_fill_gradient2(low="#2E406E",high="#C12039",midpoint=-log10(0.05),mid="white")+
  scale_size(breaks = c(1,1.2,1.4,1.6),range = c(3,15))+
  cowplot::theme_half_open()+
  labs(x = NULL, y =NULL,size="Enrichment score")+
  theme(axis.text.x = element_text(angle =45))+
  theme(axis.text.x = element_text(vjust = 1,hjust=1))+
  theme(axis.title = element_text(size = 15),axis.text.y  = element_text(size = 15),
        axis.text.x = element_text(size = 14),legend.text = element_text(size = 14)) +
  theme(panel.background = element_rect(colour = "black",size=1))
ggsave("/home/zhzhang/PG/TFBS/Homo_sapiens/Human.lncRNA.type.TFBS_enrich.pdf", 
       pp1,width = 6, height = 10,dpi=1200, units = "in", device='pdf',bg = "transparent")





```



