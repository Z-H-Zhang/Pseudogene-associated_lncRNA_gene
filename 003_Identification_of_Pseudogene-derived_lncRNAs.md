### 三.Pseudogene-derived lncRNA 鉴定
```r
###至少有一个转录本与PG显著交集(全部外显子与PG交集>=200bp 在相同链上)的lncRNA gene定义为Pseudogene-derived，任意转录本无交集为Non-pseudogene-derived。有交集但不显著的lncRNA需要在后续分析中剔除（1.同链交集＜200bp，2.异链交集）
#斑马鱼
#提取ref lncRNA基因外显子bed(GTF注释文件中，exon所属的基因biotype有的在22列，有的在24列.lncRNA gene类型包括："antisense" "lincRNA" "sense_intronic" "sense_overlapping")【chr start end gene_id transcript_id strand】
#查看ref lncRNA类型：cat "/home/zhzhang/PG/HPG/GTF_PG/Danio_rerio.GRCz11.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}' | awk '$21=="gene_biotype" {print $22} $23=="gene_biotype" {print $24}'|sort|uniq
cat "/home/zhzhang/PG/HPG/GTF_PG/Danio_rerio.GRCz11.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}'|sed "s/\"//g;s/\;//g"| awk '$22=="antisense" || $22=="lincRNA" || $22=="sense_intronic" || $22=="sense_overlapping" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $24=="antisense" || $24=="lincRNA" || $24=="sense_intronic" || $24=="sense_overlapping" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}' > /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.ref_lncRNA_exon.bed
#提取二代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】
cat /home/zhzhang/PG/2G_lnc/Danio_rerio.2G_lncRNA.gtf |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.2Gnovel_lncRNA_exon.bed
#提取三代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】
cat "/home/zhzhang/PG/ISOseqdata/Danio_rerio/novel_lnc/Danio_rerio.3G_novel_lncRNA.gtf"|awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.3Gnovel_lncRNA_exon.bed
#合并全部lncRNA gene的外显子注释
cat /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.ref_lncRNA_exon.bed /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.3Gnovel_lncRNA_exon.bed /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.2Gnovel_lncRNA_exon.bed > /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.all_lncRNA_exon.bed
#获取全部 lncRNA外显子与假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/lncRNA_class/）
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Danio_rerio/pgenes/Danio_rerio_hpg.bed -wo > /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.alllncRNAexon_intersect_pg.bed



#鸡
#提取ref lncRNA基因外显子bed(GTF注释文件中，exon所属的基因biotype有的在22列，有的在24列.lncRNA gene类型包括："lncRNA" )【chr start end gene_id transcript_id strand】
#查看ref lncRNA类型：cat "/home/zhzhang/PG/HPG/GTF_PG/Gallus_gallus.GRCg7b.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}' | awk '$21=="gene_biotype" {print $22} $23=="gene_biotype" {print $24}'|sort|uniq
cat "/home/zhzhang/PG/HPG/GTF_PG/Gallus_gallus.GRCg7b.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}'|sed "s/\"//g;s/\;//g"| awk '$22=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $24=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}' > /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.ref_lncRNA_exon.bed
#提取三代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】（人记得改列号，gtf格式不同）
cat "/home/zhzhang/PG/ISOseqdata/Gallus_gallus/novel_lnc/Gallus_gallus.3G_novel_lncRNA.gtf"|awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.3Gnovel_lncRNA_exon.bed
#提取二代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】
cat /home/zhzhang/PG/2G_lnc/Gallus_gallus.2G_lncRNA.gtf |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.2Gnovel_lncRNA_exon.bed
#合并全部lncRNA gene的外显子注释
cp /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.ref_lncRNA_exon.bed /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.all_lncRNA_exon.bed
cat /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.3Gnovel_lncRNA_exon.bed >> /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.all_lncRNA_exon.bed
cat /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.2Gnovel_lncRNA_exon.bed >> /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.all_lncRNA_exon.bed
#获取全部 lncRNA外显子与假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/lncRNA_class/）
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Gallus_gallus/pgenes/Gallus_gallus_hpg.bed -wo > /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.alllncRNAexon_intersect_pg.bed



#小鼠
#提取ref lncRNA基因外显子bed(GTF注释文件中，exon所属的基因biotype有的在22列，有的在24列.lncRNA gene类型包括："lncRNA" )【chr start end gene_id transcript_id strand】
#查看ref lncRNA类型：cat "/home/zhzhang/PG/HPG/GTF_PG/Mus_musculus.GRCm39.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}' | awk '$21=="gene_biotype" {print $22} $23=="gene_biotype" {print $24}'|sort|uniq
cat "/home/zhzhang/PG/HPG/GTF_PG/Mus_musculus.GRCm39.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}'|sed "s/\"//g;s/\;//g"| awk '$22=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $24=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}' > /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.ref_lncRNA_exon.bed
#提取三代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】（人记得改列号，gtf格式不同）
cat "/home/zhzhang/PG/ISOseqdata/Mus_musculus/novel_lnc/Mus_musculus.3G_novel_lncRNA.gtf"|awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$12"\t"$10"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.3Gnovel_lncRNA_exon.bed
#提取二代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】
cat /home/zhzhang/PG/2G_lnc/Mus_musculus.2G_lncRNA.gtf |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.2Gnovel_lncRNA_exon.bed
#合并全部lncRNA gene的外显子注释
cp /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.ref_lncRNA_exon.bed /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed
cat /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.3Gnovel_lncRNA_exon.bed >> /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed
cat /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.2Gnovel_lncRNA_exon.bed >> /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed
#获取全部 lncRNA外显子与假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/lncRNA_class/）
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Mus_musculus/pgenes/Mus_musculus_hpg.bed -wo > /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.alllncRNAexon_intersect_pg.bed

   


#人Homo_sapiens
#提取ref lncRNA基因外显子bed(GTF注释文件中，exon所属的基因biotype有的在22列，有的在24列.lncRNA gene类型包括："lncRNA" )【chr start end gene_id transcript_id strand】
#查看ref lncRNA类型：cat "/home/zhzhang/PG/HPG/GTF_PG/Homo_sapiens.GRCh38.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}' | awk '$21=="gene_biotype" {print $22} $23=="gene_biotype" {print $24}'|sort|uniq
cat "/home/zhzhang/PG/HPG/GTF_PG/Homo_sapiens.GRCh38.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}'|sed "s/\"//g;s/\;//g"| awk '$22=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $24=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}' > /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.ref_lncRNA_exon.bed
#提取三代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】（人记得改列号，gtf格式不同）
cat "/home/zhzhang/PG/ISOseqdata/Homo_sapiens/novel_lnc/Homo_sapiens.3G_novel_lncRNA.gtf"|awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.3Gnovel_lncRNA_exon.bed
#提取二代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】
cat /home/zhzhang/PG/2G_lnc/Homo_sapiens.2G_lncRNA.gtf |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.2Gnovel_lncRNA_exon.bed
#合并全部lncRNA gene的外显子注释
cp /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.ref_lncRNA_exon.bed /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.all_lncRNA_exon.bed
cat /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.3Gnovel_lncRNA_exon.bed >> /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.all_lncRNA_exon.bed
cat /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.2Gnovel_lncRNA_exon.bed >> /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.all_lncRNA_exon.bed
#获取全部lncRNA外显子与假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/lncRNA_class/）
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Homo_sapiens/pgenes/Homo_sapiens_hpg.bed -wo > /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.alllncRNAexon_intersect_pg.bed



#猕猴
#提取ref lncRNA基因外显子bed(GTF注释文件中，exon所属的基因biotype有的在22列，有的在24列.lncRNA gene类型包括："lncRNA" )【chr start end gene_id transcript_id strand】
#查看ref lncRNA类型：cat "/home/zhzhang/PG/HPG/GTF_PG/Macaca_mulatta.Mmul_10.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}' | awk '$21=="gene_biotype" {print $22} $23=="gene_biotype" {print $24}'|sort|uniq
cat "/home/zhzhang/PG/HPG/GTF_PG/Macaca_mulatta.Mmul_10.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}'|sed "s/\"//g;s/\;//g"| awk '$22=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $24=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}' > /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.ref_lncRNA_exon.bed
#提取二代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】
cat /home/zhzhang/PG/2G_lnc/Macaca_mulatta.2G_lncRNA.gtf |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.2Gnovel_lncRNA_exon.bed
#合并全部lncRNA gene的外显子注释
cat "/home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.ref_lncRNA_exon.bed" "/home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.2Gnovel_lncRNA_exon.bed" > /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.all_lncRNA_exon.bed
#获取全部 lncRNA外显子与假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/lncRNA_class/）
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Macaca_mulatta/pgenes/Macaca_mulatta_hpg.bed -wo > /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.alllncRNAexon_intersect_pg.bed



#大鼠
#提取ref lncRNA基因外显子bed(GTF注释文件中，exon所属的基因biotype有的在22列，有的在24列.lncRNA gene类型包括："lncRNA" )【chr start end gene_id transcript_id strand】
#查看ref lncRNA类型：cat "/home/zhzhang/PG/HPG/GTF_PG/Rattus_norvegicus.mRatBN7.2.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}' | awk '$21=="gene_biotype" {print $22} $23=="gene_biotype" {print $24}'|sort|uniq
cat "/home/zhzhang/PG/HPG/GTF_PG/Rattus_norvegicus.mRatBN7.2.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}'|sed "s/\"//g;s/\;//g"| awk '$22=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $24=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}' > /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.ref_lncRNA_exon.bed
#提取二代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】
cat /home/zhzhang/PG/2G_lnc/Rattus_norvegicus.2G_lncRNA.gtf |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.2Gnovel_lncRNA_exon.bed
#合并全部lncRNA gene的外显子注释
cat "/home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.ref_lncRNA_exon.bed" "/home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.2Gnovel_lncRNA_exon.bed" > /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.all_lncRNA_exon.bed
#获取全部 lncRNA外显子与假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/lncRNA_class/）
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Rattus_norvegicus/pgenes/Rattus_norvegicus_hpg.bed -wo > /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.alllncRNAexon_intersect_pg.bed



#兔
#提取ref lncRNA基因外显子bed(GTF注释文件中，exon所属的基因biotype有的在22列，有的在24列.lncRNA gene类型包括："lncRNA" )【chr start end gene_id transcript_id strand】
#查看ref lncRNA类型：cat "/home/zhzhang/PG/HPG/GTF_PG/Oryctolagus_cuniculus.OryCun2.0.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}' | awk '$21=="gene_biotype" {print $22} $23=="gene_biotype" {print $24}'|sort|uniq
cat "/home/zhzhang/PG/HPG/GTF_PG/Oryctolagus_cuniculus.OryCun2.0.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}'|sed "s/\"//g;s/\;//g"| awk '$22=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $24=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}' > /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.ref_lncRNA_exon.bed
#提取二代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】
cat /home/zhzhang/PG/2G_lnc/Oryctolagus_cuniculus.2G_lncRNA.gtf |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.2Gnovel_lncRNA_exon.bed
#合并全部lncRNA gene的外显子注释
cat "/home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.ref_lncRNA_exon.bed" "/home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.2Gnovel_lncRNA_exon.bed" > /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.all_lncRNA_exon.bed
#获取全部 lncRNA外显子与假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/lncRNA_class/）
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Oryctolagus_cuniculus/pgenes/Oryctolagus_cuniculus_hpg.bed -wo > /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.alllncRNAexon_intersect_pg.bed



#负鼠
#提取ref lncRNA基因外显子bed(GTF注释文件中，exon所属的基因biotype有的在22列，有的在24列.lncRNA gene类型包括："lncRNA" )【chr start end gene_id transcript_id strand】
#查看ref lncRNA类型：cat "/home/zhzhang/PG/HPG/GTF_PG/Monodelphis_domestica.ASM229v1.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}' | awk '$21=="gene_biotype" {print $22} $23=="gene_biotype" {print $24}'|sort|uniq
cat "/home/zhzhang/PG/HPG/GTF_PG/Monodelphis_domestica.ASM229v1.108.chr.rmpg.gtf"|awk '$1!="MT" && $3=="exon" {print $0}'|sed "s/\"//g;s/\;//g"| awk '$22=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7} $24=="lncRNA" {print $1"\t"$4"\t"$5"\t"$10"\t"$14"\t"$7}' > /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.ref_lncRNA_exon.bed
#提取二代novel lncRNA基因外显子bed【chr start end gene_id transcript_id strand】
cat /home/zhzhang/PG/2G_lnc/Monodelphis_domestica.2G_lncRNA.gtf |awk '$3=="exon" {print $1"\t"$4"\t"$5"\t"$10"\t"$12"\t"$7}'|sed "s/\"//g;s/\;//g" > /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.2Gnovel_lncRNA_exon.bed
#合并全部lncRNA gene的外显子注释
cat "/home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.ref_lncRNA_exon.bed" "/home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.2Gnovel_lncRNA_exon.bed" > /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.all_lncRNA_exon.bed
#获取全部 lncRNA外显子与假基因重叠情况以及重叠长度（传至实验室服务器/home/zhzhang/PG/lncRNA_class/）
bedtools intersect -a /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.all_lncRNA_exon.bed -b /home/zhzhang/PG/PseudoPipe_wd/ppipe_output_blastpro/Monodelphis_domestica/pgenes/Monodelphis_domestica_hpg.bed -wo > /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.alllncRNAexon_intersect_pg.bed



#生成全部lncRNA基因ID列表
cat /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.all_lncRNA_exon.bed|awk '{print $4}'|sort|uniq > /home/zhzhang/PG/lncRNA_class/Gallus_gallus/Gallus_gallus.all_lncRNA.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA_exon.bed|awk '{print $4}'|sort|uniq > /home/zhzhang/PG/lncRNA_class/Mus_musculus/Mus_musculus.all_lncRNA.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.all_lncRNA_exon.bed|awk '{print $4}'|sort|uniq > /home/zhzhang/PG/lncRNA_class/Homo_sapiens/Homo_sapiens.all_lncRNA.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.all_lncRNA_exon.bed|awk '{print $4}'|sort|uniq > /home/zhzhang/PG/lncRNA_class/Danio_rerio/Danio_rerio.all_lncRNA.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.all_lncRNA_exon.bed|awk '{print $4}'|sort|uniq > /home/zhzhang/PG/lncRNA_class/Macaca_mulatta/Macaca_mulatta.all_lncRNA.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.all_lncRNA_exon.bed|awk '{print $4}'|sort|uniq > /home/zhzhang/PG/lncRNA_class/Rattus_norvegicus/Rattus_norvegicus.all_lncRNA.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.all_lncRNA_exon.bed|awk '{print $4}'|sort|uniq > /home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus/Oryctolagus_cuniculus.all_lncRNA.geneid.txt
cat /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.all_lncRNA_exon.bed|awk '{print $4}'|sort|uniq > /home/zhzhang/PG/lncRNA_class/Monodelphis_domestica/Monodelphis_domestica.all_lncRNA.geneid.txt




#lncRNA基因分类，并确定pgdlnc对应的假基因和反式靶标（母基因）（传至集群服务器/home/zhzhang/PG/lncRNA_class/物种）
#a输入假基因与lncRNA外显子交集文件，b输入lncRNA geneid列表，c输入假基因信息
#LncRNA_class输出lnc基因分类信息,LncRNA_pg输出pglnc对应的假基因,LncRNA_transtarget输出pgdlnc的反式靶标
lncrnaclass <- function(a,b,c,LncRNA_class,LncRNA_pg,LncRNA_transtarget){
  #导入lncRNA外显子与假基因交集文件
  alllncRNAexon_intersect_pg <- read.delim(a, header=FALSE)
  colnames(alllncRNAexon_intersect_pg)[c(4,5,13,14)] <- c("lncgid","lnctid","pgid","overlap_len")
  #导入lncRNA geneid列表
  all_lncRNA_geneid <- read.table(b, quote="\"", comment.char="")
  colnames(all_lncRNA_geneid) <- c("lncgid")
  #导入假基因及其母基因信息
  pg_parent <- read.delim(c)%>%
    select(16,17)%>%
    separate(parentgene_id,c("trans_acting_target"))
  colnames(pg_parent)[1] <- "pgid"
  #计算每个lncRNA基因每个转录本的全部外显子与假基因在相同链上的总重叠长度,并筛选出总重叠长度>=200的转录本（即与假基因显著交集的转录本）
  lncRNAPG <- filter(alllncRNAexon_intersect_pg,V6==V12)%>%
    group_by(lncgid,lnctid,pgid)%>%
    summarise("total_overlap_len"=sum(overlap_len))%>%
    filter(total_overlap_len>=200)
  #获得假基因来源的lncRNA基因的id
  pgdlncRNAid <- data.frame(lncRNAPG$lncgid)%>%
    distinct(lncRNAPG.lncgid)
  colnames(pgdlncRNAid) <- "lncgid"
  pgdlncRNAid <- mutate(pgdlncRNAid,type="Pseudogene-derived lncRNA")
  #获得与假基因有交集但不显著的lncRNA基因的id（Interference）【全部交集的lncid去除显著交集的lncid即为不显著的lncRNA基因的id（Interference）】
  InterferencelncRNAid <- distinct(alllncRNAexon_intersect_pg,lncgid)%>%
    dplyr::setdiff(select(pgdlncRNAid,lncgid))%>%
    mutate(type="Interference lncRNA")
  #lncRNA基因分类文件
  lncRNA_class <- left_join(all_lncRNA_geneid,pgdlncRNAid,by="lncgid")%>%
    left_join(InterferencelncRNAid,by="lncgid")%>%
    mutate(type="un")
  lncRNA_class$type[lncRNA_class$type.x=="Pseudogene-derived lncRNA"] <- "Pseudogene-derived lncRNA"
  lncRNA_class$type[lncRNA_class$type.y=="Interference lncRNA"] <- "Interference lncRNA"
  lncRNA_class$type[lncRNA_class$type=="un"] <- "Non-pseudogene-derived lncRNA"
  lncRNA_class <- select(lncRNA_class,lncgid,type)
  colnames(lncRNA_class)[1] <- "geneid"
  #确定假基因来源的lncRNA基因 对应的假基因
  #并根据假基因-母基因，确定假基因来源的lncRNA基因潜在的反式调控靶标基因
  pgdlncRNA_PG <- select(data.frame(lncRNAPG),-2) %>%
    arrange(lncgid,pgid,desc(total_overlap_len))%>%
    distinct(lncgid,pgid,.keep_all = T)%>%
    select(-3)%>%
    left_join(pg_parent,by="pgid")
  pgdlncRNA_Pse <- select(pgdlncRNA_PG,-3)%>%
    distinct(lncgid,pgid)
  pgdlncRNA_parent <- select(pgdlncRNA_PG,-2)%>%
    distinct(lncgid,trans_acting_target)
  #输出lncRNA基因分类文件 he 假基因来源的lncRNA基因潜在的反式调控靶标基因
  data.table::fwrite(lncRNA_class,file =LncRNA_class,sep = '\t',row.names = F,quote = F,col.names = T)
  data.table::fwrite(pgdlncRNA_Pse,file =LncRNA_pg,sep = '\t',row.names = F,quote = F,col.names = T)
  data.table::fwrite(pgdlncRNA_parent,file =LncRNA_transtarget,sep = '\t',row.names = F,quote = F,col.names = T)
}

#鸡
lncrnaclass("~/PG/lncRNA_class/Gallus_gallus.alllncRNAexon_intersect_pg.bed",
            "~/PG/lncRNA_class/Gallus_gallus.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Gallus_gallus_hpg.txt",
            "/home/zhzhang/PG/lncRNA_class/Gallus_gallus.lncRNA_class.txt",
            "/home/zhzhang/PG/lncRNA_class/Gallus_gallus.pgdlncRNA_pg.txt",
            "/home/zhzhang/PG/lncRNA_class/Gallus_gallus.pgdlncRNA_transtarget.txt")
#小鼠
lncrnaclass("~/PG/lncRNA_class/Mus_musculus.alllncRNAexon_intersect_pg.bed",
            "~/PG/lncRNA_class/Mus_musculus.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Mus_musculus_hpg.txt",
            "/home/zhzhang/PG/lncRNA_class/Mus_musculus.lncRNA_class.txt",
            "/home/zhzhang/PG/lncRNA_class/Mus_musculus.pgdlncRNA_pg.txt",
            "/home/zhzhang/PG/lncRNA_class/Mus_musculus.pgdlncRNA_transtarget.txt")
#人
lncrnaclass("~/PG/lncRNA_class/Homo_sapiens.alllncRNAexon_intersect_pg.bed",
            "~/PG/lncRNA_class/Homo_sapiens.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Homo_sapiens_hpg.txt",
            "/home/zhzhang/PG/lncRNA_class/Homo_sapiens.lncRNA_class.txt",
            "/home/zhzhang/PG/lncRNA_class/Homo_sapiens.pgdlncRNA_pg.txt",
            "/home/zhzhang/PG/lncRNA_class/Homo_sapiens.pgdlncRNA_transtarget.txt")
#斑马鱼
lncrnaclass("~/PG/lncRNA_class/Danio_rerio.alllncRNAexon_intersect_pg.bed",
            "~/PG/lncRNA_class/Danio_rerio.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Danio_rerio_hpg.txt",
            "/home/zhzhang/PG/lncRNA_class/Danio_rerio.lncRNA_class.txt",
            "/home/zhzhang/PG/lncRNA_class/Danio_rerio.pgdlncRNA_pg.txt",
            "/home/zhzhang/PG/lncRNA_class/Danio_rerio.pgdlncRNA_transtarget.txt")
#猕猴
lncrnaclass("~/PG/lncRNA_class/Macaca_mulatta.alllncRNAexon_intersect_pg.bed",
            "~/PG/lncRNA_class/Macaca_mulatta.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Macaca_mulatta_hpg.txt",
            "/home/zhzhang/PG/lncRNA_class/Macaca_mulatta.lncRNA_class.txt",
            "/home/zhzhang/PG/lncRNA_class/Macaca_mulatta.pgdlncRNA_pg.txt",
            "/home/zhzhang/PG/lncRNA_class/Macaca_mulatta.pgdlncRNA_transtarget.txt")
#大鼠
lncrnaclass("~/PG/lncRNA_class/Rattus_norvegicus.alllncRNAexon_intersect_pg.bed",
            "~/PG/lncRNA_class/Rattus_norvegicus.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Rattus_norvegicus_hpg.txt",
            "/home/zhzhang/PG/lncRNA_class/Rattus_norvegicus.lncRNA_class.txt",
            "/home/zhzhang/PG/lncRNA_class/Rattus_norvegicus.pgdlncRNA_pg.txt",
            "/home/zhzhang/PG/lncRNA_class/Rattus_norvegicus.pgdlncRNA_transtarget.txt")
#兔子
lncrnaclass("~/PG/lncRNA_class/Oryctolagus_cuniculus.alllncRNAexon_intersect_pg.bed",
            "~/PG/lncRNA_class/Oryctolagus_cuniculus.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Oryctolagus_cuniculus_hpg.txt",
            "/home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus.lncRNA_class.txt",
            "/home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus.pgdlncRNA_pg.txt",
            "/home/zhzhang/PG/lncRNA_class/Oryctolagus_cuniculus.pgdlncRNA_transtarget.txt")
#负鼠
lncrnaclass("~/PG/lncRNA_class/Monodelphis_domestica.alllncRNAexon_intersect_pg.bed",
            "~/PG/lncRNA_class/Monodelphis_domestica.all_lncRNA.geneid.txt",
            "~/PG/pg_message/Monodelphis_domestica_hpg.txt",
            "/home/zhzhang/PG/lncRNA_class/Monodelphis_domestica.lncRNA_class.txt",
            "/home/zhzhang/PG/lncRNA_class/Monodelphis_domestica.pgdlncRNA_pg.txt",
            "/home/zhzhang/PG/lncRNA_class/Monodelphis_domestica.pgdlncRNA_transtarget.txt")






```


