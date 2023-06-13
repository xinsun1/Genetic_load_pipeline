# genetic load analysis for black rhino
### with snpEff

## 1.  build snpEff database  
***http://pcingola.github.io/SnpEff/se_build_db/***

``` bash
# prepare config 
cp /projects/mjolnir1/people/gnr216/a-software/snpEff/snpEff.config ./
echo "# BlackRhino, version 19Dec2016_S9zPn
black_rhino.genome : BlackRhino" > black_rhino.ref.config

cat black_rhino.ref.config >> snpEff.config

# prepare fasta
mkdir data && cd data
mkdir genomes

cp /projects/mjolnir1/data/user_owned_folders/Mick_sharing/Black_rhinos/ref_genome/black_rhino_19Dec2016_S9zPn.fasta black_rhino.fa

# prepare gff
# change name
# folder ./data/black_rhino/sequences.fa
# ./data/black_rhino/genes.fa
# change name
zcat /projects/mjolnir1/data/user_owned_folders/Mick_sharing/Black_rhinos/ref_genome/GCA_013634535.1_ASM1363453v1_genomic.fna.gz | grep "^>" | awk -F "," '{print $1}' | awk '{print substr($1,2),$NF}' OFS='\t' > list.name.JAC_S9z 
zcat /projects/mjolnir1/data/user_owned_folders/Mick_sharing/Black_rhinos/ref_genome/GCA_013634535.1_ASM1363453v1_genomic.gff.gz | awk 'NR==FNR {a[$1]=$2;next}{if($1 in a){$1=a[$1]};print $0}' OFS='\t' list.name.JAC_S9z -  > genes.gff

# build
java -jar /projects/mjolnir1/people/gnr216/a-software/snpEff/snpEff.jar build -gff3 -v black_rhino -noCheckCds -noCheckProtein -quiet

```

## 2. genotype calling
***samtools***

```bash
# get bam dp


# get scaffold used 
awk '$4>0 {print $1,$3}' OFS='\t' ./2-qc/by_id/cov.mbr_TG0271 | grep -v "#" > size.chr14mb_micky

# generate 10Mb window
bedtools makewindows -g size.chr14mb_micky -w 10000000 | awk '{$2+=1; print $1":"$2"-"$3}' OFS='\t' > size.chr14mb_micky.10M.sam

# run samtools by window

# merge results
bcftools concat -f list.vcf_by_region -O z -o raw.vcf.gz
bcftools index raw.vcf.gz

# filter 
# only keep biallelic SNPs, mask GT called with DP<5, mask het GT with min(REF,ALT)<3, then remove site with less than2 MAF COUNT or missing above 0.2
bcftools view -m2 -M2 -v snps raw.vcf.gz | bcftools filter -S . -e '(GT="0/1" & AD[0:]<3) | (GT="0/1" & AD[:1]<3) | FMT/DP<5' | bcftools filter -e ' AC<2 | AC>(AN-2) | AN < 137 ' | bcftools view -Oz -o dp5.het3.mac2.mis20.vcf.gz
bcftools index dp5.het3.mac2.mis20.vcf.gz 

```

## 3. snpEff annotation

### 3.1 run snpEff


```bash
java -Xmx64g -jar /projects/mjolnir1/people/gnr216/a-software/snpEff/snpEff.jar -c $WDIR/1-snpEff_db/snpEff.config -v black_rhino $WDIR/3-gt_sam/dp5.het3.mac2.mis20.vcf.gz | bgzip > ann.dp5.het3.mac2.mis20.vcf.gz

```

### 3.2 result summary
Some description about the extremely hard to read awk commands.
 - 2. extract variant type  
 extract snp information by variant annotation groups, including synonymous, nonsynonymous and loss of function.  
 - 3. count by sample
count number of derived alleles per individual per variant group
 - 4. only transversion
 repeat 3. except using only transverion
 - 5. exclude maf> 0.5
 repeat 3. except remove SNPs with MAF > 0.5. This is a stupid decision. Should include these SNPs by polarize the ancestral state as the the major allele.


``` bash
# 1. convert vcf to gt 0129
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' ann.dp5.het3.mac2.mis20.vcf.gz | awk 'BEGIN{nbin=0;g["0/0"]=0;g["0/1"]=1;g["1/1"]=2;g["./."]=9;g["0|0"]=0;g["0|1"]=1;g["1|1"]=2;g[".|."]=9;a[0]="0";} {for(i=5;i<=NF;i++){$i=g[$i]}; print $0 }' OFS='\t' > gt.012.dp5.het3.mac2.mis20 &

# 2. extract variant type
/usr/lib/jvm/jre-19-openjdk-19.0.0.0.36-2.rolling.el8.x86_64/bin/java -jar /projects/mjolnir1/people/gnr216/a-software/snpEff/SnpSift.jar filter "(ANN[*].EFFECT = 'missense_variant') "  ann.dp5.het3.mac2.mis20.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snp.ann.dp5.het3.mac2.mis20.nonsyn  &
/usr/lib/jvm/jre-19-openjdk-19.0.0.0.36-2.rolling.el8.x86_64/bin/java -jar /projects/mjolnir1/people/gnr216/a-software/snpEff/SnpSift.jar filter "(ANN[*].EFFECT = 'stop_gained') | (ANN[*].EFFECT = 'frameshift_variant') | (ANN[*].EFFECT = 'splice_acceptor_variant') | (ANN[*].EFFECT = 'splice_donor_variant') " ann.dp5.het3.mac2.mis20.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snp.ann.dp5.het3.mac2.mis20.lof &
/usr/lib/jvm/jre-19-openjdk-19.0.0.0.36-2.rolling.el8.x86_64/bin/java -jar /projects/mjolnir1/people/gnr216/a-software/snpEff/SnpSift.jar filter "(ANN[*].EFFECT = 'synonymous_variant') | (ANN[*].EFFECT = 'start_retained') | (ANN[*].EFFECT = 'stop_retained_variant') " ann.dp5.het3.mac2.mis20.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snp.ann.dp5.het3.mac2.mis20.syn &

# 3. count by sample
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==2){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.syn gt.012.dp5.het3.mac2.mis20 | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"syn"}' ../list.id.dp5 - > ann.sum
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==2){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.lof gt.012.dp5.het3.mac2.mis20 | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"lof"}' ../list.id.dp5 - >> ann.sum
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==2){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.nonsyn gt.012.dp5.het3.mac2.mis20 | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"nonsyn"}' ../list.id.dp5 - >> ann.sum
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==1){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.nonsyn gt.012.dp5.het3.mac2.mis20 | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"nonsyn_het"}' ../list.id.dp5 - >> ann.sum
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==1){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.lof gt.012.dp5.het3.mac2.mis20 | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"lof_het"}' ../list.id.dp5 - >> ann.sum

awk '{a[$1][$3]=$2}END{for(i in a){print i,a[i]["lof"]/a[i]["syn"], a[i]["nonsyn"]/a[i]["syn"], a[i]["lof_het"]/a[i]["syn"], a[i]["nonsyn_het"]/a[i]["syn"]}}' ann.sum  | awk 'NR==FNR {a[$1]=$3;next}{$6=a[$1];print $0}' ../2-qc/dp_sum - > ann.sum.load

# 4. only transversion

awk '{if(($3=="A" && $4=="G") || ($3=="T" && $4=="C") || ($3=="C" && $4=="T") || ($3=="G" && $4=="A")){}else{print $0}}' F| OFS='\t'  gt.012.dp5.het3.mac2.mis20  > gt.012.dp5.het3.mac2.mis20.tv &

awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==2){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.syn gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"syn"}' ../list.id.dp5 - > ann.sum.tv
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==2){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.lof gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"lof"}' ../list.id.dp5 - >> ann.sum.tv
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==2){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.nonsyn gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"nonsyn"}' ../list.id.dp5 - >> ann.sum.tv
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==1){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.nonsyn gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"nonsyn_het"}' ../list.id.dp5 - >> ann.sum.tv
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==1){c[i]++;}}}}END{for(i in c){print i,c[i]}}' snp.ann.dp5.het3.mac2.mis20.lof gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"lof_het"}' ../list.id.dp5 - >> ann.sum.tv

awk '{a[$1][$3]=$2}END{for(i in a){print i,a[i]["lof"]/a[i]["syn"], a[i]["nonsyn"]/a[i]["syn"], a[i]["lof_het"]/a[i]["syn"], a[i]["nonsyn_het"]/a[i]["syn"]}}' ann.sum.tv  | awk 'NR==FNR {a[$1]=$3;next}{$6=a[$1];print $0}' ../2-qc/dp_sum - > ann.sum.tv.load

# 5. exclude maf> 0.5
bcftools filter -e 'AC/AN > 0.5' ann.dp5.het3.mac2.mis20.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snp.ann.dp5.het3.mac2.mis20.maf50

OUT="ann.sum.tv.maf50"

awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){print $0}}' snp.ann.dp5.het3.mac2.mis20.maf50 snp.ann.dp5.het3.mac2.mis20.syn | awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==2){c[i]++;}}}}END{for(i in c){print i,c[i]}}' - gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"syn"}' ../list.id.dp5 - > $OUT
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){print $0}}' snp.ann.dp5.het3.mac2.mis20.maf50 snp.ann.dp5.het3.mac2.mis20.lof | awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==2){c[i]++;}}}}END{for(i in c){print i,c[i]}}' - gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"lof"}' ../list.id.dp5 - >> $OUT
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){print $0}}' snp.ann.dp5.het3.mac2.mis20.maf50 snp.ann.dp5.het3.mac2.mis20.nonsyn | awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==2){c[i]++;}}}}END{for(i in c){print i,c[i]}}' - gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"nonsyn"}' ../list.id.dp5 - >> $OUT

awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){print $0}}' snp.ann.dp5.het3.mac2.mis20.maf50 snp.ann.dp5.het3.mac2.mis20.lof | awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==1){c[i]++;}}}}END{for(i in c){print i,c[i]}}' - gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"lof_het"}' ../list.id.dp5 - >> $OUT
awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){print $0}}' snp.ann.dp5.het3.mac2.mis20.maf50 snp.ann.dp5.het3.mac2.mis20.nonsyn | awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){for(i=1;i<=86;i++){if($(i+4)==1){c[i]++;}}}}END{for(i in c){print i,c[i]}}' - gt.012.dp5.het3.mac2.mis20.tv | awk 'NR==FNR {a[NR]=$1;next}{print a[$1],$2,"nonsyn_het"}' ../list.id.dp5 - >> $OUT


awk '{a[$1][$3]=$2}END{for(i in a){print i,a[i]["lof"]/a[i]["syn"], a[i]["nonsyn"]/a[i]["syn"], a[i]["lof_het"]/a[i]["syn"], a[i]["nonsyn_het"]/a[i]["syn"]}}' ann.sum.tv.maf50  | awk 'NR==FNR {a[$1]=$3;next}{$6=a[$1];print $0}' ../2-qc/dp_sum - > ann.sum.tv.maf50.load

```


## 4. Rscript for plot
Just an example for the final plot.


``` {r}
#### result summary for blackrhino genetic load ####

#### 0. set env ----
setwd('~/Documents/Projects/Rhino_black/1.gload/')

library(tidyverse)


#### 1. load result ----
library(googlesheets4)

gs4_deauth()

raw = read_sheet("https://docs.google.com/spreadsheets/d/1igQyKkNVh5SeKrCRToC24qiCriG3fBLQyBOUok65S-A/edit#gid=0",
                         sheet="MAF50")
rownames(raw) = raw$Sample_ID


#### 2. plot result ----
# 2.0 plot sequencing depth
raw %>%
    ggplot() + 
    geom_beeswarm(aes(x=Age,
                      y=Sequencing_depth, color=Age), cex = 2.5, method = "center") 

raw %>%
    # filter(Age %in% c("Modern")) %>%
    filter(Age %in% c("Historical")) %>%
    
    ggplot() + 
    geom_point(aes(x=Sequencing_depth,
                   y=`Nonsynonymous_Homozygous/Synonymous_Homozygous`))

## 2.1 1) All Historical vs Modern
library(ggbeeswarm)

raw %>%
    # exclude sample without information
    filter(!is.na(Age)) %>%
    filter(Sequencing_depth > 8) %>%
    # filter(Country %in% c("Kenya")) %>%
    
    # merge homo and het 
    mutate(all_lof = 2 *`Lof_Homozygous/Synonymous_Homozygous` +
               `Lof_Heterozygous/Synonymous_Homozygous`,
           all_nonsyn = 2 * `Nonsynonymous_Homozygous/Synonymous_Homozygous` +
               `Nonsynonymous_Heterozygous/Synonymous_Homozygous`) %>%
    
    pivot_longer(!c("Sample_ID", "Country",	"Age", "Given population", "Sequencing_depth"),
                 names_to = "category",
                 values_to = "gload") %>%
    
    ggplot() + 
    geom_beeswarm(aes(x=Age,
                      y=gload, color=Age), cex = 2.5, method = "center") +
    facet_grid(category ~., scales = "free") + 
    
    theme(strip.text.y = element_text(size=5))

## 2.2 2) Historical vs Modern based on given populations
## 2.3 3) Values between populations of only Historical.

```



