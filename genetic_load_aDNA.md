# Genetic load analysis

**NB!!!** Please check this https://github.com/NBISweden/GenErode before you decided to spend time on the following content.  

Genetic load calculation using various proxies for the real genetic fitness, e.g. mutation effect annnotation (snpEFF), loci(SNP level) conservation score (phyloP, phastCon), loci conservation score at AA level (???, should be better for cross species?), heterozygosity, ROH.  

This is my summary pipline for such analysis. Please let me know if anything is incorrect or should be added.

## General things (assumptions)  

- Ancestral state  
The whole analysis is based on the assumption of ancestral state (allele) known. Thus, the effect of a new mutation (derived allele) is compared to the ancestral state to infer the possible effect of this mutation (beneifical or deleterious).  
Normally, ancestral state is assumed/inferred in two ways, the simple (lazy :)) way and the correct way :). 
  - The simple way is to use one or more of the outgroups in your dataset and assume them to have the ancestral state OR use the major allele if you don't have outgroup genotype. 
  - The correct way is to infer the ancestral state based on some other assumptions, e.g. a phylogeny.  

Before doing any dirty work, please check annotations of your species' reference genome. In most cases, the ancestral state is avaiable in the annotation hub.

- Accurate genotype information  
<del> For ancient DNA folks, this is your nightmare. </del> It is important to have a reliable genotype called. It is better to have some good reliable SNPs than more SNPs with messy genotypes (low to medium).  
There are also ways to calculate genetic load based on genotype likelihood (not in this pipeline).

## Nerdy words
- Polarize  
flip the REF/ALT of your GenoType. E.g. GenoType=AA, REF=A, in vcf file your GT is 0/0 (0|0 if phased). When your ancestral state is known ANC=G, then GT should be flipped as 1/1.


## 1. snpEff
### 1.0 set env
``` bash
WDIR=/maps/projects/mjolnir1/people/gnr216/9-others/1-robin_gload
# your vcf file should be ready
```

### 1.1 use snpEff to annotate your vcf file

``` bash

```


### 1.2 filter 
I normally do the following filter
 - Remove scafffold taht is too short
 - Mask GT (mask to missing) with low depth


Remove scafffold taht is too short

``` bash
awk '$2>50000 {print $1}' ./emily_sharing/Copsychus.sechellarum.final.assembly.fasta.fai > list.chr50kb
```

Remove transversion (later)
```
```
Mask GT with low depth
N=202

```
bcftools view -m2 -M2 -v snps ./emily_sharing/snpeff.vcf | bcftools filter -S . -e '(GT="0/1" & AD[0:]<3) | (GT="0/1" & AD[:1]<3) | FMT/DP<5' | bcftools filter -e 'AC<2 | AC>(AN-2)' | bcftools view  -Oz -o snpeff.dp5.het3.mac2.vcf.gz &
bcftools index snpeff.dp5.het3.mac2.vcf.gz 
```
Maf50
```
bcftools filter -e 'AC/AN > 0.5' snpeff.dp5.het3.mac2.vcf.gz | bcftools query -f '%CHROM\t%POS\n' >  list.snp_f_maf50
```


```
# filter
# only keep biallelic SNPs, mask GT called with DP<5, mask het GT with min(REF,ALT)<3, then remove site with less than2 MAF COUNT or missing above 0.2
bcftools view -m2 -M2 -v snps raw.vcf.gz | bcftools filter -S . -e '(GT="0/1" & AD[0:]<3) | (GT="0/1" & AD[:1]<3) | FMT/DP<5' | bcftool
s filter -e ' AC<2 | AC>(AN-2) | AN < 137 ' | bcftools view -Oz -o dp5.het3.mac2.mis20.vcf.gz
bcftools index dp5.het3.mac2.mis20.vcf.gz
```


### 1.3 get ancestral state

Get ancestral state from outgroup
OUT_GROUPS : SAXI, SCRUB, SEYCH, TAEGUT

```
awk 'NR==FNR {a[$1]=1;next}{if($1 in a){print $1":"$2"-"$2}}' ../list.chr50kb ../snpeff.al2.tsv > list.snp.al2.chr50kb.samtools
samtools faidx ../emily_sharing/og_xin/outgroup_SCRUB.fa -r list.snp.al2.chr50kb.samtools | awk 'NR%2==0 {print $0}' > anc.SCRUB
wc -l * 
```

Only Kept loci where ANC == REF
```
paste list.snp.al2.chr50kb.samtools anc.SCRUB | awk 'NR==FNR {a[$1][$2]=$3;next}{split($1,b,":");split(b[2],c,"-");if($2==a[b[1]][c[1]]){print b[1],c[1]}}' ../snpeff.al2.tsv - > list.snp.anc_ref
```

OR ANC == Major Allele


### 1.4 get genetic load

Get syn, nonsyn, lof

``` bash
# 1. convert vcf to gt 0129
bcftools view -m2 -M2 -v snps snpeff.dp5.het3.mac2.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' | awk 'BEGIN{nbin=0;g["0/0"]=0;g["0/1"]=1;g["1/1"]=2;g["./."]=9;g["0|0"]=0;g["0|1"]=1;g["1|1"]=2;g[".|."]=9;a[0]="0";} {for(i=5;i<=NF;i++){$i=g[$i]}; print $0 }' OFS='\t' > gt012.snpeff.dp5.het3.mac2.al2 &

# 2. extract variant type
/usr/lib/jvm/java-19-openjdk-19.0.1.0.10-1.rolling.el8.x86_64/bin/java -jar /projects/mjolnir1/people/gnr216/a-software/snpEff/SnpSift.jar filter "(ANN[*].EFFECT = 'missense_variant') " snpeff.dp5.het3.mac2.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snpeff.dp5.het3
.mac2.nonsyn &
/usr/lib/jvm/java-19-openjdk-19.0.1.0.10-1.rolling.el8.x86_64/bin/java -jar /projects/mjolnir1/people/gnr216/a-software/snpEff/SnpSift.jar filter "(ANN[*].EFFECT = 'stop_gained') | (ANN[*].EFFECT = 'frameshift_variant') | (ANN[*].EFFECT = 'splice_acceptor_variant') | (ANN[*].EFFECT = 'splice_donor_variant')" snpeff.dp5.het3.mac2.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snpeff.dp5.het3.mac2.lof &
/usr/lib/jvm/java-19-openjdk-19.0.1.0.10-1.rolling.el8.x86_64/bin/java -jar /projects/mjolnir1/people/gnr216/a-software/snpEff/SnpSift.jar filter "(ANN[*].EFFECT = 'synonymous_variant') | (ANN[*].EFFECT = 'start_retained') | (ANN[*].EFFECT = 'stop_retained_variant') " snpeff.dp5.het3.mac2.vcf.gz | bcftools query -f '%CHROM\t%POS\n' > snpeff.dp5.het3.mac2.syn &

# 3. summary
## 3.1 major allele as anc, chr50kb
OUT="maj_50kb_all"; for type in lof syn nonsyn; do awk 'NR==FNR {a[$1]=1;next}{if($1 in a){print $0}}' list.chr50kb gt012.snpeff.dp5.het3.mac2.al2 | awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){}else{print $0}}' list.snp_f_maf50 - | awk 'NR==FNR {a[$1][$2]=1;next}{if($1 in a){if($2 in a[$1]){for(i=1;i<=202;i++){if($(i+4)==2){chom[i]++;};if($(i+4)==1){chet[i]++;};}}}}END{for(i in chom){print i,"hom",chom[i], "'${type}'"}; for(i in chet){print i,"het",chet[i],"'${type}'"};}' snpeff.dp5.het3.mac2.${type} - >> sum.$OUT ; done

## 3.2 major allele as anc, chr50kb + tv
OUT="maj_50kb_tv"; for type in lof syn nonsyn; do awk 'NR==FNR {a[$1]=1;next}{if($1 in a){print $0}}' list.chr50kb gt012.snpeff.dp5.het3.mac2.al2 | awk '{if(($3=="A" && $4=="G") || ($3=="T" && $4=="C") || ($3=="C" && $4=="T") || ($3=="G" && $4=="A")){}else{print $0}}' | awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){}else{print $0}}' list.snp_f_maf50 - | awk 'NR==FNR {a[$1][$2]=1;next}{if($1 in a){if($2 in a[$1]){for(i=1;i<=202;i++){if($(i+4)==2){chom[i]++;};if($(i+4)==1){chet[i]++;};}}}}END{for(i in chom){print i,"hom",chom[i], "'${type}'"}; for(i in chet){print i,"het",chet[i],"'${type}'"};}' snpeff.dp5.het3.mac2.${type} - >> sum.$OUT ; done


## 3.3 anc, chr50kb
# anc = REF
OUT="anc_SCRUB_50kb_all"; for type in lof syn nonsyn; do awk 'NR==FNR {a[$1]=1;next}{if($1 in a){print $0}}' list.chr50kb gt012.snpeff.dp5.het3.mac2.al2 | awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){}else{print $0}}' ./anc/list.snp.anc_ref - | awk 'NR==FNR {a[$1][$2]=1;next}{if($1 in a){if($2 in a[$1]){for(i=1;i<=202;i++){if($(i+4)==2){chom[i]++;};if($(i+4)==1){chet[i]++;};}}}}END{for(i in chom){print i,"hom",chom[i], "'${type}'"}; for(i in chet){print i,"het",chet[i],"'${type}'"};}' snpeff.dp5.het3.mac2.${type} - >> sum.$OUT ; done

## 3.4 anc, chr50kb + tv
OUT="anc_SCRUB_50kb_tv"; for type in lof syn nonsyn; do awk 'NR==FNR {a[$1]=1;next}{if($1 in a){print $0}}' list.chr50kb gt012.snpeff.dp5.het3.mac2.al2 | awk '{if(($3=="A" && $4=="G") || ($3=="T" && $4=="C") || ($3=="C" && $4=="T") || ($3=="G" && $4=="A")){}else{print $0}}' |  awk 'NR==FNR {a[$1][$2]=1;next}{if($2 in a[$1]){}else{print $0}}' ./anc/list.snp.anc_ref - | awk 'NR==FNR {a[$1][$2]=1;next}{if($1 in a){if($2 in a[$1]){for(i=1;i<=202;i++){if($(i+4)==2){chom[i]++;};if($(i+4)==1){chet[i]++;};}}}}END{for(i in chom){print i,"hom",chom[i], "'${type}'"}; for(i in chet){print i,"het",chet[i],"'${type}'"};}' snpeff.dp5.het3.mac2.${type} - >> sum.$OUT ; done

## 3.5 merge all sum

OUT="sum_merge";for f in anc_SCRUB_50kb_all anc_SCRUB_50kb_tv maj_50kb_all maj_50kb_tv; do awk '{print "'${f}'",$0}' sum.$f >> $OUT ;done 

## 3.6 anc, new, 363-bird
OUT="anc161_50kb_tv"; for type in lof syn nonsyn; do awk 'NR==FNR {a[$1]=1;next}{if($1 in a){print $0}}' list.chr50kb gt012.snpeff.dp5.het3.mac2.al2 | awk '{if(($3=="A" && $4=="G") || ($3=="T" && $4=="C") || ($3=="C" && $4=="T") || ($3=="G" && $4=="A")){}else{print $0}}' | awk 'NR==FNR {b[$4][$6]=$7;next}{if($1 in b){if($2 in b[$1]){$207=b[$1][$2];print $0}}}' ./alignment_hal/snpeff.al2.birdAnc161.c7_anc.bed - | awk '{if($207==$3){print $0}}' | awk 'NR==FNR {a[$1][$2]=1;next}{if($1 in a){if($2 in a[$1]){for(i=1;i<=202;i++){if($(i+4)==2){chom[i]++;};if($(i+4)==1){chet[i]++;};}}}}END{for(i in chom){print i,"hom",chom[i], "'${type}'"}; for(i in chet){print i,"het",chet[i],"'${type}'"};}' snpeff.dp5.het3.mac2.${type} - >> sum.$OUT ; done


```

About new ancestral state 
file is `snpeff.al2.birdAnc161.c7_anc.bed` in `emily_share`

```
$ head snpeff.al2.birdAnc161.c7_anc.bed
birdAnc161refChr340     271429  271430  scaffold_0      19      20      G
birdAnc161refChr340     271428  271429  scaffold_0      20      21      T
birdAnc161refChr340     271423  271424  scaffold_0      25      26      G
birdAnc161refChr340     271422  271423  scaffold_0      26      27      T
birdAnc161refChr340     271009  271010  scaffold_0      290     291     G
birdAnc161refChr340     270940  270941  scaffold_0      359     360     C
birdAnc161refChr340     270638  270639  scaffold_0      661     662     A
birdAnc161refChr340     270461  270462  scaffold_0      856     857     G
birdAnc161refChr340     270172  270173  scaffold_0      1141    1142    G
birdAnc161refChr340     269945  269946  scaffold_0      1369    1370    G
```

Columns are chromosme_in_ancestral_sequence Start_in_ANC END_in_ANC SMR_chromosome Start_in_SMR END_in_SMR Ancestral_Alelle  
One SNP per line  
Snpeff load were calculated in the previous section `#3.6`



### 1.5 Rscript for plot
new ancestral state

``` {r}
raw=read_table("./sum.anc161_50kb_tv",
               col_names = c("idx", "gt_type", "n", "snp_type"))

raw_meta = raw %>% 
    mutate(idx_vcf=idx +9) %>%
    mutate(id=meta_sample[as.character(idx_vcf),]$ind...22,
           dp=meta_sample[as.character(idx_vcf),]$COVERAGE,
           age=meta_sample[as.character(idx_vcf),]$age...20,
           pop=meta_sample[as.character(idx_vcf),]$island...21,
    ) %>% 
    
    pivot_wider(names_from = c(snp_type, gt_type),
                values_from = n
    ) %>% 
    # get gload
    mutate(
        rload_nonsyn= nonsyn_hom / (syn_hom *2 + syn_het),
        mload_nonsyn= nonsyn_het / (syn_hom *2 + syn_het),
        aload_nonsyn=(nonsyn_hom * 2 + nonsyn_het) / (syn_hom *2 + syn_het),
        rload_lof= lof_hom / (syn_hom *2 + syn_het),
        mload_lof= lof_het / (syn_hom *2 + syn_het),
        aload_lof=(lof_hom * 2 + lof_het) / (syn_hom *2 + syn_het),
    ) 

save(raw_meta, file='snpeff_newANC.rdata')

raw_meta %>%
    # view()
    
    mutate(
        nonsyn_hom=NULL, nonsyn_het=NULL,
        syn_hom=NULL, syn_het=NULL,
        lof_hom=NULL, lof_het=NULL,
        idx=NULL, idx_vcf=NULL,
    ) %>% 
    
    pivot_longer(
        cols = !c(id, age, dp, pop),
        names_to = "gload_type",
        values_to = "load"
    ) %>% 
    
    # !!!NB Change this filter when necessary
    filter(dp>4) %>%  
    
    # view()
    
    ggplot() +
    geom_beeswarm(aes(x=pop, color=age, y=load)) +
    
    facet_grid(gload_type~., scales = "free_y") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90)
    )


ggsave(filename = "gload_smr.anc161.png", width = 16, height = 9, dpi = 500)
```


# 2. phyloP

## 2.1 set env and prep files
list of files get from (363-way avian hub)[https://cglgenomics.ucsc.edu/data/cactus/]
NB!!! remove 363-avian-2020.hal ASAP. File too large
```
363-avian-2020.hal
363-avian-2020.hal.md5
ancRep_separate_models_rev.bw
ancRep_separate_models_rev.bw.conserved.bed
```

Things to consider:
 - same filter strategy
 - lift robin to chicken
 - remove transversion
 - remove loci with low DP
 - only include CDS?
 - extract anc from alignment?
 - 

``` bash
awk '{print $1, $2-1, $2}' OFS='\t' snpeff.al2.tsv > snpeff.al2.col3.bed


# extract phyloP
/maps/projects/mjolnir1/people/gnr216/a-software/ucsc/bigWigToBedGraph ancRep_separate_models_rev.bw ancRep_separate_models_rev.bedGraph &
```

## 2.2 liftover

``` bash
module load hal/2.2
# get stat from alignment file
halStats 363-avian-2020.hal > stat.363-avian-2020.hal

# liftover robin to chicken, around 2mio SNP left
halLiftover --bedType 3 363-avian-2020.hal Copsychus_sechellarum ../snpeff.al2.col6.bed Gallus_gallus snpeff.al2.gallus.bed

# liftover robin to anc_node, ~ 3.5 mio SNP left 
halLiftover --bedType 3 363-avian-2020.hal Copsychus_sechellarum ../snpeff.al2.col6.bed birdAnc161 snpeff.al2.birdAnc161.bed
halLiftover --bedType 3 363-avian-2020.hal Copsychus_sechellarum ../snpeff.al2.col6.bed birdAnc160 snpeff.al2.birdAnc160.bed
```


## 2.3 Phylop

only include significant conserved loci

NB!!! bedtools assumes files sorted

``` bash
# sort -k1,1 -k2,2n in.bed > in.sorted.bed.
```

``` bash
module load bedtools/2.30.0

# Get phyloP score, large memory used
sort -k1,1 -k2,2n snpeff.al2.gallus.bed > snpeff.al2.gallus.sorted.bed &
bedtools intersect -a snpeff.al2.gallus.sorted.bed -b ancRep_separate_models_rev.bedGraph -wb > snpeff.al2.gallus.sorted.phylop.bed &

# Get anc state for each SNP
bedtools getfasta -fi birdAnc161.fa  -bed snpeff.al2.birdAnc161.bed -bedOut | awk '{$7=toupper($7); print $0}' OFS='\t' > snpeff.al2.birdAnc161.c7_anc.bed
bedtools getfasta -fi birdAnc160.fa  -bed snpeff.al2.birdAnc160.bed -bedOut | awk '{$7=toupper($7); print $0}' OFS='\t' > snpeff.al2.birdAnc160.c7_anc.bed


# Intersect SMR_SNP_filtered with Conserved_gallus_phyloP_SNP
bedtools intersect -a snpeff.al2.gallus.sorted.phylop.bed -b ancRep_separate_models_rev.bw.conserved.bed -wb > snpeff.al2.gallus.sort
ed.phylop.conserved.bed &


# filter and get load
# filter order: only keep highly conserved phyloP in gallus; remove transversion; only keep SNP with ANC state; get load
# remove loci if no ANC state available

BATCH="anc161.tv.conserved.full"; awk 'NR==FNR {a[$4][$6]=$10;next}{if($1 in a){if($2 in a[$1]){$207=a[$1][$2];print $0}}}' snpeff.al2.gallus.sorted.phylop.conserved.bed ../gt012.snpeff.dp5.het3.mac2.al2 | awk '{if(($3=="A" && $4=="G") || ($3=="T" && $4=="C") || ($3=="C" && $4=="T") || ($3=="G" && $4=="A")){}else{print $0}}' | awk 'NR==FNR {b[$4][$6]=$7;next}{if($2 in b[$1]){$208=b[$1][$2];print $0}}' snpeff.al2.birdAnc161.c7_anc.bed - | awk '{for(i=1;i<=202;i++){if($(i+4)!=9){  if($4==$208){a[i]+=($207*(2-$(i+4)));n[i]++;m[i]+=($207*($(i+4)==1));r[i]+=$207*2*($(i+4)==0); }else{if($3==$208){a[i]+=($207*$(i+4));n[i]++; m[i]+=($207*($(i+4)==1));r[i]+=$207*2*($(i+4)==2);}}} }}END{for(i in n){print i,n[i],a[i],m[i],r[i],"'${BATCH}'"}}' > load.phylop.${BATCH}

BATCH="anc160.tv.conserved.full"; awk 'NR==FNR {a[$4][$6]=$10;next}{if($1 in a){if($2 in a[$1]){$207=a[$1][$2];print $0}}}' snpeff.al2.gallus.sorted.phylop.conserved.bed ../gt012.snpeff.dp5.het3.mac2.al2 | awk '{if(($3=="A" && $4=="G") || ($3=="T" && $4=="C") || ($3=="C" && $4=="T") || ($3=="G" && $4=="A")){}else{print $0}}' | awk 'NR==FNR {b[$4][$6]=$7;next}{if($2 in b[$1]){$208=b[$1][$2];print $0}}' snpeff.al2.birdAnc160.c7_anc.bed - | awk '{for(i=1;i<=202;i++){if($(i+4)!=9){  if($4==$208){a[i]+=($207*(2-$(i+4)));n[i]++;m[i]+=($207*($(i+4)==1));r[i]+=$207*2*($(i+4)==0); }else{if($3==$208){a[i]+=($207*$(i+4));n[i]++; m[i]+=($207*($(i+4)==1));r[i]+=$207*2*($(i+4)==2);}}} }}END{for(i in n){print i,n[i],a[i],m[i],r[i],"'${BATCH}'"}}' > load.phylop.${BATCH}

```



