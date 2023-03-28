# relate-clues
9375953
## Input files

### 1000g meta info

    export CHROM=chrX
    export RAW_VCF_FILE=/home/kmt/simons/faststorage/data/1000Genomes/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
    export METAINFO_FILE=/home/kmt/simons/faststorage/data/1000Genomes//metainfo/sample_info.csv

(`sample_info.csv` is converted to csv from the metainfo exel file that goes with the 1000 genomes data set)

Extract sample names from vcf and write files with sample names divided by population and sex (also writes pop_names.tsv). 

    mkdir -p steps/metainfo_${CHROM}
    zcat < $RAW_VCF_FILE | head -n 10000 | grep CHROM | perl -pe 's/\s+/\n/g' > steps/metainfo_${CHROM}/sample_names.txt
    python scripts/write_1000gen_meta_info.py steps/metainfo_${CHROM}/sample_names.txt $METAINFO_FILE steps/metainfo_${CHROM}

### VCF file with index

Create an indexed VCF file with only biallelic SNPs in males

    export POP=CEU
    export VCF_FILE=steps/input_vcf/chrX_males_1000g.vcf.gz
    export MALE_SAMPLES_FILE=steps/metainfo_${CHROM}/${POP}_male.txt


    # make dir
    mkdir -p `dirname $VCF_FILE`
    # extract male samples (retaining only sites variable in sample)
    vcftools --gzvcf $RAW_VCF_FILE --keep $MALE_SAMPLES_FILE --non-ref-ac-any 1 --recode --recode-INFO-all --stdout > tmp1.vcf
    # remove all duplicate sites https://www.internationalgenome.org/faq/why-are-there-duplicate-calls-in-the-phase-3-call-set
    grep "#" tmp1.vcf > tmp2.vcf
    grep -v "#" tmp1.vcf | sort -k2,2 -n -u >> tmp2.vcf
    # keep only biallelic snps
    bcftools view -m2 -M2 -v snps -O z tmp2.vcf > tmp3.vcf.gz
    # remove any snps that are no longer variable
    vcftools --gzvcf tmp3.vcf.gz --non-ref-ac-any 1 --recode --recode-INFO-all --stdout > tmp4.vcf
    # compress
    bcftools view -O z tmp4.vcf > $VCF_FILE
    # index
    bcftools index -t $VCF_FILE
    rm -f tmp1.vcf tmp.vcf.gz

#########


    # make dir
    mkdir -p `dirname $VCF_FILE`
    # extract male samples (retaining only sites variable in sample)
    vcftools --gzvcf $RAW_VCF_FILE --keep $MALE_SAMPLES_FILE --non-ref-ac-any 1 --recode --recode-INFO-all --stdout > $VCF_FILE
    # remove all duplicate sites https://www.internationalgenome.org/faq/why-are-there-duplicate-calls-in-the-phase-3-call-set
    grep "#" $VCF_FILE > tmp.vcf
    grep -v "#" $VCF_FILE | sort -k2,2 -n -u >> tmp.vcf
    bcftools view -m2 -M2 -v snps -O z tmp.vcf > $VCF_FILE
    bcftools index -t $VCF_FILE
    rm -f tmp.vcf


###############




    https://www.internationalgenome.org/faq/why-are-there-duplicate-calls-in-the-phase-3-call-set
    #zcat $RAW_VCF_FILE | grep -v '#' | cut -f 2 | sort -n | uniq -d | awk '{ print "'\\t"$0"\\t'"}' > duplicate_sites_re.txt
    zcat $RAW_VCF_FILE | grep -v '#' | cut -f 2 | sort -n | uniq -d > duplicate_sites.txt

    # Check that regexps in duplicate_sites_re.txt do not match anything else than the duplicates:
    wc -l duplicate_sites_re.txt
    zcat < data/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz | grep -f duplicate_sites.txt | wc -l




    mkdir -p `dirname $VCF_FILE`
    bcftools view -m2 -M2 -v snps -O z -S $MALE_SAMPLES_FILE $RAW_VCF_FILE > $VCF_FILE
    # remove all duplicate sites
    zcat < $VCF_FILE | grep "#"  > test.vcf
    zcat < $VCF_FILE | grep -v "#" | sort -k2,2 -n -u >> test.vcf
    bgzip test.vcf
    mv test.vcf.gz $VCF_FILE
    bcftools index -t $VCF_FILE




    | bcftools sort -O z - | bcftools norm -d none > tmp.vcf $VCF_FILE

    bcftools view -m2 -M2 -v snps -O v -S $MALE_SAMPLES_FILE $RAW_VCF_FILE | grep -v -f duplicate_sites.txt | bcftools view -O z - > $VCF_FILE





### HDF5 summaries of SNPs in VCF files

Make a file with the VCFs you have produced:

    ls steps/input_vcf/*.vcf.gz | perl -ne '/.*(chr\S).*/; print "$1 $_"' > steps/input_vcf/vcf_list.txt

This workflow reads `steps/input_vcf/vcf_list.txt` and generates hdf5 summary info on SNPs for each one. These are used in `workflow.py` to quickly select snps for analysis.

    gwf -f workflow_1000g_derived_info.py

### Genetic map

got `decode_hg38_sexavg_per_gen_lifted_tohg19_chrX.tsv` from simons. Produced a map for relate like this:

    cut -f 2,4,5 data/decode_hg38_sexavg_per_gen_lifted_tohg19_chrX.tsv | sed -e 's/\t/ /g' > data/chrX_genetic_map.txt

and added the propor header by hand:

    position COMBINED_rate.cM.Mb. Genetic_Map.cM.

### Demography and selection inference

Use tennesen as starting point for relates demography estimation:

    python scripts/tennesen2coal_chrX.py data/tennessen_popsize_fine.txt > data/tennessen_chrX.coal

But I estimate a demography on the samples themselves using relate:

    gwf run




