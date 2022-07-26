# relate-clues

## Input files

### VCF file with index

    bcftools index -t VCF_FILE

### HDF5 summary of SNPs in VCF file

Extract sample names from vcf and write files with sample names divided by population and sex (also writes pop_names.tsv). `sample_info.csv` is part of the repository. It is converted to csv from the metainfo exel file that goes with the 1000 genomes data set.

    zcat < VCF_FILE | head -n 10000 | grep CHROM | perl -pe 's/\s+/\n/g' > data/metainfo/sample_names.txt
    python scripts/write_1000gen_meta_info.py steps/metainfo/sample_names.txt  data/metainfo/sample_info.csv  steps/metainfo

    gwf -f workflow_1000g_derived_info.py



### Demography



### Genetic map

got `decode_hg38_sexavg_per_gen_lifted_tohg19_chrX.tsv` from simons 

    cut -f 2,4,5 data/decode_hg38_sexavg_per_gen_lifted_tohg19_chrX.tsv | sed -e 's/\t/ /g' > data/chrX_genetic_map.txt

and added header

    position COMBINED_rate.cM.Mb. Genetic_Map.cM.

