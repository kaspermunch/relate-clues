import os
from gwf import Workflow, AnonymousTarget
from subprocess import PIPE, Popen

gwf = Workflow(
    defaults={'account': 'simons'}
    )


def modpath(p, parent=None, base=None, suffix=None):
    """
    Modifies dir, basename, or suffix of path.
    """

    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base

    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path


###############################################################################
# Functions for retrieving snps
###############################################################################

def execute(cmd, stdin=None):
    process = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate(stdin)
    if not process.returncode == 0:
        print(cmd)
        print(stderr.decode())
    return stdout, stderr


def read_snp_info(snp_file):
    snp_list = list()
    with open('snps.txt', 'r') as snp_file:
        for line in snp_file:
            chrom, snp_pos, ancestral_allele, derived_allele, derived_freq = line.split()
            snp_pos = int(snp_pos)
            derived_freq = float(derived_freq)
            snp_list.append((chrom, snp_pos, ancestral_allele, derived_allele, derived_freq))
    return snp_list


def get_single_snp(freq_data_file, chrom, pop, snp_pos):
    snp_file_name = 'snps.txt'
    if os.path.exists(snp_file_name):
        os.remove(snp_file_name)  
    execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom.replace('chr', '')} {pop} {snp_file_name} --snppos {snp_pos}")
    snp_list = read_snp_info(snp_file_name)
    return snp_list


def get_snps(freq_data_file, chrom, pop, window_start, window_end, min_freq, nr_snps):
    snp_file_name = 'snps.txt'
    if os.path.exists(snp_file_name):
        os.remove(snp_file_name)  
    cmd = f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom.replace('chr', '')} {pop} {snp_file_name} --start {window_start} --end {window_end} --minfreq {min_freq} --maxsnps {nr_snps}"
    execute(cmd)
    snp_list = read_snp_info(snp_file_name)
    return snp_list


###############################################################################
# Templates
###############################################################################

def relate(snp_pos, chrom, pop, vcf_file_name, male_samples_file_name):
    """
    Runs RELATE sampling and CLUES on a SNP
    """

    chrom = chrom.replace('chr', '')
    mutation_rate = 1.52e-08 # 5.25e-10 * 29

    inputs = [vcf_file_name, male_samples_file_name]
    output_base_name = f'steps/clues/{chrom}_{snp_pos}_{pop}/{chrom}_{snp_pos}_{pop}'
    outputs = [output_base_name + suffix for suffix in ['.anc', '.dist', '.mut', '.timeb']]


    options = {'memory': '8g',
               'walltime': '4:00:00',
               'cores': 1,
              } 

    pwd = os.getcwd()
    def rel(a, b):
        "Computes path to a relative b"
        return os.path.relpath(a, os.path.dirname(b))

    spec = f"""
    # conda environment
    source ./scripts/conda_init.sh
    conda activate relate-clues

    # dir for output
    mkdir -p `dirname {output_base_name}`

    # extract vcf for samples in the region around snp
    bcftools view \
        -O z \
        -r {chrom}:{snp_pos-1000000}-{snp_pos+1000000} \
        -S {male_samples_file_name} {vcf_file_name} \
        > {output_base_name}.vcf.gz

    # turn vcf into clues input format
    {relate_path}/bin/RelateFileFormats \
        --mode ConvertFromVcf \
        --haps {output_base_name}.haps \
        --sample {output_base_name}.sample \
        -i {output_base_name}

    # move to output dir (relate only outputs to current dir...)
    cd {os.path.dirname(output_base_name)}
    rm -rf relate
    {rel(relate_path, output_base_name)}/bin/Relate \
        --mode All --haps {os.path.basename(output_base_name)}.haps \
        --sample {os.path.basename(output_base_name)}.sample \
        --map {rel(genetic_map_file, output_base_name)} \
        --coal {rel(demography_file_name, output_base_name)} \
        -m {mutation_rate} \
        -o relate

    # resample coalescence times at snp pos
    {rel(relate_path, output_base_name)}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
        -i relate -o {os.path.basename(output_base_name)} \
        -m {mutation_rate} \
        --coal {rel(demography_file_name, output_base_name)} \
        --format b \
        --num_samples 11000 \
        --first_bp {snp_pos} \
        --last_bp {snp_pos}

    # back to orig dir
    cd '{pwd}'
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
    

def clues(snp_pos, chrom, pop):
    """
    Runs RELATE sampling and CLUES on a SNP
    """
    chrom = chrom.replace('chr', '')
    output_base_name = f'steps/clues/{chrom}_{snp_pos}_{pop}/{chrom}_{snp_pos}_{pop}'
    inputs = [output_base_name + suffix for suffix in ['.anc', '.dist', '.mut', '.timeb']]
    outputs = [output_base_name + suffix for suffix in ['.epochs.npy', '.freqs.npy', '.post.npy']]

    options = {'memory': '8g',
               'walltime': '4:00:00',
               'cores': 1,
              } 

    spec = f"""
    # conda environment
    source ./scripts/conda_init.sh
    conda activate relate-clues

    # run clues
    python {clues_path}/inference.py \
        --times {output_base_name} \
        --tCutoff 2100 \
        --burnin 1000 --thin 100 \
        --out {output_base_name} \
        --coal {demography_file_name} 
        # \
        # --timeBins {time_bin_file_name}

    python {clues_path}/plot_traj.py {output_base_name} {output_base_name}        
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
        



###############################################################################
# Parameters
###############################################################################


# TODO: IT MAKES MOST SENSE TO PUT IN THE README THE COMMAND AND EXPLANATION TO PRODUCE THE INPUT FILES: REC RATE, DEMOG, *MALE* VCF ETC.
# AND SHOULD ALSO INCLUDE HOW TO BUILD THE FREQUENCY H5 FILE OF SNPS IN THE MALE VCF
# THE WORKFLOW SHOULD THEN START FROM THERE


# TODO: BCFTOOLS SHOULD NOT FILTER FOR MALES, THAT SHOULD BE DONE IN THE VCF FILE BEFORE STARTING THE WORKFLOW...

# location of software
relate_path = '../relate_v1.1.9_MacOSX_Intel'
clues_path = '../clues'

# assumes one VCF file containing one chromosome:
vcf_file_name = './data/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
chrom = 'chrX'
populations = ['CEU']

time_bin_file_name = './data/test_timeBins.txt'
male_samples_file_names = dict([(pop, f'./steps/metainfo_{chrom}/{pop}_male.txt') for pop in populations])
# TODO: make a proper recombination map from decode
genetic_map_file = './data/chrX_genetic_map.txt'

# got decode_hg38_sexavg_per_gen_lifted_tohg19_chrX.tsv from simons 
#cut -f 2,4,5 data/decode_hg38_sexavg_per_gen_lifted_tohg19_chrX.tsv | sed -e 's/\t/ /g' > data/chrX_genetic_map.txt
# and added header
# position COMBINED_rate.cM.Mb. Genetic_Map.cM.

# TODO:
# I WOULD ADD THE workflow_1000g.py and scripts from argweaver-clues to generate this file:
freq_data_file = 'derived_pop_freqs.h5'

# TODO: Iran this beforehand: SHOULD RUN AS TARGET
# python scripts/tennesen2coal.py data/tennessen_popsize_fine.txt > data/tennessen.coal
demography_file_name = 'data/tennessen.coal'

# candidate sweep centers:
# window_centers = [19800000, 21200000]
window_centers = [21200000]
flank = 250000 # relate flanks are 1000000 so it should be much smaller than that
# windows to sample snps from
clues_windows = [(int(pos - flank), int(pos + flank)) for pos in window_centers]
min_freq = 0.5
nr_snps = 10

###############################################################################
# Index VCF file
###############################################################################
# index vcf file and 
gwf.target(name=f'tabix_{chrom}',
     inputs=[vcf_file_name], 
     outputs=[vcf_file_name + '.tbi'], 
     walltime='03:00:00', 
     memory='8g') << f"""
bcftools index -t {vcf_file_name}
"""

###############################################################################
# Generate meta info for VCF files
###############################################################################
# extract sample names from vcf and write files with sample names devided 
# by population and sex (also writes pop_names.tsv)
# 
# sample_info.csv is part of the repository. It is converted to csv from 
# the metainfo exel file that goes with the 1000 genomes data set.
gwf.target(name=f'metainf_{chrom}',
     inputs=[vcf_file_name], 
     outputs=male_samples_file_names, 
     walltime='03:00:00', 
     memory='8g') << f"""
mkdir -p steps/metainfo_{chrom}
zcat < {vcf_file_name} | head -n 10000 | grep CHROM | perl -pe 's/\s+/\n/g' > steps/metainfo_{chrom}/sample_names.txt
python scripts/write_1000gen_meta_info.py steps/metainfo_{chrom}/sample_names.txt  data/metainfo/sample_info.csv  steps/metainfo_{chrom}
"""

###############################################################################
# Run Relate and Clues
###############################################################################
for pop, male_samples_file_name in male_samples_file_names.items():
    for window_start, window_end in clues_windows:
        snp_info = get_snps(freq_data_file, chrom, pop, window_start, window_end, min_freq, nr_snps)
        for chrom_nr, snp_pos, ancestral_allele, derived_allele, derived_freq in snp_info:
            gwf.target_from_template(f"relate_{chrom_nr}_{snp_pos}", 
                relate(snp_pos, chrom, pop, vcf_file_name, male_samples_file_name))
            gwf.target_from_template(f"clues_{chrom_nr}_{snp_pos}", 
                clues(snp_pos, chrom, pop))
