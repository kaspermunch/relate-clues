import os
from gwf import Workflow, AnonymousTarget
from subprocess import PIPE, Popen

gwf = Workflow(
    defaults={'account': 'simons'}
    )

###############################################################################
# Utility functions
###############################################################################

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

def relate_demography(chrom_start, chrom_end, chrom, pop, vcf_file_name, male_samples_file_name, demography_file_name):
    """
    Runs RELATE to infer demography
    """

    chrom = chrom.replace('chr', '')
    mutation_rate = 1.52e-08 # 5.25e-10 * 29

    inputs = [vcf_file_name]
    output_base_name = f'steps/relate_demog/{chrom}_{pop}/{chrom}_{pop}'
    outputs = {'coal': output_base_name + '.coal',
               'dist': output_base_name + '.dist',
               'mut': output_base_name + '.mut',
               'pairwise.bin': output_base_name + '.pairwise.bin',
               'pairwise.coal': output_base_name + '.pairwise.coal',
               'avg.rate': output_base_name + '_avg.rate'}

    cores = 10
    options = {'memory': f'{8*cores}g',
               'walltime': '6-00:00:00',
               'cores': cores,
              } 

    pwd = os.getcwd()
    def rel(a, b):
        "Computes path to a relative b"
        return os.path.relpath(a, os.path.dirname(b))

    spec = f"""
    # conda environment
    source ./scripts/conda_init.sh
    conda activate relate-clues

    # dir for outputll
    mkdir -p `dirname {output_base_name}`

    # extract vcf for samples in the region around snp
    bcftools view \
        -m2 -M2 -v snps \
        -O z \
        -r {chrom}:{chrom_start}-{chrom_end} \
        {vcf_file_name} > {output_base_name}.vcf.gz

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
        -N 30000 \
        -m {mutation_rate} \
        --seed 1 \
        -o relate

#        --coal {rel(demography_file_name, output_base_name)} 

    cat {rel(male_samples_file_name, output_base_name)} | perl -ne '/(\S+)[\n]+/ ; print $1 . " CEU WORLD\n"' > {os.path.basename(output_base_name)}.poplabels

    {rel(relate_path, output_base_name)}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
        -i relate \
        --num_iter 20 \
        -m {mutation_rate} \
        --years_per_gen 29 \
        --poplabels {os.path.basename(output_base_name)}.poplabels \
        --seed 1 \
        -o {os.path.basename(output_base_name)}
        --threads cores

    # back to orig dir
    cd '{pwd}'
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
    

def relate_samples(snp_pos, chrom, pop, vcf_file_name, male_samples_file_name, demography_file_name):
    """
    Runs RELATE sampling and CLUES on a SNP
    """

    chrom = chrom.replace('chr', '')
    mutation_rate = 1.52e-08 # 5.25e-10 * 29

    inputs = [vcf_file_name, male_samples_file_name, demography_file_name]
    output_base_name = f'steps/clues/{chrom}_{snp_pos}_{pop}/{chrom}_{snp_pos}_{pop}'
    outputs = dict((suffix, output_base_name + '.' + suffix) for suffix in ['anc', 'dist', 'mut', 'timeb'])

    options = {'memory': '1g',
               'walltime': '00:10:00',
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
        -m2 -M2 -v snps \
        -O z \
        -r {chrom}:{snp_pos-1000000}-{snp_pos+1000000} \
        {vcf_file_name} > {output_base_name}.vcf.gz

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
        --seed 1 \
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
    

def clues(snp_pos, chrom, pop, demography_file_name, timebins=[]):
    """
    Runs RELATE sampling and CLUES on a SNP
    """
    chrom = chrom.replace('chr', '')


    if timebins:
        input_base_name = f'steps/clues/{chrom}_{snp_pos}_{pop}/{chrom}_{snp_pos}_{pop}'
        output_base_name = f'steps/clues/{chrom}_{snp_pos}_{pop}/{chrom}_{snp_pos}_{pop}_{"_".join(map(str, timebins))}'
        timebins_option =  f'--timeBins {output_base_name}_timebins.txt'
    else:
        input_base_name = f'steps/clues/{chrom}_{snp_pos}_{pop}/{chrom}_{snp_pos}_{pop}'
        output_base_name = f'steps/clues/{chrom}_{snp_pos}_{pop}/{chrom}_{snp_pos}_{pop}'
        timebins_option = ''

    inputs = [input_base_name + suffix for suffix in ['.anc', '.dist', '.mut', '.timeb']] + [demography_file_name]
    outputs = dict((suffix, output_base_name + '.' + suffix + '.npy') for suffix in ['epochs', 'freqs', 'post'])
    outputs['out'] = output_base_name + '.out'

    options = {'memory': '8g',
               'walltime': '4:00:00',
               'cores': 1,
              } 

    spec = f"""
    # conda environment
    source ./scripts/conda_init.sh
    conda activate relate-clues

    rm -f {output_base_name}_timebins.txt
    for T in {' '.join(map(str, timebins))}; do echo $T >> {output_base_name}_timebins.txt ; done

    # run clues
    python {clues_path}/inference.py {timebins_option} \
        --times {input_base_name} \
        --tCutoff 2100 \
        --burnin 1000 --thin 100 \
        --out {output_base_name} \
        --output {output_base_name}.txt \
        --coal {demography_file_name} > {output_base_name}.out

    python {clues_path}/plot_traj.py {output_base_name} {output_base_name}        
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


###############################################################################
# Parameters
###############################################################################

relate_path = '../relate_v1.1.9_x86_64_static'
clues_path = '../clues'

#vcf_file_name = './data/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
vcf_file_name = 'steps/input_vcf/chrX_males_1000g.vcf.gz'
chrom = 'chrX'
populations = ['CEU']

time_bin_file_name = './data/timeBins.txt'
male_samples_file_names = dict([(pop, f'./steps/metainfo_{chrom}/{pop}_male.txt') for pop in populations])
genetic_map_file = './data/chrX_genetic_map.txt'
freq_data_file = f'steps/freq_data/derived_pop_freqs_{chrom}.h5'

tennesen_demography_file_name = 'data/tennessen_chrX.coal'

# candidate sweep centers:
window_centers = [19800000, 21200000]
#window_centers = [21200000]
flank = 250000 # relate flanks are 1000000 so it should be much smaller than that
# windows to sample snps from
clues_windows = [(int(pos - flank), int(pos + flank)) for pos in window_centers]
# min_freq = 0.5
min_freq = 0.1
nr_snps = 10



###############################################################################
# Run Relate and Clues
###############################################################################
for pop, male_samples_file_name in male_samples_file_names.items():

    target = gwf.target_from_template(f"relate_demog", 
                    relate_demography(3000000, 154000000, 
                        chrom, pop, vcf_file_name, male_samples_file_name, tennesen_demography_file_name))

    relate_demography_file_name = target.outputs['coal']

    clues_output_files = []
    for window_start, window_end in clues_windows:
        snp_info = get_snps(freq_data_file, chrom, pop, window_start, window_end, min_freq, nr_snps)
        for chrom_nr, snp_pos, ancestral_allele, derived_allele, derived_freq in snp_info:
            print(derived_freq)            
            gwf.target_from_template(f"relate_{chrom_nr}_{snp_pos}", 
                relate_samples(snp_pos, chrom, pop, vcf_file_name, male_samples_file_name, 
                # relate_demography_file_name
                tennesen_demography_file_name
                ))
            target = gwf.target_from_template(f"clues_{chrom_nr}_{snp_pos}", 
                clues(snp_pos, chrom, pop, 
                # relate_demography_file_name
                tennesen_demography_file_name
                ))
            target = gwf.target_from_template(f"clues_{chrom_nr}_{snp_pos}_epoque", 
                clues(snp_pos, chrom, pop, 
                # relate_demography_file_name
                tennesen_demography_file_name,
                timebins=[1379.0, 2068.0]
                ))                
            # clues_output_files.append(target.outputs['out'])

#     clues_output_summary_file = f'steps/clues/{chrom}_{pop}.h5'
#     gwf.target('gather_clues_output', 
#         inputs=clues_output_files, 
#         outputs=[hdf_file_name], 
#         walltime='10:00:00', memory='36g') << f"""


#     """
# re.search('logLR: (\S+).*selection\n(\S+)\s+(\S+)', s, re.MULTILINE | re.DOTALL).groups()