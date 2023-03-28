
from gwf import Workflow, AnonymousTarget
import os, re, sys
import numpy as np

gwf = Workflow()

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

def extract_pop_frequencies(vcf_file_name, freq_file_name):

    inputs = [vcf_file_name]
    outputs = [freq_file_name + '.frq']
    options = {
        'cores': 1,
        'memory': '16g'
    }
    spec = '''
    mkdir -p steps/freq_data

    vcftools --gzvcf {} --freq --out {}
    '''.format(vcf_file_name, freq_file_name)

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

vcf_input_list_file = 'steps/input_vcf/vcf_list.txt'

included_populations = ['CEU', 'CHB', 'YRI']

with open(vcf_input_list_file) as f:
    for line in f:
        chrom, vcf_file_name = line.split()

        # read population identifiers
        pop_info_file_name = f'steps/metainfo_{chrom}/pop_names.tsv'
        populations = []
        with open(pop_info_file_name) as f:
            for line in f:
                pop, _ = line.split(maxsplit=1)
                if pop not in included_populations:
                    continue
                populations.append(pop)

        freq_task_list = list()
        for pop in populations:
            freq_file_name = f'steps/freq_data/{pop}_{chrom}'
            target = gwf.target_from_template(f"freqs_{pop}_{chrom}", 
                extract_pop_frequencies(vcf_file_name, freq_file_name))

            freq_task_list.append(target)


        freq_files = [output for task in freq_task_list for output in task.outputs]

        # Build data set:
        hdf_file_name = f'steps/freq_data/derived_pop_freqs_{chrom}.h5'

        gwf.target('compile_freq_data', 
            inputs=freq_files, 
            outputs=[hdf_file_name], 
            walltime='10:00:00', memory='36g') << f"""

        python scripts/build_derived_freq_data.py {hdf_file_name} {chrom} {','.join(included_populations)}
        """
