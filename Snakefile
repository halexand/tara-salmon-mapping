shell.executable("/bin/bash")
#shell.prefix("source $HOME/.bashrc; ")
import glob
import json
import os
import pandas as pd

IDS, = glob_wildcards("transcripts/{id}.fasta")
ERR_df = pd.read_csv('PRJEB6609_metaT_ERR.list')
EG_df = pd.read_csv('PRJEB4352_metaG_ERR.list')
metaT_dir = 'PRJEB6609'
metaG_dir = 'PRJEB4352'
ERR_list = list(ERR_df.run_accession)
GERR_list= list(EG_df.run_accession)
#ERR_list=ERR_list[300:]
trimmed_fasta_dir = '/vortexfs1/omics/alexander/data/TARA/PRJEB4352-snakmake-output/trimmed/PRJEB6609'
G_trimmed_fasta_dir = '/vortexfs1/omics/alexander/data/TARA/PRJEB4352-snakmake-output/trimmed/PRJEB4352'
rule all:
    input: 
        directory(expand('transcripts/{id}{ext}', id=IDS, ext=['.index'])),  
        expand(os.path.join('salmon', metaG_dir, '{id}', '{err}', 'quant.sf'), id=IDS, err=GERR_list), 
        expand(os.path.join('salmon', metaT_dir, '{id}', '{gerr}', 'quant.sf'), id=IDS, gerr=ERR_list), 
        expand(os.path.join('salmon', metaG_dir, '{id}'+ '.merged.tpm'), id = IDS), 
        expand(os.path.join('salmon', metaT_dir, '{id}' + '.merged.tpm'), id = IDS),
#        expand(os.path.join('salmon', metaT_dir, '{id}' + '.mapping.stats.csv'), id = IDS), 
#        expand(os.path.join('salmon', metaG_dir, '{id}' + '.mapping.stats.csv'), id = IDS)

localrules: summarize_mappingG, summarize_mapping
rule salmon_index:
    input: 'transcripts/{id}.fasta', 
    output: directory('transcripts/{id}.index')
    conda: 'envs/salmon.yaml'
    shell:
        '''
        salmon index -t {input} -i {output}
        '''
rule salmon_quant:
    input: index = directory('transcripts/{id}.index'), r1 = os.path.join(trimmed_fasta_dir, '{err}'+'_1.trimmed.fastq.gz'), r2 = os.path.join(trimmed_fasta_dir, '{err}'+'_2.trimmed.fastq.gz')
    output: os.path.join('salmon', metaT_dir, '{id}','{err}', 'quant.sf')
    params: outdir = os.path.join('salmon', metaT_dir , '{id}','{err}')
    conda: 'envs/salmon.yaml'
    shell:
        '''
        mkdir -p {params.outdir}
        salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -o {params.outdir} -p 16 || true
        touch {output}
        '''
rule salmon_quantG:
    input: index = directory('transcripts/{id}.index'), r1 = os.path.join(G_trimmed_fasta_dir, '{gerr}'+'_1.trimmed.fastq.gz'), r2 = os.path.join(G_trimmed_fasta_dir, '{gerr}'+'_2.trimmed.fastq.gz')
    output: os.path.join('salmon', metaG_dir, '{id}','{gerr}', 'quant.sf')
    params: outdir = os.path.join('salmon', metaG_dir, '{id}','{gerr}')
    conda: 'envs/salmon.yaml'
    shell:
        '''
        mkdir -p {params.outdir}
        salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -o {params.outdir} -p 16 --meta|| true
        touch {output}
        '''

rule salmon_merge:
    input: expand(os.path.join('salmon', metaT_dir, '{id}','{err}', 'quant.sf'), id = IDS, err=ERR_list)
    output: tpm = os.path.join('salmon', metaT_dir, '{id}'+'.merged.tpm'), numreads = os.path.join('salmon', metaT_dir, '{id}'+'.merged.numreads')
    params: indir = os.path.join('salmon', metaT_dir, '{id}', 'ERR*')
    conda: 'envs/salmon.yaml'
    shell:
        '''
        salmon quantmerge --column NumReads --quants {params.indir} -o {output.numreads}
        salmon quantmerge --column TPM --quants {params.indir} -o {output.tpm}

        '''

rule salmon_mergeG:
    input: expand(os.path.join('salmon', metaG_dir, '{id}','{gerr}', 'quant.sf'), id = IDS, gerr=GERR_list)
    output: tpm = os.path.join('salmon', metaG_dir, '{id}'+'.merged.tpm'), numreads = os.path.join('salmon', metaG_dir, '{id}'+'.merged.numreads')
    params: indir = os.path.join('salmon', metaG_dir, '{id}', 'ERR*')
    conda: 'envs/salmon.yaml'
    shell:
        '''
        salmon quantmerge --column NumReads --quants {params.indir} -o {output.numreads}
        salmon quantmerge --column TPM --quants {params.indir} -o {output.tpm}
        '''

#rule summarize_mapping:
#    input: os.path.join('salmon', metaT_dir, '{id}'+'.merged.tpm')
#    output: os.path.join('salmon', metaT_dir, '{id}' + '.mapping.stats.csv')
#    params: indir = os.path.join('salmon', metaT_dir, '{id}')
#    run: 
#        columns = ['num_mapped', 'num_processed', 'percent_mapped']
#        df = pd.DataFrame(columns = columns) 
#        for err in glob.glob(params.indir):
#            name = os.path.basename(err)
#            with open(os.path.join(err,'aux_info','meta_info.json')) as f:
#                data = json.load(f)
#                for col in columns: 
#                    df.loc[name, col] = data[col]
#        df.to_csv(output[0])
#
#rule summarize_mappingG:
#    input: os.path.join('salmon', metaG_dir, '{id}'+'.merged.tpm')
#    output: os.path.join('salmon', metaG_dir, '{id}' + '.mapping.stats.csv')
#    params: indir = os.path.join('salmon', metaG_dir, '{id}')
#    run:
#        columns = ['num_mapped', 'num_processed', 'percent_mapped']
#        df = pd.DataFrame(columns = columns)
#        for err in glob.glob(params.indir):
#            name = os.path.basename(err)
#            with open(os.path.join(err,'aux_info','meta_info.json')) as f:
#                data = json.load(f)
#                for col in columns:
#                    df.loc[name, col] = data[col]
#        df.to_csv(output[0])
#
