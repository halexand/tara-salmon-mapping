shell.executable("/bin/bash")
#shell.prefix("source $HOME/.bashrc; ")
import glob
import json
import os
import pandas as pd

configfile: "config.yaml"

# Load general variables and tables (ERR sample lists etc.) 
T_ERR_df = pd.read_csv(config['metaT_ERR'])
G_ERR_df = pd.read_csv(config['metaG_ERR'])
metaT_dir = 'PRJEB6609'
metaG_dir = 'PRJEB4352'
T_ERR_list = list(T_ERR_df.run_accession)
G_ERR_list= list(G_ERR_df.run_accession)

T_trimmed_fasta_dir = config['metaT_dir']
G_trimmed_fasta_dir = config['metaG_dir']


IDS, = glob_wildcards("transcripts/{id}.fasta")

rule all:
    input: 
        directory(expand('transcripts/{id}{ext}', id=IDS, ext=['.index'])),  
        expand(os.path.join('salmon', metaG_dir, '{id}', '{gerr}', 'quant.sf'), id=IDS, gerr=G_ERR_list), 
        expand(os.path.join('salmon', metaT_dir, '{id}', '{terr}', 'quant.sf'), id=IDS, terr=T_ERR_list), 
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
    input: index = directory('transcripts/{id}.index'), r1 = os.path.join(T_trimmed_fasta_dir, '{terr}'+'_1.trimmed.fastq.gz'), r2 = os.path.join(T_trimmed_fasta_dir, '{terr}'+'_2.trimmed.fastq.gz')
    output: os.path.join('salmon', metaT_dir, '{id}','{terr}', 'quant.sf')
    params: outdir = os.path.join('salmon', metaT_dir , '{id}','{terr}')
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
    input: expand(os.path.join('salmon', metaT_dir, '{id}','{err}', 'quant.sf'), id = IDS, err=T_ERR_list)
    output: tpm = os.path.join('salmon', metaT_dir, '{id}'+'.merged.tpm'), numreads = os.path.join('salmon', metaT_dir, '{id}'+'.merged.numreads')
    params: indir = os.path.join('salmon', metaT_dir, '{id}', 'ERR*')
    conda: 'envs/salmon.yaml'
    shell:
        '''
        salmon quantmerge --column NumReads --quants {params.indir} -o {output.numreads}
        salmon quantmerge --column TPM --quants {params.indir} -o {output.tpm}

        '''

rule salmon_mergeG:
    input: expand(os.path.join('salmon', metaG_dir, '{id}','{gerr}', 'quant.sf'), id = IDS, gerr=G_ERR_list)
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
