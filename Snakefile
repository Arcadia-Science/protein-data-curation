import pandas as pd
        
samples_df = pd.read_table('inputs/samples.tsv').set_index(["sample","genus_species"],drop=False)
SAMPLES = list(samples_df['sample'])
GENUSSP = list(samples_df['genus_species'])          
          
rule all:
    input:
        expand("outputs/annotated/{genussp}_combined.tsv",genussp=GENUSSP)

rule download_proteome:
    """
    This rule downloads each transcriptome/proteome file specified in samples_df
    """
    output: "outputs/proteomes/{sample}.gz"
    conda:"envs/curl.yml"
    params:
        # get the download link from the "link" column in the samples data frame
        download_link = lambda wildcards: samples_df.loc[(wildcards.sample.rstrip(".gz"), slice(None)), "link"].values[0]
    shell: """
        curl -L {params.download_link} -o {output}
        """

rule proteome_combine:
    """
    This rule combines input files for species that have multiple transcriptomes/proteomes and moves all files to a new folder with one file per species.
    If there's only one file in the input table for the species, it copies it to the combinedproteomes folder.
    """
    output: "outputs/combinedproteomes/{genussp}.gz"
    input:
        samples = lambda wildcards: [f"outputs/proteomes/{sample}.gz" for sample in samples_df.xs(key=wildcards.genussp, level='genus_species').index]
    run:
        import shutil

        # Collect all sample files for the current genussp
        sample_files = input.samples

        # If there's more than one sample, combine them
        if len(sample_files) > 1:
            with open(output[0], 'wb') as wfd:
                for f in sample_files:
                    with open(f, 'rb') as fd:
                        shutil.copyfileobj(fd, wfd)
        # If there's only one sample, copy it to the combinedproteomes folder
        else:
            shutil.copy(sample_files[0], output[0])


rule transdecoder_predict:
    """
    This rule runs transdecoder, which will predict transcripts in an assembled transcriptome, on files that have 'yes' in the transdecoder column of the input sample sheet. If the sample is a proteome
    and does not need transdecoder, it will copy the proteome file into the 'readytopreprocess' folder
    """
    input:
        "outputs/combinedproteomes/{genussp}.gz"
    output:
        "outputs/readytopreprocess/{genussp}.faa"
    conda:"envs/transdecoder.yml"
    params: 
        db=lambda wildcards: samples_df.xs(wildcards.genussp, level="genus_species")["transdecoder"].iloc[0]
    shell: """
    if [ "{params.db}" == "yes" ]; then
        if ! file {input} | grep -q "gzip compressed"; then
            gzip -c {input} > "{input}.tmp"
            mv "{input}.tmp" {input}
        fi
        TransDecoder.LongOrfs -t {input} --output_dir outputs/transdecoder/
        TransDecoder.Predict -t {wildcards.genussp} --output_dir outputs/transdecoder/
        rm {wildcards.genussp}
        cp outputs/transdecoder/{wildcards.genussp}.transdecoder.pep outputs/readytopreprocess/{wildcards.genussp}.faa
    elif [ "{params.db}" == "no" ]; then
        gunzip -c outputs/combinedproteomes/{wildcards.genussp}.gz > outputs/combinedproteomes/{wildcards.genussp} 2>/dev/null || cat outputs/combinedproteomes/{wildcards.genussp}.gz > outputs/combinedproteomes/{wildcards.genussp}
        cp outputs/combinedproteomes/{wildcards.genussp} outputs/readytopreprocess/{wildcards.genussp}.faa
    fi
    """

rule preprocess_one:
    """
    This rule prepares data to be compatible with the NovelTree workflow. This rule will remove stop codons at the end of sequences and replace any others with ambiguous AAs, and will
    remove any proteins less than 25 aa.
    """
    input:
        "outputs/readytopreprocess/{genussp}.faa"
    output:
        tmp1=temp("outputs/tmp/{genussp}temp1.pep"),
        tmp2=temp("outputs/tmp/{genussp}temp2.pep"),
        tmp3=temp("outputs/tmp/{genussp}temp3.pep")
    conda:
        "envs/protein_preprocess_env.yml"
    threads: 10
    shell:
        """
        # Remove any stop codons (encoded as asterisks) from the end of sequences and replace any others with an ambiguous AA (X), ignoring the sequence name
        sed "s/*$//g" {input} > {output.tmp1}
        sed -E '/>/!s/\*/X/g' {output.tmp1} > {output.tmp2}
        # Remove sequences shorter than 25 AA.
        seqkit seq --threads {threads} --min-len 25 {output.tmp2} > {output.tmp3}
        """

rule preprocess_two:
    """
    This rule prepares data to be compatible with the NovelTree workflow. 
    This runs a python script that will calculate isoform protein lengths for entries that have 'yes' 
    in the isoform column of the input sheet and only keep the longest isoform.
    """
    input:
        fasta="outputs/tmp/{genussp}temp3.pep"
    output:
        filtered=temp("outputs/tmp/{genussp}temp4.pep")
    params:
        isoform=lambda wildcards: samples_df.xs(wildcards.genussp, level="genus_species")["isoform"].iloc[0],
        transdecoder=lambda wildcards: samples_df.xs(wildcards.genussp, level="genus_species")["isoform"].iloc[0]
    threads: 10
    shell: """
    if [ "{params.isoform}" == "yes" ] || [ "{params.transdecoder}" == "yes" ]; then
        python scripts/filter_fasta.py {input.fasta} {output.filtered}
    else
        cp {input.fasta} {output.filtered}
    fi
    """


rule preprocess_three:
    """
    This rule prepares data to be compatible with the NovelTree workflow. It will will cluster the sequences at 90% for transcriptomes/non-uniprot proteomes and 95% for uniprot proteomes, and will rename protein entries with the species name plus protein id.
    """
    input:
        fastain="outputs/tmp/{genussp}temp4.pep"
    output:
        tmp4clstr=temp("outputs/tmp/{genussp}temp4.pep.clstr"),
        fastaout="outputs/processedprots/{genussp}.fasta"
    conda:
        "envs/protein_preprocess_env.yml"
    params:
        genus=lambda wildcards: samples_df.loc[(slice(None), wildcards.genussp), "genus_species"].values[0],
        uniprot=lambda wildcards: samples_df.xs(wildcards.genussp, level="genus_species")["uniprot"].iloc[0]
    threads: 10
    shell:
        """
        # Use CD-HIT to reduce proteome redundancy, clustered at 90% for transcriptome and 95% for proteome
        cd-hit -T {threads} -c $(if [ "{params.uniprot}" == "yes" ]; then echo "0.95"; else echo "0.90"; fi) -i {input.fastain} -o {output.tmp4clstr}
        # Rename protein sequences in the format of "genus_species-proteinid"
        awk -v genus={params.genus} '/^>/ {{
          num = $0;
          sub(/^>/, "", num);      # Remove the ">" at the beginning
          gsub(/_/, "-", num);     # Replace all underscores with hyphens
          print ">"genus"_"num;
          next;
         }} 
        1' {output.tmp4clstr} > {output.fastaout}       
        """

rule download_kofamscan_ko_list:
    """
    This rule downloads the kofamscan KEGG list file.
    """
    output: "outputs/databases/kofamscandb/ko_list"
    shell:'''
    curl -JLo {output}.gz ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz && gunzip -c {output}.gz > {output}
    '''

rule download_kofamscan_profiles:
    """
    This rule downloads the kofamscan KEGG hmm profiles.
    """
    output: "outputs/databases/kofamscandb/profiles/eukaryote.hal"
    params: outdir = "outputs/databases/kofamscandb/"
    shell:'''
    curl -JLo {params.outdir}/profiles.tar.gz ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz && tar xf {params.outdir}/profiles.tar.gz -C {params.outdir}
    '''

rule kofamscan_proteome:
    """
    This rule uses the kofamscan to perform KEGG ortholog annotation on the proteomes. 
    """
    input:
        fa="outputs/processedprots/{genussp}.fasta",
        kolist="outputs/databases/kofamscandb/ko_list",
        profiles="outputs/databases/kofamscandb/profiles/eukaryote.hal"
    output: "outputs/annotated/kofamscan/{genussp}_kofamscan.tsv"
    conda: "envs/kofamscan.yml"
    params: profilesdir = "outputs/databases/kofamscandb/profiles"
    threads: 8
    shell:'''
    exec_annotation --format detail-tsv --ko-list {input.kolist} --profile {params.profilesdir} --cpu {threads} -o {output} {input.fa}
    '''

rule download_eggnog_db:
    """
    This rule downloads the eggnog annotation database.
    The script download_eggnog_data.py is exported by the eggnog mapper tool.
    """
    output: "outputs/databases/eggnog_db/eggnog.db"
    conda: "envs/eggnog.yml"
    shell:'''
    download_eggnog_data.py -H -d 2 -y --data_dir outputs/databases/eggnog_db
    '''

rule eggnog_proteome:
    '''
    This rule uses the EggNOG database to functionally annotate the HGT candidate genes. 
    It runs the EggNOG-Mapper tool on the translated candidate gene sequences, generating a file with the annotations (NOG, KEGG, PFAM, CAZys, EC numbers) for each gene. 
    The script emapper.py is exported by the eggnog mapper tool.
    '''
    input:
        db="outputs/databases/eggnog_db/eggnog.db",
        fa="outputs/processedprots/{genussp}.fasta"
    output: "outputs/annotated/eggnog/{genussp}.emapper.annotations",
    conda: "envs/eggnog.yml"
    params:
        outdir="outputs/annotated/eggnog/",
        dbdir = "outputs/databases/eggnog_db/" 
    threads: 8
    shell:'''
    mkdir -p tmp
    emapper.py --cpu {threads} -i {input.fa} --output {wildcards.genussp} \
       --output_dir {params.outdir} -m diamond --tax_scope none \
       --seed_ortholog_score 60 --override --temp_dir tmp/ \
       --data_dir {params.dbdir}
    '''

rule deepsig_proteome:
    '''
    This rule uses deepsig to predict signal peptides in proteins using deep learning.
    '''
    input:
        fa="outputs/processedprots/{genussp}.fasta"
    output: "outputs/annotated/deepsig/{genussp}_deepsig"
    conda: "envs/deepsig.yml"
    threads: 8
    shell:'''
    deepsig -f {input} -o {output} -k euk
    '''    

rule calc_protlengths:
    '''
    This rule calculates the lengths of proteins.
    '''
    input:
        "outputs/processedprots/{genussp}.fasta"
    output:
        "outputs/processedprots/lengths/{genussp}length.tsv"
    run:
        def fasta_to_tsv_length(fasta_file, tsv_file):
            with open(fasta_file, 'r') as f, open(tsv_file, 'w') as out:
                out.write("gene_name\tLength\n")  # Writing the headers

                sequence = ""
                protein_name = ""

                for line in f:
                    line = line.strip()
                    if line.startswith('>'):  # A header line in FASTA format
                         if protein_name:  # This checks if it's not the first sequence
                            out.write(f"{protein_name}\t{len(sequence)}\n")
                            sequence = ""
                         # Extract only the desired part of the protein name
                         protein_name = line[1:].split()[0] 
                    else:
                        sequence += line

                # Write the last sequence's length
                if protein_name and sequence:
                    out.write(f"{protein_name}\t{len(sequence)}\n")

        fasta_to_tsv_length(input[0], output[0])

rule combine_annotations:
    '''
    This rule combines the EggNOG, DeepSig, and KofamScan, and length annotations into a single output.
    '''
    input:
        egg="outputs/annotated/eggnog/{genussp}.emapper.annotations",
        ko="outputs/annotated/kofamscan/{genussp}_kofamscan.tsv",
        deep="outputs/annotated/deepsig/{genussp}_deepsig",
        length="outputs/processedprots/lengths/{genussp}length.tsv"
    output:
        "outputs/annotated/{genussp}_combined.tsv"
    run:
                # Read TSVs into pandas DataFrames
        df1 = pd.read_csv(input.egg, sep='\t',skiprows=4,comment=None)
        df2 = pd.read_csv(input.ko, sep='\t',comment=None)
        df3 = pd.read_csv(input.deep, sep='\t',header=None)
        df4 = pd.read_csv(input.length, sep='\t',comment=None)
        df4['Length'] = df4['Length'].astype(str)

        #filter kofamscan annotations to keep max top five scores 
        df2 = df2.drop(df2.index[0])
        df2['score'] = df2['score'].astype(float)
        df2 = df2.groupby('gene name').head(5).reset_index(drop=True)      
        df2['score'] = df2['score'].astype(str)
    
        # Assign column names to df1 and df2
        df2.columns = ['KO_pass','gene_name', 'KO', 'KO_thrshld', 'KO_score', 'KO_E-value', 'KO_definition']
        df1.columns = ['gene_name', "egg_seed_ortholog", "egg_evalue", "egg_score", "eggNOG_OGs","egg_max_annot_lvl", "egg_COG_category", "egg_Description", "egg_Preferred_name","egg_GOs", "egg_EC",
 "egg_KEGG_ko", "egg_KEGG_Pathway", "egg_KEGG_Module", "egg_KEGG_Reaction","egg_KEGG_rclass", "egg_BRITE", "egg_KEGG_TC", "egg_CAZy", "egg_BiGG_Reaction","egg_PFAMs"]        

        # make one row per gene for kofam table
        # Define a function to aggregate columns
        def concatenate_unique(x):
            # join rows for same gene with ";"
            values = ['nan' if pd.isna(val) else str(val) for val in x]
            return ";".join(values)
        # Group by 'gene_name' and aggregate other columns
        df2 = df2.groupby('gene_name').agg(concatenate_unique).reset_index()
         
        # Drop unnecessary columns and assign column names to deepsig df
        df3 = df3.drop(df3.columns[[1, 6, 7]], axis=1)
        df3 = df3.astype(str)
        df3.columns = ['gene_name', 'deepsig_feature', 'deepsig_start','deepsig_end','deepsig_sp_score','deepsig_sp_evidence']
        
        #merge data together, keep a row for every protein 
        base_df = df4[['gene_name']].copy()
        unique_gene_names_count = base_df['gene_name'].nunique()
        merged_df = df4.merge(df3, on='gene_name', how='outer').merge(df2, on='gene_name', how='outer').merge(df1, on='gene_name', how='outer')

        # Select desired columns
        columns_to_keep = df1.columns.tolist() + ['KO_pass','KO', 'KO_thrshld', 'KO_score', 'KO_E-value', 'KO_definition','deepsig_feature', 'deepsig_start','deepsig_end','deepsig_sp_score','deepsig_sp_evidence','Length']
        result = merged_df[columns_to_keep]

        # Handle duplicates introduced by deepsig. First get names that are duplicated
        duplicates = result[result['gene_name'].duplicated(keep=False)]['gene_name'].unique()
        # For those duplicates, keep only the ones with 'Signal peptide'
        filtered_rows = result[(result['gene_name'].isin(duplicates)) & (result['deepsig_feature'] == 'Signal peptide')]
        # Get the non-duplicated rows
        non_duplicates = result[~result['gene_name'].isin(duplicates)]
        # Combine the two to get the result
        resultx = pd.concat([filtered_rows, non_duplicates], ignore_index=True)
        resultx = resultx.iloc[:-3]
        # Save the processed dataframe as a TSV
        resultx.to_csv(output[0], sep='\t', index=False)
