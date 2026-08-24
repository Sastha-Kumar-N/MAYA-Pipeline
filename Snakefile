
# Define your samples here
# FASTQ: Paired-end reads (will look for sample_1.fastq.gz and sample_2.fastq.gz in data/)
FASTQ_SAMPLES = ["fastq_SRR15006059"]
FASTQ_INPUT_SAMPLES = {"fastq_SRR15006059": "SRR15006059"}

# FASTA: Pre-assembled contigs (will look for sample.fasta in data/)
FASTA_SAMPLES = []

# Combine samples for analysis
ANALYSIS_SAMPLES = FASTQ_SAMPLES + FASTA_SAMPLES

# Build wildcard constraints dynamically
FASTQ_PATTERN = "|".join(FASTQ_SAMPLES) if FASTQ_SAMPLES else "(?!.*)"
FASTA_PATTERN = "|".join(FASTA_SAMPLES) if FASTA_SAMPLES else "(?!.*)"

############################################
# RULE ALL
############################################

rule all:
    input:
        # FASTQ only (quality control and assembly)
        *(expand("results/fastp/{sample}/{sample}_1.trim.fastq.gz", sample=FASTQ_SAMPLES) if FASTQ_SAMPLES else []),
        *(expand("results/fastp/{sample}/{sample}_2.trim.fastq.gz", sample=FASTQ_SAMPLES) if FASTQ_SAMPLES else []),
        *(expand("results/spades/{sample}/contigs.fasta", sample=FASTQ_SAMPLES) if FASTQ_SAMPLES else []),
        *(expand("results/fastqc_trimmed/{sample}/{sample}_1.trim_fastqc.html", sample=FASTQ_SAMPLES) if FASTQ_SAMPLES else []),
        *(expand("results/fastqc_trimmed/{sample}/{sample}_2.trim_fastqc.html", sample=FASTQ_SAMPLES) if FASTQ_SAMPLES else []),
        *(expand("results/jellyfish/{sample}/{sample}_histo.txt", sample=FASTQ_SAMPLES) if FASTQ_SAMPLES else []),
        *(expand("results/multiqc/{sample}/multiqc_report.html", sample=FASTQ_SAMPLES) if FASTQ_SAMPLES else []),
        
        # FASTA samples (copy to standard location)
        *(expand("results/spades/{sample}/contigs.fasta", sample=FASTA_SAMPLES) if FASTA_SAMPLES else []),

        # Analysis on contigs (from either FASTQ or FASTA)
        *(expand("results/prokka/{sample}/{sample}.gff", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/abricate/{sample}/{sample}.tab", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/mlst/{sample}/{sample}.tsv", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/diamond/{sample}/{sample}.tsv", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/hmmer/{sample}/{sample}_pfam.tsv", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/antismash/{sample}", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/busco/{sample}", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/trnascan/{sample}/{sample}_trna.tsv", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/barrnap/{sample}/{sample}.gff", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/trf/{sample}/{sample}.dat", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/kofam/{sample}/{sample}_kofam.txt", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/checkm/{sample}/qa_summary.tsv", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/quast/{sample}/report.html", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),
        *(expand("results/mob_recon/{sample}", sample=ANALYSIS_SAMPLES) if ANALYSIS_SAMPLES else []),

############################################
# FASTP
############################################
rule fastp:
    input:
        r1=lambda wildcards: f"data/{FASTQ_INPUT_SAMPLES[wildcards.sample]}_1.fastq",
        r2=lambda wildcards: f"data/{FASTQ_INPUT_SAMPLES[wildcards.sample]}_2.fastq"

    output:
        r1="results/fastp/{sample}/{sample}_1.trim.fastq.gz",
        r2="results/fastp/{sample}/{sample}_2.trim.fastq.gz"

    conda:
        "envs/fastp.yaml"

    threads: 2
    wildcard_constraints:
        sample=FASTQ_PATTERN

    shell:
        """
        mkdir -p results/fastp/{wildcards.sample}

        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2}
        """
############################################
# FASTQC
############################################
rule fastqc_trimmed:
    input:
        r1="results/fastp/{sample}/{sample}_1.trim.fastq.gz",
        r2="results/fastp/{sample}/{sample}_2.trim.fastq.gz"
    output:
        html1="results/fastqc_trimmed/{sample}/{sample}_1.trim_fastqc.html",
        html2="results/fastqc_trimmed/{sample}/{sample}_2.trim_fastqc.html"
    conda:
        "envs/fastqc.yaml"
    wildcard_constraints:
        sample=FASTQ_PATTERN
    shell:
        """
        mkdir -p results/fastqc_trimmed/{wildcards.sample}
        fastqc {input.r1} {input.r2} -o results/fastqc_trimmed/{wildcards.sample}
        """
############################################
# MULTIQC
############################################
rule multiqc:
    input:
        expand("results/fastqc_trimmed/{sample}/{sample}_1.trim_fastqc.html", sample=FASTQ_SAMPLES),
        expand("results/fastqc_trimmed/{sample}/{sample}_2.trim_fastqc.html", sample=FASTQ_SAMPLES)
    output:
        "results/multiqc/fastq_SRR15006059/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        mkdir -p results/multiqc
        mkdir -p results/multiqc/fastq_SRR15006059
        multiqc results/fastqc_trimmed/fastq_SRR15006059 -o results/multiqc/fastq_SRR15006059 --force
        """
####################################
# SPADES FOR FASTQ
####################################

rule spades_fastq:
    input:
        r1="results/fastp/{sample}/{sample}_1.trim.fastq.gz",
        r2="results/fastp/{sample}/{sample}_2.trim.fastq.gz"

    output:
        "results/spades/{sample}/contigs.fasta"

    conda:
        "envs/spades.yaml"

    threads: 4
    wildcard_constraints:
        sample=FASTQ_PATTERN

    shell:
        """
        mkdir -p results/spades/{wildcards.sample}

        spades.py --careful \
            -1 {input.r1} \
            -2 {input.r2} \
            -o results/spades/{wildcards.sample}
        """


####################################
# FASTA INPUT
####################################

rule fasta_copy:
    input:
        "data/{sample}.fasta"
    wildcard_constraints:
        sample=FASTA_PATTERN

    output:
        "results/spades/{sample}/contigs.fasta"

    shell:
        """
        mkdir -p results/spades/{wildcards.sample}

        cp {input} {output}
        """
        
           
       
############################################
# QUAST
############################################

rule quast:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        "results/quast/{sample}/report.html"
    conda:
        "envs/quast.yaml"
    shell:
        """
        mkdir -p results/quast/{wildcards.sample}
        quast {input} -o results/quast/{wildcards.sample}
        """

############################################
# ABRICATE
############################################

rule abricate:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        "results/abricate/{sample}/{sample}.tab"
    conda:
        "envs/abricate.yaml"
    shell:
        """
        mkdir -p results/abricate/{wildcards.sample}
        abricate {input} > {output}
        """

############################################
# MLST
############################################

rule mlst:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        "results/mlst/{sample}/{sample}.tsv"
    conda:
        "envs/mlst.yaml"
    shell:
        """
        mkdir -p results/mlst/{wildcards.sample}
        mlst {input} > {output}
        """

############################################
# DIAMOND
############################################

rule diamond:
    input:
        faa="results/prokka/{sample}/{sample}.faa",
        db="databases/swissprot.dmnd"
    output:
        "results/diamond/{sample}/{sample}.tsv"
    conda:
        "envs/diamond.yaml"
    threads: 4
    shell:
        """
        mkdir -p results/diamond/{wildcards.sample}
        diamond blastp \
            -d {input.db} \
            -q {input.faa} \
            -o {output} \
            -f 6 \
            -k 1 \
            --threads {threads}
        """

############################################
# HMMER
############################################

rule hmmer:
    input:
        faa="results/prokka/{sample}/{sample}.faa",
        db="databases/Pfam-A.hmm"
    output:
        "results/hmmer/{sample}/{sample}_pfam.tsv"
    conda:
        "envs/hmmer.yaml"
    threads: 4
    shell:
        """
        mkdir -p results/hmmer/{wildcards.sample}
        hmmscan --cpu {threads} \
            --domtblout {output} \
            {input.db} {input.faa}
        """

############################################
# CHECKM
############################################

rule checkm:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        "results/checkm/{sample}/qa_summary.tsv"
    conda:
        "envs/checkm.yaml"
    threads: 1
    shell:
        """
        mkdir -p results/checkm/{wildcards.sample}

        checkm taxonomy_wf domain Bacteria \
            -t 1 \
            -x fasta \
            results/spades/{wildcards.sample} \
            results/checkm/{wildcards.sample}

        checkm qa \
            results/checkm/{wildcards.sample}/Bacteria.ms \
            results/checkm/{wildcards.sample} \
            -o 2 > {output}
        """
############################################
# PROKKA
############################################
rule prokka:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        gff="results/prokka/{sample}/{sample}.gff",
        faa="results/prokka/{sample}/{sample}.faa"
    conda:
        "envs/prokka.yaml"
    threads: 2
    shell:
        """
        prokka \
            --outdir results/prokka/{wildcards.sample} \
            --prefix {wildcards.sample} \
            --force \
            {input}
        """
############################################
# KOFAMSCAN
############################################

rule kofamscan:
    input:
        faa="results/prokka/{sample}/{sample}.faa"
    output:
        "results/kofam/{sample}/{sample}_kofam.txt"
    conda:
        "envs/kofamscan.yaml"
    threads: 4
    shell:
        """
        mkdir -p results/kofam/{wildcards.sample}
        exec_annotation --cpu {threads} \
            -o {output} \
            --profile kofam_db/profiles \
            --ko-list kofam_db/ko_list \
            {input.faa}
        """

############################################
# ANTISMASH
############################################
rule antismash:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        directory("results/antismash/{sample}")
    conda:
        "envs/antismash.yaml"
    threads: 1
    shell:
        """
        mkdir -p {output}

        antismash \
            --genefinding-tool prodigal \
            --databases /home/shabari/.conda/envs/test_antismash/lib/python3.11/site-packages/antismash/databases \
            --output-dir {output} \
            {input}
        """

############################################
# BUSCO
############################################

rule busco:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        directory("results/busco/{sample}")
    conda:
        "envs/busco.yaml"
    threads: 4
    shell:
        """
        busco \
            -i {input} \
            -o {wildcards.sample} \
            -l bacteria_odb10 \
            -m genome \
            -c {threads} \
            --out_path results/busco
        """

############################################
# TRNASCAN
############################################

rule trnascan:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        "results/trnascan/{sample}/{sample}_trna.tsv"
    conda:
        "envs/trnascan.yaml"
    shell:
        """
        mkdir -p results/trnascan/{wildcards.sample}
        tRNAscan-SE -B {input} -o {output}
        """

############################################
# BARRNAP
############################################

rule barrnap:
    input:
        contigs="results/spades/{sample}/contigs.fasta"
    output:
        "results/barrnap/{sample}/{sample}.gff"
    conda:
        "envs/barrnap.yaml"
    shell:
        """
        mkdir -p results/barrnap/{wildcards.sample}
        barrnap {input.contigs} > {output}
        """
############################################
# TRF
############################################

rule trf:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        "results/trf/{sample}/{sample}.dat"
    conda:
        "envs/trf.yaml"
    shell:
        """
        mkdir -p results/trf/{wildcards.sample}
        cp {input} {wildcards.sample}.temp.fasta
        trf {wildcards.sample}.temp.fasta 2 7 7 80 10 50 500 -d -h

        if ls {wildcards.sample}.temp.fasta.*.dat 1> /dev/null 2>&1; then
            mv {wildcards.sample}.temp.fasta.*.dat {output}
        else
            touch {output}
        fi

        rm -f {wildcards.sample}.temp.fasta*
        """

############################################
# JELLYFISH
############################################
rule jellyfish:
    input:
        r1="results/fastp/{sample}/{sample}_1.trim.fastq.gz",
        r2="results/fastp/{sample}/{sample}_2.trim.fastq.gz"

    output:
        "results/jellyfish/{sample}/{sample}_histo.txt"

    conda:
        "envs/jellyfish.yaml"
    
    wildcard_constraints:
        sample=FASTQ_PATTERN

    shell:
        """
        mkdir -p results/jellyfish/{wildcards.sample}

        zcat {input.r1} > r1.fastq
        zcat {input.r2} > r2.fastq

        jellyfish count -m 21 -s 100M -t 2 \
        r1.fastq r2.fastq -o jf.tmp

        rm r1.fastq r2.fastq

        jellyfish histo jf.tmp > {output}

        rm jf.tmp
        """

############################################
# MOB-SUITE
############################################
rule mob_recon:
    input:
        "results/spades/{sample}/contigs.fasta"
    output:
        directory("results/mob_recon/{sample}")
    conda:
        "envs/mob_suite.yaml"
    threads: 2
    shell:
        """
        rm -rf {output}
        mkdir -p {output}

        mob_recon \
            -i {input} \
            -o {output} \
            -n {threads} \
            -f
        """
     


           

 

        

  

        
  
        
        
