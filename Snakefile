
############################################
# SAMPLES
############################################
import os

def get_input(sample):
    if os.path.exists(f"data/{sample}_1.fastq.gz"):
        return {
            "type": "fastq",
            "r1": f"data/{sample}_1.fastq.gz",
            "r2": f"data/{sample}_2.fastq.gz"
        }
    elif os.path.exists(f"data/{sample}.fasta"):
        return {
            "type": "fasta",
            "fasta": f"data/{sample}.fasta"
        }
    else:
        raise ValueError(f"No input found for {sample}")

FASTQ_SAMPLES = ["ERR3771946"]
FASTA_SAMPLES = ["Bsubtilis"]
SAMPLES = FASTQ_SAMPLES + FASTA_SAMPLES
############################################
# RULE ALL
############################################

rule all:
    input:
        # FASTQ only
        expand("results/fastp/{sample}_1.trim.fastq.gz", sample=FASTQ_SAMPLES),
        expand("results/fastp/{sample}_2.trim.fastq.gz", sample=FASTQ_SAMPLES),
        expand("results/spades/{sample}/contigs.fasta", sample=FASTQ_SAMPLES),
        expand("results/fastqc_trimmed/{sample}_1.trim_fastqc.html", sample=FASTQ_SAMPLES),
        expand("results/fastqc_trimmed/{sample}_2.trim_fastqc.html", sample=FASTQ_SAMPLES),

        # BOTH
        expand("results/prokka/{sample}/{sample}.gff", sample=SAMPLES),
        "results/multiqc/multiqc_report.html",
        expand("results/jellyfish/{sample}_histo.txt", sample=SAMPLES),
        expand("results/abricate/{sample}/{sample}.tab", sample=SAMPLES),
        expand("results/mlst/{sample}/{sample}.tsv", sample=SAMPLES),
        expand("results/diamond/{sample}.tsv", sample=SAMPLES),
        expand("results/hmmer/{sample}_pfam.tsv", sample=SAMPLES),
        expand("results/antismash/{sample}", sample=SAMPLES),
        expand("results/busco/{sample}", sample=SAMPLES),
        expand("results/trnascan/{sample}_trna.tsv", sample=SAMPLES),
        expand("results/barrnap/{sample}/{sample}.gff", sample=SAMPLES),
        expand("results/trf/{sample}/{sample}.dat", sample=SAMPLES),
        expand("results/kofam/{sample}_kofam.txt", sample=SAMPLES),
        expand("results/checkm/{sample}/qa_summary.tsv", sample=SAMPLES),
        expand("results/quast/{sample}/report.html", sample=SAMPLES),
############################################
# FASTP
############################################
rule fastp:
    output:
        r1="results/fastp/{sample}_1.trim.fastq.gz",
        r2="results/fastp/{sample}_2.trim.fastq.gz"
    conda:
        "envs/fastp.yaml"
    threads: 2
    wildcard_constraints:
        sample="|".join(FASTQ_SAMPLES)
    run:
        inp = get_input(wildcards.sample)

        shell(f"""
        mkdir -p results/fastp
        fastp \
            -i {inp["r1"]} \
            -I {inp["r2"]} \
            -o {output.r1} \
            -O {output.r2}
        """)
############################################
# FASTQC
############################################
rule fastqc_trimmed:
    input:
        r1="results/fastp/{sample}_1.trim.fastq.gz",
        r2="results/fastp/{sample}_2.trim.fastq.gz"
    output:
        html1="results/fastqc_trimmed/{sample}_1.trim_fastqc.html",
        html2="results/fastqc_trimmed/{sample}_2.trim_fastqc.html"
    conda:
        "envs/fastqc.yaml"
    wildcard_constraints:
        sample="|".join(FASTQ_SAMPLES)   # ✅ KEY FIX
    shell:
        """
        mkdir -p results/fastqc_trimmed
        fastqc {input.r1} {input.r2} -o results/fastqc_trimmed
        """
############################################
# MULTIQC
############################################
rule multiqc:
    input:
        expand("results/fastqc_trimmed/{sample}_1.trim_fastqc.html", sample=SAMPLES),
        expand("results/fastqc_trimmed/{sample}_2.trim_fastqc.html", sample=SAMPLES)
    output:
        "results/multiqc/multiqc_report.html"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        mkdir -p results/multiqc
        multiqc results/fastqc_trimmed -o results/multiqc --force
        """
############################################
# SPADES
############################################
rule spades:
    input: []
    output:
        "results/spades/{sample}/contigs.fasta"
    conda:
        "envs/spades.yaml"
    threads: 4
    run:
        inp = get_input(wildcards.sample)

        if inp["type"] == "fastq":
            shell(f"""
                mkdir -p results/spades/{wildcards.sample}
                spades.py \
                    -1 {inp["r1"]} \
                    -2 {inp["r2"]} \
                    -o results/spades/{wildcards.sample}
            """)
        else:
            shell(f"""
                mkdir -p results/spades/{wildcards.sample}
                cp {inp["fasta"]} {output}
            """)
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
        "results/diamond/{sample}.tsv"
    conda:
        "envs/diamond.yaml"
    threads: 4
    shell:
        """
        mkdir -p results/diamond
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
        "results/hmmer/{sample}_pfam.tsv"
    conda:
        "envs/hmmer.yaml"
    threads: 4
    shell:
        """
        mkdir -p results/hmmer
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
        "results/kofam/{sample}_kofam.txt"
    conda:
        "envs/kofamscan.yaml"
    threads: 4
    shell:
        """
        mkdir -p results/kofam
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
        "results/trnascan/{sample}_trna.tsv"
    conda:
        "envs/trnascan.yaml"
    shell:
        """
        mkdir -p results/trnascan
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
        r1="results/fastp/{sample}_1.trim.fastq.gz",
        r2="results/fastp/{sample}_2.trim.fastq.gz"
    output:
        "results/jellyfish/{sample}_histo.txt"
    conda:
        "envs/jellyfish.yaml"
    shell:
        """
        mkdir -p results/jellyfish
        jellyfish count -m 21 -s 100M -t 2 \
            <(zcat {input.r1}) <(zcat {input.r2}) -o jf.tmp
        jellyfish histo jf.tmp > {output}
        rm jf.tmp
        """
