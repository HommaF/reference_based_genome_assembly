configfile: "config.json"

rule all:
	input:
		expand("bwa_aligned/{reference_genome}_{sample}_sorted.bam", sample=config["samples"], reference_genome=config["reference_genome"])

rule trimmomatic:
	input:
		fwd="genomic_reads/{sample}_R1.fastq.gz",
		rev="genomic_reads/{sample}_R2.fastq.gz"
	output:
		fwd_paired="trimmomatic/{sample}_R1_paired.fastq.gz",
		fwd_unpaired="trimmomatic/{sample}_R1_unpaired.fastq.gz",
		rev_paired="trimmomatic/{sample}_R2_paired.fastq.gz",
		rev_unpaired="trimmomatic/{sample}_R2_unpaired.fastq.gz"

	conda:
		"envs/dn_assembly.yaml"

	threads: 6
	shell:
		"trimmomatic PE -threads 12 {input.fwd} {input.rev} {output.fwd_paired} {output.fwd_unpaired} {output.rev_paired} {output.rev_unpaired} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule fastQC:
	input:
		fwd_paired="trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="trimmomatic/{sample}_R2_paired.fastq.gz"
	output:
		directory("fastqc/{sample}")
	threads: 2
	conda:
		"envs/dn_assembly.yaml"
	shell:
		"mkdir {output} && fastqc --quiet --threads 2 -o {output} -f fastq {input.fwd_paired} {input.rev_paired}"


rule filter_scaffold:
	input:
		"reference_genome/{reference_genome}.fasta"
	output:
		filtered = "reference_genome/{reference_genome}_filtered.fasta",
		index = "reference_genome/{reference_genome}_filtered.fasta.fai"
	threads: 2
	conda:
		"envs/dn_assembly.yaml"

	shell:
		"seqkit seq -m 10000 {input} > {output.filtered}; "
		"samtools faidx {output.filtered}"
rule bwa_align:
	input:
		fwd_paired="trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="trimmomatic/{sample}_R2_paired.fastq.gz",
		genome="reference_genome/{reference_genome}_filtered.fasta",
		fastqc="fastqc/{sample}/"

	output:
		"bwa_aligned/{reference_genome}_{sample}_sorted.bam"

	params:
		rg="@RG\\tID:{sample}\\tSM:{sample}"


	threads: 6
	conda:
		"envs/dn_assembly.yaml"

	shell:
		"bwa-mem2 -R '{params.rg}' -t 6 {input.genome} {input.fwd_paired} {input.rev_paired} | samtools view -bS - | samtools sort - -o {output}"
