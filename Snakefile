configfile: "config.json"

rule all:
	input:
		expand("bwa_aligned/{reference_genome}_{sample}_done.txt", sample=config["samples"], reference_genome=config["reference_genome"])

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
		index = "reference_genome/{reference_genome}_filtered.bwt"
	threads: 2
	conda:
		"envs/dn_assembly.yaml"

	shell:
		"seqkit seq -m 10000 {input} > {output.filtered}; "
		"bwa index {output.filtered}"
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


	threads: 10
	conda:
		"envs/dn_assembly.yaml"

	shell:
		"bwa mem -R '{params.rg}' -t 10 {input.genome} {input.fwd_paired} {input.rev_paired} | samtools view -bS - | samtools sort - -o {output}"

rule index_mapped:
	input:
                "bwa_aligned/{reference_genome}_{sample}_sorted.bam"
	output:
                "bwa_aligned/{reference_genome}_{sample}_sorted.bam.bai"

	threads: 1
	conda:
		"envs/dn_assembly.yaml"
	shell:
		"samtools index {input}"

rule extract_mapped:
	input:
                bwa_in = "bwa_aligned/{reference_genome}_{sample}_sorted.bam",
                bai_in = "bwa_aligned/{reference_genome}_{sample}_sorted.bam.bai"
	output:
		bwa_out = "bwa_aligned/{reference_genome}_{sample}_mapped.bam",
		bai_out = "bwa_aligned/{reference_genome}_{sample}_mapped.bam.bai"
	threads: 4
	conda:
		"envs/dn_assembly.yaml"
	shell:
		"samtools view -b -@ 3 -F 4 {input.bwa_in} > {output.bwa_out}; "
		"samtools index {output.bwa_out}"

rule filter_mapped:
	input:
		"bwa_aligned/{reference_genome}_{sample}_mapped.bam"
	output:
		"bwa_aligned/{reference_genome}_{sample}_mapped_filtered.bam"
	threads: 4
	conda:
		"envs/dn_assembly.yaml"
	shell:
		"samtools view -b -@ 3 -q 10 {input} > {output}"

rule extract_unmapped:
	input:
		bwa_in = "bwa_aligned/{reference_genome}_{sample}_sorted.bam",
		bai_in = "bwa_aligned/{reference_genome}_{sample}_sorted.bam.bai",
		fwd_paired="trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="trimmomatic/{sample}_R2_paired.fastq.gz"
	output:
		unm_fwd = "bwa_aligned/{reference_genome}_{sample}_unmapped_R1.fq.gz",
		unm_rev = "bwa_aligned/{reference_genome}_{sample}_unmapped_R2.fq.gz"

	
	params:
		unmapped_ids = "bwa_aligned/{reference_genome}_{sample}_unmapped_ids.txt"

	threads: 1
	conda:
		"envs/dn_assembly.yaml"
	shell:
		"samtools view -b -f 4 {input.bwa_in} | samtools view -f 9 - | cut -f1 > {params.unmapped_ids}; "
		"printf '1\n2\n' | xargs -I@ -n 1 -P 2 bash -c 'scripts/extract_fastq.sh @ {params.unmapped_ids} {input.fwd_paired} {input.rev_paired} {output.unm_fwd} {output.unm_rev}'; "
		"rm {params.unmapped_ids}"


rule combine:
	input:
		unm_rev = "bwa_aligned/{reference_genome}_{sample}_unmapped_R2.fq.gz",
		bam = "bwa_aligned/{reference_genome}_{sample}_mapped_filtered.bam"

	output:
		"bwa_aligned/{reference_genome}_{sample}_done.txt"
	threads: 1

	shell:
		"touch {output}"
