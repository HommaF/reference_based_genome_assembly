configfile: "config.json"

rule all:
	input:
		expand("results/bwa_aligned/{reference_genome}_{sample}_done.txt", sample=config["samples"], reference_genome=config["reference_genome"])

rule trimmomatic:
	input:
		fwd="genomic_reads/{sample}_R1.fastq.gz",
		rev="genomic_reads/{sample}_R2.fastq.gz"
	output:
		fwd_paired="results/trimmomatic/{sample}_R1_paired.fastq.gz",
		fwd_unpaired="results/trimmomatic/{sample}_R1_unpaired.fastq.gz",
		rev_paired="results/trimmomatic/{sample}_R2_paired.fastq.gz",
		rev_unpaired="results/trimmomatic/{sample}_R2_unpaired.fastq.gz"

	conda:
		"envs/dn_assembly.yaml"

	threads: 6
	shell:
		"trimmomatic PE -threads 12 {input.fwd} {input.rev} {output.fwd_paired} {output.fwd_unpaired} {output.rev_paired} {output.rev_unpaired} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule fastQC:
	input:
		fwd_paired="results/trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="results/trimmomatic/{sample}_R2_paired.fastq.gz"
	output:
		directory("results/fastqc/{sample}")
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
		fwd_paired="results/trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="results/trimmomatic/{sample}_R2_paired.fastq.gz",
		genome="reference_genome/{reference_genome}_filtered.fasta",
		fastqc="results/fastqc/{sample}/"

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
                "results/bwa_aligned/{reference_genome}_{sample}_sorted.bam"
	output:
                "results/bwa_aligned/{reference_genome}_{sample}_sorted.bam.bai"

	threads: 1
	conda:
		"envs/dn_assembly.yaml"
	shell:
		"samtools index {input}"

rule extract_mapped:
	input:
                bwa_in = "results/bwa_aligned/{reference_genome}_{sample}_sorted.bam",
                bai_in = "results/bwa_aligned/{reference_genome}_{sample}_sorted.bam.bai"
	output:
		bwa_out = "results/bwa_aligned/{reference_genome}_{sample}_mapped.bam",
		bai_out = "results/bwa_aligned/{reference_genome}_{sample}_mapped.bam.bai"
	threads: 4
	conda:
		"envs/dn_assembly.yaml"
	shell:
		"samtools view -b -@ 3 -F 4 -F 8 {input.bwa_in} > {output.bwa_out}; "
		"samtools index {output.bwa_out}"

rule filter_mapped:
	input:
		mapped = "results/bwa_aligned/{reference_genome}_{sample}_mapped.bam",
                bwa_in = "results/bwa_aligned/{reference_genome}_{sample}_sorted.bam"
	output:
		mapped = "results/bwa_aligned/{reference_genome}_{sample}_mapped_filtered.bam",
		bwa_in = "results/bwa_aligned/{reference_genome}_{sample}_sorted_filtered.bam",
		bwa_in_index = "results/bwa_aligned/{reference_genome}_{sample}_sorted_filtered.bam.bai"

	threads: 4
	conda:
		"envs/dn_assembly.yaml"
	shell:
		"samtools view -b -@ 3 -q 10 {input.mapped} > {output.mapped};"
		"samtools view -b -@ 3 -q 10 {input.bwa_in} > {output.bwa_in};"
		"samtools index -@ 3 {output.bwa_in_index}"

rule mapped_genomecov:
	input:
		"results/bwa_aligned/{reference_genome}_{sample}_mapped_filtered.bam"
	output:
		"results/bwa_aligned/{reference_genome}_{sample}_mapped_filtered_histogram.tsv.gz"

	threads: 1
	conda: 
		"envs/dn_assembly.yaml"

	shell:
		"bedtools genomecov -ibam {input} -bga | bgzip > {output}"

rule propper_mapped_genomecov:
	input:
		"results/bwa_aligned/{reference_genome}_{sample}_mapped_filtered.bam"
	output:
		"results/bwa_aligned/{reference_genome}_{sample}_proper_mapped_filtered_histogram.tsv.gz"

	params: ref_genome_index = "reference_genome/{reference_genome}.fasta.fai"

	threads: 1

	conda:
		"envs/dn_assembly.yaml"
	
	shell:
		"samtools view -bf 0x2 {input} | samtools sort - -n | bedtools bamtobed -i - -bedpe | awk '$1 == $4' | cut -f1,2,6 | sort -k 1,1 | bedtools genomecov -i - -bga -g {params.ref_genome_index} | bgzip > {output}"
	

rule define_superblocks:
	input:
		"results/bwa_aligned/{reference_genome}_{sample}_proper_mapped_filtered_histogram.tsv.gz"
	output:
		tmp_blocks = "results/blocks/{reference_genome}_{sample}_tmp_blocks.tsv",
		blocks = "results/blocks/{reference_genome}_{sample}_blocks.tsv",
		superblocks = "results/blocks/{reference_genome}_{sample}_superblocks.tsv"

	threads: 12
	conda:
		"envs/dn_assembly.yaml"

	shell:
		"scripts/gen_superblocks.sh {input} {output.tmp_blocks};"
		"scripts/gen_superblocks.py {output.tmp_blocks} {output.blocks} {output.superblocks};"
		"scripts/superblock_reads.sh {output.superblocks}"
	
rule extract_unmapped:
	input:
		bwa_in = "results/bwa_aligned/{reference_genome}_{sample}_sorted.bam",
		bai_in = "results/bwa_aligned/{reference_genome}_{sample}_sorted.bam.bai",
		fwd_paired="results/trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="results/trimmomatic/{sample}_R2_paired.fastq.gz"
	output:
		unm_fwd = "results/bwa_aligned/{reference_genome}_{sample}_unmapped_R1.fq.gz",
		unm_rev = "results/bwa_aligned/{reference_genome}_{sample}_unmapped_R2.fq.gz"

	
	params:
		unmapped_ids = "results/bwa_aligned/{reference_genome}_{sample}_unmapped_ids.txt"

	threads: 1
	conda:
		"envs/dn_assembly.yaml"
	shell:
		"samtools view -b -f 4 {input.bwa_in} | samtools view -f 9 - | cut -f1 > {params.unmapped_ids}; "
		"printf '1\n2\n' | xargs -I@ -n 1 -P 2 bash -c 'scripts/extract_fastq.sh @ {params.unmapped_ids} {input.fwd_paired} {input.rev_paired} {output.unm_fwd} {output.unm_rev}'; "
		"rm {params.unmapped_ids}"


rule combine:
	input:
		unm_rev = "results/bwa_aligned/{reference_genome}_{sample}_unmapped_R2.fq.gz",
		bam = "results/bwa_aligned/{reference_genome}_{sample}_mapped_filtered.bam",
		mapped_filtered = "results/bwa_aligned/{reference_genome}_{sample}_mapped_filtered_histogram.tsv.gz",
		proper_mapped_filtered = "results/bwa_aligned/{reference_genome}_{sample}_proper_mapped_filtered_histogram.tsv.gz",
		superblocks = "results/blocks/{reference_genome}_{sample}_superblocks.tsv"

	output:
		"results/bwa_aligned/{reference_genome}_{sample}_done.txt"
	threads: 1

	shell:
		"touch {output}"
