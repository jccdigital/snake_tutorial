SAMPLES, = glob_wildcards("fastq/{sample}.fastq.gz")

rule all:
	input:
		"analyses/results/counts.txt"
rule trimse:
	input:
		"fastq/{sample}.fastq.gz" 
	output:
		"analyses/results/{sample}.trimmed.fastq.gz" 
	log:
		"analyses/logs/{sample}.trimse" 
	benchmark:
		"analyses/benchmarks/{sample}.trimse"
	params:
		qualtype="sanger" 
	shell:
		"sickle se -g -t {params.qualtype} -f {input} -o {output}"
		" 2> {log}"
rule makeidx: 
	input:
		fasta = "data/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa" 
	output:
 		touch("data/makeidx.done")	 
	log:
		"analyses/logs/makeidx.log" 
	benchmark:
		"analyses/benchmarks/makeidx.txt" 
	conda:
		"envs/map.yaml" 
	shell:
		"bwa index {input.fasta} 2> {log}"

rule map:
	input:
		reads = "analyses/results/{sample}.trimmed.fastq.gz",
 		idxdone = "data/makeidx.done"	 
	output:
		"analyses/results/{sample}.bam" 
	log:
		"analyses/logs/{sample}.map" 
	benchmark:
		"analyses/benchmarks/{sample}.map"
	threads: 8	 
	conda:
		"envs/map.yaml" 
	params:
		idx = "data/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa" 
	shell:
		"bwa mem -t {threads} {params.idx} {input.reads} | samtools view -Sbh > {output} 2> {log}"

rule featurecount: 
	input:
		gtf = "data/Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz",
 		bams = expand("analyses/results/{sample}.bam", sample=SAMPLES)	 
	output:
		"analyses/results/counts.txt"
	log:
		"analyses/logs/featurecount.log" 
	benchmark:
		"analyses/benchmarks/featurecount.txt" 
	conda:
		"envs/subread.yaml" 
	threads: 4
	shell:
		"featureCounts -T {threads} -t exon -g gene_id -a {input.gtf} -o {output} {input.bams}" 
		" 2> {log}"





