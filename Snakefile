configfile: "cfg/config.yaml"

# final output
rule all:
  input:
    "out/msi.summary",
    "out/msi.measure.summary",
    "out/msi.common",
    "out/mmr.common",

##### msi #####
include: 'Snakefile.msi'
include: 'Snakefile.mmr'

### overall stats

rule count_mutations:
  input:
    expand("in/{sample}.{caller}.vcf", sample=config['samples'], caller=config['callers'])
  output:
    "out/mutations.summary"
  log:
    stderr="log/count_mutations.stderr",
  shell:
    "src/count_mutations.py {input} 1>{output} 2>{log.stderr}"

### general setup

# generates bed regions to filter vcfs on based on the provided genes, and refseq exons
rule make_mmr_regions:
  input: 
    config["genes_mmr"],
    config["refseq"],
    config["genome_lengths"]
  output:
    "out/regions.mmr.bed"
  shell:
    "module load bedtools-intel/2.27.1; "
    "src/make_bed.py {input[0]} < {input[1]} | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct | bedtools slop -b {config[slop]} -i - -g {input[2]} > {output}"

# builds a genome length file
rule make_genome_lengths:
  input:
    config["genome"]
  output:
    config["genome_lengths"]
  shell:
    "src/make_lengths.py < {input} > {output}"
