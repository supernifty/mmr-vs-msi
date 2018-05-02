configfile: "cfg/config.yaml"

# final output
rule all:
  input:
    "out/msi.summary",
    "out/msi.measure.summary",
    "out/msi.cluster",
    "out/msi.dendrogram.png",
    expand("out/msi.{caller}.common", caller=config['msi_callers']),
    expand("out/mmr.{caller}.common", caller=config['callers']),

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

# builds a genome length file
rule make_genome_lengths:
  input:
    config["genome"]
  output:
    config["genome_lengths"]
  shell:
    "src/make_lengths.py < {input} > {output}"
