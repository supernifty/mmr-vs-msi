configfile: "cfg/config.yaml"

# final output
rule all:
  input:
    "out/msi.summary"

##### msi #####

# overall likelihood of each sample to be affected by msi
rule msi_assess:
  input:
    "out/mmr.summary",
    expand("out/{sample}.msi.{caller}.annotated.vcf", sample=config['samples'], caller=config['msi_callers'])
  output:
    "out/msi.summary"
  log:
    stderr="log/msi.summary.stderr",
  shell:
    "src/msi_stats.py {input} 1>{output} 2>{log.stderr}"

# msi output is list of samples and the degree of MSI found
rule annotate_msi:
  input:
    "out/regions.msi.final.bed",
    "in/{sample}.{caller}.vcf" # each unfiltered msi_caller
  output:
    "out/{sample}.msi.{caller}.annotated.vcf" # makes a filtered vcf
  log:
    stderr="log/msi_filter.{sample}.{caller}.stderr",
  shell:
    "src/annotate_vcf.py {input[0]} < {input[1]} 1>{output} 2>{log.stderr}"

# add oncogene annotations to msi list
rule annotate_msi_oncogene_regions:
  input:
    "out/regions.oncogenes.bed", # list of oncogene regions 
    "out/regions.msi.bed", # list of msi regions
  output:
    "out/regions.msi.final.bed" # annotated regions
  log:
    stderr="log/annotate_bed.stderr",
  shell:
    "src/annotate_bed.py oncogene {input[0]} < {input[1]} > {output} 2> {log.stderr}"

# add oncogene annotations to msi list
rule make_oncogene_regions:
  input:
    "out/oncogenes.msi", # list of genes 
    config["refseq"], # list of refseq regions
  output:
    "out/regions.oncogenes.bed" # oncogene regions
  shell:
    "module load bedtools-intel/2.27.1; "
    "src/make_bed.py {input[0]} < {input[1]} | sort -k1,1 -k2,2n > {output}" 

# make a list of genes relevant to msi for the cancer being studied
rule make_msi_oncogenes:
  input:
    config['cosmic']
  output:
    "out/oncogenes.msi"
  shell:
    "cut -f1,8,9,10,11 < {input} | grep -i {config[cancer_type]} | cut -f1 | sort > {output}"

# limits msi regions by length and repeat type
rule make_msi_regions:
  input:
    "reference/msi.candidates.bed"
  output:
    "out/regions.msi.bed"
  log:
    stderr="log/make_msi_regions.stderr",
  shell:
    "src/filter_msi.py --minlen 10 --maxlen 100000 --minrepeat 1 --maxrepeat 4 < {input} 1>{output} 2>{log.stderr}"


##### mmr #####

# mmr analysis output is a list of candidates likely to be dMMR
rule mmr_stats:
  input:
    expand("out/{sample}.mmr.{caller}.filter.vcf", sample=config['samples'], caller=config['callers'])
  output:
    "out/mmr.summary"
  log:
    stderr="log/mmr.summary.stderr",
  shell:
    "src/mmr_stats.py --all --verbose {input} 1>{output} 2>{log.stderr}"

# applies filtering to the annotated vcf files
rule mmr_filter:
  input:
    "out/{sample}.mmr.{caller}.vep.vcf",
  output:
    "out/{sample}.mmr.{caller}.filter.vcf"
  log:
    stderr="log/filter.{sample}.{caller}.stderr",
  shell:
    "src/filter_vcf.py {wildcards.caller} < {input} 1>{output} 2>{log.stderr}"

# annotates vcf files with vep
rule mmr_annotate:
  input:
    "out/regions.mmr.bed",
    "in/{sample}.{caller}.vcf"
  output:
    "out/{sample}.mmr.{caller}.vep.vcf",
  log:
    stderr="log/annotate.{sample}.{caller}.stderr",
    stdout="log/annotate.{sample}.{caller}.stdout"
  shell:
    "src/annotate.sh {config[genome]} {input} {output} 1>{log.stdout} 2>{log.stderr}"

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
