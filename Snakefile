configfile: "cfg/config.yaml"

# final output
rule all:
  input:
    "out/mmr.summary"

# generate mmr candidates
rule mmr_stats:
  input:
    expand("out/{sample}.mmr.varscan.filter.vcf", sample=config['samples'])
  output:
    "out/mmr.summary"
  log:
    stderr="log/mmr.summary.stderr",
  shell:
    "src/mmr_stats.py {input} 1>{output} 2>{log.stderr}"

# applies filtering to the annotated vcf files
rule filter:
  input:
    "out/{sample}.mmr.varscan.vep.vcf"
  output:
    "out/{sample}.mmr.varscan.filter.vcf"
  log:
    stderr="log/filter.{sample}.stderr",
  shell:
    "src/filter_vcf.py < {input} 1>{output} 2>{log.stderr}"

# annotates vcf files with vep
rule annotate:
  input:
    "in/{sample}.varscan.vcf"
  output:
    "out/{sample}.mmr.varscan.vep.vcf"
  log:
    stderr="log/annotate.{sample}.stderr",
    stdout="log/annotate.{sample}.stdout"
  shell:
    "src/annotate.sh {wildcards.sample} {config[genome]} out/regions.bed 1>{log.stdout} 2>{log.stderr}"

# generates bed regions to filter vcfs on based on the provided genes, and refseq exons
rule make_regions:
  input: 
    config["genes"],
    config["refseq"],
    config["genome_lengths"]
  output:
    "out/regions.bed"
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
