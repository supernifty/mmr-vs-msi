configfile: "cfg/config.yaml"

rule all:
  input:
    "out/mmr.summary"

rule mmr_stats:
  input:
    expand("out/{sample}.mmr.varscan.filter.vcf", sample=config['samples'])
  output:
    "out/mmr.summary"
  log:
    stderr="log/mmr.summary.stderr",
  shell:
    "src/mmr_stats.py {input} 1>{output} 2>{log.stderr}"

rule filter:
  input:
    "out/{sample}.mmr.varscan.vep.vcf"
  output:
    "out/{sample}.mmr.varscan.filter.vcf"
  log:
    stderr="log/filter.{sample}.stderr",
  shell:
    "src/filter_vcf.py < {input} 1>{output} 2>{log.stderr}"

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

rule make_genome_lengths:
  input:
    config["genome"]
  output:
    config["genome_lengths"]
  shell:
    "src/make_lengths.py < {input} > {output}"
