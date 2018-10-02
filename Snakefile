configfile: "cfg/config.yaml"

def germline_samples():
  samples = set(config['samples'])
  tumours = set(config['tumours'])
  return list(samples.difference(tumours))

# final output
rule all:
  input:
    "out/msi.summary.tsv",
    "out/msi.summary.rotated.transcribed.tsv",
    "out/msi.summary.mini.tsv",
    "out/msi.measure.summary",
    "out/msi.cluster.tumours",
    "out/msi.primaries.cluster",
    "out/msi.blood.cluster",
    "out/msi.tumours.png",
    "out/msi.samples.png",
    expand("out/msi.cluster.{range}.tsv", range=config['indel_length_cluster']),
    expand("out/msi.cluster.{range}.png", range=config['indel_length_cluster']),
    expand("out/msi.{caller}.common", caller=config['msi_callers']),
    expand("out/mmr.{caller}.common", caller=config['callers']),
    "out/mmr.cluster.cosmic.png",
    "out/mmr.cluster.mmr.png",
    "out/mmr.cluster.cancertype.png",
    "out/regions.msi.stats",
    "out/msi.repeat_context.all.heatmap.png",
    "out/msi.indel_context.all.heatmap.png",
    "out/msi.repeat_indel_context.all.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.heatmap.png",
    "out/msi.repeat_context_emast.all.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.raw.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.threshold_0.3.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.threshold_0.3.transcribed.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.threshold_0.3.log.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.representative.transcribed.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.representative.transcribed.sample_norm.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.representative.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.representative.complemented.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.sample_norm.heatmap.png",
    "out/msi.repeat_context.all.rotated.transcribed.heatmap.png",
    "out/msi.repeat_context.all.rotated.transcribed.sample_norm.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.normalized.heatmap.png",
    "out/msi.repeat_indel_context_rotated.all.representative.normalized.heatmap.png",
    "out/msi.unique.tsv",
    "out/ddr.summary",
    "out/mutational_signatures.tsv", # somatic mutational signatures
    "out/mutational_signatures_germline.tsv", # germline mutational signatures
    "out/aggregated/msi.affected_proportion.png"

##### msi #####
include: 'Snakefile.msi'
include: 'Snakefile.mmr'
include: 'Snakefile.differential'

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

#rule exon_mutation_rate_somatic:
#  input:
#    vcf="in/{tumour}.{caller}.vcf", sample=config['tumours'], caller=config['callers']),
#  output:
#    "out/{tumour}.{caller}.somatic.vcf"
#  log:
#    stderr="log/exon_mutation_rate_somatic.{sample}.{caller}.err"
#  shell:
#    "src/subtract_germline.sh {input.vcf} >{output} 2>{log.stderr}"

#rule exon_mutation_rate:
#  input:
#    vcfs=expand("out/{sample}.{caller}.somatic.vcf", sample=config['samples'], caller=config['callers']),
#    bed="out/hg19.exons"
#  output:
#    "out/mutation_rate.exons.summary"
#  params:
#    min_dp=config['min_dp']
#  shell:
#    "src/mutation_rate.py --vcfs {input} --bed {input.bed} --min_dp {params.min_dp} 1>{output} 2>{log.stderr}"
#
#

# mutational signatures
rule mutational_signature:
  input:
    reference=config["genome"],
    vcf=expand("in/{{tumour}}.{caller}.vcf", caller=config["snv_caller"])
  output:
    "out/{tumour}.mutational_signature.tumour.exposures"
  log:
    stderr="log/{tumour}.mutational_signature.stderr"
  params:
    af=config["snv_af"],
    dp=config["snv_dp"]
  shell:
    "(module load htslib-intel/1.5 && "
    "src/filter_af.py --sample TUMOR --af {params.af} --dp {params.dp} < {input.vcf} | bgzip > tmp/{wildcards.tumour}.filter.vcf.gz && "
    "python software/mutational_signature-0.1/count.py --genome {input.reference} --vcf tmp/{wildcards.tumour}.filter.vcf.gz > out/{wildcards.tumour}.mutational_signature.counts && "
    "python software/mutational_signature-0.1/plot.py out/{wildcards.tumour}.mutational_signature.png < out/{wildcards.tumour}.mutational_signature.counts && "
    "python software/mutational_signature-0.1/decompose.py --signatures software/mutational_signature-0.1/signatures30.txt --counts out/{wildcards.tumour}.mutational_signature.counts > {output}) 2>{log.stderr}"

rule combine_mutational_signatures:
  input:
    expand("out/{tumour}.mutational_signature.tumour.exposures", tumour=config['tumours']),
  output:
    "out/mutational_signatures.tsv"
  shell:
    "src/combine_tsv.py {input} | sed 's/^out\\/\\(.*\)\\.mutational_signature\\.tumour\\.exposures/\\1/' > {output}"

# mutational signatures - germline (experiment)
rule mutational_signature_germline:
  input:
    reference=config["genome"],
    vcf=expand("in/{{germline}}.{caller}.vcf", caller=config["snv_caller"])
  output:
    "out/{germline}.mutational_signature.germline.exposures"
  log:
    stderr="log/{germline}.mutational_signature.stderr"
  params:
    af=config["snv_af"],
    dp=config["snv_dp"]
  shell:
    "(module load htslib-intel/1.5 && "
    "src/filter_af.py --sample {wildcards.germline} --af {params.af} --dp {params.dp} < {input.vcf} | bgzip > tmp/{wildcards.germline}.filter.vcf.gz && "
    "python software/mutational_signature-0.1/count.py --genome {input.reference} --vcf tmp/{wildcards.germline}.filter.vcf.gz > out/{wildcards.germline}.mutational_signature.counts && "
    "python software/mutational_signature-0.1/plot.py out/{wildcards.germline}.mutational_signature.png < out/{wildcards.germline}.mutational_signature.counts && "
    "python software/mutational_signature-0.1/decompose.py --signatures software/mutational_signature-0.1/signatures30.txt --counts out/{wildcards.germline}.mutational_signature.counts > {output}) 2>{log.stderr}"

rule combine_mutational_signatures_germline:
  input:
    expand("out/{germline}.mutational_signature.germline.exposures", germline=germline_samples()),
  output:
    "out/mutational_signatures_germline.tsv"
  shell:
    "src/combine_tsv.py {input} | sed 's/^out\\/\\(.*\)\\.mutational_signature\\.germline\\.exposures/\\1/' > {output}"

### general setup

# builds a genome length file
rule make_genome_lengths:
  input:
    config["genome"]
  output:
    config["genome_lengths"]
  shell:
    "src/make_lengths.py < {input} > {output}"
