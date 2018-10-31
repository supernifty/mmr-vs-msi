configfile: "cfg/config.yaml"

def germline_samples():
  samples = set(config['samples'])
  tumours = set(config['tumours'])
  return list(samples.difference(tumours))

# final output
rule all:
  input:
    "out/aggregate/msi.summary.tsv",
    "out/aggregate/msi.summary.rotated.transcribed.tsv",
    "out/aggregate/msi.summary.mini.tsv",
    "out/aggregate/msi.measure.summary",
    "out/aggregate/msi.cluster.tumours",
    "out/aggregate/msi.primaries.cluster",
    "out/aggregate/msi.blood.cluster",
    "out/aggregate/msi.tumours.png",
    "out/aggregate/msi.samples.png",
    expand("out/aggregate/msi.cluster.{range}.tsv", range=config['indel_length_cluster']),
    expand("out/aggregate/msi.cluster.{range}.png", range=config['indel_length_cluster']),
    expand("out/aggregate/msi.{caller}.common", caller=config['msi_callers']),
#    expand("out/aggregate/mmr.{caller}.common", caller=config['callers']),
#    "out/aggregate/mmr.cluster.cosmic.png",
#    "out/aggregate/mmr.cluster.mmr.png",
#    "out/aggregate/mmr.cluster.cancertype.png",
    "out/ref/regions.msi.stats",
    "out/aggregate/msi.repeat_context.all.heatmap.png",
    "out/aggregate/msi.indel_context.all.heatmap.png",
    "out/aggregate/msi.repeat_indel_context.all.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.heatmap.png",
    "out/aggregate/msi.repeat_context_emast.all.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.raw.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.threshold_0.3.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.threshold_0.3.transcribed.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.threshold_0.3.log.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.representative.transcribed.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.representative.transcribed.sample_norm.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.representative.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.representative.complemented.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.sample_norm.heatmap.png",
    "out/aggregate/msi.repeat_context.all.rotated.transcribed.heatmap.png",
    "out/aggregate/msi.repeat_context.all.rotated.transcribed.sample_norm.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.normalized.heatmap.png",
    "out/aggregate/msi.repeat_indel_context_rotated.all.representative.normalized.heatmap.png",
    "out/aggregate/msi.unique.tsv",
    "out/aggregate/ddr.summary.tsv",
    "out/aggregate/mutational_signatures.tsv", # somatic mutational signatures
    "out/aggregate/mutational_signatures_germline.tsv", # germline mutational signatures
    "out/aggregate/msi.affected_proportion.png", # comparing groups
    "out/aggregate/mutations.summary.tsv",
    "out/aggregate/mutational_signatures.indel.counts.tsv"

##### msi #####
include: 'Snakefile.msi'
include: 'Snakefile.mmr'
include: 'Snakefile.differential'

### overall stats

rule count_mutations:
  input:
    expand("in/{sample}.{caller}.vcf", sample=config['samples'], caller=config['callers'])
  output:
    "out/aggregate/mutations.summary.tsv"
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
#    bed="out/refseq.exons"
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
    "({config[module_htslib]} && "
    "src/filter_af.py --sample TUMOR --af {params.af} --dp {params.dp} < {input.vcf} | bgzip > tmp/{wildcards.tumour}.filter.vcf.gz && "
    "python software/mutational_signature-0.1/count.py --genome {input.reference} --vcf tmp/{wildcards.tumour}.filter.vcf.gz > out/{wildcards.tumour}.mutational_signature.counts && "
    "python software/mutational_signature-0.1/plot.py out/{wildcards.tumour}.mutational_signature.png < out/{wildcards.tumour}.mutational_signature.counts && "
    "python software/mutational_signature-0.1/decompose.py --signatures software/mutational_signature-0.1/signatures30.txt --counts out/{wildcards.tumour}.mutational_signature.counts > {output}) 2>{log.stderr}"

rule mutational_signature_indels:
  input:
    regions="out/ref/regions.msi.final.bed",
    vcf=expand("in/{{tumour}}.{caller}.vcf", caller=config["snv_caller"])
  output:
    "out/{tumour}.mutational_signature.tumour.indel.counts"
  log:
    stderr="log/{tumour}.mutational_signature_indel.stderr"
  params:
    af=config["snv_af"],
    dp=config["snv_dp"]
  shell:
    "({config[module_htslib]} && "
    "src/filter_af.py --sample TUMOR --af {params.af} --dp {params.dp} < {input.vcf} | bgzip > tmp/{wildcards.tumour}.filter.vcf.gz && "
    "software/mutational_signature-0.1/count_indels.py --annotation {input.regions} --rotate --complement < tmp/{wildcards.tumour}.filter.vcf.gz > {output}) 2>{log.stderr}"

rule combine_mutational_indels:
  input:
    expand("out/{tumour}.mutational_signature.tumour.indel.counts", tumour=config['tumours']),
  output:
    "out/aggregate/mutational_signatures.indel.counts.tsv"
  shell:
    "src/combine_tsv.py {input} | sed 's/^out\\/\\(.*\)\\.mutational_signature\\.tumour\\.indel.counts/\\1/' > {output}"


rule combine_mutational_signatures:
  input:
    expand("out/{tumour}.mutational_signature.tumour.exposures", tumour=config['tumours']),
  output:
    "out/aggregate/mutational_signatures.tsv"
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
    "({config[module_htslib]} && "
    "src/filter_af.py --sample NORMAL --af {params.af} --dp {params.dp} < {input.vcf} | bgzip > tmp/{wildcards.germline}.filter.vcf.gz && "
    "python software/mutational_signature-0.1/count.py --genome {input.reference} --vcf tmp/{wildcards.germline}.filter.vcf.gz > out/{wildcards.germline}.mutational_signature.counts && "
    "python software/mutational_signature-0.1/plot.py out/{wildcards.germline}.mutational_signature.png < out/{wildcards.germline}.mutational_signature.counts && "
    "python software/mutational_signature-0.1/decompose.py --signatures software/mutational_signature-0.1/signatures30.txt --counts out/{wildcards.germline}.mutational_signature.counts > {output}) 2>{log.stderr}"

rule combine_mutational_signatures_germline:
  input:
    expand("out/{germline}.mutational_signature.germline.exposures", germline=germline_samples()),
  output:
    "out/aggregate/mutational_signatures_germline.tsv"
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

rule make_msi_candidates:
  input:
    config["genome"]
  output:
    tb1="tmp/tb1.out",
    tb2="tmp/tb2.out",
    tb3="tmp/tb3.out",
    tb4="tmp/tb4.out",
    tb5="tmp/tb5.out" 
  shell:
    "software/tandem_repeats --repeat:1 --min:6 < {config[genome]} >{output.tb1} && "
    "software/tandem_repeats --repeat:2 --min:6 < {config[genome]} >{output.tb2} && "
    "software/tandem_repeats --repeat:3 --min:6 < {config[genome]} >{output.tb3} && "
    "software/tandem_repeats --repeat:4 --min:6 < {config[genome]} >{output.tb4} && "
    "software/tandem_repeats --repeat:5 --min:6 < {config[genome]} >{output.tb5}"

rule annotate_msi_candidates_bed:
  input:
    tb1="tmp/tb1.out",
    tb2="tmp/tb2.out",
    tb3="tmp/tb3.out",
    tb4="tmp/tb4.out",
    tb5="tmp/tb5.out",
    exons="out/refseq.exons" 
  output:
    config["msi_loci"]
  shell:
    "src/annotate.py --name exon {input.tb1} {input.tb2} {input.tb3} {input.tb4} {input.tb5} < {input.exons} > tmp/msi.exons && "
    "src/annotate.py --name bethesda tmp/msi.exons < reference/bethesda.bed | sort -k1,1 -k2,2n > reference/msi.candidates.bed "

