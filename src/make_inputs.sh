
for sample in $(cat cfg/samples); do
  ln -s /scratch/VR0211/pan-prostate/out/${sample}.varscan.vcf in/
done
