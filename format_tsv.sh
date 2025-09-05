cat $1 | \
awk -v OFS='\t' '{split($1, a, ","); print a[2],a[1],$0}' | cut -f 1,2,4- | \
sed 's/hprc-v1.1-mc-sampled32d/HPRC-v1-Sampled32d/g' | \
sed 's/hprc-v2.0-mc-eval-d46/HPRC-d46/g' | \
sed 's/hprc-v2.0-mc-eval.ec1M-sampled16o/HPRC-sampled16/g' | \
sed 's/hprc-v2.0-minigraph-eval/HPRC-minigraph/g' | \
sed 's/primary/CHM13/g' | \
sed 's/giraffe/Giraffe/g' | \
sed 's/bwa/BWA-MEM/g' | \
sed 's/mini/Mini/g' | \
sed 's/graphaligner/GraphAligner/g' | \
sed 's/pbmm2/PBMM2/g' | \
sed 's/winnowmap/Winnowmap/g' | \
sed 's/HPRC-d46\tBWA-MEM/CHM13\tBWA-MEM/g' | \
sed 's/HPRC-d46\tMinimap2/CHM13\tMinimap2/g' | \
sed 's/HPRC-d46\tWinnowmap/CHM13\tWinnowmap/g' | \
sed 's/HPRC-d46\tPBMM2/CHM13\tPBMM2/g' | \
sed 's/\t/\t\& /g'  | sed 's/$/ \\\\/g' | \
sort
