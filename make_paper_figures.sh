set -ex

ROOT_DIR=/private/groups/patenlab/project-lrg
OUT_DIR=./plots
mkdir -p "${OUT_DIR}"


# This can be used to get all the input files
# I left the --dry-run in because it keeps trying to re-run anything
# Asking for all_paper_figures also works
#(umask 0002 ; snakemake -p --dry-run --keep-going --rerun-triggers mtime --resources full_cluster_nodes=2 threads=256 --slurm --rerun-incomplete --latency-wait 120 -j128 --use-singularity --singularity-args '-B /private/home/jmonlong/workspace/lreval/data/:/private/home/jmonlong/workspace/lreval/data/ -B /private/groups/patenlab/project-lrg/:/private/groups/patenlab/project-lrg/ -B /private/groups/patenlab/anovak/projects/hprc/lr-giraffe/:/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/' \
#/private/groups/patenlab/project-lrg/experiments/hifi_sim_1m_headline/results/compared.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10y2025_sim_1m_headline/results/compared.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_sim_1m_headline/results/mapping_stats_sim.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10y2025_sim_1m_headline/results/mapping_stats_sim.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_real_full_headline/results/clipped_or_unmapped_percent.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10_real_full_headline/results/clipped_or_unmapped_percent.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_real_full_headline/results/run_and_sampling_time_from_benchmark.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10_real_full_headline/results/run_and_sampling_time_from_benchmark.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_real_full_headline/results/index_load_time.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_real_full_headline/results/memory_from_benchmark.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10_real_full_headline/results/memory_from_benchmark.tsv)

######################################### Define colors

PURPLE=#9B8BF4
LIGHT_PURPLE=#C49AB6
BLUE=#53AEDF
LIGHT_BLUE=#B3C7F7
TEAL=#6396B0
GREEN=#27A974
YELLOW=#F2C300
DARK_ORANGE=#C44601
ORANGE=#F57600
PINK=#D066A3
LIGHT_GREEN=#9CDA4D


BWA_COLOR=${LIGHT_GREEN}
MINIMAP_COLOR=${PURPLE}
WINNOWMAP_COLOR=${LIGHT_PURPLE}
GIRAFFE_PRIMARY_COLOR=${BLUE}
PBMM_COLOR=${LIGHT_BLUE}
MINIGRAPH_MINIGRAPH_COLOR=${TEAL}
GRAPHALIGNER_MINIGRAPH_COLOR=${GREEN}
GRAPHALIGNER_DEFAULT_COLOR=${YELLOW}
GRAPHALIGNER_FAST_COLOR=${PINK}
GIRAFFE_COLOR=${DARK_ORANGE}
GIRAFFE_SAMPLED_COLOR=${ORANGE}

######################################### Define conditions
MINIMAP_HIFI=minimap2-map-hifi,hprc-v2.0-mc-eval-d46
MINIMAP_R10=minimap2-lr:hqae,hprc-v2.0-mc-eval-d46
MINIMAP_SR=minimap2-sr,hprc-v2.0-mc-eval-d46
WINNOWMAP=winnowmap,hprc-v2.0-mc-eval-d46
GIRAFFE_PRIMARY_HIFI=giraffe,primary
GIRAFFE_PRIMARY_R10=giraffe,primary
GIRAFFE_PRIMARY_SR=giraffe,primary
PBMM=pbmm2,hprc-v2.0-mc-eval-d46
BWA=bwa,hprc-v2.0-mc-eval-d46

MINIGRAPH_MINIGRAPH=minigraph,hprc-v2.0-minigraph-eval
GRAPHALIGNER_MINIGRAPH=graphaligner-default,hprc-v2.0-minigraph-eval

GRAPHALIGNER_DEFAULT=graphaligner-default,hprc-v2.0-mc-eval-d46
GRAPHALIGNER_FAST=graphaligner-fast,hprc-v2.0-mc-eval-d46
GIRAFFE_HIFI=giraffe,hprc-v2.0-mc-eval-d46
GIRAFFE_R10=giraffe,hprc-v2.0-mc-eval-d46
GIRAFFE_SR=giraffe,hprc-v2.0-mc-eval-d46

GIRAFFE_SAMPLED_HIFI=giraffe,hprc-v2.0-mc-eval.ec1M-sampled16o
GIRAFFE_SAMPLED_R10=giraffe,hprc-v2.0-mc-eval.ec1M-sampled16o
GIRAFFE_SAMPLED_SR=giraffe,hprc-v2.0-mc-eval.ec1M-sampled16o

###################################### Define conditions and colors for each experiment

HIFI_SIM_1M_HEADLINE_CATEGORIES="${MINIMAP_HIFI} ${WINNOWMAP} ${PBMM} ${GIRAFFE_PRIMARY_HIFI} ${MINIGRAPH_MINIGRAPH} ${GRAPHALIGNER_MINIGRAPH} ${GRAPHALIGNER_DEFAULT} ${GIRAFFE_HIFI} ${GIRAFFE_SAMPLED_HIFI}"
HIFI_SIM_1M_HEADLINE_COLORS="${MINIMAP_COLOR} ${WINNOWMAP_COLOR} ${PBMM_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${MINIGRAPH_MINIGRAPH_COLOR} ${GRAPHALIGNER_MINIGRAPH_COLOR} ${GRAPHALIGNER_DEFAULT_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"

R10_SIM_1M_HEADLINE_CATEGORIES="${MINIMAP_R10} ${WINNOWMAP} ${GIRAFFE_PRIMARY_R10} ${MINIGRAPH_MINIGRAPH} ${GRAPHALIGNER_MINIGRAPH} ${GRAPHALIGNER_DEFAULT} ${GIRAFFE_R10} ${GIRAFFE_SAMPLED_R10}"
R10_SIM_1M_HEADLINE_COLORS="${MINIMAP_COLOR} ${WINNOWMAP_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${MINIGRAPH_MINIGRAPH_COLOR} ${GRAPHALIGNER_MINIGRAPH_COLOR} ${GRAPHALIGNER_DEFAULT_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"

SR_SIM_1M_CATEGORIES="${MINIMAP_SR} ${BWA} ${GIRAFFE_PRIMARY_SR} ${GIRAFFE_SR} ${GIRAFFE_SAMPLED_SR}"
SR_SIM_1M_COLORS="${MINIMAP_COLOR} ${BWA_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"

# graphaligner-default on the M/C graph and graphaligner-default on minigraph timed out after 7 days
HIFI_REAL_FULL_HEADLINE_CATEGORIES="${MINIMAP_HIFI} ${WINNOWMAP} ${PBMM} ${GIRAFFE_PRIMARY_HIFI} ${MINIGRAPH_MINIGRAPH} ${GRAPHALIGNER_FAST} ${GIRAFFE_HIFI} ${GIRAFFE_SAMPLED_HIFI}"
HIFI_REAL_FULL_HEADLINE_COLORS="${MINIMAP_COLOR} ${WINNOWMAP_COLOR} ${PBMM_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${MINIGRAPH_MINIGRAPH_COLOR} ${GRAPHALIGNER_FAST_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"

# graphalinger-default on the M/C graph timed out after 7 days
R10_REAL_FULL_HEADLINE_CATEGORIES="${MINIMAP_R10} ${WINNOWMAP} ${GIRAFFE_PRIMARY_R10} ${MINIGRAPH_MINIGRAPH} ${GRAPHALIGNER_MINIGRAPH} ${GRAPHALIGNER_FAST} ${GIRAFFE_R10} ${GIRAFFE_SAMPLED_R10}"
R10_REAL_FULL_HEADLINE_COLORS="${MINIMAP_COLOR} ${WINNOWMAP_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${MINIGRAPH_MINIGRAPH_COLOR} ${GRAPHALIGNER_MINIGRAPH_COLOR} ${GRAPHALIGNER_FAST_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"

SR_REAL_FULL_HEADLINE_CATEGORIES="${MINIMAP_SR} ${BWA} ${GIRAFFE_PRIMARY_SR} ${GIRAFFE_SR} ${GIRAFFE_SAMPLED_SR}"
SR_REAL_FULL_HEADLINE_COLORS="${MINIMAP_COLOR} ${BWA_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"

########################################### ROCs and QQs

# Zoomed in headline roc plots
HIFI_ACCURACY=${ROOT_DIR}/experiments/hifi_sim_1m/results/compared.tsv
R10_ACCURACY=${ROOT_DIR}/experiments/r10y2025_sim_1m/results/compared.tsv
ILLUMINA_ACCURACY=${ROOT_DIR}/experiments/illumina_sim_1m/results/compared.tsv
ELEMENT_ACCURACY=${ROOT_DIR}/experiments/element_sim_1m/results/compared.tsv

Rscript plot-roc.R ${HIFI_ACCURACY} ${OUT_DIR}/roc_hifi_headline_zoomed.pdf $(echo $HIFI_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "HIFI ROC" $(echo $HIFI_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' ) '0.0001,0.03,0.9545,0.995'

Rscript plot-roc.R ${R10_ACCURACY} ${OUT_DIR}/roc_r10_headline_zoomed.pdf $(echo $R10_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "R10 ROC" $(echo $R10_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' ) '0.001,0.030,0.965,0.994'

# Zoomed out headline roc plots
Rscript plot-roc.R ${HIFI_ACCURACY} ${OUT_DIR}/roc_hifi_headline.pdf $(echo $HIFI_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "HIFI ROC" $(echo $HIFI_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' )

Rscript plot-roc.R ${R10_ACCURACY} ${OUT_DIR}/roc_r10_headline.pdf $(echo $R10_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "R10 ROC" $(echo $R10_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' )

Rscript plot-roc.R ${ILLUMINA_ACCURACY} ${OUT_DIR}/roc_illumina_headline.pdf $(echo $SR_SIM_1M_CATEGORIES | sed 's/ /;/g' ) "Illumina ROC" $(echo $SR_SIM_1M_COLORS | sed 's/ /,/g' ) '0.0001,0.03,0.918,0.985'

Rscript plot-roc.R ${ELEMENT_ACCURACY} ${OUT_DIR}/roc_element_headline.pdf $(echo $SR_SIM_1M_CATEGORIES | sed 's/ /;/g' ) "Element ROC" $(echo $SR_SIM_1M_COLORS | sed 's/ /,/g' )

# headline qq plots
Rscript plot-qq.R ${HIFI_ACCURACY} ${OUT_DIR}/qq_hifi_headline.pdf $(echo $HIFI_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "HIFI QQ" $(echo $HIFI_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' )

Rscript plot-qq.R ${R10_ACCURACY} ${OUT_DIR}/qq_r10_headline.pdf $(echo $R10_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "R10 QQ" $(echo $R10_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' )

########################################## Incorrect counts

HIFI_INCORRECT=${ROOT_DIR}/experiments/hifi_sim_1m/results/mapping_stats_sim.tsv
R10_INCORRECT=${ROOT_DIR}/experiments/r10y2025_sim_1m/results/mapping_stats_sim.tsv

python3 barchart.py <(tail -n +2 ${HIFI_INCORRECT} | cut -f1,5) --divisions <(tail -n +2 ${HIFI_INCORRECT} | cut -f1,4)  --title 'HiFi incorrect read count' --y_label 'Count' --x_label 'Condition' --x_sideways --no_n --categories ${HIFI_SIM_1M_HEADLINE_CATEGORIES}  --colors ${HIFI_SIM_1M_HEADLINE_COLORS} --save ${OUT_DIR}/incorrectness_hifi.pdf

python3 barchart.py <(tail -n +2 ${R10_INCORRECT} | cut -f1,5) --divisions <(tail -n +2 ${R10_INCORRECT} | cut -f1,4)  --title 'R10 incorrect read count' --y_label 'Count' --x_label 'Condition' --x_sideways --no_n --categories ${R10_SIM_1M_HEADLINE_CATEGORIES}  --colors ${R10_SIM_1M_HEADLINE_COLORS}  --save ${OUT_DIR}/incorrectness_r10.pdf

########################################## Clips

HIFI_REAL=${ROOT_DIR}/experiments/hifi_real_full/results/mapping_stats_real.tsv
R10_REAL=${ROOT_DIR}/experiments/r10y2025_real_full/results/mapping_stats_real.tsv
ILLUMINA_REAL=${ROOT_DIR}/experiments/illumina_real_full/results/mapping_stats_real.tsv
ELEMENT_REAL=${ROOT_DIR}/experiments/element_real_full/results/mapping_stats_real.tsv

HIFI_UNMAPPED=${OUT_DIR}/hifi_unmapped.tsv 
R10_UNMAPPED=${OUT_DIR}/r10_unmapped.tsv
HIFI_CLIPS=${OUT_DIR}/hifi_clips.tsv 
R10_CLIPS=${OUT_DIR}/r10_clips.tsv

grep -E $(echo ${HIFI_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${HIFI_REAL} | awk -v OFS='\t' '{{print $1,100*$9/$10}}' > ${HIFI_CLIPS}
grep -E $(echo ${R10_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${R10_REAL} | awk -v OFS='\t' '{{print $1,100*$9/$10}}' > ${R10_CLIPS}


grep -E $(echo ${HIFI_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${HIFI_REAL} | awk -v OFS='\t' '{{print $1,100*$8/$10}}' > ${HIFI_UNMAPPED}
grep -E $(echo ${R10_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${R10_REAL} | awk -v OFS='\t' '{{print $1,100*$8/$10}}' > ${R10_UNMAPPED}
cat $HIFI_CLIPS
cat $HIFI_UNMAPPED

#hifi softclips low
LIMIT=$(python3 get_outlier_limit.py ${HIFI_CLIPS} small)
python3 barchart.py ${HIFI_CLIPS} --divisions ${HIFI_UNMAPPED} --max ${LIMIT}  --title 'HiFi Clipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/softclips_hifi_headline_low.pdf

#hifi softclips high
LIMIT=$(python3 get_outlier_limit.py ${HIFI_CLIPS} big)
python3 barchart.py ${HIFI_CLIPS} --divisions ${HIFI_UNMAPPED} --min ${LIMIT}  --title 'HiFi Clipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/softclips_hifi_headline_high.pdf

#r10 softclips
python3 barchart.py ${R10_CLIPS} --divisions ${R10_UNMAPPED} --title 'R10 Clipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories ${R10_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${R10_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/softclips_r10_headline_high.pdf

rm $HIFI_UNMAPPED
rm $R10_UNMAPPED
rm $HIFI_CLIPS
rm $R10_CLIPS

########################################### Runtime

HIFI_RUNTIME=${OUT_DIR}/hifi_runtimes.tsv
R10_RUNTIME=${OUT_DIR}/r10_runtimes.tsv
ILLUMINA_RUNTIME=${OUT_DIR}/illumina_runtimes.tsv
ELEMENT_RUNTIME=${OUT_DIR}/element_runtimes.tsv
grep -E $(echo ${HIFI_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${HIFI_REAL} > ${HIFI_RUNTIME}
grep -E $(echo ${R10_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${R10_REAL} > ${R10_RUNTIME}
grep -E $(echo ${SR_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${ILLUMINA_REAL} > ${ILLUMINA_RUNTIME}
grep -E $(echo ${SR_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${ELEMENT_REAL} > ${ELEMENT_RUNTIME}



# Higher hifi runtimes
LIMIT=$(python3 get_outlier_limit.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) big)
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,($3+$4)/60}}' ${HIFI_RUNTIME}) --min ${LIMIT}  --title 'HiFi Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_hifi_headline_slow.pdf

# Lower hifi runtimes
LIMIT=$(python3 get_outlier_limit.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) small)
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,($3+$4)/60}}' ${HIFI_RUNTIME}) --max ${LIMIT}  --title 'HiFi Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_hifi_headline_fast.pdf

# R10 runtimes as one plot
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${R10_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,($3+$4)/60}}' ${R10_RUNTIME}) --title 'R10 Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${R10_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${R10_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_r10_headline.pdf

# Higher r10 runtimes
LIMIT=$(python3 get_outlier_limit.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${R10_RUNTIME}) big)
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${R10_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,($3+$4)/60}}' ${R10_RUNTIME}) --min ${LIMIT}  --title 'R10 Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${R10_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${R10_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_r10_headline_slow.pdf

# Lower R10 runtimes
LIMIT=$(python3 get_outlier_limit.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${R10_RUNTIME}) small)
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${R10_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,($3+$4)/60}}' ${R10_RUNTIME}) --max ${LIMIT}  --title 'R10 Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${R10_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${R10_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_r10_headline_fast.pdf

# Illumina runtimes as one plot
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${ILLUMINA_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,($3+$4)/60}}' ${ILLUMINA_RUNTIME}) --title 'Illumina Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${SR_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${SR_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_illumina_headline.pdf

# Element runtimes as one plot
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${ELEMENT_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,($3+$4)/60}}' ${ELEMENT_RUNTIME}) --title 'Element Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${SR_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${SR_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_element_headline.pdf

rm ${HIFI_RUNTIME}
rm ${R10_RUNTIME}
rm ${ILLUMINA_RUNTIME}
rm ${ELEMENT_RUNTIME}
####################################### Memory
R10_MEMORY=${OUT_DIR}/r10_memory.tsv
ILLUMINA_MEMORY=${OUT_DIR}/illumina_memory.tsv
ELEMENT_MEMORY=${OUT_DIR}/element_memory.tsv

grep -E $(echo ${HIFI_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${HIFI_REAL} | awk -v OFS='\t' '{{print $1,$5}}' > ${HIFI_MEMORY}
grep -E $(echo ${R10_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${R10_REAL} | awk -v OFS='\t' '{{print $1,$5}}' > ${R10_MEMORY}
grep -E $(echo ${SR_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${ILLUMINA_REAL} | awk -v OFS='\t' '{{print $1,$5}}' > ${ILLUMINA_MEMORY}
grep -E $(echo ${SR_REAL_FULL_HEADLINE_CATEGORIES} | sed 's/ /|/g') ${ELEMENT_REAL} | awk -v OFS='\t' '{{print $1,$5}}' > ${ELEMENT_MEMORY}

python3 barchart.py ${HIFI_MEMORY} --title 'HiFi Memory From Benchmark' --y_label 'Memory (GB)' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/memory_hifi.pdf

python3 barchart.py ${R10_MEMORY} --title 'R10 Memory From Benchmark' --y_label 'Memory (GB)' --x_label 'Mapper' --x_sideways --no_n --categories ${R10_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${R10_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/memory_r10.pdf

python3 barchart.py ${ILLUMINA_MEMORY} --title 'Illumina Memory From Benchmark' --y_label 'Memory (GB)' --x_label 'Mapper' --x_sideways --no_n --categories ${SR_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${SR_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/memory_illumina.pdf

python3 barchart.py ${ELEMENT_MEMORY} --title 'Element Memory From Benchmark' --y_label 'Memory (GB)' --x_label 'Mapper' --x_sideways --no_n --categories ${SR_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${SR_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/memory_element.pdf

rm $HIFI_MEMORY
rm $R10_MEMORY
rm $ILLUMINA_MEMORY
rm $ELEMENT_MEMORY

####################################### SV Calling

GIRAFFE_VGCALL_HIFI=chm13,chm13,vgcall,hifi,giraffe-hifi-mmcs100,eval-d46
GIRAFFE_SAMPLED_VGCALL_HIFI=chm13,chm13,vgcall,hifi,giraffe-hifi-mmcs100,eval.ec1M-sampled16o
MINIMAP_SNIFFLES_HIFI=chm13,chm13,sniffles,hifi,minimap2-map-hifi,eval-d46

GIRAFFE_VGCALL_R10=chm13,chm13,vgcall,r10y2025,giraffe-r10-noflags,eval-d46
GIRAFFE_SAMPLED_VGCALL_R10=chm13,chm13,vgcall,r10y2025,giraffe-r10-noflags,eval.ec1M-sampled16o
MINIMAP_SNIFFLES_R10=chm13,chm13,sniffles,r10y2025,minimap2-lr:hqae,eval-d46

SV_CALLING_CONDITIONS_HIFI="${GIRAFFE_VGCALL_HIFI};${GIRAFFE_SAMPLED_VGCALL_HIFI};${MINIMAP_SNIFFLES_HIFI}"
SV_CALLING_CONDITIONS_R10="${GIRAFFE_VGCALL_R10};${GIRAFFE_SAMPLED_VGCALL_R10};${MINIMAP_SNIFFLES_R10}"

Rscript plot-calling-results.R ${ROOT_DIR}/experiments/sv_calling/results/sv_summary.tsv ${OUT_DIR}/sv_summary_hifi.pdf ${SV_CALLING_CONDITIONS_HIFI}
Rscript plot-calling-results.R ${ROOT_DIR}/experiments/sv_calling/results/sv_summary.tsv ${OUT_DIR}/sv_summary_r10.pdf ${SV_CALLING_CONDITIONS_R10}

############################################# DV Calling

GIRAFFE_DV_HIFI=hifi,giraffe-hifi-mmcs100,eval-d46,.model2025-03-26noinfo
GIRAFFE_SAMPLED_DV_HIFI=hifi,giraffe-hifi-mmcs100,eval.ec1M-sampled16o,.model2025-03-26noinfo
MINIMAP_DV_HIFI=hifi,minimap2-map-hifi,eval-d46,.nomodel

DV_CALLING_CONDITIONS_HIFI="${GIRAFFE_DV_HIFI};${GIRAFFE_SAMPLED_DV_HIFI};${MINIMAP_DV_HIFI}"

GIRAFFE_DV_R10=r10y2025,giraffe-r10-noflags,eval-d46,.nomodel
GIRAFFE_SAMPLED_DV_R10=r10y2025,giraffe-r10-noflags,eval.ec1M-sampled16o,.nomodel
MINIMAP_DV_R10=r10y2025,minimap2-lr:hqae,eval-d46,.nomodel

DV_CALLING_CONDITIONS_R10="${GIRAFFE_DV_R10};${GIRAFFE_SAMPLED_DV_R10};${MINIMAP_DV_R10}"


Rscript plot-calling-results.R ${ROOT_DIR}/experiments/dv_calling/results/dv_snp_summary.tsv ${OUT_DIR}/dv_snp_summary_hifi.pdf ${DV_CALLING_CONDITIONS_HIFI}
Rscript plot-calling-results.R ${ROOT_DIR}/experiments/dv_calling/results/dv_indel_summary.tsv ${OUT_DIR}/dv_indel_summary_hifi.pdf ${DV_CALLING_CONDITIONS_HIFI}

Rscript plot-calling-results.R ${ROOT_DIR}/experiments/dv_calling/results/dv_snp_summary.tsv ${OUT_DIR}/dv_snp_summary_r10.pdf ${DV_CALLING_CONDITIONS_R10}
#Rscript plot-calling-results.R ${ROOT_DIR}/experiments/dv_calling/results/dv_indel_summary.tsv ${OUT_DIR}/dv_indel_summary_r10.pdf ${DV_CALLING_CONDITIONS_R10}
