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
#/private/groups/patenlab/project-lrg/experiments/hifi_real_full_headline/results/softclipped_or_unmapped_percent.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10_real_full_headline/results/softclipped_or_unmapped_percent.tsv \
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
WINNOWMAP=winnowmap,hprc-v2.0-mc-eval-d46
GIRAFFE_PRIMARY_HIFI=giraffe-noflags,primary
GIRAFFE_PRIMARY_R10=giraffe-mcspb0.00,primary
PBMM=pbmm2,hprc-v2.0-mc-eval-d46

MINIGRAPH_MINIGRAPH=minigraph,hprc-v2.0-minigraph-eval
GRAPHALIGNER_MINIGRAPH=graphaligner-default,hprc-v2.0-minigraph-eval

GRAPHALIGNER_DEFAULT=graphaligner-default,hprc-v2.0-mc-eval-d46
GRAPHALIGNER_FAST=graphaligner-fast,hprc-v2.0-mc-eval-d46
GIRAFFE_HIFI=giraffe-noflags,hprc-v2.0-mc-eval-d46
GIRAFFE_R10=giraffe-mcspb0.00,hprc-v2.0-mc-eval-d46

GIRAFFE_SAMPLED_HIFI=giraffe-noflags,hprc-v2.0-mc-eval.ec1M-sampled16o
GIRAFFE_SAMPLED_R10=giraffe-mcspb0.00,hprc-v2.0-mc-eval.ec1M-sampled16o

###################################### Define conditions and colors for each experiment

HIFI_SIM_1M_HEADLINE_CATEGORIES="${MINIMAP_HIFI} ${WINNOWMAP} ${PBMM} ${GIRAFFE_PRIMARY_HIFI} ${MINIGRAPH_MINIGRAPH} ${GRAPHALIGNER_MINIGRAPH} ${GRAPHALIGNER_DEFAULT} ${GIRAFFE_HIFI} ${GIRAFFE_SAMPLED_HIFI}"
HIFI_SIM_1M_HEADLINE_COLORS="${MINIMAP_COLOR} ${WINNOWMAP_COLOR} ${PBMM_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${MINIGRAPH_MINIGRAPH_COLOR} ${GRAPHALIGNER_MINIGRAPH_COLOR} ${GRAPHALIGNER_DEFAULT_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"

R10_SIM_1M_HEADLINE_CATEGORIES="${MINIMAP_R10} ${WINNOWMAP} ${GIRAFFE_PRIMARY_R10} ${MINIGRAPH_MINIGRAPH} ${GRAPHALIGNER_MINIGRAPH} ${GRAPHALIGNER_DEFAULT} ${GIRAFFE_R10} ${GIRAFFE_SAMPLED_R10}"
R10_SIM_1M_HEADLINE_COLORS="${MINIMAP_COLOR} ${WINNOWMAP_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${MINIGRAPH_MINIGRAPH_COLOR} ${GRAPHALIGNER_MINIGRAPH_COLOR} ${GRAPHALIGNER_DEFAULT_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"

# graphaligner-default on the M/C graph and graphaligner-default on minigraph timed out after 7 days
HIFI_REAL_FULL_HEADLINE_CATEGORIES="${MINIMAP_HIFI} ${WINNOWMAP} ${PBMM} ${GIRAFFE_PRIMARY_HIFI} ${MINIGRAPH_MINIGRAPH} ${GRAPHALIGNER_FAST} ${GIRAFFE_HIFI} ${GIRAFFE_SAMPLED_HIFI}"
HIFI_REAL_FULL_HEADLINE_COLORS="${MINIMAP_COLOR} ${WINNOWMAP_COLOR} ${PBMM_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${MINIGRAPH_MINIGRAPH_COLOR} ${GRAPHALIGNER_FAST_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"

# graphalinger-default on the M/C graph timed out after 7 days
R10_REAL_FULL_HEADLINE_CATEGORIES="${MINIMAP_R10} ${WINNOWMAP} ${GIRAFFE_PRIMARY_R10} ${MINIGRAPH_MINIGRAPH} ${GRAPHALIGNER_MINIGRAPH} ${GRAPHALIGNER_FAST} ${GIRAFFE_R10} ${GIRAFFE_SAMPLED_R10}"
R10_REAL_FULL_HEADLINE_COLORS="${MINIMAP_COLOR} ${WINNOWMAP_COLOR} ${GIRAFFE_PRIMARY_COLOR} ${MINIGRAPH_MINIGRAPH_COLOR} ${GRAPHALIGNER_MINIGRAPH_COLOR} ${GRAPHALIGNER_FAST_COLOR} ${GIRAFFE_COLOR} ${GIRAFFE_SAMPLED_COLOR}"


########################################## ROCs and QQs

# Zoomed in headline roc plots
HIFI_ACCURACY=${ROOT_DIR}/experiments/hifi_sim_1m_headline/results/compared.tsv
R10_ACCURACY=${ROOT_DIR}/experiments/r10y2025_sim_1m_headline/results/compared.tsv

Rscript plot-roc.R ${HIFI_ACCURACY} ${OUT_DIR}/roc_hifi_headline_zoomed.pdf $(echo $HIFI_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "HIFI ROC" $(echo $HIFI_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' ) '0.005,0.0085,0.9905,0.9945'

Rscript plot-roc.R ${R10_ACCURACY} ${OUT_DIR}/roc_r10_headline_zoomed.pdf $(echo $R10_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "R10 ROC" $(echo $R10_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' ) '0.004,0.010,0.990,0.994'

# Zoomed out headline roc plots
Rscript plot-roc.R ${HIFI_ACCURACY} ${OUT_DIR}/roc_hifi_headline.pdf $(echo $HIFI_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "HIFI ROC" $(echo $HIFI_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' )

Rscript plot-roc.R ${R10_ACCURACY} ${OUT_DIR}/roc_r10_headline.pdf $(echo $R10_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "R10 ROC" $(echo $R10_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' )

# headline qq plots
Rscript plot-qq.R ${HIFI_ACCURACY} ${OUT_DIR}/qq_hifi_headline.pdf $(echo $HIFI_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "HIFI QQ" $(echo $HIFI_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' )

Rscript plot-qq.R ${R10_ACCURACY} ${OUT_DIR}/qq_r10_headline.pdf $(echo $R10_SIM_1M_HEADLINE_CATEGORIES | sed 's/ /;/g' ) "R10 QQ" $(echo $R10_SIM_1M_HEADLINE_COLORS | sed 's/ /,/g' )

########################################## Incorrect counts

HIFI_INCORRECT=${ROOT_DIR}/experiments/hifi_sim_1m_headline/results/mapping_stats_sim.tsv
R10_INCORRECT=${ROOT_DIR}/experiments/r10y2025_sim_1m_headline/results/mapping_stats_sim.tsv

python3 barchart.py <(tail -n +2 ${HIFI_INCORRECT} | cut -f1,5) --divisions <(tail -n +2 ${HIFI_INCORRECT} | cut -f1,4)  --title 'HiFi incorrect read count' --y_label 'Count' --x_label 'Condition' --x_sideways --no_n --categories ${HIFI_SIM_1M_HEADLINE_CATEGORIES}  --colors ${HIFI_SIM_1M_HEADLINE_COLORS} --save ${OUT_DIR}/incorrectness_hifi.pdf

python3 barchart.py <(tail -n +2 ${R10_INCORRECT} | cut -f1,5) --divisions <(tail -n +2 ${R10_INCORRECT} | cut -f1,4)  --title 'R10 incorrect read count' --y_label 'Count' --x_label 'Condition' --x_sideways --no_n --categories ${R10_SIM_1M_HEADLINE_CATEGORIES}  --colors ${R10_SIM_1M_HEADLINE_COLORS}  --save ${OUT_DIR}/incorrectness_r10.pdf

########################################## Softclips

HIFI_SOFTCLIPS=${ROOT_DIR}/experiments/hifi_real_full_headline/results/softclipped_or_unmapped_percent.tsv
R10_SOFTCLIPS=${ROOT_DIR}/experiments/r10y2025_real_full_headline/results/softclipped_or_unmapped_percent.tsv

#hifi softclips low
LIMIT=$(python3 get_outlier_limit.py ${HIFI_SOFTCLIPS} small)
python3 barchart.py ${HIFI_SOFTCLIPS} --max ${LIMIT}  --title 'HiFi Softclipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/softclips_hifi_headline_low.pdf

#hifi softclips high
LIMIT=$(python3 get_outlier_limit.py ${HIFI_SOFTCLIPS} big)
python3 barchart.py ${HIFI_SOFTCLIPS} --min ${LIMIT}  --title 'HiFi Softclipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/softclips_hifi_headline_high.pdf

#r10 softclips
python3 barchart.py ${R10_SOFTCLIPS} --title 'R10 Softclipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories ${R10_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${R10_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/softclips_r10_headline_high.pdf


############################################ Runtime

HIFI_RUNTIME=${ROOT_DIR}/experiments/hifi_real_full_headline/results/run_and_sampling_time_from_benchmark.tsv
R10_RUNTIME=${ROOT_DIR}/experiments/r10y2025_real_full_headline/results/run_and_sampling_time_from_benchmark.tsv
HIFI_INDEX_TIME=${ROOT_DIR}/experiments/hifi_real_full_headline/results/index_load_time.tsv
R10_INDEX_TIME=${ROOT_DIR}/experiments/r10y2025_real_full_headline/results/index_load_time.tsv

# Higher hifi runtimes
LIMIT=$(python3 get_outlier_limit.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) big)
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_INDEX_TIME}) --min ${LIMIT}  --title 'HiFi Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_hifi_headline_slow.pdf

# Lower hifi runtimes
LIMIT=$(python3 get_outlier_limit.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) small)
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,$2/60}}' ${HIFI_INDEX_TIME}) --max ${LIMIT}  --title 'HiFi Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_hifi_headline_fast.pdf

# R10 runtimes as one plot
python3 barchart.py <(awk -v OFS='\t' '{{print $1,$2/60}}' ${R10_RUNTIME}) --divisions <(awk -v OFS='\t' '{{print $1,$2/60}}' ${R10_INDEX_TIME}) --title 'R10 Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories ${R10_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${R10_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/runtime_r10_headline.pdf

###################################### Memory

HIFI_MEMORY=${ROOT_DIR}/experiments/hifi_real_full_headline/results/memory_from_benchmark.tsv
R10_MEMORY=${ROOT_DIR}/experiments/r10y2025_real_full_headline/results/memory_from_benchmark.tsv

python3 barchart.py ${HIFI_MEMORY} --title 'HiFi Memory From Benchmark' --y_label 'Memory (GB)' --x_label 'Mapper' --x_sideways --no_n --categories ${HIFI_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${HIFI_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/memory_hifi.pdf

python3 barchart.py ${R10_MEMORY} --title 'R10 Memory From Benchmark' --y_label 'Memory (GB)' --x_label 'Mapper' --x_sideways --no_n --categories ${R10_REAL_FULL_HEADLINE_CATEGORIES}  --colors ${R10_REAL_FULL_HEADLINE_COLORS} --save ${OUT_DIR}/memory_r10.pdf

