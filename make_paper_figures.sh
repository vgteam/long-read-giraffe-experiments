set -ex

ROOT_DIR="/private/groups/patenlab/project-lrg"
OUT_DIR="./plots"


# This can be used to get all the input files
# I left the --dry-run in because it keeps trying to re-run anything
# Asking for all_paper_figures also works
#(umask 0002 ; snakemake -p --dry-run --keep-going --rerun-triggers mtime --resources full_cluster_nodes=2 threads=256 --slurm --rerun-incomplete --latency-wait 120 -j128 --use-singularity --singularity-args '-B /private/home/jmonlong/workspace/lreval/data/:/private/home/jmonlong/workspace/lreval/data/ -B /private/groups/patenlab/project-lrg/:/private/groups/patenlab/project-lrg/ -B /private/groups/patenlab/anovak/projects/hprc/lr-giraffe/:/private/groups/patenlab/anovak/projects/hprc/lr-giraffe/' \
#/private/groups/patenlab/project-lrg/experiments/hifi_sim_1m_headline/results/compared.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10_sim_1m_headline/results/compared.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_sim_1m_headline/results/mapping_stats_sim.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10_sim_1m_headline/results/mapping_stats_sim.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_real_full_headline/results/softclipped_or_unmapped_percent.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10_real_full_headline/results/softclipped_or_unmapped_percent.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_real_full_headline/results/run_and_sampling_time_from_benchmark.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10_real_full_headline/results/run_and_sampling_time_from_benchmark.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_real_full_headline/results/index_load_time.tsv \
#/private/groups/patenlab/project-lrg/experiments/hifi_real_full_headline/results/memory_from_benchmark.tsv \
#/private/groups/patenlab/project-lrg/experiments/r10_real_full_headline/results/memory_from_benchmark.tsv)



########################################## ROCs and QQs

# Zoomed in headline roc plots
HIFI_ACCURACY=${ROOT_DIR}/experiments/hifi_sim_1m_headline/results/compared.tsv
R10_ACCURACY=${ROOT_DIR}/experiments/r10_sim_1m_headline/results/compared.tsv

Rscript plot-roc-hifi.R ${HIFI_ACCURACY} ${OUT_DIR}/roc_hifi_headline_zoomed.pdf

Rscript plot-roc-r10.R ${R10_ACCURACY} ${OUT_DIR}/roc_r10_headline_zoomed.pdf

# Zoomed out headline roc plots
Rscript plot-roc.R ${HIFI_ACCURACY} ${OUT_DIR}/roc_hifi_headline.pdf

Rscript plot-roc.R ${R10_ACCURACY} ${OUT_DIR}/roc_r10_headline.pdf

# headline qq plots
Rscript plot-qq.R ${HIFI_ACCURACY} ${OUT_DIR}/qq_hifi_headline.pdf

Rscript plot-qq.R ${R10_ACCURACY} ${OUT_DIR}/qq_r10_headline.pdf

######################################### Incorrect counts

HIFI_INCORRECT=${ROOT_DIR}/experiments/hifi_sim_1m_headline/results/mapping_stats_sim.tsv
R10_INCORRECT=${ROOT_DIR}/experiments/r10_sim_1m_headline/results/mapping_stats_sim.tsv

python3 barchart.py <(tail -n +2 ${HIFI_INCORRECT} | cut -f1,5) --divisions <(tail -n +2 ${HIFI_INCORRECT} | cut -f1,4)  --title 'HiFi incorrect read count' --y_label 'Count' --x_label 'Condition' --x_sideways --no_n --save ${OUT_DIR}/incorrectness_hifi.pdf

python3 barchart.py <(tail -n +2 ${R10_INCORRECT} | cut -f1,5) --divisions <(tail -n +2 ${R10_INCORRECT} | cut -f1,4)  --title 'R10 incorrect read count' --y_label 'Count' --x_label 'Condition' --x_sideways --no_n --save ${OUT_DIR}/incorrectness_r10.pdf

########################################## Softclips

HIFI_SOFTCLIPS=${ROOT_DIR}/experiments/hifi_real_full_headline/results/softclipped_or_unmapped_percent.tsv
R10_SOFTCLIPS=${ROOT_DIR}/experiments/r10_real_full_headline/results/softclipped_or_unmapped_percent.tsv

#hifi softclips low
LIMIT=$(python3 get_outlier_limit.py ${HIFI_SOFTCLIPS} small)
python3 barchart.py ${HIFI_SOFTCLIPS} --max ${LIMIT}  --title 'HiFi Softclipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories minimap2-map-hifi,hprc-v1.1-mc-d9 winnowmap,hprc-v1.1-mc-d9 giraffe,primary minigraph,hprc-v1.0-minigraph graphaligner-default, hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc giraffe,hprc-v1.1-mc-sampled16  --colors '#F57600' '#C44601' '#3C96B0' '#F2C300' '#C49AB6' '#53AEDF' '#27A974' '#9B8BF4' --save ${OUT_DIR}/softclips_hifi_headline_low.pdf

#hifi softclips high
LIMIT=$(python3 get_outlier_limit.py ${HIFI_SOFTCLIPS} big)
python3 barchart.py ${HIFI_SOFTCLIPS} --min ${LIMIT}  --title 'HiFi Softclipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories minimap2-map-hifi,hprc-v1.1-mc-d9 winnowmap,hprc-v1.1-mc-d9 giraffe,primary minigraph,hprc-v1.0-minigraph graphaligner-default, hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc giraffe,hprc-v1.1-mc-sampled16  --colors '#F57600' '#C44601' '#3C96B0' '#F2C300' '#C49AB6' '#53AEDF' '#27A974' '#9B8BF4' --save ${OUT_DIR}/softclips_hifi_headline_high.pdf

#r10 softclips
python3 barchart.py ${R10_SOFTCLIPS} --title 'R10 Softclipped or Unmapped Bases' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories minimap2-lr:hqae,hprc-v1.1-mc-d9 winnowmap,hprc-v1.1-mc-d9 giraffe,primary minigraph,hprc-v1.0-minigraph graphaligner-default, hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc giraffe,hprc-v1.1-mc-sampled16  --colors '#F57600' '#C44601' '#3C96B0' '#F2C300' '#C49AB6' '#53AEDF' '#27A974' '#9B8BF4' --save ${OUT_DIR}/softclips_r10_headline_high.pdf


############################################ Runtime

HIFI_RUNTIME=${ROOT_DIR}/experiments/hifi_real_full_headline/results/run_and_sampling_time_from_benchmark.tsv
R10_RUNTIME=${ROOT_DIR}/experiments/r10_real_full_headline/results/run_and_sampling_time_from_benchmark.tsv
HIFI_INDEX_TIME=${ROOT_DIR}/experiments/hifi_real_full_headline/results/index_load_time.tsv
R10_INDEX_TIME=${ROOT_DIR}/experiments/r10_real_full_headline/results/index_load_time.tsv

# Higher hifi runtimes
LIMIT=$(python3 get_outlier_limit.py ${HIFI_RUNTIME} big)
python3 barchart.py <(awk -v OFS='\\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) --divisions <(awk -v OFS='\\t' '{{print $1,$2/60}}' ${HIFI_INDEX_TIME}) --min ${LIMIT}  --title 'HiFi Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories minimap2-map-hifi,hprc-v1.1-mc-d9 winnowmap,hprc-v1.1-mc-d9 giraffe,primary minigraph,hprc-v1.0-minigraph graphaligner-default, hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc giraffe,hprc-v1.1-mc-sampled16  --colors '#F57600' '#C44601' '#3C96B0' '#F2C300' '#C49AB6' '#53AEDF' '#27A974' '#9B8BF4' --save ${OUT_DIR}/runtime_hifi_headline_slow.pdf

# Lower hifi runtimes
LIMIT=$(python3 get_outlier_limit.py ${HIFI_RUNTIME} small)
python3 barchart.py <(awk -v OFS='\\t' '{{print $1,$2/60}}' ${HIFI_RUNTIME}) --divisions <(awk -v OFS='\\t' '{{print $1,$2/60}}' ${HIFI_INDEX_TIME}) --max ${LIMIT}  --title 'HiFi Runtime' --y_label 'Time (hours)' --x_label 'Mapper' --x_sideways --no_n --categories minimap2-map-hifi,hprc-v1.1-mc-d9 winnowmap,hprc-v1.1-mc-d9 giraffe,primary minigraph,hprc-v1.0-minigraph graphaligner-default, hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc giraffe,hprc-v1.1-mc-sampled16  --colors '#F57600' '#C44601' '#3C96B0' '#F2C300' '#C49AB6' '#53AEDF' '#27A974' '#9B8BF4' --save ${OUT_DIR}/runtime_hifi_headline_fast.pdf

# R10 runtimes as one plot
python3 barchart.py <(awk -v OFS='\\t' '{{print $1,$2/60}}' ${R10_RUNTIME}) --divisions <(awk -v OFS='\\t' '{{print $1,$2/60}}' ${R10_INDEX_TIME}) --title 'R10 Runtime' --y_label 'Percent of bases' --x_label 'Mapper' --x_sideways --no_n --categories minimap2-lr:hqae,hprc-v1.1-mc-d9 winnowmap,hprc-v1.1-mc-d9 giraffe,primary minigraph,hprc-v1.0-minigraph graphaligner-default, hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc giraffe,hprc-v1.1-mc-sampled16  --colors '#F57600' '#C44601' '#3C96B0' '#F2C300' '#C49AB6' '#53AEDF' '#27A974' '#9B8BF4' --save ${OUT_DIR}/runtime_r10_headline.pdf

###################################### Memory

HIFI_MEMORY=${ROOT_DIR}/experiments/hifi_real_full_headline/results/memory_from_benchmark.tsv
R10_MEMORY=${ROOT_DIR}/experiments/r10_real_full_headline/results/memory_from_benchmark.tsv

python3 barchart.py ${HIFI_MEMORY} --title 'HiFi Memory From Benchmark' --y_label 'Memory (GB)' --x_label 'Mapper' --x_sideways --no_n --categories minimap2-map-hifi,hprc-v1.1-mc-d9 winnowmap,hprc-v1.1-mc-d9 giraffe,primary minigraph,hprc-v1.0-minigraph graphaligner-default, hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc giraffe,hprc-v1.1-mc-sampled16  --colors '#F57600' '#C44601' '#3C96B0' '#F2C300' '#C49AB6' '#53AEDF' '#27A974' '#9B8BF4'  --save ${OUT_DIR}/memory_hifi.pdf

python3 barchart.py ${HIFI_MEMORY} --title 'HiFi Memory From Benchmark' --y_label 'Memory (GB)' --x_label 'Mapper' --x_sideways --no_n --categories minimap2-lr:hqae,hprc-v1.1-mc-d9 winnowmap,hprc-v1.1-mc-d9 giraffe,primary minigraph,hprc-v1.0-minigraph graphaligner-default, hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc-d9 giraffe,hprc-v1.1-mc giraffe,hprc-v1.1-mc-sampled16  --colors '#F57600' '#C44601' '#3C96B0' '#F2C300' '#C49AB6' '#53AEDF' '#27A974' '#9B8BF4' --save ${OUT_DIR}/memory_r10.pdf

