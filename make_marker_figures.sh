#!/usr/bin/env bash
# make_marker_figures.sh: Make Long Read Giraffe panel for HPRC R2 marker paper
# Assuming we've made the total summary for the experiment, give its conditions readable names.
cat output/experiments/hprc_y2_marker/results/dv_total_summary.tsv | sed -e 's/hifi,giraffe-hifi,v1.1-sampled4o,chm13v1,.model2025-03-26noinfo/Giraffe-HPRCv1.1-HiFi/g' -e 's/hifi,giraffe-hifi,v2.1-eval-sampled16o,chm13,.model2025-03-26noinfo/Giraffe-HPRCv2.1-HiFi/g' -e 's/hifi,minimap2-map-hifi,v2.1-eval-sampled16o,chm13,.nomodel/Minimap2-CHM13-HiFi/g' -e 's/r10y2025,giraffe-r10,v1.1-sampled4o,chm13v1,.nomodel/Giraffe-HPRCv1.1-R10/g' -e 's/r10y2025,giraffe-r10,v2.1-eval-sampled16o,chm13,.nomodel/Giraffe-HPRCv2.1-R10/g' -e 's/r10y2025,minimap2-lr!hqae,v2.1-eval-sampled16o,chm13,.nomodel/Minimap2-CHM13-R10/g' >renamed.tsv
# Actually plot the plot
Rscript plot-calling-results.R renamed.tsv ./sv_total_summary.pdf "$(cat renamed.tsv | tail -n +1 | cut -f1 | tr '\n' ';')" "Long Read Giraffe + DeepVariant"
