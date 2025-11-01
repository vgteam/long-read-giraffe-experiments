set -ex

ROOT_DIR=/private/groups/patenlab/project-lrg
OUT_DIR=./hprc_plots
mkdir -p "${OUT_DIR}"

GIRAFFE_V1_HIFI=hifi,giraffe-hifi-mmcs100,v1.1-sampled32d,.model2025-03-26noinfo,chm13v1
GIRAFFE_V2_HIFI=hifi,giraffe-hifi-mmcs100,v2.0-eval.ec1M-sampled16o,.model2025-03-26noinfo,chm13
MINIMAP_HIFI=hifi,minimap2-map-hifi,v2.0-eval-d46,.nomodel,chm13

GIRAFFE_V1_R10=r10y2025,giraffe-r10-noflags,v1.1-sampled32d,.nomodel,chm13v1
GIRAFFE_V2_R10=r10y2025,giraffe-r10-noflags,v2.0-eval.ec1M-sampled16o,.nomodel,chm13
MINIMAP_R10=r10y2025,minimap2-lr-hqae,v2.0-eval-d46,.nomodel,chm13

#Rename conditions for output
GIRAFFE_V1_HIFI_NEW=Giraffe-HPRCv1.1-HiFi
GIRAFFE_V2_HIFI_NEW=Giraffe-HPRCv2-HiFi
MINIMAP_HIFI_NEW=Minimap2-CHM13-HiFi

GIRAFFE_V1_R10_NEW=Giraffe-HPRCv1.1-R10
GIRAFFE_V2_R10_NEW=Giraffe-HPRCv2-R10
MINIMAP_R10_NEW=Minimap2-CHM13-R10

DV_CALLING_CONDITIONS="${GIRAFFE_V1_HIFI_NEW};${GIRAFFE_V2_HIFI_NEW};${MINIMAP_HIFI_NEW};${GIRAFFE_V1_R10_NEW};${GIRAFFE_V2_R10_NEW};${MINIMAP_R10_NEW}"

SNP_SUMMARY=${ROOT_DIR}/experiments/hprc_dv_calling/results/dv_snp_summary.tsv

sed "s/${GIRAFFE_V1_HIFI}/${GIRAFFE_V1_HIFI_NEW}/g" ${SNP_SUMMARY} | sed "s/${GIRAFFE_V2_HIFI}/${GIRAFFE_V2_HIFI_NEW}/g" | sed "s/${MINIMAP_HIFI}/${MINIMAP_HIFI_NEW}/g" | sed "s/${GIRAFFE_V1_R10}/${GIRAFFE_V1_R10_NEW}/g"  | sed "s/${GIRAFFE_V2_R10}/${GIRAFFE_V2_R10_NEW}/g" | sed "s/${MINIMAP_R10}/${MINIMAP_R10_NEW}/g" > ${OUT_DIR}/renamed_stats.tsv
Rscript plot-calling-results.R ${OUT_DIR}/renamed_stats.tsv ${OUT_DIR}/dv_snp_summary.pdf ${DV_CALLING_CONDITIONS}
rm ${OUT_DIR}/renamed_stats.tsv

