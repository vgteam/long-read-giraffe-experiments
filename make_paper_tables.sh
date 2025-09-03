set -ex

ROOT_DIR=/private/groups/patenlab/project-lrg
OUT_DIR=./tables
mkdir -p "${OUT_DIR}"

############################################### All input files


HIFI_SIM=${ROOT_DIR}/experiments/hifi_sim_1m/results/mapping_stats_sim.tsv
R10_SIM=${ROOT_DIR}/experiments/r10y2025_sim_1m/results/mapping_stats_sim.tsv
ILLUMINA_SIM=${ROOT_DIR}/experiments/illumina_sim_1m/results/mapping_stats_sim.tsv
ELEMENT_SIM=${ROOT_DIR}/experiments/element_sim_1m/results/mapping_stats_sim.tsv

HIFI_REAL=${ROOT_DIR}/experiments/hifi_real_full/results/mapping_stats_real.tsv
R10_REAL=${ROOT_DIR}/experiments/r10y2025_real_full/results/mapping_stats_real.tsv
ILLUMINA_REAL=${ROOT_DIR}/experiments/illumina_real_full/results/mapping_stats_real.tsv
ELEMENT_REAL=${ROOT_DIR}/experiments/element_real_full/results/mapping_stats_real.tsv

SV_SUMMARY=${ROOT_DIR}/experiments/sv_calling/results/sv_summary.tsv
SNP_SUMMARY=${ROOT_DIR}/experiments/dv_calling/results/dv_snp_summary.tsv
INDEL_SUMMARY=${ROOT_DIR}/experiments/dv_calling/results/dv_indel_summary.tsv

################################### Simulated read tables

# Assumes 1 mil reads
for SIM_FILE in $HIFI_SIM $R10_SIM $ILLUMINA_SIM $ELEMENT_SIM ;
do

    CURRENT_OUTFILE=${OUT_DIR}/trash.txt
    if [ $SIM_FILE == $HIFI_SIM ]; then
        CURRENT_OUTFILE=${OUT_DIR}/simulated_hifi.tsv
    elif [ $SIM_FILE == $R10_SIM ]; then
        CURRENT_OUTFILE=${OUT_DIR}/simulated_r10.tsv
    elif [ $SIM_FILE == $ILLUMINA_SIM ]; then
        CURRENT_OUTFILE=${OUT_DIR}/simulated_illumina.tsv
    elif [ $SIM_FILE == $ELEMENT_SIM ]; then
        CURRENT_OUTFILE=${OUT_DIR}/simulated_element.tsv
    fi    
    printf "Reference\t& Mapper\t& Correct\t& Mapq60\t& Wrong \\\\\\ \n" >$CURRENT_OUTFILE
    printf "\t& \t& (\\%%) \t& (\\%%) \t& Mapq60 (\\%%) \\\\\\ \n" >>$CURRENT_OUTFILE
    printf "\\hline\n" >> $CURRENT_OUTFILE
    tail -n +2 ${SIM_FILE} | awk -v OFS='\t' '{print $1,$2/10000,$3/10000,$4/10000}' | ./format_tsv.sh >> $CURRENT_OUTFILE 
done

################################### Real read tables


for REAL_FILE in $HIFI_REAL $R10_REAL ;
do

    CURRENT_OUTFILE=${OUT_DIR}/trash.txt
    if [ $REAL_FILE == $HIFI_REAL ]; then
        CURRENT_OUTFILE=${OUT_DIR}/real_hifi.tsv
    elif [ $REAL_FILE == $R10_REAL ]; then
        CURRENT_OUTFILE=${OUT_DIR}/real_r10.tsv
    fi    
    printf "Reference\t& Mapper\t& Runtime (min)\t& Index time (min)\t& Memory (GB)\t& Soft- or Hard-&\t Unmapped \\\\\\ \n" >$CURRENT_OUTFILE
    printf "\t& \t& \t& (+ Sampling (min)) \t& \t& clipped Bases (\\%%)\t& Bases (\\%%) \\\\\\ \n" >>$CURRENT_OUTFILE
    printf "\\hline\n" >> $CURRENT_OUTFILE
    tail -n +2 ${REAL_FILE} | awk -v OFS='\t' '{printf "%0s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",$1,$2,$3,$5,(100*($6+$7))/$10,(100*$8)/$10}' | ./format_tsv.sh >> $CURRENT_OUTFILE 
    printf "SAMPLING TIMES\n" >>$CURRENT_OUTFILE
    tail -n +2 ${REAL_FILE} | awk -v OFS='\t' '$4 != 0 {print $1,$4}'  >> $CURRENT_OUTFILE 
done

################################# DV Table

for TECH in HIFI R10 ;
do
    for VAR in SNP INDEL ;
        do
        VAR_FILE=${SNP_SUMMARY}
        if [ $VAR == SNP ]; then
            VAR_FILE=${SNP_SUMMARY}
        elif [ $VAR == INDEL ]; then
            VAR_FILE=${INDEL_SUMMARY} 
        elif [ $VAR == SV ]; then
            VAR_FILE=${SV_SUMMARY} 
        fi 
        CURRENT_OUTFILE=${OUT_DIR}/${VAR}_${TECH}.tsv
        printf "Reference\t& Mapper\t& True\t& False\t& False\t& Recall\t& Precision\t& F1 \\\\\\ \n" >$CURRENT_OUTFILE
        printf "\t& \t& Positive\t& Negative\t& Positive\t& \t& \t&  \\\\\\ \n" >>$CURRENT_OUTFILE
        printf "\\hline\n" >> $CURRENT_OUTFILE
        TEMP_FILE=${OUT_DIR}/tmp.txt
        if [ $TECH == HIFI ]; then
            tail -n +2 ${VAR_FILE} | grep hifi | sed 's/^hifi,//g' >$TEMP_FILE 
        elif [ $TECH == R10 ]; then
            tail -n +2 ${VAR_FILE} | grep r10 | sed 's/^r10y2025,//g' >$TEMP_FILE 
        fi
        cat ${TEMP_FILE} | sed -r 's/,.(model2025-03-26noinfo|nomodel)//' | sed 's/eval/HPRC/g' | sed -r 's/giraffe-(r10|hifi)-(noflags|mmcs100)/Giraffe/g' | awk -v OFS='\t' '{printf "%0s\t%0i\t%0i\t%0i\t%.4f\t%.4f\t%.4f\n",$1,$2,$3,$4,$5,$6,$7}' | ./format_tsv.sh >> $CURRENT_OUTFILE
        rm $TEMP_FILE
    done
done
