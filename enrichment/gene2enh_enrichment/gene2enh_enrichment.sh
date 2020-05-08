
step=$1

CODE=/ahg/regevdata/projects/Cell2CellCommunication/code/genome_utils/enrichment/gene2enh_enrichment
source /broad/software/scripts/useuse
use .bedtools-2.26.0

my_rscript=/ahg/regevdata/users/oursu/software/anaconda3/bin/Rscript
LOLA_script=/ahg/regevdata/projects/Cell2CellCommunication/code/genome_utils/enrichment/LOLA/run_LOLA_basic.R

TSS_FILE=$2
REG_FILE=$3
GENE_REG_WINDOW=$4
WINDOW_ENRICHMENT_DATASET=$5
OUT=$6
PREF=$7
GENE_SETS_METADATA=$8
ROI_SETS_METADATA=$9
GOI_COL=${10}

mkdir -p ${OUT}

if [[ ${step} == "genes2reg" ]];
then
    zcat -f ${TSS_FILE} | cut -f1-4 | gzip > ${OUT}/${PREF}.genes2reg.tss.gz
    zcat -f ${REG_FILE} | cut -f1-3 | gzip > ${OUT}/${PREF}.genes2reg.reg.gz
    
    bedtools window -w ${GENE_REG_WINDOW} -a ${OUT}/${PREF}.genes2reg.tss.gz -b ${OUT}/${PREF}.genes2reg.reg.gz | gzip > ${OUT}/${PREF}.genes2reg.gz
    zcat -f ${OUT}/${PREF}.genes2reg.gz | cut -f5,6,7 | sort | uniq | gzip > ${OUT}/${PREF}.universe.gz
    echo "here you go: ${OUT}/${PREF}.genes2reg.gz"
    echo "universe: "$(zcat -f ${OUT}/${PREF}.universe.gz | wc -l)
    rm ${OUT}/${PREF}.genes2reg.tss.gz ${OUT}/${PREF}.genes2reg.reg.gz 
fi

if [[ ${step} == "enrich" ]];
then

    outfile=${OUT}/${PREF}.enrichment.csv
    echo "Gene_set,Region_set,universe_size,region_size,gene_size,overlap_size,odds_ratio,p" | sed 's/,/\t/g' > ${outfile}
    universe=${OUT}/${PREF}.universe.gz
    while IFS= read -r line #roi sets                                                      
    do
	line_roi=${line}
	roi_set=$(echo "${line_roi}" | cut -f1)
	roi_set_file=$(echo "${line_roi}" | cut -f2)
	echo "roi set: ${roi_set}"
	
	#annotate universe with roi
	universe_in_roi=${OUT}/${PREF}.ROI${roi_set}.std.gz
	bedtools window -w ${WINDOW_ENRICHMENT_DATASET} -a ${universe} -b ${roi_set_file} | cut -f1-3 | sort | uniq | gzip > ${universe_in_roi}
    
	while IFS= read -r line #gene sets                                                                   
	do
	    line_gene=${line}
            gene_set=$(echo "${line_gene}" | cut -f1)
            gene_set_file=$(echo "${line_gene}" | cut -f2)

	    #annotate universe with the genes, using the gene2reg file
	    universe_in_goi=${OUT}/${PREF}.GOI${gene_set}.std.gz
	    join -2 ${GOI_COL} -1 4 <(zcat -f ${OUT}/${PREF}.genes2reg.gz | sort -k4) <(zcat -f ${gene_set_file} | sort -k${GOI_COL} | uniq) | sed 's/ /\t/g' | cut -f5,6,7 | sort | uniq | gzip > ${universe_in_goi}

	    #overlap the gene elts with roi elts
	    overlap=${OUT}/${PREF}.GOI${gene_set}.ROI${roi_set}.overlap.gz
	    bedtools window -w 0 -a ${universe_in_roi} -b ${universe_in_goi} | gzip > ${overlap}

	    size_universe=$(zcat -f ${universe} | sort | uniq | wc -l)
	    size_roi=$(zcat -f ${universe_in_roi} | sort |uniq | wc -l)
	    size_goi=$(zcat -f ${universe_in_goi} | sort |uniq | wc -l)
	    size_overlap=$(zcat -f ${overlap} | sort | uniq | wc -l)
	    test_out=$(python ${CODE}/enrich.py --universe ${size_universe} --set1 ${size_roi} --set2 ${size_goi} --overlap ${size_overlap})
	    odds=$(echo "${test_out}" | cut -d " " -f1)
	    p=$(echo "${test_out}" | cut -d " " -f2)
	    echo "${gene_set},${roi_set},${size_universe},${size_roi},${size_goi},${size_overlap},${odds},${p}" | sed 's/,/\t/g' >> ${outfile}
	    rm ${universe_in_goi} ${overlap}
	done < "${GENE_SETS_METADATA}"
	rm ${universe_in_roi}
    done < "${ROI_SETS_METADATA}"
    echo "this will make your day :) ${outfile}"

    rm ${universe} ${OUT}/${PREF}.genes2reg.gz
fi

if [[ ${step} == "LOLA" ]];
then

    outfile=${OUT}/${PREF}.LOLA.csv
    universe=${OUT}/${PREF}.universe.gz
    
    while IFS= read -r line #gene sets                                                                   
    do
	line_gene=${line}
        gene_set=$(echo "${line_gene}" | cut -f1)
        gene_set_file=$(echo "${line_gene}" | cut -f2)
	
	#annotate universe with the genes, using the gene2reg file
	universe_in_goi=${OUT}/${PREF}.GOI${gene_set}.std.gz
	join -2 ${GOI_COL} -1 4 <(zcat -f ${OUT}/${PREF}.genes2reg.gz | sort -k4) <(zcat -f ${gene_set_file} | sort -k${GOI_COL} | uniq) | sed 's/ /\t/g' | cut -f5,6,7 | sort | uniq | gzip > ${universe_in_goi}

	outpref=${OUT}/${PREF}.GOI${gene_set}
	
	direction=enrichment
	echo "${my_rscript} ${LOLA_script} --set ${universe_in_goi} --outpref ${outpref} --direction ${direction} --universe ${universe}"
    done < "${GENE_SETS_METADATA}"
fi
