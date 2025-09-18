#!/usr/bin/bash

DEST_DIR="CRISPResso2_indels_plots"
mkdir -p "$DEST_DIR" 

for DIR in *CRISPResso2_analysis
do
	SAMPLE="${DIR/_CRISPResso2_analysis/}"
	PLOT_FILE="${SAMPLE}_indels_plot.png"
	cp ${DIR}/CRISPResso_on*/3a.Indel_size_distribution.png "${DEST_DIR}/$PLOT_FILE"
done
	
	
	
