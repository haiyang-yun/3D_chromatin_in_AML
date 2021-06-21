# histone modifications - data mapping and peak calling

## reads mapping
get_data.sh
run_getData.sh

## process reads
process_aligned_reads.sh
run_process_aligned_reads.sh

## peakcalling for H3K4me1 and H3K4me3 (using HOMER)
findpeaks.sh
run_findpeaks.sh
filter_peaks.sh
filter_overlapped_intervals.py
filterpeaks_based_on_index.py

## define H3K4me1-marked enhancers
enhancer_calling.R

## peakcalling for H3K27ac (using MACS2)
peakCalling_H3K27ac.sh
run_peakCalling_H3K27ac.sh
run_MACS.NoModel.pl






