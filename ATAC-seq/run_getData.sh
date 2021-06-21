#!/bin/bash
get_data.sh -g LN.ATAC -m LN.ATAC.R1 -i LN.ATAC.R1.fastq -x mm10
get_data.sh -g LN.ATAC -m LN.ATAC.R2 -i LN.ATAC.R2.fastq -x mm10
get_data.sh -g DM.ATAC -m DM.ATAC.R1 -i DM.ATAC.R1 .fastq -x mm10
get_data.sh -g DM.ATAC -m DM.ATAC.R2 -i DM.ATAC.R2.fastq -x mm10
get_data.sh -g FLT3.ATAC -m FLT3.ATAC.R1 -i FLT3.ATAC.R1.fastq -x mm10
get_data.sh -g FLT3.ATAC -m FLT3.ATAC.R2 -i FLT3.ATAC.R2.fastq -x mm10
get_data.sh -g NPM1.ATAC -m NPM1.ATAC.R1 -i NPM1.ATAC.R1.fastq -x mm10
get_data.sh -g NPM1.ATAC -m NPM1.ATAC.R2 -i NPM1.ATAC.R2.fastq -x mm10
