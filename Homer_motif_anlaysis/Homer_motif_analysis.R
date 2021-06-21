## motif predition for regions with mutation-induced altered accessiblity
findMotifsGenome.pl ATAC.DM.down.bed mm10 ATAC.DM.down/ -size given -mask -p 12 -S 10
findMotifsGenome.pl ATAC.DM.up.bed mm10 ATAC.DM.up/ -size given -mask -p 12 -S 10
findMotifsGenome.pl ATAC.FLT3.down.bed mm10 ATAC.FLT3.down/ -size given -mask -p 12 -S 10
findMotifsGenome.pl ATAC.FLT3.up.bed mm10 ATAC.FLT3.up/ -size given -mask -p 12 -S 10
findMotifsGenome.pl ATAC.NPM1.down.bed mm10 ATAC.NPM1.down/ -size given -mask -p 12 -S 10
findMotifsGenome.pl ATAC.NPM1.up.bed mm10 ATAC.NPM1.up/ -size given -mask -p 12 -S 10

## motif prediction for regions with DM gain or loss enhancer marks (for identifying leukemia-specific changes)
findMotifsGenome.pl DM.gain-1.bed mm10 DM.gain-1/ -size given -mask -p 12 -S 10
findMotifsGenome.pl DM.gain-2.bed mm10 DM.gain-2/ -size given -mask -p 12 -S 10
findMotifsGenome.pl DM.loss-1.bed mm10 DM.loss-1/ -size given -mask -p 12 -S 10
findMotifsGenome.pl DM.loss-2.bed mm10 DM.loss-2/ -size given -mask -p 12 -S 10

## motif predition for regions from Seurat-guided clustering analysis
findMotifsGenome.pl tsne.cluster-6.bed mm10 tsne.cluster-6/ -size given -mask -p 12 -S 10
findMotifsGenome.pl tsne.cluster-7.bed mm10 tsne.cluster-7/ -size given -mask -p 12 -S 10
findMotifsGenome.pl tsne.cluster-8.bed mm10 tsne.cluster-8/ -size given -mask -p 12 -S 10
findMotifsGenome.pl tsne.cluster-9.bed mm10 tsne.cluster-9/ -size given -mask -p 12 -S 10
findMotifsGenome.pl tsne.cluster-10.bed mm10 tsne.cluster-10/ -size given -mask -p 12 -S 10