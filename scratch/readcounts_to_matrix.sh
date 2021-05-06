### this is the clunky way I turned output from STAR into a matrix of read counts for edgeR

ls *ReadsPerGene.out.tab | paste -s - > samplenames.txt

for i in *ReadsPerGene.out.tab; do awk '{print $2}' $i > $i.unstranded; done

cat samplenames.txt unstranded_counts.txt > unstranded_counts_matrix.txt

awk '{print $1}' X2654PLUS_S14_L002_ReadsPerGene.out.tab | paste - unstranded_counts_matrix.txt > readcountsmatrix.txt
