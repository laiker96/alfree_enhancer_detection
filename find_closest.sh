#! /bin/bash

BED_FILE=$1
GENE_FILE=$2

sort-bed --max-mem 8G "$GENE_FILE" > tmp && mv tmp "$GENE_FILE"
sort-bed --max-mem 8G "$BED_FILE" > tmp && mv tmp "$BED_FILE"

bedtools closest -io -iu -D ref -t first \
	-a $BED_FILE -b $GENE_FILE > downstreamGenes_$BED_FILE


bedtools closest -io -id -D ref -t first \
	-a $BED_FILE -b $GENE_FILE > upstreamGenes_$BED_FILE


bedmap --echo-map --exact downstreamGenes_$BED_FILE upstreamGenes_$BED_FILE > tmp1
bedmap --echo-map --exact upstreamGenes_$BED_FILE downstreamGenes_$BED_FILE > tmp2

sort-bed --max-mem 8G tmp1 > tmp1_s && mv tmp1_s tmp1
sort-bed --max-mem 8G tmp2 > tmp2_s && mv tmp2_s tmp2
paste tmp1 tmp2 > closest_genes_"$BED_FILE"
rm tmp1 tmp2 downstreamGenes_$BED_FILE upstreamGenes_$BED_FILE

sort-bed --max-mem 8G closest_genes_"$BED_FILE" > tmp && mv tmp closest_genes_"$BED_FILE"
awk '$10 == $26' closest_genes_"$BED_FILE" \
	| awk '$10 != "." && $26 != "." {print}' \
	| awk 'BEGIN{FS="\t";OFS="\t"} {if($11 < $27) print $1,$2,$3,$10,$12,$27; else print $1,$2,$3,$10,$28,$11}'  \
	| awk '$6 > $5' | awk '$3 - $2 > 200' | awk '$6 - $5 < 600000' \
	| awk 'BEGIN{FS="\t";OFS="\t"} i+=1 {print $1,$2,$3,"enhancer_"i"_"$4":"$5"-"$6}' > annot_"$BED_FILE"


grep -oP "chr[[:alnum:]]+:[0-9]+-[0-9]+" annot_"$BED_FILE" \
	| sed 's/:/\t/' | sed 's/-/\t/' \
	| awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$3,$1":"$2"-"$3}' \
	| sort-bed --max-mem 8G - | uniq > rregions.bed 

#bedtools getfasta -nameOnly -fi hg19.fa -bed annot_FANTOM5_enhancers.bed > annot_FANTOM5_enhancers.fa
#bedtools getfasta -nameOnly -fi mm10.fa -bed rregions.bed > rregions.fa
