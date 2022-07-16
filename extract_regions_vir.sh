#! /bin/bash
GFF=$1
FASTA_FILE=$2


function extractCDS() {

	ID=$(basename "$1" .gtf.gz)
	zgrep -v '^#' "$1" \
	  | sed 's/gene_id /gene_id=/' | sed 's/transcript_id /transcript_id=/' \
	  | grep -P 'transcript_id="NM_[0-9]+"' \
	  | sed  's/\"//g' \
	  | awk 'BEGIN{FS="\t";OFS="\t"}($3=="CDS")' \
	  | cut -f 1,4,5,7,9 \
	  | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2-1,$3,$5,"1000",$4}' \
	  | sort -u | sort -k1,1 -k2,2n > "$ID"_CDS.bed
	
	cut -f 4 "$ID"_CDS.bed \
		| grep -owP "transcript_id=NM_[0-9]+" > column4_tmp

	paste <(cut -f 1-3 "$ID"_CDS.bed) column4_tmp <(cut -f 1-4 --complement "$ID"_CDS.bed) > tmp
	mv tmp "$ID"_CDS.bed
  rm column4_tmp

	awk -F "\t" '$6 == "+" {print}' "$ID"_CDS.bed > "$ID"_CDS_pluss.bed
	awk -F "\t" '$6 == "-" {print}' "$ID"_CDS.bed > "$ID"_CDS_negs.bed

}

function joindCDS() {

	cut -f 2 -d ";" "$1" | sort -g | uniq | grep -oP "NM.*\.[0-9]+" | sort -u > ids.txt
	while read j; do
		
		grep -Fw "$j" "$1" > tmp.bed
		if [[ $2 == "+" ]]; then
			bedtools getfasta -nameOnly -s -fi "$3" -bed tmp.bed \
				| sed 's/(.*)//' > tmp.fa
		else
			tac tmp.bed | bedtools getfasta -nameOnly -s -fi "$3" -bed - \
				| sed 's/(.*)//' > tmp.fa

		fi
		
		cat tmp.fa | sed -e '1!{/^>.*/d;}' \
			| sed  ':a;N;$!ba;s/\n//2g' >> "$4".fa
		
	done < ids.txt
	rm ids.txt tmp.bed tmp.fa


}

function extractGeneModels() {

	ID=$(basename "$1" .gtf.gz)
	zgrep -v '^#' "$1" \
	  | sed 's/gene_id /gene_id=/' | sed 's/transcript_id /transcript_id=/' \
	  | grep -P 'transcript_id="NM_[0-9]+\.[0-9]+"' \
	  | sed  's/\"//g' \
	  | awk 'BEGIN{FS="\t";OFS="\t"}($3=="exon")' \
	  | cut -f 1,4,5,7,9 \
	  | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2-1,$3,$5,"1000",$4}' \
	  | sort -u | sort -k1,1 -k2,2n > "$ID"_exons.bed
	
	cut -f 4 "$ID"_exons.bed \
		| grep -owP "transcript_id=NM_[0-9]+\.[0-9]+" > column4_tmp

	paste <(cut -f 1-3 "$ID"_exons.bed) column4_tmp <(cut -f 1-4 --complement "$ID"_exons.bed) > tmp
	mv tmp "$ID"_exons.bed
	sort -g column4_tmp | uniq > tmp && mv tmp column4_tmp

	while read j; do
		grep -Fw "$j" "$ID"_exons.bed | bedtools merge -i - -c 4,5,6 \
			-o distinct,distinct,distinct \
			-d 10000000000000 >> "$ID"_geneModels.bed
	done < column4_tmp
  rm column4_tmp

}



ID=$(basename "$GFF" .gtf.gz)
extractCDS "$GFF"
joindCDS "$ID"_CDS_pluss.bed "+" $2 "$ID"
joindCDS "$ID"_CDS_negs.bed "-" $2 "$ID"
source activate ML
python translate_sequence.py "$ID".fa
extractGeneModels "$GFF"
