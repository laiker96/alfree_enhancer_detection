#! /bin/bash
MM10_FILE=$1
HG38_FILE=$2

makeblastdb -dbtype prot -in $2 -out blast_db/hg19/hg19
makeblastdb -dbtype prot -in $1 -out blast_db/mm10/mm10


blastp -outfmt 6 -num_threads 12 -max_target_seqs 10 -query $2 -db blast_db/mm10/mm10 > hg19_to_mm10_blastp.txt
blastp -outfmt 6 -num_threads 12 -max_target_seqs 10 -query $1 -db blast_db/hg19/hg19 > mm10_to_hg19_blastp.txt

sort -k1,1 -k12,12gr -k11,11g -k3,3gr hg19_to_mm10_blastp.txt | sort -u -k1,1 --merge > bestHits_hg19_to_mm10_blastp.txt
sort -k1,1 -k12,12gr -k11,11g -k3,3gr mm10_to_hg19_blastp.txt | sort -u -k1,1 --merge > bestHits_mm10_to_hg19_blastp.txt



cat <(awk '{print $2"\t"$1}' bestHits_mm10_to_hg19_blastp.txt) \
	<(awk '{print $1"\t"$2}' bestHits_hg19_to_mm10_blastp.txt) \
	| sort -k 1 | uniq -c | awk '$1 == 2{print$2"\t"$3}' \
	| sort -k 1 | cat <(echo 'hg10ID\tmm10ID') - > ortho_pairs.txt
