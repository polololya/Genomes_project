rm -r output
mkdir output || true


awk 'BEGIN {OFS="\t"}
{
    if ($7 == "+") {print $1,$9,$3,$4,$7}
	else if ($7 == "-") {print $1,$9,$3,$5-1,$7}
}' data/TAIR10_GFF3_genes.gff | grep 'mRNA' | awk 'BEGIN {OFS="\t"}{print $1,$2,$4,$5}' > output/tss_arabidopsis.txt



