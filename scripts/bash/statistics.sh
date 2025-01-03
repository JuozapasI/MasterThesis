#!/bin/bash

bam=$1
unassigned_bam_10x=$2
intersecting_bam_10x=$3
intergenic_bam_10x=$4
AT_seq_bam=$5
intergenic_good_bed=$6
intergenic_trash_bed=$7

unassigned_bam_gencode=$8
intersecting_bam_gencode=$9
intergenic_bam_gencode=${10}

unassigned_bam_ncbi=${11}
intersecting_bam_ncbi=${12}
intergenic_bam_ncbi=${13}

unassigned_bam_lnc=${14}
intersecting_bam_lnc=${15}
intergenic_bam_lnc=${16}


total_reads=$(samtools flagstat $bam | head -2 | tail -1 | cut -d ' ' -f1)
unassigned_reads=$(samtools flagstat $unassigned_bam_10x | head -2 | tail -1 | cut -d ' ' -f1)
intersecting_10x=$(samtools flagstat $intersecting_bam_10x | head -2 | tail -1 | cut -d ' ' -f1)
intergenic_10x=$(samtools flagstat $intergenic_bam_10x | head -2 | tail -1 | cut -d ' ' -f1)
AT_intergenic_10x=$(samtools flagstat $AT_seq_bam | head -2 | tail -1 | cut -d ' ' -f1)
intergenic_good_lnc=$(awk -F'\'t '{sum += $5} END {print sum}' $intergenic_good_bed)
intergenic_trash_lnc=$(awk -F'\'t '{sum += $5} END {print sum}' $intergenic_trash_bed)
unassigned_gencode=$(samtools flagstat $unassigned_bam_gencode | head -2 | tail -1 | cut -d ' ' -f1)
intersecting_gencode=$(samtools flagstat $intersecting_bam_gencode | head -2 | tail -1 | cut -d ' ' -f1)
intergenic_gencode=$(samtools flagstat $intergenic_bam_gencode | head -2 | tail -1 | cut -d ' ' -f1)
unassigned_ncbi=$(samtools flagstat $unassigned_bam_ncbi | head -2 | tail -1 | cut -d ' ' -f1)
intersecting_ncbi=$(samtools flagstat $intersecting_bam_ncbi | head -2 | tail -1 | cut -d ' ' -f1)
intergenic_ncbi=$(samtools flagstat $intergenic_bam_ncbi | head -2 | tail -1 | cut -d ' ' -f1)
unassigned_lnc=$(samtools flagstat $unassigned_bam_lnc | head -2 | tail -1 | cut -d ' ' -f1)
intersecting_lnc=$(samtools flagstat $intersecting_bam_lnc | head -2 | tail -1 | cut -d ' ' -f1)
intergenic_lnc=$(samtools flagstat $intergenic_bam_lnc | head -2 | tail -1 | cut -d ' ' -f1)

print_percentage() {
    value=$1
    total=$2
    label=$3
    percentage=$(echo "scale=7; $value / $total * 100" | bc)
    echo "$label,$value,$percentage"
}

# Print mapped to 10x reference
echo "Mapped to 10x reference,,"
print_percentage $total_reads $total_reads "Total reads:"
print_percentage $unassigned_reads $total_reads "Unassigned reads:"
print_percentage $intersecting_10x $total_reads "Unassigned intersecting 10x reference:"
print_percentage $AT_intergenic_10x $total_reads "AT rich (intergenic) reads (excluded from below stats):"
print_percentage $intergenic_10x $total_reads "Intergenic:"

# Print mapped to gencode reference
echo "Intergenic mapped to gencode reference,,"
print_percentage $unassigned_gencode $total_reads "Unassigned reads:"
print_percentage $intersecting_gencode $total_reads "Unassigned intersecting gencode reference:"
print_percentage $intergenic_gencode $total_reads "Intergenic:"

# Print mapped to ncbi reference
echo "Intergenic mapped to ncbi reference,,"
print_percentage $unassigned_ncbi $total_reads "Unassigned reads:"
print_percentage $intersecting_ncbi $total_reads "Unassigned intersecting ncbi reference:"
print_percentage $intergenic_ncbi $total_reads "Intergenic:"

# Print mapped to lncRNA reference
echo "Intergenic mapped to lncRNA reference,,"
print_percentage $unassigned_lnc $total_reads "Unassigned reads:"
print_percentage $intersecting_lnc $total_reads "Unassigned intersecting lncpedia reference:"
print_percentage $intergenic_lnc $total_reads "Intergenic:"
print_percentage $intergenic_trash_lnc $total_reads "    Intergenic small clusters:"
print_percentage $intergenic_good_lnc $total_reads "    Intergenic big clusters:"


