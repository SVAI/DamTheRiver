




python2.7 /hive/users/anjowall/bin/vcf2fasta.py /hive/users/anjowall/genomes/hg19/STAR_genome_hg19/hg19.fa bla_snp rarekidneycancer_patient_0%2FF18FTSUSAT0015_HUMaasR%2FWBA%2Fresult_variation%2Fsnp%2FWBA.snp.vcf
python2.7 /hive/users/anjowall/bin/vcf2fasta.py /hive/users/anjowall/genomes/hg19/STAR_genome_hg19/hg19.fa t1_1_snp rarekidneycancer_patient_0%2FF18FTSUSAT0015_HUMaasR%2FT1_1A%2Fresult_variation%2Fsnp%2FT1_1A.snp.vcf


awk '$1 ~ /chr/' gencode.v27lift37.annotation.gtf > gencode.v27lift37.annotation_standard_chrom.gtf
gffread gencode.v27lift37.annotation_standard_chrom.gtf -g bla_snp.fasta -w gencode.v27lift37.annotation_standard_chrom.bla_snp.fa &
gffread gencode.v27lift37.annotation_standard_chrom.gtf -g t1_1_snp.fasta -w gencode.v27lift37.annotation_standard_chrom.t1_1_snp.fa



python2.7 cds_insertion.py --transcript_gtf gencode.v27lift37.annotation_standard_chrom.gtf --transcript_fasta gencode.v27lift37.annotation_standard_chrom.bla_snp.fa --annotation_gtf gencode.v27lift37.annotation.gtf --CCDS --outdir bla_snp_cds_insertion/ &
python2.7 cds_insertion.py --transcript_gtf gencode.v27lift37.annotation_standard_chrom.gtf --transcript_fasta gencode.v27lift37.annotation_standard_chrom.t1_1_snp.fa --annotation_gtf gencode.v27lift37.annotation.gtf --CCDS --outdir t1_1_snp_cds_insertion/ &


python2.7 misc.py --input_tsv bla_snp_cds_insertion/transcript_aa_seq.tsv --outdir bla_snp_cds_insertion/

python2.7 misc.py --input_tsv t1_1_snp_cds_insertion/transcript_aa_seq.tsv --outdir t1_1_snp_cds_insertion/

sort -k1,1 kmers_by_gene.tsv > kmers_by_gene.sorted.tsv
sort -k1,1 kmers_by_gene.tsv > kmers_by_gene.sorted.tsv

join -v 2 -1 1 -2 1 bla_snp_cds_insertion/kmers_by_gene.sorted.tsv t1_1_snp_cds_insertion/kmers_by_gene.sorted.tsv > tumor_kmers_by_gene.tsv

awk '{ print $1}' tumor_kmers_by_gene.tsv | sort | uniq > tumor_kmers.txt









awk '!/#/ {OFS="\t"; print $1,$2 - 1, $2, $4"_"$5"_WBA"}' rarekidneycancer_patient_0%2FF18FTSUSAT0015_HUMaasR%2FWBA%2Fresult_variation%2Fsnp%2FWBA.snp.vcf > wba_var.bed
awk '!/#/ {OFS="\t"; print $1,$2 - 1, $2, $4"_"$5"_T1_1"}' rarekidneycancer_patient_0%2FF18FTSUSAT0015_HUMaasR%2FT1_1A%2Fresult_variation%2Fsnp%2FT1_1A.snp.vcf > t1_1_var.bed
cat t1_1_var.bed wba_var.bed > all_var.bed

python2.7 find_source_var.py --peptides tumor_kmers_by_gene.sep_fixed.tsv --full_proteins t1_1_snp_cds_insertion/t1_1_snp_transcript_aa_seq.tsv --outdir .
bedtools intersect -a peptide_positions.bed -b all_var.bed | head
python2.7 find_source_var.py --peptides tumor_kmers_by_gene.sep_fixed.tsv --full_proteins t1_1_snp_cds_insertion/t1_1_snp_transcript_aa_seq.tsv --outdir .

bedtools intersect -a peptide_positions.bed -b all_var.bed -wa -wb > peptides_ol_var.tsv

bedtools intersect -a peptide_positions.bed -b all_var.bed -wa -wb | awk '{ print $7,$8,$9,$10}' | sort | uniq | awk -F"_" '{OFS="\t"; print $1"_"$2}' | sort | uniq -d > duplicated_vars.bed

bedtools intersect -a <(bedtools intersect -a peptide_positions.bed -b all_var.bed -wa) -b duplicated_vars.bed -v | sort | uniq > peptide_positions_candidates_shared_var_filt.bed

awk -F"," 'NR > 1 {OFS="\t"; print $2,$3}' data_sorted_peptide_scores.csv > data_sorted_peptide_scores.tsv

join <(sort -k1,1 data_sorted_peptide_scores_shared_var_filt.tsv ) <(sort -k1,1 tumor_kmers_by_gene.sep_fixed.tsv) | awk '{OFS="\t"; print $1,$2,$3}' | sort -n -k2,2 | uniq > data_sorted_peptide_scores_shared_var_filt_symbol.tsv