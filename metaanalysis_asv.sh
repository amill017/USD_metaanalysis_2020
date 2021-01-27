#!/bin/bash
##bash analysis.sh forward.fastq reverse.fastq barcode.fastq mapping.txt directory_of_reads 16s/shotgun
##If reads are not demultiplexed, provide the files for forward.fastq reverse.fastq barcode.fastq and put the date in form MMDDYY for directory_of_reads, where demultiplexed files will go
##If reads are demultiplexed, put "NA" for forward.fastq reverse.fastq barcode.fastq and provide the directory_of_reads in the form of the date MMDDYY
##Raw, zipped sequence files will be added to the parent directory for the central server
##Unzip files, remove zipped files
##Differentiate into 16S or shotgun
# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT
##If the files are not demultiplexed, do so below
if [[ $1 = "NA" ]]
then
	mkdir $5\_asv
	cp -r $5/* $5\_asv
else
	mkdir $5\_asv
	cp $1 $2 $3 $5\_asv
	bash scripts/demultiplex.sh $1 $2 $3 $5\_asv $4
fi
cp -r $5\_asv data/downloads
if [[ $6 = "shotgun" ]]
then
	mv $1 $5\_asv
	mv $4 $5\_asv
	bash scripts/shotgun.sh $1 $4 $5\_asv
	exit 1
else [[ $6 = "16s" ]]
	##Provides asv assigned taxonomy/count, pd_alpha, weighted unifrac dm, rep_set.tre
  set +e
	rm $5\_asv*R2*
  set -e
	Rscript scripts/dada2_paired2.r $5\_asv

fi
#Single study analysis
##Raw OTU table
tr ',' '\t' < $5\_asv/otu_table_raw.csv > $5\_asv/otu_table_raw.txt
tr ',' '\t' < $5\_asv/taxonomy_table_raw.csv > $5\_asv/taxonomy_table_raw.txt
sed '1s/^/\#OTU_ID/' $5\_asv/otu_table_raw.txt | sed 's/\"//g' > $5\_asv/otu_table_raw_head.txt
sed '1s/^/\#OTU_ID/' $5\_asv/taxonomy_table_raw.txt | sed 's/\"//g' > $5\_asv/taxonomy_raw_head.txt
/Users/millera25/Desktop/server/new_ref_db/new_pipeline/stats-master/merge -k -e 0 $5\_asv/otu_table_raw_head.txt $5\_asv/taxonomy_raw_head.txt > $5\_asv/otu_table_raw_wtax.txt
sed 's/\#OTU_ID/\#OTU\ ID/g' $5\_asv/otu_table_raw_wtax.txt > $5\_asv/otu_table_raw_wtax2.txt
rm $5\_asv/otu_table_raw_head.txt
rm $5\_asv/taxonomy_raw_head.txt
rm $5\_asv/otu_table_raw_wtax.txt
mv $5\_asv/otu_table_raw_wtax2.txt $5\_asv/otu_table_raw_wtax.txt
sed 's/\"//g' $5\_asv/otu_table_raw_wtax.txt > $5\_asv/otu_table_raw_wtax2.txt
awk 'BEGIN{print""}1' $5\_asv/otu_table_raw_wtax2.txt | sed '1s/^/\#Constructed\ from\ biom\ file/' > $5\_asv/otu_table_raw_wtax3.txt
mv $5\_asv/otu_table_raw_wtax3.txt $5\_asv/otu_table_raw_wtax.txt
rm $5\_asv/otu_table_raw.txt
rm $5\_asv/otu_table_raw_wtax2.txt
rm $5\_asv/otu_table_raw.csv
rm $5\_asv/taxonomy_table_raw.csv
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_raw_wtax.txt > $5\_asv/otu_table_raw_wtax2.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_raw_wtax2.txt > $5\_asv/otu_table_raw_wtax3.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_raw_wtax3.txt > $5\_asv/otu_table_raw_wtax4.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_raw_wtax4.txt > $5\_asv/otu_table_raw_wtax5.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_raw_wtax5.txt > $5\_asv/otu_table_raw_wtax6.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_raw_wtax6.txt > $5\_asv/otu_table_raw_wtax7.txt
perl -pe 's/\x{0}//g' $5\_asv/otu_table_raw_wtax7.txt > $5\_asv/otu_table_raw_wtax8.txt
sed 's/Kingdom\;\ Phylum\;\ Class\;\ Order\;\ Family\;\ Genus\;\ Species/taxonomy/g' $5\_asv/otu_table_raw_wtax8.txt > $5\_asv/otu_table_raw_wtax9.txt
mv $5\_asv/otu_table_raw_wtax9.txt $5\_asv/otu_table_raw_wtax.txt
rm $5\_asv/otu_table_raw_wtax2.txt
rm $5\_asv/otu_table_raw_wtax3.txt
rm $5\_asv/otu_table_raw_wtax4.txt
rm $5\_asv/otu_table_raw_wtax5.txt
rm $5\_asv/otu_table_raw_wtax6.txt
rm $5\_asv/otu_table_raw_wtax7.txt
rm $5\_asv/otu_table_raw_wtax8.txt
biom convert -i $5\_asv/otu_table_raw_wtax.txt -o $5\_asv/otu_table_raw_wtax.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy
##Normal OTU table
tr ',' '\t' < $5\_asv/otu_table_normal.csv > $5\_asv/otu_table_normal.txt
tr ',' '\t' < $5\_asv/taxonomy_table_normal.csv > $5\_asv/taxonomy_table_normal.txt
sed '1s/^/\#OTU_ID/' $5\_asv/otu_table_normal.txt | sed 's/\"//g' > $5\_asv/otu_table_normal_head.txt
sed '1s/^/\#OTU_ID/' $5\_asv/taxonomy_table_normal.txt | sed 's/\"//g' > $5\_asv/taxonomy_normal_head.txt
/Users/millera25/Desktop/server/new_ref_db/new_pipeline/stats-master/merge -k -e 0 $5\_asv/otu_table_normal_head.txt $5\_asv/taxonomy_normal_head.txt > $5\_asv/otu_table_normal_wtax.txt
sed 's/\#OTU_ID/\#OTU\ ID/g' $5\_asv/otu_table_normal_wtax.txt > $5\_asv/otu_table_normal_wtax2.txt
rm $5\_asv/otu_table_normal_head.txt
rm $5\_asv/taxonomy_normal_head.txt
rm $5\_asv/otu_table_normal_wtax.txt
mv $5\_asv/otu_table_normal_wtax2.txt $5\_asv/otu_table_normal_wtax.txt
sed 's/\"//g' $5\_asv/otu_table_normal_wtax.txt > $5\_asv/otu_table_normal_wtax2.txt
awk 'BEGIN{print""}1' $5\_asv/otu_table_normal_wtax2.txt | sed '1s/^/\#Constructed\ from\ biom\ file/' > $5\_asv/otu_table_normal_wtax3.txt
mv $5\_asv/otu_table_normal_wtax3.txt $5\_asv/otu_table_normal_wtax.txt
rm $5\_asv/otu_table_normal.txt
rm $5\_asv/otu_table_normal_wtax2.txt
rm $5\_asv/otu_table_normal.csv
rm $5\_asv/taxonomy_table_normal.csv
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_normal_wtax.txt > $5\_asv/otu_table_normal_wtax2.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_normal_wtax2.txt > $5\_asv/otu_table_normal_wtax3.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_normal_wtax3.txt > $5\_asv/otu_table_normal_wtax4.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_normal_wtax4.txt > $5\_asv/otu_table_normal_wtax5.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_normal_wtax5.txt > $5\_asv/otu_table_normal_wtax6.txt
sed 's/\(.*\)	/\1; /' $5\_asv/otu_table_normal_wtax6.txt > $5\_asv/otu_table_normal_wtax7.txt
perl -pe 's/\x{0}//g' $5\_asv/otu_table_normal_wtax7.txt > $5\_asv/otu_table_normal_wtax8.txt
sed 's/Kingdom\;\ Phylum\;\ Class\;\ Order\;\ Family\;\ Genus\;\ Species/taxonomy/g' $5\_asv/otu_table_normal_wtax8.txt > $5\_asv/otu_table_normal_wtax9.txt
mv $5\_asv/otu_table_normal_wtax9.txt $5\_asv/otu_table_normal_wtax.txt
rm $5\_asv/otu_table_normal_wtax2.txt
rm $5\_asv/otu_table_normal_wtax3.txt
rm $5\_asv/otu_table_normal_wtax4.txt
rm $5\_asv/otu_table_normal_wtax5.txt
rm $5\_asv/otu_table_normal_wtax6.txt
rm $5\_asv/otu_table_normal_wtax7.txt
rm $5\_asv/otu_table_normal_wtax8.txt
biom convert -i $5\_asv/otu_table_normal_wtax.txt -o $5\_asv/otu_table_normal_wtax.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy
####sort alpha & map
awk 'NR<2{print $0;next}{print $0| "sort -r"}' $5\_asv/current_alpha.txt | sed 's/\"//g' > $5\_asv/current_alpha_sorted.txt
sed 's/PD/Sample	PD/g' $5\_asv/current_alpha_sorted.txt | sed 's/\ /	/g' > $5\_asv/current_alpha_sorted2.txt
awk 'NR<2{print $0;next}{print $0| "sort -r"}' $4 | sed 's/\"//g' > $5\_asv/map_sorted.txt
set +e
Rscript scripts/alpha_singlestudy.R $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt
mv data/downloads/alpha_group.pdf data/downloads/$5\_asv/alpha_group_sampletype.pdf
Rscript scripts/alpha_singlestudy_grouponly.R $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool
mv data/downloads/alpha_group.pdf data/downloads/$5\_asv/alpha_group_stool.pdf
Rscript scripts/alpha_singlestudy_grouponly.R $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine
mv data/downloads/alpha_group.pdf data/downloads/$5\_asv/alpha_group_urine.pdf
Rscript scripts/alpha_singlestudy_grouponly.R $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stone
mv data/downloads/alpha_group.pdf data/downloads/$5\_asv/alpha_group_stone.pdf
set -e
#####beta diversity
sed 's/\"//g' $5\_asv/wunifrac.txt > $5\_asv/wunifrac2.txt
sed '1s/^/\ /' $5\_asv/wunifrac2.txt > $5\_asv/wunifrac22.txt
tr ' ' '\t' < $5\_asv/wunifrac22.txt > $5\_asv/wunifrac23.txt
sort $4 > $5\_asv/map_sorted.txt
(head -n 1 $5\_asv/wunifrac23.txt && tail -n +2 $5\_asv/wunifrac23.txt | sort) > $5\_asv/wunifrac32.txt
datamash transpose < $5\_asv/wunifrac32.txt > $5\_asv/wunifrac32t.txt
(head -n 1 $5\_asv/wunifrac32t.txt && tail -n +2 $5\_asv/wunifrac32t.txt | sort) > $5\_asv/wunifrac3.txt
Rscript scripts/beta_singlestudy.R $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt
set +e
Rscript scripts/beta_singlestudy_group_sampletype.R $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt
mv data/downloads/beta_group.pdf data/downloads/$5\_asv/beta_group.pdf
mv data/downloads/beta_group_sampletype.pdf data/downloads/$5\_asv/beta_group_sampletype.pdf
set -e
####Differential abundance
head -1 $4 > $5\_asv/map_header.txt
tail -n +2 $5\_asv/map_sorted.txt > $5\_asv/map_nohead.txt
set +e
grep 'stool' $5\_asv/map_nohead.txt > $5\_asv/stool.txt
grep 'urine' $5\_asv/map_nohead.txt > $5\_asv/urine.txt
grep 'stone' $5\_asv/map_nohead.txt > $5\_asv/stone.txt
cat $5\_asv/map_header.txt $5\_asv/stool.txt > $5\_asv/map_stool.txt
cat $5\_asv/map_header.txt $5\_asv/urine.txt > $5\_asv/map_urine.txt
cat $5\_asv/map_header.txt $5\_asv/stone.txt > $5\_asv/map_stone.txt
filter_samples_from_otu_table.py -i $5\_asv/otu_table_raw_wtax.biom -o $5\_asv/otu_table_raw_wtax_stool.biom --sample_id_fp $5\_asv/map_stool.txt
filter_samples_from_otu_table.py -i $5\_asv/otu_table_raw_wtax.biom -o $5\_asv/otu_table_raw_wtax_urine.biom --sample_id_fp $5\_asv/map_urine.txt
filter_samples_from_otu_table.py -i $5\_asv/otu_table_raw_wtax.biom -o $5\_asv/otu_table_raw_wtax_stone.biom --sample_id_fp $5\_asv/map_stone.txt
awk 'NR<2{print $0;next}{print $0| "sort -r"}' $5\_asv/map_stool.txt | sed 's/\"//g' > $5\_asv/map_stool_sorted.txt
awk 'NR<2{print $0;next}{print $0| "sort -r"}' $5\_asv/map_urine.txt | sed 's/\"//g' > $5\_asv/map_urine_sorted.txt
awk 'NR<2{print $0;next}{print $0| "sort -r"}' $5\_asv/map_stone.txt | sed 's/\"//g' > $5\_asv/map_stone_sorted.txt
##Stool diff_abun use MetagenomeSeq since larger sample size
filter_otus_from_otu_table.py -i $5\_asv/otu_table_raw_wtax_stool.biom -o $5\_asv/otu_table_raw_wtax_stool_present.biom -n 1
filter_otus_from_otu_table.py -i $5\_asv/otu_table_raw_wtax_urine.biom -o $5\_asv/otu_table_raw_wtax_urine_present.biom -n 1
biom convert -i $5\_asv/otu_table_raw_wtax_stool_present.biom -o $5\_asv/otu_table_raw_wtax_stool_present.txt --to-tsv --header-key taxonomy
biom convert -i $5\_asv/otu_table_raw_wtax_urine_present.biom -o $5\_asv/otu_table_raw_wtax_urine_present.txt --to-tsv --header-key taxonomy
sed 's/\;\ /\;_/g' $5\_asv/otu_table_raw_wtax_stool_present.txt | awk '{print $NF}' | tail -n +2 > $5\_asv/stool_taxa.txt
sed 's/\;\ /\;_/g' $5\_asv/otu_table_raw_wtax_urine_present.txt | awk '{print $NF}' | tail -n +2 > $5\_asv/urine_taxa.txt
biom convert -i $5\_asv/otu_table_raw_wtax_stool_present.biom -o $5\_asv/otu_table_raw_stool_present.txt --to-tsv
biom convert -i $5\_asv/otu_table_raw_wtax_urine_present.biom -o $5\_asv/otu_table_raw_urine_present.txt --to-tsv
tail -n +2 $5\_asv/otu_table_raw_stool_present.txt | sed 's/\#OTU\ ID/OTU_ID/g' > $5\_asv/otu_table_raw_stool_present_r.txt
tail -n +2 $5\_asv/otu_table_raw_urine_present.txt | sed 's/\#OTU\ ID/OTU_ID/g' > $5\_asv/otu_table_raw_urine_present_r.txt
###Natively run deseq2 algorithm
Rscript scripts/diff_abun_stool.R $5\_asv/otu_table_raw_stool_present_r.txt $5\_asv/map_stool_sorted.txt $5\_asv/otu_table_raw_urine_present_r.txt $5\_asv/map_urine_sorted.txt $5\_asv/
Rscript scripts/diff_abun_urine.R $5\_asv/otu_table_raw_stool_present_r.txt $5\_asv/map_stool_sorted.txt $5\_asv/otu_table_raw_urine_present_r.txt $5\_asv/map_urine_sorted.txt $5\_asv/
sed 's/\"//g' $5\_asv/diff_abun_stool.txt | sed 's/	baseMean/ASV	baseMean/g' > $5\_asv/diff_abun_stool2.txt
mv $5\_asv/diff_abun_stool2.txt $5\_asv/diff_abun_stool.txt
sed 's/\"//g' $5\_asv/diff_abun_urine.txt | sed 's/	baseMean/ASV	baseMean/g' > $5\_asv/diff_abun_urine2.txt
mv $5\_asv/diff_abun_urine2.txt $5\_asv/diff_abun_urine.txt
awk -v OFS="\t" '$1=$1' $5\_asv/diff_abun_urine.txt | sed '1s/^/\	/' > $5\_asv/diff_abun_urine2.txt
python3 scripts/merge_files.py $5\_asv/diff_abun_urine2.txt $5\_asv/urine_taxa.txt $5\_asv/diff_abun_urine_taxa.txt
cut -d "	" -f2- $5\_asv/diff_abun_urine_taxa.txt | sed 's/\;_/\;\ /g' | sed 's/Unnamed\:\ 0/ASV/g' > $5\_asv/diff_abun_urine_taxa2.txt
mv $5\_asv/diff_abun_urine_taxa2.txt $5\_asv/diff_abun_urine_taxa.txt
awk -v OFS="\t" '$1=$1' $5\_asv/diff_abun_stool.txt | sed '1s/^/\	/' > $5\_asv/diff_abun_stool2.txt
python3 scripts/merge_files.py $5\_asv/diff_abun_stool2.txt $5\_asv/stool_taxa.txt $5\_asv/diff_abun_stool_taxa.txt
cut -d "	" -f2- $5\_asv/diff_abun_stool_taxa.txt | sed 's/\;_/\;\ /g' | sed 's/Unnamed\:\ 0/ASV/g' > $5\_asv/diff_abun_stool_taxa2.txt
mv $5\_asv/diff_abun_stool_taxa2.txt $5\_asv/diff_abun_stool_taxa.txt
cp $5\_asv/diff_abun_stool_taxa.txt data/downloads/$5\_asv
#####Urine diff_abun use DESeq2 since smaller sample size - can switch when larger sample size is possible
cp $5\_asv/diff_abun_urine_taxa.txt data/downloads/$5\_asv
############Format diff abun files for heatmap
##############Get significant urine
awk '{ if ($6 < 0.05) { print }}' $5\_asv/diff_abun_urine_taxa.txt > $5\_asv/diff_abun_urine_signif.txt
awk '{ if ($3 < 0) { print }}' $5\_asv/diff_abun_urine_signif.txt > $5\_asv/diff_abun_urine_signif_control.txt
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_urine_signif_control.txt > $5\_asv/diff_abun_urine_signif_control2.txt
mv $5\_asv/diff_abun_urine_signif_control2.txt $5\_asv/diff_abun_urine_signif_control.txt
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_urine_signif_control.txt | sed 's/NA\;//g' | awk -v OFS='\t' '{print $NF}' | sed 's/\;//g' | awk '! /[0-9]+$/' > $5\_asv/diff_abun_urine_signif_control_genus.txt
sed 's/\[//g' $5\_asv/diff_abun_urine_signif_control_genus.txt | sed 's/\]//g' | sort | uniq -c > $5\_asv/diff_abun_urine_signif_control_genus_formatted.txt
awk 'BEGIN{print "control_urine Taxa"}; {print}' $5\_asv/diff_abun_urine_signif_control_genus_formatted.txt > $5\_asv/diff_abun_urine_forheatmap_control.txt
awk '{$1=$1};1' $5\_asv/diff_abun_urine_forheatmap_control.txt > $5\_asv/diff_abun_urine_forheatmap_control2.txt
tr ' ' \\t < $5\_asv/diff_abun_urine_forheatmap_control2.txt > $5\_asv/diff_abun_urine_forheatmap_control3.txt
rm $5\_asv/diff_abun_urine_signif_control_genus.txt $5\_asv/diff_abun_urine_signif_control_genus_formatted.txt $5\_asv/diff_abun_urine_forheatmap_control.txt $5\_asv/diff_abun_urine_forheatmap_control2.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_asv/diff_abun_urine_forheatmap_control3.txt | awk '!/Bacteria/' > $5\_asv/diff_abun_urine_forheatmap_control4.txt
awk '{ if ($3 > 0) { print }}' $5\_asv/diff_abun_urine_signif.txt > $5\_asv/diff_abun_urine_signif_USD.txt
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_urine_signif_USD.txt > $5\_asv/diff_abun_urine_signif_USD2.txt
mv $5\_asv/diff_abun_urine_signif_USD2.txt $5\_asv/diff_abun_urine_signif_USD.txt
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_urine_signif_USD.txt | sed 's/NA\;//g' | awk -v OFS='\t' '{print $NF}' | sed 's/\;//g' | awk '! /[0-9]+$/' > $5\_asv/diff_abun_urine_signif_USD_genus.txt
sed 's/\[//g' $5\_asv/diff_abun_urine_signif_USD_genus.txt | sed 's/\]//g' | sort | uniq -c > $5\_asv/diff_abun_urine_signif_USD_genus_formatted.txt
awk 'BEGIN{print "USD_urine Taxa"}; {print}' $5\_asv/diff_abun_urine_signif_USD_genus_formatted.txt > $5\_asv/diff_abun_urine_forheatmap_USD.txt
awk '{$1=$1};1' $5\_asv/diff_abun_urine_forheatmap_USD.txt > $5\_asv/diff_abun_urine_forheatmap_USD2.txt
tr ' ' \\t < $5\_asv/diff_abun_urine_forheatmap_USD2.txt > $5\_asv/diff_abun_urine_forheatmap_USD3.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_asv/diff_abun_urine_forheatmap_USD3.txt | awk '!/Bacteria/' > $5\_asv/diff_abun_urine_forheatmap_USD4.txt
####Get significant stool
awk '{ if ($7 < 0.05) { print }}' $5\_asv/diff_abun_stool_taxa.txt > $5\_asv/diff_abun_stool_signif.txt
awk '{ if ($3 < 0) { print }}' $5\_asv/diff_abun_stool_signif.txt > $5\_asv/diff_abun_stool_signif_control.txt
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_stool_signif_control.txt > $5\_asv/diff_abun_stool_signif_control2.txt
mv $5\_asv/diff_abun_stool_signif_control2.txt $5\_asv/diff_abun_stool_signif_control.txt
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_stool_signif_control.txt | sed 's/NA\;//g' | awk -v OFS='\t' '{print $NF}' | sed 's/\;//g' | awk '! /[0-9]+$/' > $5\_asv/diff_abun_stool_signif_control_genus.txt
sed 's/\[//g' $5\_asv/diff_abun_stool_signif_control_genus.txt | sed 's/\]//g' | sort | uniq -c > $5\_asv/diff_abun_stool_signif_control_genus_formatted.txt
awk 'BEGIN{print "control_stool Taxa"}; {print}' $5\_asv/diff_abun_stool_signif_control_genus_formatted.txt > $5\_asv/diff_abun_stool_forheatmap_control.txt
awk '{$1=$1};1' $5\_asv/diff_abun_stool_forheatmap_control.txt > $5\_asv/diff_abun_stool_forheatmap_control2.txt
tr ' ' \\t < $5\_asv/diff_abun_stool_forheatmap_control2.txt > $5\_asv/diff_abun_stool_forheatmap_control3.txt
rm $5\_asv/diff_abun_stool_signif_control_genus.txt $5\_asv/diff_abun_stool_signif_control_genus_formatted.txt $5\_asv/diff_abun_stool_forheatmap_control.txt $5\_asv/diff_abun_stool_forheatmap_control2.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_asv/diff_abun_stool_forheatmap_control3.txt | awk '!/Bacteria/' > $5\_asv/diff_abun_stool_forheatmap_control4.txt
awk '{ if ($3 > 0) { print }}' $5\_asv/diff_abun_stool_signif.txt > $5\_asv/diff_abun_stool_signif_USD.txt
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_stool_signif_USD.txt > $5\_asv/diff_abun_stool_signif_USD2.txt
mv $5\_asv/diff_abun_stool_signif_USD2.txt $5\_asv/diff_abun_stool_signif_USD.txt
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_stool_signif_USD.txt | sed 's/NA\;//g' | awk -v OFS='\t' '{print $NF}' | sed 's/\;//g' | awk '! /[0-9]+$/' > $5\_asv/diff_abun_stool_signif_USD_genus.txt
sed 's/\[//g' $5\_asv/diff_abun_stool_signif_USD_genus.txt | sed 's/\]//g' | sort | uniq -c > $5\_asv/diff_abun_stool_signif_USD_genus_formatted.txt
awk 'BEGIN{print "USD_stool Taxa"}; {print}' $5\_asv/diff_abun_stool_signif_USD_genus_formatted.txt > $5\_asv/diff_abun_stool_forheatmap_USD.txt
awk '{$1=$1};1' $5\_asv/diff_abun_stool_forheatmap_USD.txt > $5\_asv/diff_abun_stool_forheatmap_USD2.txt
tr ' ' \\t < $5\_asv/diff_abun_stool_forheatmap_USD2.txt > $5\_asv/diff_abun_stool_forheatmap_USD3.txt
rm $5\_asv/diff_abun_stool_signif_USD_genus.txt $5\_asv/diff_abun_stool_signif_USD_genus_formatted.txt $5\_asv/diff_abun_stool_forheatmap_USD.txt $5\_asv/diff_abun_stool_forheatmap_USD2.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_asv/diff_abun_stool_forheatmap_USD3.txt | awk '!/Bacteria/' > $5\_asv/diff_abun_stool_forheatmap_USD4.txt
##Get counts per lowest taxonomy
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_urine_taxa.txt > $5\_asv/diff_abun_urine_genusfamily_count.txt
sed 's/NA\;//g' $5\_asv/diff_abun_urine_genusfamily_count.txt | awk -v OFS='\t' '{print $NF}' | sed 's/\;//g' | awk '! /[0-9]+$/' | awk '!/pvalue/' > $5\_asv/diff_abun_urine_genusfamily_count2.txt
sed 's/\[//g' $5\_asv/diff_abun_urine_genusfamily_count2.txt | sed 's/\]//g' | sort | uniq -c > $5\_asv/diff_abun_urine_genusfamily_count3.txt
awk 'BEGIN{print "Count Taxa"}; {print}' $5\_asv/diff_abun_urine_genusfamily_count3.txt > $5\_asv/diff_abun_urine_genusfamily_count4.txt
awk '{$1=$1};1' $5\_asv/diff_abun_urine_genusfamily_count4.txt > $5\_asv/diff_abun_urine_genusfamily_count5.txt
tr ' ' \\t < $5\_asv/diff_abun_urine_genusfamily_count5.txt > $5\_asv/diff_abun_urine_genusfamily_count6.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_asv/diff_abun_urine_genusfamily_count6.txt > $5\_asv/diff_abun_urine_genusfamily_count7.txt
awk -v OFS='\t' 'NF{NF-=1};1' $5\_asv/diff_abun_stool_taxa.txt > $5\_asv/diff_abun_stool_genusfamily_count.txt
sed 's/NA\;//g' $5\_asv/diff_abun_stool_genusfamily_count.txt | awk -v OFS='\t' '{print $NF}' | sed 's/\;//g' | awk '! /[0-9]+$/' | awk '!/pvalue/' > $5\_asv/diff_abun_stool_genusfamily_count2.txt
sed 's/\[//g' $5\_asv/diff_abun_stool_genusfamily_count2.txt | sed 's/\]//g' | sort | uniq -c > $5\_asv/diff_abun_stool_genusfamily_count3.txt
awk 'BEGIN{print "Count Taxa"}; {print}' $5\_asv/diff_abun_stool_genusfamily_count3.txt > $5\_asv/diff_abun_stool_genusfamily_count4.txt
awk '{$1=$1};1' $5\_asv/diff_abun_stool_genusfamily_count4.txt > $5\_asv/diff_abun_stool_genusfamily_count5.txt
tr ' ' \\t < $5\_asv/diff_abun_stool_genusfamily_count5.txt > $5\_asv/diff_abun_stool_genusfamily_count6.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_asv/diff_abun_stool_genusfamily_count6.txt > $5\_asv/diff_abun_stool_genusfamily_count7.txt
rm $5\_asv/diff_abun_urine_genusfamily_count.txt $5\_asv/diff_abun_urine_genusfamily_count2.txt $5\_asv/diff_abun_urine_genusfamily_count3.txt $5\_asv/diff_abun_urine_genusfamily_count4.txt $5\_asv/diff_abun_urine_genusfamily_count5.txt $5\_asv/diff_abun_urine_genusfamily_count6.txt $5\_asv/diff_abun_stool_genusfamily_count.txt $5\_asv/diff_abun_stool_genusfamily_count2.txt $5\_asv/diff_abun_stool_genusfamily_count3.txt $5\_asv/diff_abun_stool_genusfamily_count4.txt $5\_asv/diff_abun_stool_genusfamily_count5.txt $5\_asv/diff_abun_stool_genusfamily_count6.txt
##Scale diff abun analyses
stats-master/merge -k -e 0 $5\_asv/diff_abun_stool_forheatmap_USD4.txt $5\_asv/diff_abun_stool_genusfamily_count7.txt > $5\_asv/diff_abun_stool_forheatmap_USD_toscale.txt
stats-master/merge -k -e 0 $5\_asv/diff_abun_stool_forheatmap_control4.txt $5\_asv/diff_abun_stool_genusfamily_count7.txt > $5\_asv/diff_abun_stool_forheatmap_control_toscale.txt
stats-master/merge -k -e 0 $5\_asv/diff_abun_urine_forheatmap_USD4.txt $5\_asv/diff_abun_urine_genusfamily_count7.txt > $5\_asv/diff_abun_urine_forheatmap_USD_toscale.txt
stats-master/merge -k -e 0 $5\_asv/diff_abun_urine_forheatmap_control4.txt $5\_asv/diff_abun_urine_genusfamily_count7.txt > $5\_asv/diff_abun_urine_forheatmap_control_toscale.txt
Rscript scripts/scale_dysbiosis.R $5\_asv/diff_abun_urine_forheatmap_control_toscale.txt $5\_asv/diff_abun_urine_forheatmap_USD_toscale.txt $5\_asv/cont_urin_dev.txt $5\_asv/usd_urin_dev.txt $5\_asv/diff_abun_stool_forheatmap_control_toscale.txt $5\_asv/cont_stool_dev.txt $5\_asv/diff_abun_stool_forheatmap_USD_toscale.txt $5\_asv/usd_stool_dev.txt
awk -v OFS='\t' '{print $2,$3,$4,$5,$6}' $5\_asv/cont_urin_dev.txt > $5\_asv/cont_urin_dev2.txt
sed 's/\"//g' $5\_asv/cont_urin_dev2.txt > $5\_asv/cont_urin_dev3.txt
awk -v OFS='\t' '{print $2,$3,$4,$5,$6}' $5\_asv/usd_urin_dev.txt > $5\_asv/usd_urin_dev2.txt
sed 's/\"//g' $5\_asv/usd_urin_dev2.txt > $5\_asv/usd_urin_dev3.txt
awk -v OFS='\t' '{print $2,$3,$4,$5,$6}' $5\_asv/cont_stool_dev.txt > $5\_asv/cont_stool_dev2.txt
sed 's/\"//g' $5\_asv/cont_stool_dev2.txt > $5\_asv/cont_stool_dev3.txt
awk -v OFS='\t' '{print $2,$3,$4,$5,$6}' $5\_asv/usd_stool_dev.txt > $5\_asv/usd_stool_dev2.txt
sed 's/\"//g' $5\_asv/usd_stool_dev2.txt > $5\_asv/usd_stool_dev3.txt
awk -v OFS='\t' '{print $1,$5}' $5\_asv/usd_stool_dev3.txt > $5\_asv/usd_stool_dev4.txt
awk -v OFS='\t' '{print $1,$5}' $5\_asv/cont_stool_dev3.txt > $5\_asv/cont_stool_dev4.txt
awk -v OFS='\t' '{print $1,$5}' $5\_asv/usd_urin_dev3.txt > $5\_asv/usd_urin_dev4.txt
awk -v OFS='\t' '{print $1,$5}' $5\_asv/cont_urin_dev3.txt > $5\_asv/cont_urin_dev4.txt
##Merge scaled files
stats-master/merge -k -e 0 $5\_asv/cont_urin_dev4.txt $5\_asv/usd_urin_dev4.txt > $5\_asv/scaled_difference_urine.txt
stats-master/merge -k -e 0 $5\_asv/cont_stool_dev4.txt $5\_asv/usd_stool_dev4.txt > $5\_asv/scaled_difference_stool.txt
sed '1s/.*/Taxon	Control_urine	USD_urine/g'  $5\_asv/scaled_difference_urine.txt | grep -v "taxonomy" > $5\_asv/scaled_difference_urine2.txt
sed '1s/.*/Taxon	Control_stool	USD_stool/g'  $5\_asv/scaled_difference_stool.txt | grep -v "taxonomy" > $5\_asv/scaled_difference_stool2.txt
rm $5\_asv/*dev.txt
rm $5\_asv/*dev2.txt
rm $5\_asv/*dev3.txt
rm $5\_asv/*dev4.txt
rm $5\_asv/scaled_difference_urine.txt
rm $5\_asv/scaled_difference_stool.txt
awk '{ if ($2 > 30 || $3 > 30) { print }}' $5\_asv/scaled_difference_urine2.txt > $5\_asv/scaled_difference_urine3.txt
awk '{ if ($2 > 20 || $3 > 20) { print }}' $5\_asv/scaled_difference_stool2.txt > $5\_asv/scaled_difference_stool3.txt
awk -v OFS='\t' 'NR==1; NR > 1 {print $1, $2, -$3}' $5\_asv/scaled_difference_urine3.txt > $5\_asv/scaled_difference_urine4.txt
awk -v OFS='\t' 'NR==1; NR > 1 {print $1, $2, -$3}' $5\_asv/scaled_difference_stool3.txt > $5\_asv/scaled_difference_stool4.txt
awk -v OFS='\t' 'NR==1; NR > 1 {print $1, $2, $3, $2+$3}' $5\_asv/scaled_difference_urine4.txt > $5\_asv/scaled_difference_urine5.txt
awk -v OFS='\t' 'NR==1; NR > 1 {print $1, $2, $3, $2+$3}' $5\_asv/scaled_difference_stool4.txt > $5\_asv/scaled_difference_stool5.txt
awk 'NR==1; NR > 1 {print $0 | "sort -n -r -k4"}' $5\_asv/scaled_difference_urine5.txt > $5\_asv/scaled_difference_urine5_sorted.txt
awk 'NR==1; NR > 1 {print $0 | "sort -n -r -k4"}' $5\_asv/scaled_difference_stool5.txt > $5\_asv/scaled_difference_stool5_sorted.txt
sed '1s/.*/Taxon	Control_urine	USD_urine	Descend/g' $5\_asv/scaled_difference_urine5_sorted.txt | grep -v "taxonomy" > $5\_asv/scaled_difference_urine5_sorted2.txt
sed '1s/.*/Taxon	Control_stool	USD_stool	Descend/g' $5\_asv/scaled_difference_stool5_sorted.txt | grep -v "taxonomy" > $5\_asv/scaled_difference_stool5_sorted2.txt
awk '{print $1, $2, $3}' $5\_asv/scaled_difference_urine5_sorted2.txt > $5\_asv/scaled_difference_urine5_sorted3.txt
awk '{print $1, $2, $3}' $5\_asv/scaled_difference_stool5_sorted2.txt > $5\_asv/scaled_difference_stool5_sorted3.txt
##Extract stone uclust_open, get average relative abundance, and plot heatmap
filter_otus_from_otu_table.py -i $5\_asv/otu_table_raw_wtax_stone.biom -o $5\_asv/otu_table_raw_wtax_stone_present.biom -n 2
biom convert -i $5\_asv/otu_table_raw_wtax_stone_present.biom -o $5\_asv/otu_table_raw_wtax_stone_present.txt --to-tsv --header-key taxonomy
tail -n +2 $5\_asv/otu_table_raw_wtax_stone_present.txt  > $5\_asv/otu_table_over_10_filtered_stone_over2_2.txt
sed 's/\#OTU\ ID/OTU_ID/g' $5\_asv/otu_table_over_10_filtered_stone_over2_2.txt > $5\_asv/otu_table_over_10_filtered_stone_over2_3.txt
sed 's/\ /_/g' $5\_asv/otu_table_over_10_filtered_stone_over2_3.txt > $5\_asv/otu_table_over_10_filtered_stone_over2_4.txt
Rscript scripts/rel_abun_stone.R $5\_asv/otu_table_over_10_filtered_stone_over2_4.txt $5\_asv/stone_relabun.txt
awk -v OFS='\t' '{print $NF}' $5\_asv/otu_table_over_10_filtered_stone_over2_4.txt > $5\_asv/stone_taxonomy.txt
sed 's/\;_/\;\ /g' $5\_asv/stone_taxonomy.txt | awk 'NF{NF-=1};1' | sed 's/NA\;//g' | sed 's/NA//g' | awk -v OFS='\t' '{print $NF}' > $5\_asv/stone_taxonomy2.txt
awk -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; print sum}' $5\_asv/stone_relabun.txt > $5\_asv/stone_relabun2.txt
sed '1s/.*/Average/g' $5\_asv/stone_relabun2.txt > $5\_asv/stone_relabun3.txt
paste $5\_asv/stone_relabun3.txt $5\_asv/stone_taxonomy2.txt > $5\_asv/stone_relabun_taxonomy.txt
awk -v OFS='\t' '{print $2,$1}' $5\_asv/stone_relabun_taxonomy.txt | sed 's/\;//g' > $5\_asv/stone_relabun_heatmap.txt
awk '{ if ($2 > 0.01) { print }}' $5\_asv/stone_relabun_heatmap.txt > $5\_asv/stone_relabun_heatmap2.txt
sort -n -r -k2  $5\_asv/stone_relabun_heatmap2.txt | awk '!/Bacteria/' | awk 'BEGIN{print "Taxon    Stone"}1' > $5\_asv/stone_relabun_heatmap3.txt
set +e
Rscript scripts/heatmap_singlestudy.R $5\_asv/scaled_difference_urine5_sorted3.txt $5\_asv/scaled_difference_stool5_sorted3.txt $5\_asv/stone_relabun_heatmap3.txt
Rscript scripts/alpha_by_choice_singlestudy_stonecomp.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Stonecomp
Rscript scripts/alpha_by_choice_singlestudy_stonecomp.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Stonecomp
Rscript scripts/alpha_by_choice_singlestudy_stonecomp.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stone Stonecomp
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt All Group Sample_type
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group Sample_type ######one way
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group Study_location
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group age_group
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group Sex
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group Weight_group
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group antibiotic_12m
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group antibiotic_30d
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group gout
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group diabetes
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group HTN
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group diet_type
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group water
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group dessert
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group meat
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group fruit
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group veggie
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt stool Group bread
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group Ston_comp  #####one way
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group Study_location
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group age_group
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group Sex
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group Weight_group
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group antibiotic_12m
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group antibiotic_30d
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group gout
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group diabetes
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group HTN
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group diet
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group water
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group dessert
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group meat
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group fruit
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group veggie
Rscript scripts/alpha_by_choice_singlestudy.r $5\_asv/current_alpha_sorted2.txt $5\_asv/map_sorted.txt urine Group bread
set +e
Rscript scripts/beta_by_choice_singlestudy_stonecomp.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt stone StoneComp
Rscript scripts/beta_by_choice_singlestudy_stonecomp.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt stool StoneComp
Rscript scripts/beta_by_choice_singlestudy_stonecomp.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt urine StoneComp
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  all Group Sample_type
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group Study_location
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group age_group
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group Sex
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group Weight_group
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group antibiotic_12m
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group antibiotic_30d
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group gout
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group diabetes
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group HTN
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group diet
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group water
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group dessert
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group meat
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group fruit
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group veggie
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  stool Group bread
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group age_group
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group Sex
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group Weight_group
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group antibiotic_12m
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group antibiotic_30d
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group gout
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group diabetes
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group HTN
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group diet
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group water
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group dessert
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group meat
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group fruit
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group veggie
Rscript scripts/beta_by_choice_singlestudy.r $5\_asv/wunifrac3.txt $5\_asv/map_sorted.txt  urine Group bread
set -e
cp $5\_asv/map_sorted.txt $5\_asv/otu_table_raw_wtax.biom $5\_asv/rep_set.tre data/downloads/$5\_asv
mv data/downloads/*.pdf data/downloads/$5\_asv
gzip data/downloads/$5\_asv/*.fastq
