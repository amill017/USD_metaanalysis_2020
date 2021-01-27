#!/bin/bash
##bash analysis.sh forward.fastq reverse.fastq barcode.fastq mapping.txt directory_of_reads 16s/shotgun
##Raw, zipped sequence files will be added to the parent directory for the central server
##Unzip files, remove zipped files
##Differentiate into 16S or shotgun
# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT
if [[ $6 = "shotgun" ]]
then
	mkdir $5\_otu
	cp $1 $5\_otu
	cp $4 $5\_otu
	bash scripts/shotgun.sh $1 $4 $5\_otu
	exit 1
elif [[ $6 = "16s" ]] && [[ $2 = "NA" ]]
then
	mkdir $5\_otu
	multiple_split_libraries_fastq.py -i $5 -o $5\_otu --mapping_indicator map -p parameter.txt

else
#join pairs
	mkdir $5\_otu
	join_paired_ends.py -f $1 -r $2 -b $3 -o joined
####### #Demultiplex
	split_libraries_fastq.py -i joined/fastqjoin.join.fastq -b joined/fastqjoin.join_barcodes.fastq --rev_comp_mapping_barcodes --rev_comp_barcode --barcode_type 12 -m $4 -o $5\_otu
fi
##########Single study analysis
##############open reference pick otus
mv $5\_otu/seqs.fna $5\_otu/"$5".fna
pick_open_reference_otus.py -i $5\_otu/"$5".fna -o $5\_otu/uclust_open -r new_refseqs_asv.fna -s 1 -n $5 --min_otu_size 10
mkdir data/downloads/$5\_otu
cp $5\_otu/"$5".fna data/downloads/$5\_otu
cp $5\_otu/"$5".fna data/downloads
rm $5\_otu/"$5".fna
#########Filter mitochondria and chloroplasts -verfied
filter_taxa_from_otu_table.py -i $5\_otu/uclust_open/otu_table_mc10_w_tax_no_pynast_failures.biom -o $5\_otu/uclust_open/otu_table_over_10_filtered.biom -n c__chloroplast,f__mitochondria
##############Normalize data -verfied
normalize_table.py -i $5\_otu/uclust_open/otu_table_over_10_filtered.biom -o $5\_otu/uclust_open/otu_table_over_10_filtered_DESeq2.biom -a DESeq2 -z
##########Add taxonomy back to file -verfied
biom convert -i $5\_otu/uclust_open/otu_table_over_10_filtered.biom -o $5\_otu/uclust_open/"$5".txt --to-tsv --header-key taxonomy
mv $5\_otu/uclust_open/"$5".txt $5\_otu/uclust_open/current.txt
sed 's/\;\ /\;\_/g' $5\_otu/uclust_open/current.txt | awk '{print $NF}' | sed 's/\;\_/\;\ /g' | sed 's/file//g' | sed 's/taxonomy//g' > $5\_otu/uclust_open/taxonomy.txt
biom convert -i $5\_otu/uclust_open/otu_table_over_10_filtered_DESeq2.biom -o $5\_otu/uclust_open/"$5"_deseq2.txt --to-tsv --header-key taxonomy
mv $5\_otu/uclust_open/"$5"_deseq2.txt $5\_otu/uclust_open/current_deseq2.txt
paste $5\_otu/uclust_open/current_deseq2.txt $5\_otu/uclust_open/taxonomy.txt > $5\_otu/uclust_open/current_deseq2_taxonomy.txt
biom convert -i $5\_otu/uclust_open/current_deseq2_taxonomy.txt -o $5\_otu/uclust_open/"$5"_deseq2_taxonomy.biom --table-type="OTU table" --to-json --process-obs-metadata taxonomy
mv $5\_otu/uclust_open/"$5"_deseq2_taxonomy.biom $5\_otu/uclust_open/current_deseq2_taxonomy.biom
#####alpha diversity
alpha_diversity.py -i $5\_otu/uclust_open/current_deseq2_taxonomy.biom -o $5\_otu/uclust_open/current_alpha.txt -m PD_whole_tree -t $5\_otu/uclust_open/rep_set.tre
############sort alpha
awk 'NR<2{print $0;next}{print $0| "sort"}' $5\_otu/uclust_open/current_alpha.txt > $5\_otu/uclust_open/current_alpha_sorted.txt
sort $4 > $5\_otu/uclust_open/map_sorted.txt
sed 's/PD_whole_tree/PD/g' $5\_otu/uclust_open/current_alpha_sorted.txt > $5\_otu/uclust_open/current_alpha_sorted2.txt
set +e
Rscript scripts/alpha_singlestudy_dornbier.R $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt
mv data/downloads/alpha_group.pdf data/downloads/$5\_otu/alpha_group_sampletype.pdf
Rscript scripts/alpha_singlestudy_grouponly.R $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool
mv data/downloads/alpha_group.pdf data/downloads/$5\_otu/alpha_group_stool.pdf
Rscript scripts/alpha_singlestudy_grouponly.R $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine
mv data/downloads/alpha_group.pdf data/downloads/$5\_otu/alpha_group_urine.pdf
Rscript scripts/alpha_singlestudy_grouponly.R $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stone
mv data/downloads/alpha_group.pdf data/downloads/$5\_otu/alpha_group_stone.pdf
set -e
#############beta diversity
beta_diversity_through_plots.py -i $5\_otu/uclust_open/current_deseq2_taxonomy.biom -o $5\_otu/uclust_open/beta -m $5\_otu/uclust_open/map_sorted.txt -t $5\_otu/uclust_open/rep_set.tre
(head -n 1 $5\_otu/uclust_open/beta/weighted_unifrac_dm.txt && tail -n +2 $5\_otu/uclust_open/beta/weighted_unifrac_dm.txt | sort) > $5\_otu/uclust_open/beta/weighted_unifrac_dm2.txt
datamash transpose < $5\_otu/uclust_open/beta/weighted_unifrac_dm2.txt > $5\_otu/uclust_open/beta/weighted_unifrac_dm2t.txt
(head -n 1 $5\_otu/uclust_open/beta/weighted_unifrac_dm2t.txt && tail -n +2 $5\_otu/uclust_open/beta/weighted_unifrac_dm2t.txt | sort) > $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt
Rscript scripts/beta_singlestudy_dornbier.R $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt
set +e
Rscript scripts/beta_singlestudy_group_sampletype.R $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt
mv data/downloads/beta_group.pdf data/downloads/$5\_otu/beta_group.pdf
mv data/downloads/beta_group_sampletype.pdf data/downloads/$5\_otu/beta_group_sampletype.pdf
#######################Differential abundance
set +e
head -1 $5\_otu/uclust_open/map_sorted.txt > $5\_otu/uclust_open/map_header.txt
tail -n +2 $5\_otu/uclust_open/map_sorted.txt > $5\_otu/uclust_open/map_nohead.txt
grep 'stool' $5\_otu/uclust_open/map_nohead.txt > $5\_otu/uclust_open/stool.txt
grep 'urine' $5\_otu/uclust_open/map_nohead.txt > $5\_otu/uclust_open/urine.txt
grep 'stone' $5\_otu/uclust_open/map_nohead.txt > $5\_otu/uclust_open/stone.txt
cat $5\_otu/uclust_open/map_header.txt $5\_otu/uclust_open/stool.txt > $5\_otu/uclust_open/map_stool.txt
cat $5\_otu/uclust_open/map_header.txt $5\_otu/uclust_open/urine.txt > $5\_otu/uclust_open/map_urine.txt
cat $5\_otu/uclust_open/map_header.txt $5\_otu/uclust_open/stone.txt > $5\_otu/uclust_open/map_stone.txt
filter_samples_from_otu_table.py -i $5\_otu/uclust_open/otu_table_over_10_filtered.biom -o $5\_otu/uclust_open/otu_table_over_10_filtered_stool.biom --sample_id_fp $5\_otu/uclust_open/map_stool.txt
filter_samples_from_otu_table.py -i $5\_otu/uclust_open/otu_table_over_10_filtered.biom -o $5\_otu/uclust_open/otu_table_over_10_filtered_urine.biom --sample_id_fp $5\_otu/uclust_open/map_urine.txt
filter_samples_from_otu_table.py -i $5\_otu/uclust_open/otu_table_over_10_filtered.biom -o $5\_otu/uclust_open/otu_table_over_10_filtered_stone.biom --sample_id_fp $5\_otu/uclust_open/map_stone.txt
##Stool diff_abun use MetagenomeSeq since larger sample size
filter_otus_from_otu_table.py -i $5\_otu/uclust_open/otu_table_over_10_filtered_stool.biom -o $5\_otu/uclust_open/otu_table_over_10_filtered_stool_present.biom -n 1
filter_otus_from_otu_table.py -i $5\_otu/uclust_open/otu_table_over_10_filtered_urine.biom -o $5\_otu/uclust_open/otu_table_over_10_filtered_urine_present.biom -n 1
differential_abundance.py -i $5\_otu/uclust_open/otu_table_over_10_filtered_stool_present.biom -o $5\_otu/uclust_open/diff_abun_stool -m $5\_otu/uclust_open/map_stool.txt -a DESeq2_nbinom -d -c Group -x Control -y USD
cp $5\_otu/uclust_open/diff_abun_stool data/downloads/$5\_otu
##Urine diff_abun use DESeq2 since smaller sample size - can switch when larger sample size is possible
differential_abundance.py -i $5\_otu/uclust_open/otu_table_over_10_filtered_urine_present.biom -o $5\_otu/uclust_open/diff_abun_urine -a DESeq2_nbinom -m $5\_otu/uclust_open/map_urine.txt -c Group -x Control -y USD -d
cp $5\_otu/uclust_open/diff_abun_urine data/downloads/$5\_otu
#############Format diff abun files for heatmap
######################Get significant urine
awk '{ if ($7 < 0.05) { print }}' $5\_otu/uclust_open/diff_abun_urine > $5\_otu/uclust_open/diff_abun_urine_signif.txt
awk '{ if ($3 > 0) { print }}' $5\_otu/uclust_open/diff_abun_urine_signif.txt > $5\_otu/uclust_open/diff_abun_urine_signif_control.txt
sed 's/k__//g' $5\_otu/uclust_open/diff_abun_urine_signif_control.txt | sed 's/p__//g' | sed 's/c__//g' | sed 's/o__//g' | sed 's/f__//g' | sed 's/g__//g' | sed 's/s__.*//g' | sed 's/\;//g' | sed 's/NA//g' | awk -v OFS='\t' '{print $NF}' > $5\_otu/uclust_open/diff_abun_urine_signif_control_genus.txt
sed 's/\[//g' $5\_otu/uclust_open/diff_abun_urine_signif_control_genus.txt | sed 's/\]//g' | sort | uniq -c > $5\_otu/uclust_open/diff_abun_urine_signif_control_genus_formatted.txt
awk 'BEGIN{print "control_urine Taxa"}; {print}' $5\_otu/uclust_open/diff_abun_urine_signif_control_genus_formatted.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_control.txt
awk '{$1=$1};1' $5\_otu/uclust_open/diff_abun_urine_forheatmap_control.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_control2.txt
tr ' ' \\t < $5\_otu/uclust_open/diff_abun_urine_forheatmap_control2.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_control3.txt
rm $5\_otu/uclust_open/diff_abun_urine_signif_control_genus.txt $5\_otu/uclust_open/diff_abun_urine_signif_control_genus_formatted.txt $5\_otu/uclust_open/diff_abun_urine_forheatmap_control.txt $5\_otu/uclust_open/diff_abun_urine_forheatmap_control2.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_otu/uclust_open/diff_abun_urine_forheatmap_control3.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_control4.txt
awk '{ if ($3 < 0) { print }}' $5\_otu/uclust_open/diff_abun_urine_signif.txt > $5\_otu/uclust_open/diff_abun_urine_signif_USD.txt
sed 's/k__//g' $5\_otu/uclust_open/diff_abun_urine_signif_USD.txt | sed 's/p__//g' | sed 's/c__//g' | sed 's/o__//g' | sed 's/f__//g' | sed 's/g__//g' | sed 's/s__.*//g' | sed 's/\;//g' | sed 's/NA//g' | awk -v OFS='\t' '{print $NF}' > $5\_otu/uclust_open/diff_abun_urine_signif_USD_genus.txt
sed 's/\[//g' $5\_otu/uclust_open/diff_abun_urine_signif_USD_genus.txt | sed 's/\]//g' | sort | uniq -c > $5\_otu/uclust_open/diff_abun_urine_signif_USD_genus_formatted.txt
awk 'BEGIN{print "USD_urine Taxa"}; {print}' $5\_otu/uclust_open/diff_abun_urine_signif_USD_genus_formatted.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD.txt
awk '{$1=$1};1' $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD2.txt
tr ' ' \\t < $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD2.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD3.txt
rm $5\_otu/uclust_open/diff_abun_urine_signif_USD_genus.txt $5\_otu/uclust_open/diff_abun_urine_signif_USD_genus_formatted.txt $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD.txt $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD2.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD3.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD4.txt
##############Get significant stool
awk '{ if ($7 < 0.05) { print }}' $5\_otu/uclust_open/diff_abun_stool > $5\_otu/uclust_open/diff_abun_stool_signif.txt
awk '{ if ($3 < 0) { print }}' $5\_otu/uclust_open/diff_abun_stool_signif.txt > $5\_otu/uclust_open/diff_abun_stool_signif_control.txt
sed 's/k__//g' $5\_otu/uclust_open/diff_abun_stool_signif_control.txt | sed 's/p__//g' | sed 's/c__//g' | sed 's/o__//g' | sed 's/f__//g' | sed 's/g__//g' | sed 's/s__.*//g' | sed 's/\;//g' | sed 's/NA//g' | awk -v OFS='\t' '{print $NF}' > $5\_otu/uclust_open/diff_abun_stool_signif_control_genus.txt
sed 's/\[//g' $5\_otu/uclust_open/diff_abun_stool_signif_control_genus.txt | sed 's/\]//g' | sort | uniq -c > $5\_otu/uclust_open/diff_abun_stool_signif_control_genus_formatted.txt
awk 'BEGIN{print "control_stool Taxa"}; {print}' $5\_otu/uclust_open/diff_abun_stool_signif_control_genus_formatted.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_control.txt
awk '{$1=$1};1' $5\_otu/uclust_open/diff_abun_stool_forheatmap_control.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_control2.txt
tr ' ' \\t < $5\_otu/uclust_open/diff_abun_stool_forheatmap_control2.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_control3.txt
rm $5\_otu/uclust_open/diff_abun_stool_signif_control_genus.txt $5\_otu/uclust_open/diff_abun_stool_signif_control_genus_formatted.txt $5\_otu/uclust_open/diff_abun_stool_forheatmap_control.txt $5\_otu/uclust_open/diff_abun_stool_forheatmap_control2.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_otu/uclust_open/diff_abun_stool_forheatmap_control3.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_control4.txt
awk '{ if ($3 > 0) { print }}' $5\_otu/uclust_open/diff_abun_stool_signif.txt > $5\_otu/uclust_open/diff_abun_stool_signif_USD.txt
sed 's/\[//g' $5\_otu/uclust_open/diff_abun_stool_signif_USD_genus.txt | sed 's/\]//g' | sort | uniq -c > $5\_otu/uclust_open/diff_abun_stool_signif_USD_genus_formatted.txt
awk 'BEGIN{print "USD_stool Taxa"}; {print}' $5\_otu/uclust_open/diff_abun_stool_signif_USD_genus_formatted.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD.txt
awk '{$1=$1};1' $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD2.txt
tr ' ' \\t < $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD2.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD3.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD3.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD4.txt
##Get counts per lowest taxonomy
sed 's/k__//g' $5\_otu/uclust_open/diff_abun_urine | sed 's/p__//g' | sed 's/c__//g' | sed 's/o__//g' | sed 's/f__//g' | sed 's/g__//g' | sed 's/s__.*//g' | sed 's/\;//g' | sed 's/NA//g' | awk -v OFS='\t' '{print $NF}' | sort | uniq -c > $5\_otu/uclust_open/diff_abun_urine_genusfamily_count.txt
awk 'BEGIN{print "Count Taxa"}; {print}' $5\_otu/uclust_open/diff_abun_urine_genusfamily_count.txt > $5\_otu/uclust_open/diff_abun_urine_genusfamily_count2.txt
awk '{$1=$1};1' $5\_otu/uclust_open/diff_abun_urine_genusfamily_count2.txt > $5\_otu/uclust_open/diff_abun_urine_genusfamily_count3.txt
tr ' ' \\t < $5\_otu/uclust_open/diff_abun_urine_genusfamily_count3.txt > $5\_otu/uclust_open/diff_abun_urine_genusfamily_count4.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_otu/uclust_open/diff_abun_urine_genusfamily_count4.txt > $5\_otu/uclust_open/diff_abun_urine_genusfamily_count5.txt
sed 's/k__//g' $5\_otu/uclust_open/diff_abun_stool | sed 's/p__//g' | sed 's/c__//g' | sed 's/o__//g' | sed 's/f__//g' | sed 's/g__//g' | sed 's/s__.*//g' | sed 's/\;//g' | sed 's/NA//g' | awk -v OFS='\t' '{print $NF}' | sort | uniq -c > $5\_otu/uclust_open/diff_abun_stool_genusfamily_count.txt
awk 'BEGIN{print "Count Taxa"}; {print}' $5\_otu/uclust_open/diff_abun_stool_genusfamily_count.txt > $5\_otu/uclust_open/diff_abun_stool_genusfamily_count2.txt
awk '{$1=$1};1' $5\_otu/uclust_open/diff_abun_stool_genusfamily_count2.txt > $5\_otu/uclust_open/diff_abun_stool_genusfamily_count3.txt
tr ' ' \\t < $5\_otu/uclust_open/diff_abun_stool_genusfamily_count3.txt > $5\_otu/uclust_open/diff_abun_stool_genusfamily_count4.txt
awk -F '\t' 'BEGIN {OFS = FS} {print $2,$1}' $5\_otu/uclust_open/diff_abun_stool_genusfamily_count4.txt > $5\_otu/uclust_open/diff_abun_stool_genusfamily_count5.txt
rm $5\_otu/uclust_open/diff_abun_urine_genusfamily_count.txt $5\_otu/uclust_open/diff_abun_urine_genusfamily_count2.txt $5\_otu/uclust_open/diff_abun_urine_genusfamily_count3.txt $5\_otu/uclust_open/diff_abun_urine_genusfamily_count4.txt $5\_otu/uclust_open/diff_abun_stool_genusfamily_count.txt $5\_otu/uclust_open/diff_abun_stool_genusfamily_count2.txt $5\_otu/uclust_open/diff_abun_stool_genusfamily_count3.txt $5\_otu/uclust_open/diff_abun_stool_genusfamily_count4.txt
##Scale diff abun analyses
stats-master/merge -k -e 0 $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD4.txt $5\_otu/uclust_open/diff_abun_stool_genusfamily_count5.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD_toscale.txt
stats-master/merge -k -e 0 $5\_otu/uclust_open/diff_abun_stool_forheatmap_control4.txt $5\_otu/uclust_open/diff_abun_stool_genusfamily_count5.txt > $5\_otu/uclust_open/diff_abun_stool_forheatmap_control_toscale.txt
stats-master/merge -k -e 0 $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD4.txt $5\_otu/uclust_open/diff_abun_urine_genusfamily_count5.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD_toscale.txt
stats-master/merge -k -e 0 $5\_otu/uclust_open/diff_abun_urine_forheatmap_control4.txt $5\_otu/uclust_open/diff_abun_urine_genusfamily_count5.txt > $5\_otu/uclust_open/diff_abun_urine_forheatmap_control_toscale.txt
Rscript scripts/scale_dysbiosis.R $5\_otu/uclust_open/diff_abun_urine_forheatmap_control_toscale.txt $5\_otu/uclust_open/diff_abun_urine_forheatmap_USD_toscale.txt $5\_otu/uclust_open/cont_urin_dev.txt $5\_otu/uclust_open/usd_urin_dev.txt $5\_otu/uclust_open/diff_abun_stool_forheatmap_control_toscale.txt $5\_otu/uclust_open/cont_stool_dev.txt $5\_otu/uclust_open/diff_abun_stool_forheatmap_USD_toscale.txt $5\_otu/uclust_open/usd_stool_dev.txt
awk -v OFS='\t' '{print $2,$3,$4,$5,$6}' $5\_otu/uclust_open/cont_urin_dev.txt > $5\_otu/uclust_open/cont_urin_dev2.txt
sed 's/\"//g' $5\_otu/uclust_open/cont_urin_dev2.txt > $5\_otu/uclust_open/cont_urin_dev3.txt
awk -v OFS='\t' '{print $2,$3,$4,$5,$6}' $5\_otu/uclust_open/usd_urin_dev.txt > $5\_otu/uclust_open/usd_urin_dev2.txt
sed 's/\"//g' $5\_otu/uclust_open/usd_urin_dev2.txt > $5\_otu/uclust_open/usd_urin_dev3.txt
awk -v OFS='\t' '{print $2,$3,$4,$5,$6}' $5\_otu/uclust_open/cont_stool_dev.txt > $5\_otu/uclust_open/cont_stool_dev2.txt
sed 's/\"//g' $5\_otu/uclust_open/cont_stool_dev2.txt > $5\_otu/uclust_open/cont_stool_dev3.txt
awk -v OFS='\t' '{print $2,$3,$4,$5,$6}' $5\_otu/uclust_open/usd_stool_dev.txt > $5\_otu/uclust_open/usd_stool_dev2.txt
sed 's/\"//g' $5\_otu/uclust_open/usd_stool_dev2.txt > $5\_otu/uclust_open/usd_stool_dev3.txt
awk -v OFS='\t' '{print $1,$5}' $5\_otu/uclust_open/usd_stool_dev3.txt > $5\_otu/uclust_open/usd_stool_dev4.txt
awk -v OFS='\t' '{print $1,$5}' $5\_otu/uclust_open/cont_stool_dev3.txt > $5\_otu/uclust_open/cont_stool_dev4.txt
awk -v OFS='\t' '{print $1,$5}' $5\_otu/uclust_open/usd_urin_dev3.txt > $5\_otu/uclust_open/usd_urin_dev4.txt
awk -v OFS='\t' '{print $1,$5}' $5\_otu/uclust_open/cont_urin_dev3.txt > $5\_otu/uclust_open/cont_urin_dev4.txt
##Merge scaled files
stats-master/merge -k -e 0 $5\_otu/uclust_open/cont_urin_dev4.txt $5\_otu/uclust_open/usd_urin_dev4.txt > $5\_otu/uclust_open/scaled_difference_urine.txt
stats-master/merge -k -e 0 $5\_otu/uclust_open/cont_stool_dev4.txt $5\_otu/uclust_open/usd_stool_dev4.txt > $5\_otu/uclust_open/scaled_difference_stool.txt
sed '1s/.*/Taxon	Control_urine	USD_urine/g'  $5\_otu/uclust_open/scaled_difference_urine.txt | grep -v "taxonomy" > $5\_otu/uclust_open/scaled_difference_urine2.txt
sed '1s/.*/Taxon	Control_stool	USD_stool/g'  $5\_otu/uclust_open/scaled_difference_stool.txt | grep -v "taxonomy" > $5\_otu/uclust_open/scaled_difference_stool2.txt
rm $5\_otu/uclust_open/*dev.txt
rm $5\_otu/uclust_open/*dev2.txt
rm $5\_otu/uclust_open/*dev3.txt
rm $5\_otu/uclust_open/*dev4.txt
awk '{ if ($2 > 5 || $3 > 5) { print }}' $5\_otu/uclust_open/scaled_difference_urine2.txt > $5\_otu/uclust_open/scaled_difference_urine3.txt
awk '{ if ($2 > 5 || $3 > 5) { print }}' $5\_otu/uclust_open/scaled_difference_stool2.txt > $5\_otu/uclust_open/scaled_difference_stool3.txt
awk -v OFS='\t' 'NR==1; NR > 1 {print $1, $2, -$3}' $5\_otu/uclust_open/scaled_difference_urine3.txt > $5\_otu/uclust_open/scaled_difference_urine4.txt
awk -v OFS='\t' 'NR==1; NR > 1 {print $1, $2, -$3}' $5\_otu/uclust_open/scaled_difference_stool3.txt > $5\_otu/uclust_open/scaled_difference_stool4.txt
awk -v OFS='\t' 'NR==1; NR > 1 {print $1, $2, $3, $2+$3}' $5\_otu/uclust_open/scaled_difference_urine4.txt > $5\_otu/uclust_open/scaled_difference_urine5.txt
awk -v OFS='\t' 'NR==1; NR > 1 {print $1, $2, $3, $2+$3}' $5\_otu/uclust_open/scaled_difference_stool4.txt > $5\_otu/uclust_open/scaled_difference_stool5.txt
awk 'NR==1; NR > 1 {print $0 | "sort -n -r -k4"}' $5\_otu/uclust_open/scaled_difference_urine5.txt > $5\_otu/uclust_open/scaled_difference_urine5_sorted.txt
awk 'NR==1; NR > 1 {print $0 | "sort -n -r -k4"}' $5\_otu/uclust_open/scaled_difference_stool5.txt > $5\_otu/uclust_open/scaled_difference_stool5_sorted.txt
sed '1s/.*/Taxon	Control_urine	USD_urine	Descend/g' $5\_otu/uclust_open/scaled_difference_urine5_sorted.txt | grep -v "taxonomy" > $5\_otu/uclust_open/scaled_difference_urine5_sorted2.txt
sed '1s/.*/Taxon	Control_stool	USD_stool	Descend/g' $5\_otu/uclust_open/scaled_difference_stool5_sorted.txt | grep -v "taxonomy" > $5\_otu/uclust_open/scaled_difference_stool5_sorted2.txt
awk '{print $1, $2, $3}' $5\_otu/uclust_open/scaled_difference_urine5_sorted2.txt > $5\_otu/uclust_open/scaled_difference_urine5_sorted3.txt
awk '{print $1, $2, $3}' $5\_otu/uclust_open/scaled_difference_stool5_sorted2.txt > $5\_otu/uclust_open/scaled_difference_stool5_sorted3.txt
##Extract stone $5\_otu/uclust_open, get average relative abundance, and plot heatmap
filter_otus_from_otu_table.py -i $5\_otu/uclust_open/otu_table_over_10_filtered_stone.biom -o $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2.biom -n 2
biom convert -i $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2.biom -o $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2.txt --to-tsv --header-key taxonomy
tail -n +2 $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2.txt  > $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2_2.txt
sed 's/\#OTU\ ID/OTU_ID/g' $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2_2.txt > $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2_3.txt
sed 's/\ /_/g' $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2_3.txt > $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2_4.txt
Rscript scripts/rel_abun_stone.R $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2_4.txt $5\_otu/uclust_open/stone_relabun.txt
awk -v OFS='\t' 'NR>1 {$1=""}1' $5\_otu/uclust_open/stone_relabun.txt | awk -v OFS='\t' '{$1=$1}1' > $5\_otu/uclust_open/stone_relabun2.txt
awk -v OFS='\t' '{print $NF}' $5\_otu/uclust_open/otu_table_over_10_filtered_stone_over2_4.txt > $5\_otu/uclust_open/stone_taxonomy.txt
sed 's/\;_/\;\ /g' $5\_otu/uclust_open/stone_taxonomy.txt | sed 's/k__//g' | sed 's/p__//g' | sed 's/c__//g' | sed 's/o__//g' | sed 's/f__//g' | sed 's/g__//g' | sed 's/s__.*//g' | sed 's/\;//g' | sed 's/NA//g' | awk -v OFS='\t' '{print $NF}'> $5\_otu/uclust_open/stone_taxonomy2.txt
awk -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; print sum}' $5\_otu/uclust_open/stone_relabun2.txt > $5\_otu/uclust_open/stone_relabun3.txt
sed '1s/.*/Average/g' $5\_otu/uclust_open/stone_relabun3.txt > $5\_otu/uclust_open/stone_relabun4.txt
paste $5\_otu/uclust_open/stone_relabun4.txt $5\_otu/uclust_open/stone_taxonomy2.txt > $5\_otu/uclust_open/stone_relabun_taxonomy.txt
awk -v OFS='\t' '{print $2,$1}' $5\_otu/uclust_open/stone_relabun_taxonomy.txt > $5\_otu/uclust_open/stone_relabun_heatmap.txt
awk '{ if ($2 > 0.01) { print }}' $5\_otu/uclust_open/stone_relabun_heatmap.txt > $5\_otu/uclust_open/stone_relabun_heatmap2.txt
awk 'NR==1; NR > 1 {print $0 | "sort -n -r -k2"}' $5\_otu/uclust_open/stone_relabun_heatmap2.txt > $5\_otu/uclust_open/stone_relabun_heatmap3.txt
set +e
Rscript scripts/heatmap_singlestudy.R $5\_otu/uclust_open/scaled_difference_urine5_sorted3.txt $5\_otu/uclust_open/scaled_difference_stool5_sorted3.txt $5\_otu/uclust_open/stone_relabun_heatmap3.txt
Rscript scripts/alpha_by_choice_singlestudy_stonecomp.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Stonecomp
Rscript scripts/alpha_by_choice_singlestudy_stonecomp.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stone Stonecomp
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group age_group
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group Study_location
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group Sex
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group Weight_group
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group antibiotic_12m
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group antibiotic_30d
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group gout
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group diabetes
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group HTN
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group diet
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group water
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group dessert
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group meat
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group fruit
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group veggie
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt stool Group bread
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group Ston_comp
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group age_group
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group Study_location
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group Sex
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group Weight_group
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group antibiotic_12m
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group antibiotic_30d
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group gout
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group diabetes
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group HTN
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group diet
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group water
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group dessert
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group meat
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group fruit
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group veggie
Rscript scripts/alpha_by_choice_singlestudy.r $5\_otu/uclust_open/current_alpha_sorted2.txt $5\_otu/uclust_open/map_sorted.txt urine Group bread
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group Ston_comp
Rscript scripts/beta_by_choice_singlestudy_stonecomp.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stone StoneComp
Rscript scripts/beta_by_choice_singlestudy_stonecomp.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool StoneComp
Rscript scripts/beta_by_choice_singlestudy_stonecomp.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine StoneComp
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group age_group
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group Sex
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt  stool Group Study_location
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group Weight_group
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group antibiotic_12m
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group antibiotic_30d
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group gout
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group diabetes
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group HTN
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group diet
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group water
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group dessert
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group meat
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group fruit
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group veggie
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt stool Group bread
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group Ston_comp
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group age_group
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group Sex
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt  urine Group Study_location
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group Weight_group
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group antibiotic_12m
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group antibiotic_30d
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group gout
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group diabetes
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group HTN
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group diet
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group water
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group dessert
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group meat
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group fruit
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group veggie
Rscript scripts/beta_by_choice_singlestudy.r $5\_otu/uclust_open/beta/weighted_unifrac_dm_sorted.txt $5\_otu/uclust_open/map_sorted.txt urine Group bread
cp $5\_otu/uclust_open/*.pdf $5\_otu/uclust_open/otu_table_over_10_filtered_map.txt $5\_otu/uclust_open/otu_table_over_10_filtered.biom data/downloads/$5\_otu
mv data/downloads/*.pdf data/downloads/$5\_otu
gzip data/downloads/$5\_otu/*.fna
