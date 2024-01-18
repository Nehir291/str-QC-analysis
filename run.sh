#To run commands;
#bash run.sh list sampleId
#The "list" is a file containing the list of VariantId. One VariantId per line.

list=$1
sampleId=$2

cat $list | while read line;
do
  echo $line
  python3 get_reviewer_image_sections_new.py "$line".svg "$line"_output1 "$line"_output2
  bash count_nucleotides_from_reads_Total_v3.sh "$line".svg "$line"_output2 > "$line"_TotalCounts
  bash count_nucleotides_from_reads_ORANGE_v3.sh "$line".svg "$line"_output2 > "$line"_orangeCounts
  python3 combine_files.py "$line".metrics.tsv "$line"_output1 "$line"_TotalCounts \
  "$line"_orangeCounts \
  "$line"_AllMetrices
done

rm *_output1
rm *_output2
rm *_TotalCounts
rm *_orangeCounts
ls *_AllMetrices > image_list

python merge_csv_files.py -i image_list -p $sampleId

rm image_list
rm *_AllMetrices

#Run the Rscript to get image metrices
Rscript Collect_xgbtree_scores.R HG00514.1_json_files.variants_reduced.tsv $sampleId"_metrics.csv" $sampleId
