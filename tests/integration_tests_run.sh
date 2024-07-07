echo "Running bulk tests..."

tests_folder="strandedness_tests/"

for t in "F1R2_pair_test-single_end_mode" "F1R2_pair_test" "F2R1_end_second_in_pair_test" "same_pos_dif_reads_test" "tax1bp3_chr17_3665556_read_test" "pair_test"

do
    echo $t
    echo "Removing old files..."
    rm $tests_folder$t/* -r

    echo "Running tests..."
    bash $tests_folder/scripts/$t.sh 
   
done


echo "Running single-cell tests..."


tests_folder="singlecell_tests/"
for t in "only_5_cells_test" "long_read_sc_test"

do
    echo $t
    echo "Removing old files..."
    rm $tests_folder$t/* -r

    echo "Running old tests..."
    bash $tests_folder/scripts/$t.sh 
   
done


echo "Checking results..."
python integration_tests_auto_check.py

exitcode=$?

echo "Exit code: $exitcode"
exit $exitcode