mypython=$1

echo "Python is $mypython"


# MARINE environment variable must be set from main marine directory:
#     export MARINE=$(pwd)

echo "Running bulk tests..."

tests_folder="strandedness_tests/"

echo "Bulk tests scripts"
ls -lh $MARINE/tests/$tests_folder/scripts/

#for t in "F1R2_pair_test-single_end_mode" "F1R2_pair_test" "F2R1_end_second_in_pair_test" "same_pos_dif_reads_test" "tax1bp3_chr17_3665556_read_test" "pair_test"
for t in "F1R2_pair_test-single_end_mode"
do
    echo $t
    echo "Removing old files..."
    rm $MARINE/tests/$tests_folder$t/* -r

    echo "Running tests..."
    bash $MARINE/tests/$tests_folder/scripts/$t.sh $mypython
   
done


echo "Running single-cell tests..."


tests_folder="singlecell_tests/"

echo "SC tests scripts"
ls -lh $MARINE/tests/$tests_folder/scripts/


#for t in "only_5_cells_test" "long_read_sc_test"
#
#do
#    echo $t
#    echo "Removing old files..."
#    rm $MARINE/tests/$tests_folder$t/* -r
#
#    echo "Running tests..."
#    bash $MARINE/tests/$tests_folder/scripts/$t.sh $mypython
#   
#done


echo "Checking results..."
$mypython $MARINE/tests/integration_tests_auto_check.py

exitcode=$?

echo "Exit code: $exitcode"
exit $exitcode