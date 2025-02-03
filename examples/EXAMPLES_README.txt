Before running examples, MARINE environment variable must be set from main MARINE directory (i.e. the directory containingg the examples directory). The path to the main MARINE directory should also be appended to the PATH environment variable. 

    export MARINE=$(pwd)
    export PATH=$PATH:$(pwd)

Run the examples from the main MARINE directory like so:

    bash examples/test_sc_subset_CT.sh