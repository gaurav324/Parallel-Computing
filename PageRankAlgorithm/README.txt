# As you have already unzipped the project. Instructions for executing:

# Compile
g++ -O -Wno-deprecated  *.cpp *.h -fopenmp

# For executing the serial version.
./a.out <Name of the data set>

# For executing the parallel version.
./a.out <Name of the data set> parallel

# I have created two sample test cases. Incase, you dont't want to go and submit a job to TACC, you can simply run for these smaller test cases.

# SAMPLE COMMANDS:

# Executes "test_case" in serial.
./a.out test_case 

# Executest
./a.out test_case_medium parallel
