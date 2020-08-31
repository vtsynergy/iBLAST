#!/bin/bash
echo "Extracting command-line BLAST source"
tar -vxzf ../ncbi-blast/ncbi-blast-2.8.1+-src-iBLAST.tar.gz -C ../ncbi-blast

echo "Building NCBI BLAST command line tools"
cd ../ncbi-blast/ncbi-blast-2.8.1+-src-iBLAST/c++
./configure
cd ReleaseMT/build
make all_r


