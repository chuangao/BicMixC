# BicMixC

BicMix is a sparse matrix decomposition tool. Given a matrix **Y** with dimension of P by N, BicMix decompose it into the product of two sparse matrices

# How to install

git clone https://github.com/chuangao/BicMixC.git <br/>
cd BicMix <br/>
make <br/>

This will generate the executible BicMix <br/>
For linux, the default compilation have openmp turned on <br/>
For macOS, openmp is turned off because by default it is not supported (yet) <br/>

# Usage
### To run without using openmp
./BicMix --y input_file --out output_directory <br/>
### To run using openmp <br/>
OMP_NUM_THREADS=10 ./BicMix --y input_file --out output_directory <br/>

### Parameters accepted
**--y** The matrix that is to be decomposed <br/>
**--out** The output directory to dump the results <br/>
**--nf** The initial number of factors that the algorithm should start with, BicMix will shrink this number to the true number (assume there is) <br/>
**--itr** The maximum number of iterations to run, default to 5000 <br/>
**--interval** How many iterations that the number of nonzero elements remain unchanged (to claim that the algorithm converges), default to 500 iterations <br/>
**--write_itr** Every how many iterations that the result will be written, default to 50






