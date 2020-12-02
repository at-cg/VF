#Purpose: download and install dependencies

mainwd=$(pwd)  #project top-level directory
mkdir -p build && cd build
buildwd=$(pwd)

echo "downloading vcftools"
wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz

#compile vcftools
tar xzf vcftools-0.1.16.tar.gz && cd vcftools-0.1.16
./autogen.sh
./configure --prefix=$(pwd)
make -j -s && make install
rm -f "vcftools-0.1.16.tar.gz"
echo "vcftools download and compilation finished"

#get Gurobi
echo "downloading gurobi"
cd $buildwd 
wget  https://packages.gurobi.com/9.1/gurobi9.1.0_linux64.tar.gz
tar xzf gurobi9.1.0_linux64.tar.gz
make -j -C gurobi910/linux64/src/build #re-compile gurobi cpp files using user's c++ compiler
cp gurobi910/linux64/src/build/libgurobi_c++.a gurobi910/linux64/lib
rm -f "gurobi9.1.0_linux64.tar.gz"
echo "gurobi download and compilation finished"

#check if appropriate files exist
if [ ! -f gurobi910/linux64/include/gurobi_c++.h ]; then echo "gurobi download failed"; fi
if [ ! -f gurobi910/linux64/src/build/libgurobi_c++.a ]; then echo "gurobi compilation failed"; fi
if [ ! -f vcftools-0.1.16/bin/vcftools ]; then echo "vcftools compilation failed"; fi

echo "Looks like it went okay, now run <make>"
#Next, run make
