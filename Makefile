GUROBI_INSTALL=$(shell pwd)/build/gurobi910/linux64
VCFTOOLS_INSTALL=$(shell pwd)/build/vcftools-0.1.16/bin/vcftools
TARGET_DIR=$(shell pwd)/build
CPPFLAGS= -g -std=c++11 -DNDEBUG -O3 

all:
	mkdir -p build
	$(CXX) $(CPPFLAGS) -D VCFTOOLSPATH=$(VCFTOOLS_INSTALL) -o $(TARGET_DIR)/greedy_snp src/greedy_snp.cpp
	$(CXX) $(CPPFLAGS) -o $(TARGET_DIR)/greedy_sv_indels src/greedy_sv_indels.cpp 
	$(CXX) $(CPPFLAGS) -D VCFTOOLSPATH=$(VCFTOOLS_INSTALL) -m64 -o $(TARGET_DIR)/lp_snp -I $(GUROBI_INSTALL)/include/ -L  $(GUROBI_INSTALL)/lib/ src/lp_snp.cpp -lgurobi_c++ $(GUROBI_INSTALL)/lib/libgurobi91.so -lm
	$(CXX) $(CPPFLAGS) -D VCFTOOLSPATH=$(VCFTOOLS_INSTALL) -m64 -o $(TARGET_DIR)/ilp_sv_indels -I $(GUROBI_INSTALL)/include/ -L  $(GUROBI_INSTALL)/lib/ src/ilp_sv_indels.cpp -lgurobi_c++ $(GUROBI_INSTALL)/lib/libgurobi91.so -lm
	@echo "check executables in build directory"


clean:
	rm -f $(TARGET_DIR)/greedy_snp
	rm -f $(TARGET_DIR)/greedy_sv_indels
	rm -f $(TARGET_DIR)/lp_snp
	rm -f $(TARGET_DIR)/ilp_sv_indels
