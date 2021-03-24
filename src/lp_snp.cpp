#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector> 
#include <algorithm>
#include <ctime>
#include <chrono>
#include <numeric>
#include <cstdlib>
#include <random>
#include "common.hpp"
#include "gurobi_c++.h"

int main(int argc, char **argv) {

  //parse command line arguments
  Parameters parameters;
  parseandSave(argc, argv, parameters);

  //*********************************************************
  // Extract SNPs and allele count from VCF
  // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
  srand(time(0)); int random = rand() % 100000;  
  std::string tmp_file = ".VF." + std::to_string(random) + ".txt";
  std::string cmd = std::string(TOSTRING(VCFTOOLSPATH)) + " --vcf " + parameters.vcffile + " --chr " + parameters.chr + " --counts --remove-indels --out " + tmp_file + " 2>/dev/null";
  std::cout << "INFO, VF::main, extracting SNPs from vcf file using command = " << cmd << std::endl;
  std::system(cmd.c_str());
  std::cout << "INFO, VF::main, vcftools finished" << std::endl;

  //*********************************************************
  // Reading from temp file to store c

  std::ifstream file(tmp_file + ".frq.count");
  std::vector<int> p, c; 
  std::string line;

  while (std::getline(file, line))
  {
    std::string col1;
    int col2, col3;

    std::istringstream iss(line);
    iss >> col1;
    if (col1.compare("CHROM") != 0)  //ignore the first header line
    {
      iss >> col2 >> col3;
      p.push_back(col2);
      c.push_back(col3 - 1);
    }
  }

  //keep only one record per loci
  ignoreDuplicateSNPrecords(p, c);
  assert (std::is_sorted(p.begin(), p.end()));

  cmd = "rm -f " + tmp_file + "*";
  std::system(cmd.c_str()); //delete tmp file
  if (p.size() == 0)
  {
    std::cerr << "ERROR, VF::main, count of variants is zero, did you provide the correct vcf file and chrommosome id?" << std::endl;
    exit(1);
  }
  //*********************************************************

  std::cout<< "INFO, VF::main, count of variant containing positions = " << p.size() << "\n";
  std::cout<< "INFO, VF::main, count of variants = " << std::accumulate(c.begin(), c.end(), 0) << "\n";

  int n = p.size();

  // To track retained positions 
  std::vector<int> new_c;

  //*********************************************************
  std::vector<bool> R(n, 0);  /* R[i] = true means variant position i is retained*/ 
  GRBVar* x = 0;

  // Lp algorithm using Gurobi
  auto tStart = std::chrono::system_clock::now();
  std::cout<< "INFO, VF::main, starting timer" << "\n";

  try
  {
    //Gurobi modeling
    std::cout<< "INFO, VF::main, Gurobi solver starting" << "\n";
    GRBEnv* env = 0;
    env = new GRBEnv();
    GRBModel model = GRBModel(*env);

    //comment out this line to enable Gurobi output log
    model.set(GRB_IntParam_LogToConsole, 0);

    // Create variables
    std::vector<double> zeros (n, 0.0);
    std::vector<double> ones (n, 1.0);
    std::vector<char> type (n, GRB_CONTINUOUS);
    //std::vector<char> type (n, GRB_BINARY);
    //NOTE: I'm finding that Gurobi requires less memory for ILP than LP, not sure why
    //TODO:use ILP instead?

    x = model.addVars(zeros.data(), ones.data(), NULL, type.data(), NULL, n);

    // Set objective
    GRBLinExpr obj = 0;
    for (int i = 0; i < n; i++)
      obj += c[i]*x[i];

    //maximize c.x
    model.setObjective(obj, GRB_MAXIMIZE);

    // Add constraints
    for (int i = 0; i < n; i++)
    {
      GRBLinExpr lhs = 0;
      for (int j = i; j >= 0 && p[i]-p[j] < parameters.alpha; j--)
        lhs += x[j];

      model.addConstr(lhs , GRB_LESS_EQUAL, 1.0 * parameters.delta);
      //this adds each contraint row of A.x <= b one by one
    }

    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
      double objval = model.get(GRB_DoubleAttr_ObjVal);
      std::cout << "Optimal objective: " << objval << std::endl;
    } 

    // To store SNP positions retained
    for(int i =0; i < n; i++){
      if(x[i].get(GRB_DoubleAttr_X) < 0.5)
        R[i] = true;
    }
  } 
  catch (GRBException e) {
    std::cout << "ERROR, VF::main, Gurobi exception raised, error code = " << e.getErrorCode() << ", ";
    std::cout << e.getMessage() << std::endl;
    if (e.getErrorCode() == 10009)
    {
      std::cout << "Step1. Get your free Gurobi academic license code by registering here: https://www.gurobi.com/downloads/end-user-license-agreement-academic" << std::endl;
      std::cout << "Step2. Add your licence key by using build/gurobi910/linux64/bin/grbgetkey tool" << std::endl;
    }
    exit(1);
  } 
  catch (...) {
    std::cout << "Error during optimization" << std::endl;
    exit(1);
  }

  // End of LP
  std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
  std::cout<< "INFO, VF::main, time taken by variant selection algorithm = " << wctduration.count() << " seconds" << "\n"; 

  std::cout<< "INFO, VF::main, count of variant containing positions retained = " << std::count(R.begin(), R.end(), true) << "\n";

  for (int i = 0; i < n; i++)
    if (R[i]) new_c.push_back(c[i]);

  std::cout<< "INFO, VF::main, count of variants retained = " << std::accumulate(new_c.begin(), new_c.end(), 0) << "\n";
  printVariantGapStats (R, p);
  if (parameters.prefix.length() > 0) print_snp_vcf(R, p, parameters);

  //clear memory
  delete[] x;
  return 0;
}
