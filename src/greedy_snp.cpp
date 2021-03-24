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
#include "common.hpp"

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

  // To track retained positions 
  std::vector<int> new_c;
  int index;

  // Greedy algorithm
  int n = p.size(), i=0, j=0, count=0, event1, event2;  
  std::vector<bool> R(n, 0);

  auto tStart = std::chrono::system_clock::now();
  std::cout<< "INFO, VF::main, starting timer" << "\n";

  // Note: VCF variant positions are 1-based (i.e., they must be >=1)

  // Block of greedy
  while (i < n)
  {                                                 /* we are done when the end event of last SNP is processed */
    event1 = std::max(1, p[i]-parameters.alpha+1);  /* position of next beginning event */
    event2 = p[j]+1;                                /* position of next ending event */

    if(event2 <= event1) {                          /* we are processing ending event */
      if (!R[j]) count--;
      j++;
    }   

    if(event1 <= event2) {                          /* we are processing beginning event */
      count++;
      if (count > parameters.delta) {
        R[i] = 1;
        count--;                                    /* note the SNP position to be retained is still given by P[i] */
      }
      i++;
    }
  }       

  // End of greedy

  std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
  std::cout<< "INFO, VF::main, time taken by variant selection algorithm = " << wctduration.count() << " seconds" << "\n"; 

  std::cout<< "INFO, VF::main, count of variant containing positions retained = " << std::count(R.begin(), R.end(), true) << "\n";

  for (int i = 0; i < n; i++)
    if (R[i]) new_c.push_back(c[i]);

  std::cout<< "INFO, VF::main, count of variants retained = " << std::accumulate(new_c.begin(), new_c.end(), 0) << "\n";
  printVariantGapStats (R, p);
  if (parameters.prefix.length() > 0) print_snp_vcf(R, p, parameters);

  return 0;
}
