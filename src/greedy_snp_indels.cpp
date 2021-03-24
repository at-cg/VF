#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector> 
#include <algorithm>
#include <ctime>
#include <cstring>
#include <chrono>
#include <numeric>
#include <cassert>
#include <unordered_map>
#include "ext/prettyprint.hpp"
#include "common.hpp"

/********* Helper functions ******/

/**
 * @brief  parse VCF file to record info of insertion and deletion SVs
 */
void parseVCF_indel (const std::string &vcf_file, const std::string &chromosomeId, std::vector<int> &indelpos, std::vector<int> &indellen)
{  
  // Extract indels from VCF
  // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
  srand(time(0)); int random = rand() % 100000;  
  std::string tmp_file = ".VF." + std::to_string(random) + ".txt";
  //the following command assumes that vcf record will contain "VT=INDEL" for indel variants
  std::string cmd = " grep \"VT=INDEL\" " + vcf_file + " >  " + tmp_file;
  std::cout << "INFO, VF::main, extracting indels from vcf file using command: " << cmd << std::endl;
  std::system(cmd.c_str());
  std::cout << "INFO, VF::main, grep finished" << std::endl;

  std::ifstream file (tmp_file);
  std::string line;
  while (std::getline(file, line))
  {
    if (line[0] != '#') //ignore beginning header lines 
    {
      std::string col1, col3, col4, col5;
      int col2;

      std::istringstream iss(line);
      iss >> col1 >> col2 >> col3 >> col4 >> col5;

      if (col1.compare(chromosomeId) == 0)
      {
        //ref and alt sequences should only contain alphabetic letter
        assert( std::all_of(std::begin(col4), std::end(col4),
                [](char c){ return std::isalpha(c); }));

        assert( std::all_of(std::begin(col5), std::end(col5),
                [](char c){ return std::isalpha(c); }));

        int var_size = col5.length() - col4.length(); //alt - ref, -ve if deletion, +ve for insertion
        if (var_size != 0)
        {
          indelpos.push_back(col2);
          indellen.push_back(var_size);
        }
      }
    }
  }

  if (indelpos.size() == 0 || indellen.size() == 0)
  {
    std::cerr << "ERROR, VF::parseVCF_indel, count of indels found is zero, did you provide the correct vcf file and chrommosome id?" << std::endl;
    exit(1);
  }

  cmd = "rm -f " + tmp_file + "*";
  std::system(cmd.c_str()); //delete tmp file
}

/**
 * @brief  parse VCF file to record info of insertion and deletion SVs
 */
void parseVCF_SNP (const std::string &vcf_file, const std::string &chromosomeId, std::vector<int> &snppos, std::vector<int> &snpcount)
{
  srand(time(0)); int random = rand() % 100000;  
  std::string tmp_file = ".VF." + std::to_string(random) + ".txt";
  std::string cmd = std::string(TOSTRING(VCFTOOLSPATH)) + " --vcf " + vcf_file + " --chr " + chromosomeId + " --counts --remove-indels --out " + tmp_file + " 2>/dev/null";
  std::cout << "INFO, VF::main, extracting SNPs from vcf file using command: " << cmd << std::endl;
  std::system(cmd.c_str());
  std::cout << "INFO, VF::main, vcftools finished" << std::endl;

  //*********************************************************
  // Reading from temp file to store c

  std::ifstream file(tmp_file + ".frq.count");
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
      snppos.push_back(col2);
      snpcount.push_back(col3 - 1);
    }
  }

  //keep only one record per loci
  ignoreDuplicateSNPrecords(snppos, snpcount);
  assert (std::is_sorted(snppos.begin(), snppos.end()));

  cmd = "rm -f " + tmp_file + "*";
  std::system(cmd.c_str()); //delete tmp file
  if (snppos.size() == 0)
  {
    std::cerr << "ERROR, VF::parseVCF_SNP, count of SNPs is zero, did you provide the correct vcf file and chrommosome id?" << std::endl;
    exit(1);
  }
}

/**
 * @brief   compute left-most reachable vertex from each variant position
 *          using up to alpha-1 labeled edges     
 */
void calculateLeftMostReachable (std::vector<int> &reach, const std::vector<int> &pos_u, const std::vector<int> &indelpos, const std::vector<int> &indellen, const int &alpha)
{
  assert(alpha > 2); //the logic below requires this
  assert (reach.size() == pos_u.size());

  /* 1. build edge-labeled graph g */
  //reminder: vcf file contains 1-based position offsets
  //suppose x = last position containing variant
  //then build backbone chain of x+1 vertices
  int g_edge_count = *(std::max_element(pos_u.begin(), pos_u.end())); 
  //std::cout << "INFO, VF::calculateLeftMostReachable, max SV position recorded = " << g_edge_count << "\n";
  std::vector<int> reach_tmp (g_edge_count+1, 0); //leftmost reachable from all vertices

  //not necessary, but bool vectors may help with faster lookup
  std::vector<bool> out_edges (g_edge_count+1, 0);
  std::vector<bool> in_edges (g_edge_count+1, 0);

  //save outgoing neighbors associated with unlabeled deletion edges only
  std::unordered_map<int,std::vector<int>> g_out_d_neighbors;
  for (std::size_t i = 0; i < indelpos.size(); i++)
  {
    if (indellen[i] < 0) //only deletions
    {
      //add deletion edge from svpos[i] to svpos[i] + |svlen[i]|
      int from = indelpos[i]; int to = indelpos[i] + std::abs(indellen[i]);

      out_edges[from] = true, in_edges[to] = true;

      if (g_out_d_neighbors.find(from) != g_out_d_neighbors.end()) //search
        g_out_d_neighbors[from].push_back(to);   //add 'to'
      else
        g_out_d_neighbors[from] = std::vector<int> {to}; //initialize
    }
  }

  //initialize vector of size alpha - 1
  //first value signifies vertex reachable using 1 edge, second using 2 edges and so on
  std::vector<int> currentPos (alpha-1);
  std::vector<int>  previousPos (alpha-1);

  //each vertex sends update using its out-going (unlabeled) edge
  //these updates are maintained in a hash table, until the vertex is visited
  std::unordered_map<int, std::vector<int>> updateAhead;

  std::cout << "INFO, VF::calculateLeftMostReachable, computing window ranges...\n" << std::flush;

  for (std::size_t i = 1; i <= g_edge_count; i++)
  {
    if (i == 1) //base case
    {
      std::fill (currentPos.begin(), currentPos.end(), 1); //can only reach 1st vertex
    }
    else
    {
      //use labeled edge
      std::memmove(&currentPos[1], &currentPos[0], (alpha-2)*sizeof(int)); //right shift by 1
      currentPos[0] = i - 1;  // vertex i-1 is reachable using one edge

      //use updates sent via in-coming unlabeled edges
      if (in_edges[i])
      {
        assert (updateAhead.find(i) != updateAhead.end());

        //take pairwise minimum with vector saved at updateAhead[i]
        std::transform (updateAhead[i].begin(), updateAhead[i].end(),
            currentPos.begin(), currentPos.begin(), [&](int a, int b){return std::min(a,b);});

        updateAhead.erase(i); //no longer needed
      }
    }

    //check unlabelled outgoing edges going *out* of vertex position i
    if (out_edges[i]) //check adjacency list
    {
      assert (g_out_d_neighbors.find(i) != g_out_d_neighbors.end());

      for (auto &v: g_out_d_neighbors[i]) //iterate over all out-going vertices
      {
        assert (v>i);

        if (updateAhead.find(v) == updateAhead.end())
          updateAhead[v] = std::vector<int> (alpha-1, v); //initialize

        assert (currentPos.size() == updateAhead[v].size());

        //take pairwise minimum with vector currentPos
        std::transform (currentPos.begin(), currentPos.end(),
            updateAhead[v].begin(), updateAhead[v].begin(), [&](int a, int b){return std::min(a,b);});
      }
    }

    reach_tmp[i] = currentPos[alpha-2]; //this is the last value in currentPos vector
  }

  std::cout << "INFO, VF::calculateLeftMostReachable, done" << std::endl;

  //we need reachability info only for variant positions
  for (std::size_t i = 0; i < reach.size(); i++)
    reach[i] = reach_tmp[pos_u[i]];

}

/**
 * @brief   Compute penalty associated with dropping variants at each position.
 *          This function also computes c vector along side penalties.
 */
void calculatePenalty (std::vector<int> &penalty, std::vector<int> &c, std::vector<int> &c_snp, const std::vector<int> &pos_u, const std::vector<int> &indelpos, const std::vector<int> &indellen, const std::vector<int> &snppos, const std::vector<int> &snpcount)
{
  assert(penalty.size() == pos_u.size());

  int g_edge_count = *(std::max_element(pos_u.begin(), pos_u.end())); 
  std::vector<int> max_ins_size (g_edge_count+1, 0);
  std::vector<int> max_del_size (g_edge_count+1, 0);
  std::vector<bool> snp_present (g_edge_count+1, false);
  std::vector<int> count (g_edge_count+1, 0);
  std::vector<int> count_snp (g_edge_count+1, 0);

  for (std::size_t i = 0; i < indelpos.size(); i++)
  {
    if (indellen[i] > 0)
      max_ins_size[indelpos[i]] = std::max(max_ins_size[indelpos[i]], indellen[i]);
    else
      max_del_size[indelpos[i]] = std::max(max_del_size[indelpos[i]], std::abs(indellen[i]));
    count[indelpos[i]]++;
  }

  for (std::size_t i = 0; i < snppos.size(); i++)
  {
    snp_present[snppos[i]] = true;
    count[snppos[i]] += snpcount[i];
    count_snp[snppos[i]] += snpcount[i];
  }

  for (std::size_t i = 0; i < penalty.size(); i++)
  {
    if (max_del_size[pos_u[i]] > 0)
    {
      //deletion penalty subsumes SNP penalty
      penalty[i] = max_ins_size[pos_u[i]] + max_del_size[pos_u[i]];
    }
    else
    {
      if (snp_present[pos_u[i]] == true)
        penalty[i] = max_ins_size[pos_u[i]] + 1;
      else
        penalty[i] = max_ins_size[pos_u[i]];
    }

    c[i] = count[pos_u[i]];
    c_snp[i] = count_snp[pos_u[i]];

    assert (penalty[i] > 0);
    assert (c[i] > 0);
    if (snp_present[pos_u[i]]) assert (c_snp[i] > 0);
  }
}

int main(int argc, char **argv) {

  //parse command line arguments
  Parameters parameters;
  parseandSave(argc, argv, parameters);

  std::vector<int> indelpos, indellen; 
  parseVCF_indel (parameters.vcffile, parameters.chr, indelpos, indellen); 
  assert (indelpos.size() == indellen.size());
  assert (std::is_sorted(indelpos.begin(), indelpos.end())); //must be sorted in ascending order

  std::vector<int> snppos, snpcount; 
  parseVCF_SNP (parameters.vcffile, parameters.chr, snppos, snpcount);
  assert (snppos.size() == snpcount.size());
  assert (std::is_sorted(snppos.begin(), snppos.end())); //must be sorted in ascending order

  //variant positions (unique values)
  std::vector<int> pos_u;
  pos_u.insert (pos_u.end(), indelpos.begin(), indelpos.end());
  pos_u.insert (pos_u.end(), snppos.begin(), snppos.end());
  std::sort (pos_u.begin(), pos_u.end());
  pos_u.erase(std::unique(pos_u.begin(), pos_u.end()), pos_u.end() );

  std::cout<< "INFO, VF::main, count of variant containing positions = " << pos_u.size() << "\n";
  std::cout<< "INFO, VF::main, count of indels = " << indelpos.size() << "\n";
  std::cout<< "INFO, VF::main, count of SNP variants = " << std::accumulate(snpcount.begin(), snpcount.end(), 0) << "\n";


  // Greedy algorithm
  auto tStart = std::chrono::system_clock::now();
  std::cout<< "INFO, VF::main, starting timer" << "\n";

  int n = pos_u.size();
  //compute reachability
  std::vector<int> reach (n);
  calculateLeftMostReachable (reach, pos_u, indelpos, indellen, parameters.alpha); 

  //compute penalty of variant removal for each position
  std::vector<bool> R(n, 0);  /* R[i] = true means variant position i is retained*/
  std::vector<int> c(n, 0); //count of variants at these positions
  std::vector<int> c_snp(n, 0); //count of SNP variants at these positions
  std::vector<int> penalty (n);
  calculatePenalty (penalty, c, c_snp, pos_u, indelpos, indellen, snppos, snpcount);

  //sum of 'c_snp' values should equal sum of SNPs
  assert (std::accumulate(c_snp.begin(), c_snp.end(), 0) == std::accumulate(snpcount.begin(), snpcount.end(), 0));
  //sum of 'c' values should equal sum of indels and SNPs
  assert (std::accumulate(c.begin(), c.end(), 0) == indelpos.size() + std::accumulate(snpcount.begin(), snpcount.end(), 0));

  // Block of greedy
  std::vector<int> cumulative_penalty (n,0); 
  //cumulative_penalty[i] indicates cumulative penalty until svpos_u[i] (exclusive)
  for (std::size_t i = 0; i < n; i++)
  {
    //penalty to drop svpos_u[i]
    int pen = penalty[i];

    //check range
    auto leftMostVariantPos_it = std::upper_bound (pos_u.begin(), pos_u.end(), reach[i]);
    auto leftMostVariantPos_index = leftMostVariantPos_it - pos_u.begin();

    int penalty_already_incurred = cumulative_penalty[i] - cumulative_penalty[leftMostVariantPos_index];

    if (penalty_already_incurred + pen <= parameters.delta)
    {
      if (i < n-1)
        cumulative_penalty[i+1] = cumulative_penalty[i] + pen; //drop
    }
    else
    {
      R[i] = true; //retain
      if (i < n-1)
        cumulative_penalty[i+1] = cumulative_penalty[i];
    }
  }

  // End of greedy
   
  std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
  std::cout<< "INFO, VF::main, time taken by variant selection algorithm = " << wctduration.count() << " seconds" << "\n"; 
  std::cout<< "INFO, VF::main, count of variant containing positions retained = " << std::count(R.begin(), R.end(), true) << "\n";

  int count_variants_retained=0;
  for (std::size_t i = 0; i < n; i++) if(R[i]) count_variants_retained += c[i]; 
  std::cout<< "INFO, VF::main, count of variants retained = " << count_variants_retained << "\n";

  int count_snp_variants_retained=0;
  for (std::size_t i = 0; i < n; i++) if(R[i]) count_snp_variants_retained += c_snp[i]; 
  std::cout<< "INFO, VF::main, count of SNP variants retained = " << count_snp_variants_retained << "\n";
  std::cout<< "INFO, VF::main, count of indel variants retained = " << count_variants_retained - count_snp_variants_retained << "\n";

  printVariantGapStats (R, pos_u);
  if (parameters.prefix.length() > 0) print_snp_indel_vcf (R, pos_u, parameters);
  return 0;
}
