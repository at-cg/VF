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
#include "gurobi_c++.h"

/********* Helper functions ******/

/**
 * @brief  parse VCF file to record info of insertion and deletion SVs
 */
void parseVCF (const std::string &sv_vcf_file, const std::string &chromosomeId, std::vector<int> &svpos, std::vector<int> &svlen)
{
  std::ifstream file (sv_vcf_file);
  std::string line;
  while (std::getline(file, line))
  {
    if (line[0] != '#') //ignore beginning header lines 
    {
      int column = 0, pos, len;
      std::string chr, type;
      for (size_t i = 0; i < line.size() && line[i] != '\0' && line[i] != '\n'; i++) 
      {
        if (column == 0 && line[i] != '\t')
          chr += line[i];

        if (column == 1 && line[i-1] == '\t')
          pos = std::atoi(&line[i]);

        if (column == 7 && (std::strncmp(&line[i], "SVLEN=", 6) == 0))
          len = std::abs((int) std::atof(&line[i + 6])); //ignore sign

        if (column == 7 && strncmp(&line[i], "SVTYPE=", 7) == 0) 
          type = std::string(&line[i + 7]);

        if (line[i] == '\t') 
          column++; //next column
      }

      if (chr == chromosomeId)
      {
        //only consider INSs and DELs
        if (std::strncmp(type.c_str(), "INS", 3) == 0) 
        {
          svpos.push_back(pos);
          svlen.push_back(len);
        }
        else if (std::strncmp(type.c_str(), "DEL", 3) == 0)
        {
          svpos.push_back(pos);
          svlen.push_back(-1 * len); //negative for deletions
        }
      }
    }
  }

  if (svpos.size() == 0 || svlen.size() == 0)
  {
    std::cerr << "ERROR, VF::parseVCF, count of SVs found is zero, did you provide the correct vcf file and chrommosome id?" << std::endl;
    exit(1);
  }
}

/**
 * @brief   compute left-most reachable vertex from each variant position
 *          using up to alpha-1 labeled edges     
 */
void calculateLeftMostReachable (std::vector<int> &reach, const std::vector<int> &svpos_u, const std::vector<int> &svpos, const std::vector<int> &svlen, const int &alpha)
{
  assert(alpha > 2); //the logic below requires this
  assert (reach.size() == svpos_u.size());

  /* 1. build edge-labeled graph g */
  //reminder: vcf file contains 1-based position offsets
  //suppose x = last position containing variant
  //then build backbone chain of x+1 vertices
  int g_edge_count = *(std::max_element(svpos.begin(), svpos.end())); 
  //std::cout << "INFO, VF::calculateLeftMostReachable, max SV position recorded = " << g_edge_count << "\n";
  std::vector<int> reach_tmp (g_edge_count+1, 0); //leftmost reachable from all vertices

  //not necessary, but bool vectors may help with faster lookup
  std::vector<bool> out_edges (g_edge_count+1, 0);
  std::vector<bool> in_edges (g_edge_count+1, 0);


  //save outgoing neighbors associated with unlabeled deletion edges only
  std::unordered_map<int,std::vector<int>> g_out_d_neighbors;
  for (std::size_t i = 0; i < svpos.size(); i++)
  {
    if (svlen[i] < 0) //only deletions
    {
      //add deletion edge from svpos[i] to svpos[i] + |svlen[i]|
      int from = svpos[i]; int to = svpos[i] + std::abs(svlen[i]);

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
    reach[i] = reach_tmp[svpos_u[i]];
}

/**
 * @brief   Compute penalty associated with dropping variants at each position.
 *          This function also computes c vector along side penalties.
 */
void calculatePenalty (std::vector<int> &penalty, std::vector<int> &c, const std::vector<int> &svpos_u, const std::vector<int> &svpos, const std::vector<int> &svlen)
{
  assert(penalty.size() == svpos_u.size());

  int g_edge_count = *(std::max_element(svpos.begin(), svpos.end())); 
  std::vector<int> max_ins_size (g_edge_count+1, 0);
  std::vector<int> max_del_size (g_edge_count+1, 0);
  std::vector<int> count (g_edge_count+1, 0);

  for (std::size_t i = 0; i < svpos.size(); i++)
  {
    if (svlen[i] > 0)
      max_ins_size[svpos[i]] = std::max(max_ins_size[svpos[i]], svlen[i]);
    else
      max_del_size[svpos[i]] = std::max(max_del_size[svpos[i]], std::abs(svlen[i]));

    count[svpos[i]]++;
  }

  for (std::size_t i = 0; i < penalty.size(); i++)
  {
    penalty[i] = max_ins_size[svpos_u[i]] + max_del_size[svpos_u[i]];
    c[i] = count[svpos_u[i]];
  }
}

int main(int argc, char **argv) {

  //parse command line arguments
  Parameters parameters;
  parseandSave_ILP(argc, argv, parameters);

  //*********************************************************
  // Reading from file to store c

  std::vector<int> svpos, svlen; 
  parseVCF (parameters.vcffile, parameters.chr, svpos, svlen); 
  assert (svpos.size() == svlen.size());
  assert (std::is_sorted(svpos.begin(), svpos.end())); //must be sorted in ascending order
  //*********************************************************

  //variant positions (i.e., unique values in svpos vector)
  std::vector<int> svpos_u;
  std::unique_copy(svpos.begin(), svpos.end(), std::back_inserter(svpos_u));
  int n = svpos_u.size();
  std::vector<bool> R(n, 0);  /* R[i] = true means variant position i is retained*/
  std::vector<int> c(n, 0); //count of variants at these positions

  std::cout<< "INFO, VF::main, count of variant containing positions = " << svpos_u.size() << "\n";
  std::cout<< "INFO, VF::main, count of variants = " << svpos.size() << "\n";

  // ILP algorithm
  GRBVar* x = 0;

  auto tStart = std::chrono::system_clock::now();
  std::cout<< "INFO, VF::main, starting timer" << "\n";

  //compute reachability
  std::vector<int> reach (n);
  calculateLeftMostReachable (reach, svpos_u, svpos, svlen, parameters.alpha); 

  //compute penalty of variant removal for each position
  std::vector<int> penalty (n);
  calculatePenalty (penalty, c, svpos_u, svpos, svlen);

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
    std::vector<char> type (n, GRB_BINARY);

    x = model.addVars(zeros.data(), ones.data(), NULL, type.data(), NULL, n);

    // Set objective
    GRBLinExpr obj = 0;

    if (parameters.pos)
        std::cout << "INFO, VF::main, ILP solver will attempt to minimize variant positions" << "\n";
    else
        std::cout << "INFO, VF::main, ILP solver will attempt to minimize count of variants " << "\n";

    for (int i = 0; i < n; i++)
    {
      if (parameters.pos)
        obj += 1*x[i];
      else
        obj += c[i]*x[i];
    }

    //maximize c.x
    model.setObjective(obj, GRB_MAXIMIZE);

    // Add constraints
    for (int i = 0; i < n; i++)
    {
      GRBLinExpr lhs = 0;

      auto leftMostVariantPos_it = std::upper_bound (svpos_u.begin(), svpos_u.end(), reach[i]);
      auto leftMostVariantPos_index = leftMostVariantPos_it - svpos_u.begin();
      for (int j = i; j >= leftMostVariantPos_index; j--)
        lhs += penalty[j] * x[j];

      model.addConstr(lhs , GRB_LESS_EQUAL, 1.0 * parameters.delta);
      //this adds each contraint row of A.x <= b one by one
    }

    model.optimize();

    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
      double objval = model.get(GRB_DoubleAttr_ObjVal);
      std::cout << "Optimal objective: " << objval << std::endl;
    } 

    // To store variant positions retained
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

  std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - tStart);
  std::cout<< "INFO, VF::main, time taken by variant selection algorithm = " << wctduration.count() << " seconds" << "\n"; 

  std::cout<< "INFO, VF::main, count of variant containing positions retained = " << std::count(R.begin(), R.end(), true) << "\n";

  int count_variants_retained=0;
  for (std::size_t i = 0; i < n; i++) if(R[i]) count_variants_retained += c[i]; 
  std::cout<< "INFO, VF::main, count of variants retained = " << count_variants_retained << "\n";
  printVariantGapStats (R, svpos_u);

  return 0;
}
