#include <cassert>
#include "ext/clipp.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

struct Parameters
{
  int alpha;
  int delta;
  std::string vcffile;
  std::string chr;
  std::string prefix;
  bool pos;
};

inline bool exists (const std::string& filename) {
  std::ifstream f(filename.c_str());
  return f.good();
}

/**
 * @brief  parse and print command line arguments
 */
void parseandSave(int argc, char** argv, Parameters &param)
{
    param.pos = false; //default

  //define all arguments
  auto cli =
    (
     clipp::required("-a") & clipp::value("alpha", param.alpha).doc("path length in variation graph (e.g., 500)"),
     clipp::required("-d") & clipp::value("delta", param.delta).doc("differences allowed (e.g., 10)"),
     clipp::required("-vcf") & clipp::value("file1", param.vcffile).doc("uncompressed vcf file (something.vcf)"),
     clipp::required("-chr") & clipp::value("id", param.chr).doc("chromosome id (e.g., 1 or chr1), make it consistent with vcf file"),
     clipp::option("-prefix") & clipp::value("file2", param.prefix).doc("filename to optionally save input and output variants")
    );

  if(!clipp::parse(argc, argv, cli))
  {
    //print help page
    clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0])) << std::endl;
    exit(1);
  }

  std::cout << "INFO, VF::parseandSave, alpha = " << param.alpha << std::endl;
  std::cout << "INFO, VF::parseandSave, delta = " << param.delta << std::endl;
  std::cout << "INFO, VF::parseandSave, vcf file = " << param.vcffile << std::endl;
  std::cout << "INFO, VF::parseandSave, chromosome id = " << param.chr << std::endl;
  if (param.prefix.length() > 0) std::cout << "INFO, VF::parseandSave, prefix = " << param.prefix << std::endl;

  if (! exists(param.vcffile))
  {
    std::cerr << "ERROR, VF::parseandSave, vcf file cannot be opened" << std::endl;
    exit(1);
  }
}

/**
 * @brief  parse and print command line arguments (modified for ILP)
 */
void parseandSave_ILP(int argc, char** argv, Parameters &param)
{
    param.pos = false; //default

  //define all arguments
  auto cli =
    (
     clipp::required("-a") & clipp::value("alpha", param.alpha).doc("path length in variation graph (e.g., 500)"),
     clipp::required("-d") & clipp::value("delta", param.delta).doc("differences allowed (e.g., 10)"),
     clipp::required("-vcf") & clipp::value("file1", param.vcffile).doc("uncompressed vcf file (something.vcf)"),
     clipp::required("-chr") & clipp::value("id", param.chr).doc("chromosome id (e.g., 1 or chr1), make it consistent with vcf file"),
     clipp::option("-prefix") & clipp::value("file2", param.prefix).doc("filename to optionally save input and output variants"),
     clipp::option("--pos").set(param.pos).doc("set objective to minimize variation positions rather than variant count")
    );

  if(!clipp::parse(argc, argv, cli))
  {
    //print help page
    clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0])) << std::endl;
    exit(1);
  }

  std::cout << "INFO, VF::parseandSave, alpha = " << param.alpha << std::endl;
  std::cout << "INFO, VF::parseandSave, delta = " << param.delta << std::endl;
  std::cout << "INFO, VF::parseandSave, vcf file = " << param.vcffile << std::endl;
  std::cout << "INFO, VF::parseandSave, chromosome id = " << param.chr << std::endl;
  if (param.prefix.length() > 0) std::cout << "INFO, VF::parseandSave, prefix = " << param.prefix << std::endl;

  if (! exists(param.vcffile))
  {
    std::cerr << "ERROR, VF::parseandSave, vcf file cannot be opened" << std::endl;
    exit(1);
  }
}

/**
 * @brief   VCFtools rarely reports multiple SNP entries with same pos,
 *          here we remove the duplicate entries
 */
void ignoreDuplicateSNPrecords(std::vector<int> &p, std::vector<int> &c)
{
    assert (p.size() == c.size());
    assert (p.size() > 0);
    assert (std::is_sorted(p.begin(), p.end()));

    std::vector<int> p_new, c_new;

    p_new.push_back(p[0]);
    c_new.push_back(c[0]);

    for (int i = 1; i < p.size(); i++)
    {
        if (p[i-1] != p[i])
        {
            p_new.push_back(p[i]);
            c_new.push_back(c[i]);
        }
    }
    p = p_new; c=c_new;
}

/**
 * @brief   measure length of string between two consecutive variant loci before
 *          and after graph reduction
 */
void printVariantGapStats (const std::vector<bool> &retained, const std::vector<int> &pos)
{
    assert (retained.size() == pos.size());
    assert (pos.size() > 0);
    assert (std::is_sorted(pos.begin(), pos.end()));

    std::vector<int> gap_length_original;
    std::vector<int> gap_length_reduced;

    for (int i = 0; i < pos.size() - 1; i++)
    {
        assert(pos[i+1] != pos[i]);
        gap_length_original.push_back(pos[i+1] - pos[i] - 1);
    }

    std::vector<int> new_pos;
    for (int i = 0; i < pos.size(); i++)
        if (retained[i]) new_pos.push_back(pos[i]);

    for (int i = 0; i < new_pos.size() - 1; i++)
    {
        assert(new_pos[i+1] != new_pos[i]);
        gap_length_reduced.push_back(new_pos[i+1] - new_pos[i] - 1);
    }

    int max1 = *std::max_element(gap_length_original.begin(), gap_length_original.end());
    int min1 = *std::min_element(gap_length_original.begin(), gap_length_original.end());
    int avg1 = std::accumulate(gap_length_original.begin(), gap_length_original.end(), 0) / gap_length_original.size();

    int max2 = *std::max_element(gap_length_reduced.begin(), gap_length_reduced.end());
    int min2 = *std::min_element(gap_length_reduced.begin(), gap_length_reduced.end());
    int avg2 = std::accumulate(gap_length_reduced.begin(), gap_length_reduced.end(), 0) / gap_length_reduced.size();

    std::cout<< "INFO, VF::printVariantGapStats, before: (min, mean, max) = (" << min1 << ", " << avg1 << ", " << max1 << ")\n";
    std::cout<< "INFO, VF::printVariantGapStats, after: (min, mean, max) = (" << min2 << ", " << avg2 << ", " << max2 << ")\n";
}

void print_SV_vcf (const std::vector<bool> &retained, const std::vector<int> &pos, const Parameters &param)
{
  srand(time(0)+1); int random = rand() % 100000;  
  std::string tmp_file1 = ".VF." + std::to_string(random) + ".txt";
  random = rand() % 100000;  
  std::string tmp_file2 = ".VF." + std::to_string(random) + ".txt";

  std::ofstream myfile1 (tmp_file1);
  for (std::size_t i = 0; i < pos.size(); i++) if(retained[i]) myfile1 << pos[i] << "\n";
  myfile1.close();
  std::cout << "INFO, VF::print_SV_vcf, written retained variant loci to " << tmp_file1 << "\n";

  std::string cmd = "cat " + param.vcffile + " | grep '^#' > " + param.prefix + ".inputrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());
  cmd = "cat " + param.vcffile + " | grep -vE '^#' | grep 'INS\\|DEL' | awk -v chr=" + param.chr + " '$1 == chr {print $0}' >> " + param.prefix + ".inputrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());

  cmd = "cat " + param.prefix + ".inputrecords.vcf | grep '^#' > " + param.prefix + ".retainedrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());
  cmd = "cat " + param.prefix + ".inputrecords.vcf | grep -vE '^#' > " + tmp_file2; std::cout << cmd << "\n"; std::system(cmd.c_str());
  cmd = "awk 'NR == FNR { a[$0]; next } $2 in a' " + tmp_file1 + " " + tmp_file2 + " >> " + param.prefix + ".retainedrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());

  cmd = "rm -f " + tmp_file1 + "*"; std::cout << cmd << "\n"; std::system(cmd.c_str()); //delete tmp file
  cmd = "rm -f " + tmp_file2 + "*"; std::cout << cmd << "\n"; std::system(cmd.c_str()); //delete tmp file
}

void print_snp_vcf (const std::vector<bool> &retained, const std::vector<int> &pos, const Parameters &param)
{
  srand(time(0)+1); int random = rand() % 100000;  
  std::string tmp_file1 = ".VF." + std::to_string(random) + ".txt";
  random = rand() % 100000;  
  std::string tmp_file2 = ".VF." + std::to_string(random) + ".txt";

  std::ofstream myfile1 (tmp_file1);
  for (std::size_t i = 0; i < pos.size(); i++) if(retained[i]) myfile1 << pos[i] << "\n";
  myfile1.close();
  std::cout << "INFO, VF::print_snp_vcf, written retained variant loci to " << tmp_file1 << "\n";

  std::string cmd = "cat " + param.vcffile + " | grep '^#' > " + param.prefix + ".inputrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());
  cmd = "cat " + param.vcffile + " | grep -vE '^#' | grep 'S' | awk -v chr=" + param.chr + " '$1 == chr {print $0}' >> " + param.prefix + ".inputrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());


  cmd = "cat " + param.prefix + ".inputrecords.vcf | grep '^#' > " + param.prefix + ".retainedrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());
  cmd = "cat " + param.prefix + ".inputrecords.vcf | grep -vE '^#' > " + tmp_file2; std::cout << cmd << "\n"; std::system(cmd.c_str());
  cmd = "awk 'NR == FNR { a[$0]; next } $2 in a' " + tmp_file1 + " " + tmp_file2 + " >> " + param.prefix + ".retainedrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());

  cmd = "rm -f " + tmp_file1 + "*"; std::cout << cmd << "\n"; std::system(cmd.c_str()); //delete tmp file
  cmd = "rm -f " + tmp_file2 + "*"; std::cout << cmd << "\n"; std::system(cmd.c_str()); //delete tmp file
}

void print_snp_indel_vcf (const std::vector<bool> &retained, const std::vector<int> &pos, const Parameters &param)
{
  srand(time(0)+1); int random = rand() % 100000;  
  std::string tmp_file1 = ".VF." + std::to_string(random) + ".txt";
  random = rand() % 100000;  
  std::string tmp_file2 = ".VF." + std::to_string(random) + ".txt";

  std::ofstream myfile1 (tmp_file1);
  for (std::size_t i = 0; i < pos.size(); i++) if(retained[i]) myfile1 << pos[i] << "\n";
  myfile1.close();
  std::cout << "INFO, VF::print_snp_indel_vcf, written retained variant loci to " << tmp_file1 << "\n";

  std::string cmd = "cat " + param.vcffile + " | grep '^#' > " + param.prefix + ".inputrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());
  cmd = "cat " + param.vcffile + " | grep -vE '^#' | grep 'S\\|INS\\|DEL' | awk -v chr=" + param.chr + " '$1 == chr {print $0}' >> " + param.prefix + ".inputrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());


  cmd = "cat " + param.prefix + ".inputrecords.vcf | grep '^#' > " + param.prefix + ".retainedrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());
  cmd = "cat " + param.prefix + ".inputrecords.vcf | grep -vE '^#' > " + tmp_file2; std::cout << cmd << "\n"; std::system(cmd.c_str());
  cmd = "awk 'NR == FNR { a[$0]; next } $2 in a' " + tmp_file1 + " " + tmp_file2 + " >> " + param.prefix + ".retainedrecords.vcf"; std::cout << cmd << "\n"; std::system(cmd.c_str());

  cmd = "rm -f " + tmp_file1 + "*"; std::cout << cmd << "\n"; std::system(cmd.c_str()); //delete tmp file
  cmd = "rm -f " + tmp_file2 + "*"; std::cout << cmd << "\n"; std::system(cmd.c_str()); //delete tmp file
}
