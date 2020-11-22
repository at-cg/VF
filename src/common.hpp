#include "ext/clipp.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

struct Parameters
{
  int alpha;
  int delta;
  std::string vcffile;
  std::string chr;
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
  //define all arguments
  auto cli = 
    (
     clipp::required("-a") & clipp::value("alpha", param.alpha).doc("path length in variation graph (e.g., 500)"),
     clipp::required("-d") & clipp::value("delta", param.delta).doc("differences allowed (e.g., 10)"),
     clipp::required("-vcf") & clipp::value("file", param.vcffile).doc("uncompressed vcf file (something.vcf)"),
     clipp::required("-chr") & clipp::value("id", param.chr).doc("chromosome id (e.g., 1 or chr1), make it consistent with vcf file")
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

  if (! exists(param.vcffile))
  {
    std::cerr << "ERROR, VF::parseandSave, vcf file cannot be opened" << std::endl;
    exit(1);
  }
}

