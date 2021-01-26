Examples
========================================================================
## Download human genome variant files

Download [1KG Phase 3 VCF files](https://www.internationalgenome.org/data) from [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). Here we will use mitochondrial chromosome as one of the examples.
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
gzip -d ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
mv ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf 1KG_chrMT.vcf 
```
Also, download SV calls from [Audano et al.](https://doi.org/10.1016/j.cell.2018.12.019)
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181025_EEE_SV-Pop_1/VariantCalls_EEE_SV-Pop_1/EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz
gzip -d EEE_SV-Pop_1.ALL.sites.20181204.vcf.gz
mv EEE_SV-Pop_1.ALL.sites.20181204.vcf SV.vcf
```

## Graph reduction (only SNPs)

Assuming you have four executables `greedy_snp`, `lp_snp`, `greedy_sv_indels` and `ilp_sv_indels` in build directory, you can use them in the following ways.
* Minimise positions containing variants in variation graph
```
$ greedy_snp -a 1000 -d 100 -vcf 1KG_chrMT.vcf -chr MT
INFO, VF::parseandSave, alpha = 1000
INFO, VF::parseandSave, delta = 100
INFO, VF::parseandSave, vcf file = 1KG_chrMT.vcf
INFO, VF::parseandSave, chromosome id = MT
INFO, VF::main, extracting SNPs from vcf file using command = vcftools-0.1.16/bin/vcftools --vcf 1KG_chrMT.vcf --chr MT --counts --remove-indels --out .VF.683.txt 2>/dev/null
INFO, VF::main, vcftools finished
INFO, VF::main, count of variant containing positions = 3771
INFO, VF::main, count of variants = 3963
INFO, VF::main, starting timer
INFO, VF::main, time taken by variant selection algorithm = 5.0772e-05 seconds
INFO, VF::main, count of variant containing positions retained = 2132
INFO, VF::main, count of variants retained = 2232
INFO, VF::printVariantGapStats, before: (min, mean, max) = (0, 3, 49)
INFO, VF::printVariantGapStats, after: (min, mean, max) = (0, 6, 814)
```
* Minimise count of variants in variation graph
```
$ lp_snp -a 1000 -d 100 -vcf 1KG_chrMT.vcf -chr MT
INFO, VF::parseandSave, alpha = 1000
INFO, VF::parseandSave, delta = 100
INFO, VF::parseandSave, vcf file = 1KG_chrMT.vcf
INFO, VF::parseandSave, chromosome id = MT
INFO, VF::main, extracting SNPs from vcf file using command = vcftools-0.1.16/bin/vcftools --vcf 1KG_chrMT.vcf --chr MT --counts --remove-indels --out .VF.2546.txt 2>/dev/null
INFO, VF::main, vcftools finished
INFO, VF::main, count of variant containing positions = 3771
INFO, VF::main, count of variants = 3963
INFO, VF::main, starting timer
INFO, VF::main, Gurobi solver starting
Optimal objective: 1831
INFO, VF::main, time taken by variant selection algorithm = 0.667862 seconds
INFO, VF::main, count of variant containing positions retained = 2132
INFO, VF::main, count of variants retained = 2132
INFO, VF::printVariantGapStats, before: (min, mean, max) = (0, 3, 49)
INFO, VF::printVariantGapStats, after: (min, mean, max) = (0, 6, 350)
```

## Graph reduction (indel SVs)

Executables `greedy_sv_indels` and `ilp_sv_indels` are useful for graph reduction when variation graph has SVs. In theory, the algorithm would also work for small indels, but its VCF file IO is yet to be engineered to be able to parse them.
```
$ greedy_sv_indels -a 1000 -d 100 -vcf SV.vcf -chr chr21
INFO, VF::parseandSave, alpha = 1000
INFO, VF::parseandSave, delta = 100
INFO, VF::parseandSave, vcf file = SV.vcf
INFO, VF::parseandSave, chromosome id = chr21
INFO, VF::main, count of variant containing positions = 1724
INFO, VF::main, count of variants = 1772
INFO, VF::main, starting timer
INFO, VF::calculateLeftMostReachable, computing window ranges...
INFO, VF::calculateLeftMostReachable, done
INFO, VF::main, time taken by variant selection algorithm = 3.75185 seconds
INFO, VF::main, count of variant containing positions retained = 1310
INFO, VF::main, count of variants retained = 1356
INFO, VF::printVariantGapStats, before: (min, mean, max) = (0, 19547, 427236)
INFO, VF::printVariantGapStats, after: (min, mean, max) = (0, 25730, 541132)
```

```
$ ilp_sv_indels -a 1000 -d 100 -vcf SV.vcf -chr chr21
INFO, VF::parseandSave, alpha = 1000
INFO, VF::parseandSave, delta = 100
INFO, VF::parseandSave, vcf file = SV.vcf
INFO, VF::parseandSave, chromosome id = chr21
INFO, VF::main, count of variant containing positions = 1724
INFO, VF::main, count of variants = 1772
INFO, VF::main, starting timer
INFO, VF::calculateLeftMostReachable, computing window ranges...
INFO, VF::calculateLeftMostReachable, done
INFO, VF::main, Gurobi solver starting
INFO, VF::main, ILP solver will attempt to minimize count of variants
Optimal objective: 417
INFO, VF::main, time taken by variant selection algorithm = 3.88308 seconds
INFO, VF::main, count of variant containing positions retained = 1309
INFO, VF::main, count of variants retained = 1355
INFO, VF::printVariantGapStats, before: (min, mean, max) = (0, 19547, 427236)
INFO, VF::printVariantGapStats, after: (min, mean, max) = (0, 25749, 541132)
```

```
$ ilp_sv_indels -a 1000 -d 100 -vcf SV.vcf -chr chr21 --pos
INFO, VF::parseandSave, alpha = 1000
INFO, VF::parseandSave, delta = 100
INFO, VF::parseandSave, vcf file = SV.vcf
INFO, VF::parseandSave, chromosome id = chr21
INFO, VF::main, count of variant containing positions = 1724
INFO, VF::main, count of variants = 1772
INFO, VF::main, starting timer
INFO, VF::calculateLeftMostReachable, computing window ranges...
INFO, VF::calculateLeftMostReachable, done
INFO, VF::main, Gurobi solver starting
INFO, VF::main, ILP solver will attempt to minimize variant positions
Optimal objective: 415
INFO, VF::main, time taken by variant selection algorithm = 3.84347 seconds
INFO, VF::main, count of variant containing positions retained = 1309
INFO, VF::main, count of variants retained = 1357
INFO, VF::printVariantGapStats, before: (min, mean, max) = (0, 19547, 427236)
INFO, VF::printVariantGapStats, after: (min, mean, max) = (0, 25749, 541132)
```
