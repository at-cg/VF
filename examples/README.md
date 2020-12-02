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

Assuming you have three executables `greedy_snp`, `lp_snp`, `greedy_sv_indels` in build directory, you can use them in the following ways.
* Minimise positions containing variants in variation graph
```
$ greedy_snp -a 1000 -d 100 -vcf 1KG_chrMT.vcf -chr MT
INFO, VF::parseandSave, alpha = 1000
INFO, VF::parseandSave, delta = 100
INFO, VF::parseandSave, vcf file = 1KG_chrMT.vcf
INFO, VF::parseandSave, chromosome id = MT
INFO, VF::main, extracting SNPs from vcf file using command = VF/build/vcftools-0.1.16/bin/vcftools --vcf 1KG_chrMT.vcf --chr MT --counts --remove-indels --out .VF.7543.txt 2>/dev/null
INFO, VF::main, vcftools finished
INFO, VF::main, count of variant containing positions = 3771
INFO, VF::main, count of variants = 3963
INFO, VF::main, starting timer
INFO, VF::main, time taken by variant selection algorithm = 0.000280911 seconds
INFO, VF::main, count of variant containing positions retained = 2132
INFO, VF::main, count of variants retained = 2232
```
* Minimise count of variants in variation graph
```
$ lp_snp -a 1000 -d 100 -vcf 1KG_chrMT.vcf -chr MT
INFO, VF::parseandSave, alpha = 1000
INFO, VF::parseandSave, delta = 100
INFO, VF::parseandSave, vcf file = 1KG_chrMT.vcf
INFO, VF::parseandSave, chromosome id = MT
INFO, VF::main, extracting SNPs from vcf file using command = /global/project/projectdirs/nguest/cjain7/projects/vf/software/VF/build/vcftools-0.1.16/bin/vcftools --vcf 1KG_chrMT.vcf --chr MT --counts --remove-indels --out .VF.5127.txt 2>/dev/null
INFO, VF::main, vcftools finished
INFO, VF::main, count of variant containing positions = 3771
INFO, VF::main, count of variants = 3963
INFO, VF::main, starting timer
INFO, VF::main, Gurobi solver starting
Academic license - for non-commercial use only
Optimal objective: 1831
INFO, VF::main, time taken by variant selection algorithm = 0.661689 seconds
INFO, VF::main, count of variant containing positions retained = 2132
INFO, VF::main, count of variants retained = 2132
```

## Graph reduction (indel SVs)

Executable `greedy_sv_indels` is useful for graph reduction when variation graph has SVs. In theory, the algorithm would also work for small indels, but its VCF file IO is yet to be engineered to be able to parse them.
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
INFO, VF::main, time taken by variant selection algorithm = 6.43583 seconds
INFO, VF::main, count of variant containing positions retained = 1310
INFO, VF::main, count of variants retained = 1356
```
