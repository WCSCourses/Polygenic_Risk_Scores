## BridgePRS

### Learning Objectives
The goal of this practical is to implement the basic steps needed to implement multi-ancestry PRS prediction models successfully, using the BridgePRS software.

By the end of this session participants will understand how to:
- Set up the configuration files used as input by the software.
- Edit and interact with BridgePRS via the command line. 
- Implement the 3 PRS models integral to the functionality of the software, using the integrative "easy run" function.
- Perform a single-ancestry version of the BridgePRS method.

### Introduction to BridgePRS
BridgePRS is a Bayesian PRS method that integrates trans-ancestry GWAS summary statistics. Unlike the fine-mapping approach of PRS-CSx, BridgePRS retains all variants within loci to best tag causal variants shared across ancestries. The focus is on correctly estimating causal effect sizes, which is key when the goal is prediction, rather than on estimating their location. BridgePRS is most applicable for combining information of a well-powered GWAS performed in a (discovery) population or populations unmatched to the ancestry of the target sample, with a second GWAS of relatively limited power in a (target) population that is matched to the ancestry of the target sample.

BridgePRS is first applied to a large discovery population GWAS, by running Bayesian ridge regression, at putative loci. The BridgePRS algorithm is trained according to 3 PRS models:
1. PRS run using only the target (Non-European) dataset (prs-single)
2. PRS run using SNP-weights calculated from the European Model (prs-port)
3. PRS run using a prior effect-size distribution from the European Model (prs-prior)
The 3 models are subsequently combined to produce a weighted PRS solution

As well as applying each of the steps in sequence, models can also be run separately according to the needs of individual users. Here, we use the single-ancestry BridgePRS approach to train polygenic scores for an African target sample applying the Bayesian ridge regression approach of BridgePRS to shrink the effect size of SNPs derived from an African ancestry GWAS towards their true underlying values.
<br><br> 
### BridgePRS Scenario 1: Application of African GWAS weights to an African target group

#### Create configuration file for the target-only analysis
In this example we will run BridgePRS across chromosomes 1 - 22. The first required step is to generate the configuration file needed to run the desired single ancestry BridgePRS analysis. The following command should be run from the main directory:
```
bridgePRS check pop -o out_config-AFR-single --pop AFR --sumstats_prefix data/pop_europe/sumstats/eur.chr --genotype_prefix data/pop_africa/genotypes/afr_genotypes --phenotype_file data/pop_africa/phenotypes/afr_pheno.dat
```
(BridgePRS produces on-screen information which tells you some of the tasks the software is doing behind the scenes).

#### Questions
From the on-screen output:
1. Do you notice anything interesting in the way BridgePRS handled the phenotype file?
2. How many phenotypes was BridgePRS able to identify in the phenotype file? in what way do these phenotypes differ?
3. Which resulting files have been generated in ./out/save/ ?

#### Note
In the next code snippet the file path of the newly created .config (configuration) file has been incorporated into the run command. This information was provided as part of the on-screen output in the previous step. 


### Single ancestry BridgePRS analysis:
Now that we have our African configuration file prepared we are ready to perform the single-ancestry BridgePRS analysis.
In this step we opt to select the continuous version of the trait for analysis.
```
bridgePRS prs-single run -o out_config-AFR-single --pop africa --config_files out_config-AFR-single/save/primary.AFR.config --phenotype y
```

#### Task
- Review the contents of the output directory  out_bridge-AFR-single/prs-single_AFRICA and subfolders

#### Questions
4. What evidence can you see that the analysis was successfully executed?

<br><br> 

### BridgePRS Scenario 2:  Prediction into African target data using European and African summary statistics
BridgePRS is most commonly used to combine the power of a smaller ancestry-matched GWAS with a much larger but genetically-distant GWAS population, for the purpose of maximising PRS prediction quality for under-served target populations.

#### Create configuration file for base and target populations
```
bridgePRS check pops -o out_config-EUR-AFR-easyrun --pop AFR EUR --sumstats_prefix data/pop_africa/sumstats/afr.chr data/pop_europe/sumstats/eur.chr --sumstats_suffix .glm.linear.gz .glm.linear.gz --genotype_prefix data/pop_africa/genotypes/afr_genotypes --phenotype_file data/pop_africa/phenotypes/afr_pheno.dat
```
#### Question
5. Carefully check the information given in the on-screen output: (i) How many different .config files have been produced
   and (ii) where are they located?

#### Tasks
- Incorporate the config path information into the code template provided below (for both European and African .config files).
- Based on your-recent understanding of genetic distances between continental populations, choose a sensible    
  value of the --fst parameter to reflect the genetic distance between Africans and Europeans. This extra information will inform the prior   distribution from which posterior effect weights for the target population will be calculated. This needs to be done before attempting to
  run the code given below.

### Multi-ancestry BRIDGEPRS analysis:
Add the missing peices of information to the code below, as you enter it into your terminal.
```
bridgePRS easyrun go -o out_easyrun-EUR-AFR --config_files target.AFR.config base.EUR.config --fst --phenotype y
```

#### Tasks
- After running the above code, navigate to the output directory: ./out_config-EUR-AFR-easyrun to inspect the results.
- Open either of the 2 plots that you see in the directory.

#### Questions
6. In the summary plot, which set of values expresses the correlation between the weights calculated by BridgePRS and
   the beta weights from the initial GWASs?
8. In which output directory will you find precise values of variance explained by the prs-combined-AFR model?
9. What is the variance explained by the prs-combined-AFR model?

### Short Quiz
I have GWAS data and genotype/phenotype data for a cohort consisting of >2000 samples from a small East European population. 
The population LD structure is unique and so I would like this information to be incorporated into my PRS prediction model. I additionally have GWAS and genotype/phenotype data from the UKB biobank that I want to include. How should I formulate the relevant Config files for the BridgePRS analysis?

#### File types
```
ukr/sumstats/ukr.sumstats.out
ukr/genotypes/chr1.bed,bim,fam...chr22.bed,bin,fam,
ukr/phenotypes/ukr_test.dat, ukr/phenotypes/ukr_validation.dat
ukb/sumstats/ukb.sumstats.gz
ukb/genotypes/chr1.bed,bim,fam...chr22.bed,bin,fam
ukb/phenotypes/ukb_test.dat, ukb/phenotypes/ukb_validation.dat
```
#### Answer
```
Target config:
POP=
LDPOP=
SUMSTATS_PREFIX=
GENOTYPE_PREFIX=
PHENOTYPE_FILE=
VALIDATION_FILE=

Model config file:
POP=
LDPOP=
SUMSTATS_PREFIX=
GENOTYPE_PREFIX=
PHENOTYPE_FILE=
VALIDATION_FILE=
```
