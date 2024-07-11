## Day 3b practical
We need to move into the directory you will be working in;

```sh
cd ~/data/Data_Day4/data
```

## Introduction to Cross-Ancestry PRS computation
Before starting the practical, the following commands will need to be run from within your virtual machine. 
These commands set up an 'environment' that allows you to work outside the virtual machine, which has memory restrictions.

Set up the environment using conda:
```sh
conda create -n "PRScsx" python=3.7
conda activate PRScsx
pip install scipy
pip install h5py
```

The aim of this practical is to provide you with a basic understanding and some experience of running PRS-CSx software. After completing this practical, you should:

* Be able to perform cross-population analyses.
* Be familiar with running cross-ancestry PRS analyses using PRS-CSx.
* Understand how to evaluate linear models using Akaike’s Information Criterion.

#### 1. The 1000 Genomes datasets
The data we will be working with comes from the 1000 Genomes Project reference panel. The data relates to individuals from 26 different source populations around the world. For simplicity, the populations have been collapsed into 5 broader continental super-populations: East Asian, European, South Asian, Amerindian, African ((EAS, EUR, SAS, EUR and AFR)).

The scripts used to download and process the 1000Genomes data for the purposes of this course will be provided in the course appendix at the end of this week. 

#### 2. Cross-population allele frequency
Genetic variation is conveyed using allelic frequencies. Allele frequency is shaped by evolutionary forces and drift.  Here we compare profiles of allele frequency across the five ancestral populations. Global differences in allelic frequency has important implications for the portability of PRS across populations.

Using plink it is possible to generate allele frequency statistics for each SNP, across populations, using the annotations provided in the **file pop_info.pheno**. In _/home/manager/data/Data_Day4_:

```sh
cd ../
./software/plink_linux --bfile ./data/chr1-22 --freq --within ./data/pop_info.pheno
```
Population-stratified allele frequencies are reported in the output file **plink.frq.strat.**
For each population, print the numbers of total SNPs to screen, as follows:

**AFR**
```sh
grep -F 'AFR' plink.frq.strat | wc -l
```
From there we can print the number of SNPs with minor allele frequencies greater than 0 (and are hence potentially available for genomic analyes).

```sh
grep -F 'AFR' plink.frq.strat | awk '$6 >0' | wc -l
```

**EUR**
```sh
grep -F 'EUR' plink.frq.strat | wc -l
```
Number of SNPs with MAF > 0 in EUR.
```sh
grep -F 'EUR' plink.frq.strat | awk '$6 >0' | wc -l
```

**EAS**
```sh
grep -F 'EAS' plink.frq.strat | wc -l
```
Number of SNPs with MAF > 0 in EAS.
```sh
grep -F 'EAS' plink.frq.strat | awk '$6 >0' | wc -l
```

**SAS**
```sh
grep -F 'SAS' plink.frq.strat | wc -l
```
Number of SNPs with MAF > 0 in SAS.
```sh
grep -F 'SAS' plink.frq.strat | awk '$6 >0' | wc -l
```

**AFR**
```sh
grep -F 'AFR' plink.frq.strat | wc -l
```
Number of SNPs with MAF > 0 in AFR.
```sh
grep -F 'AFR' plink.frq.strat | awk '$6 >0' | wc -l
```

#### **Questions**
##### (i) Which population contains the most SNPs?
##### (ii) What  is the significance of the observed population order?  
&nbsp;

#### 3. Distribution of allele frequencies

In this exercise, we will analyse and compare the distribution of allele frequencies across different ancestries. 
We will use R for visualisation.

**Open a new terminal and open R**

Open a new tab in the terminal (**plus icon in the top left corner**)

In this new terminal window, make sure you are in the correct directory:

```sh
cd /home/manager/data/Data_Day4/out/
```
Now open R: 

```sh
R
```
Generate the plot:

```sh
# Install the necessary libraries
install.packages("dplyr")
install.packages("ggplot2")

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Create a function to read the files and add ancestry information
freq <-read.table("~/data/Data_Day4/plink.frq.strat", header =T)
plotDat <- freq %>%
  mutate(AlleleFrequency = cut(MAF, seq(0, 1, 0.25))) %>%
  group_by(AlleleFrequency, CLST) %>%
  summarise(FractionOfSNPs = n()/nrow(freq) * 100)

# Create a bar graph 
png('/home/manager/data/Data_Day4/out/MAF_ancestry_analysis.png', unit='px', res=300, width=3500, height=4500)

maf_ancestry <- ggplot(na.omit(plotDat),
  aes(AlleleFrequency, FractionOfSNPs, group = CLST, col = CLST)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 12)) +
  ggtitle("Distribution of allele frequency across genome")

print(maf_ancestry)

dev.off()
```

Examine the plot the MAF across each ancestry: 

```sh
# This does not work when you are in R, ensure you out of R before running:
xdg-open MAF_ancestry_analysis.png
```

#### **Questions**
##### (i) Which population has the most SNPs?
##### (ii) What  is the significance of the observed population ordering?
##### (iii) What is the reason behind these two features?


## Introduction to PRS-CSx
#### 5. Background to PRS-CSX
PRS-CSx is a Python based command line tool that integrates GWAS summary statistics and external LD reference panels from multiple populations to improve cross-population polygenic prediction. We will be using simulated trait data pertaininng to systolic blood pressure (SBP) to explore PRS performance using 2 target populations that consist of 650 individuals of African ancestry and 500 individuals of European ancestry. Please note when running PRSice that the object of the flag "--prsice" will change according to whether plink is being called within the linux-like environment of the virtual machine (PRSice_linux) or a mac (PRSice_mac). Both executables can be found in the _/home/manager/data/Data_Day4_ directory. 


#### 7. Running PRS-CSx
To model the coupling of ect sizes at individual SNPs across ancestries PRS-CSx uses an MCMC (Bayesian) sampling algorithm to determine values of the global shrinkage parameter ("phi") by Maximum likelihood. For samples of mixed or admixed genetic ancestry (which ours are not) the optimal value of the shrinkage parameter is estimated autonomously from the data. Here we use the value of phi (1e-4), which is suggested by the software authors, given that our trait is simulated to have a relatively small number (N=110) causal variants, distributed genome-wide.
  
**Step 1: Set up environment**
------------------------------

First change to the working directory with the data for this practical 
```sh
cd /home/manager/data/Data_Day4/data
```
Make a directory called **hm3_by_ancestry** within the data folder, and move a folder back out of the data folder

```sh
mkdir hm3_by_ancestry
cd ..
```

**AFR**
```
for chr in {21..22}; do \
/home/manager/data/Data_Day4/software/plink_linux \
	--bfile /home/manager/data/Data_Day4/data/AFR_1kg.hm3.only.csx \
	--chr $chr \
	--make-bed \
	--out /home/manager/data/Data_Day4/data/hm3_by_ancestry/AFR_1kg.hm3.chr${chr}_only.csx;
done
```

**EUR**
```
for chr in {21..22}; do \
/home/manager/data/Data_Day4/software/plink_linux \
	--bfile /home/manager/data/Data_Day4/data/EUR_1kg.hm3.only.csx \
	--chr $chr \
	--make-bed \
	--out /home/manager/data/Data_Day4/data/hm3_by_ancestry/EUR_1kg.hm3.chr${chr}_only.csx;
done
```

**Set up the necessary environment variables for threading and verify they are set correctly.**
```
export N_THREADS=2
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS
```
**Verify the variables are set**
```
echo $N_THREADS
echo $MKL_NUM_THREADS
echo $NUMEXPR_NUM_THREADS
echo $OMP_NUM_THREADS
```
<br>

Step 2: Run CSX. Derive new SNPs weights trained on European and African summary stats
--------------------------------------------------------------------------------------

**Generate job file containing the threaded PRScsx commands.**
First, to minimize computational resources and time, we should create a script to run the tasks in parallel. Make a script called **create_multithread.sh**

```sh
nano create_multithread.sh
```
Then copy and paste the code below into that script. After save, then close the script ctrl + x  

```sh
#!/bin/bash

# Create the script file
SCRIPT_FILE="multithread.sh"

# Write the header of the script file
echo "#!/bin/bash" > $SCRIPT_FILE



for chr in {21..22}; do
  echo "python3 /home/manager/data/Data_Day4/software/PRScsx.py \
        --ref_dir=/home/manager/data/Data_Day4/reference/csx \
        --bim_prefix=/home/manager/data/Data_Day4/data/hm3_by_ancestry/AFR_1kg.hm3.chr${chr}_only.csx \
        --sst_file=/home/manager/data/Data_Day4/data/3b/data/sumstats_by_chr/EUR-SBP-simulated.sumstats.chr${chr},/home/manager/data/Data_Day4/data/3b/data/sumstats_by_chr/AFR-SBP-simulated.sumstats.chr${chr} \
        --n_gwas=25732,4855 \
        --chrom=${chr} \
        --n_iter=1000 \
        --n_burnin=500 \
        --thin=5 \
        --pop=EUR,AFR \
        --phi=1e-4 \
        --out_dir=/home/manager/data/Data_Day4/out/csx \
        --out_name=afr.target_chr${chr}.csx" >> $SCRIPT_FILE
done
```
Run: Ctrl + O, followed by Enter '**to save**'
Run: Ctrl + X, to Exit

You will need to change permission for the script to be able to execute

```sh
chmod +x create_multithread.sh
```
Run the builder
```sh
./create_multithread.sh
```
To be able to run the next command you first install 'parallel' 

```sh
sudo apt install parallel
```
**Run the Job File with GNU Parallel:** (May take a while)
```
parallel --verbose --jobs $N_THREADS < multithread.sh
```

<br>

Step 3: Combine CSX-derived SNP weights across chromosomes (Currently Excludes Chromosome 3)
--------------------------------------------------------------------------------------------

**Load R and the necessary library**
```
R
```
Call in the package
```
library(dplyr)
```
**Define the path to the directory containing the PRS-CSX output files**
```
path <- "/home/manager/data/Data_Day4/out/csx"
```
**Define the ancestry you want to combine ("EUR" or "AFR")**
```
ancestry <- "EUR"
```
**Initialize an empty data frame to store the combined data**
```
combined_data <- data.frame()
```
**Loop through chromosomes 21 to 22, (currently excluding chromosome 3)**
```
for (chr in setdiff(21:22, 3)) {
  # Construct the file name
  file_name <- paste0("afr.target_chr", chr, ".csx_", ancestry, "_pst__a1_b0.5_phi1e-04_chr", chr, ".txt")
  file_path <- file.path(path, file_name)
  
  # Check if file exists before reading
  if (file.exists(file_path)) {
	# Read the data from the file
	data <- read.table(file_path, header = FALSE, sep = "\t", col.names = c("CHR", "rsid", "pos", "ref", "alt", "beta"))
	
	# Combine the data
	combined_data <- rbind(combined_data, data)
  } else {
	warning(paste("File not found:", file_path))
  }
}

```
**Write the combined data to a new file**
```
output_file <- file.path(path, paste0("combined_", ancestry, "_pst_.txt"))
write.table(combined_data, output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
```
### Task: Replace 'ancestry <- "EUR" ' with 'ancestry <- "AFR" ' and repeat the subsequent steps shown above
<br>

# Step 4: Merge genotype-phenotype data
-------------------------------------

**Prepare data**
The data is slow to merge unless you split the input bim file into just chr21 and chr22 (Q- why are those faster?)

```
cd /home/manager/data/Data_Day4/data/3b/data/
plink --bfile AFR_1kg.hm3.only.csx --chr 21 22 --make-bed --out AFR_1kg.hm3.only.csx_21_22
```
Start R with sudo rights to allow you to install
```
sudo R 
```

```
# Load libraries. For any unavailable package, install it with **install.packages("_package_name")**
install.packages("BiocManager")
BiocManager::install("snpStats")

library(data.table)
library(ggplot2)
library(snpStats)

# Define the path to the directory containing the PLINK files and phenotypic data
plink_path <- "/home/manager/data/Data_Day4/data/3b/data/"

# Read PLINK files and phenotype data into R
bim_file <- file.path(plink_path, "AFR_1kg.hm3.only.csx_21_22.bim")
fam_file <- file.path(plink_path, "AFR_1kg.hm3.only.csx_21_22.fam")
bed_file <- file.path(plink_path, "AFR_1kg.hm3.only.csx_21_22.bed")
pheno_file <- file.path(plink_path, "sbp_afr_1kg.sim_pheno")

# Read the genotype data using snpStats
geno_data <- read.plink(bed_file, bim_file, fam_file)

# Extract SNP IDs from geno_data$map$snp.name
snp_ids <- geno_data$map$snp.name
if (is.null(snp_ids) || !is.character(snp_ids)) {
  stop("SNP IDs are missing or not in the correct format.")
}

# Convert the genotype data to a matrix and then to a data.table
geno_matrix <- as(geno_data$genotypes, "matrix")
geno_df <- data.table(geno_matrix)
setnames(geno_df, snp_ids)

# Add IID column from the fam file
geno_df[, IID := geno_data$fam$member]

# Read the phenotypic data**
pheno_data <- fread(pheno_file, sep = " ", header = TRUE)
```

**Merge phenotype and genotype data**
```
# Do merge
combined_data <- merge(pheno_data, geno_df, by = "IID")

# Keep only genotype columns
geno <- combined_data[, !names(combined_data) %in% c("FID", "IID", "pheno100", "pheno20", "pheno33", "pheno10"), with = FALSE]
phen <- combined_data$pheno100

# Convert geno to numeric matrix**
geno <- as.matrix(geno)
mode(geno) <- "numeric"

# Convert phen to vector
phen <- as.vector(phen)
```
<br>

Step 5: Split data into validation and test sets 
------------------------------------------------
**Specify Proportion**
```
# Here we specify that 40% of all IDs will be used to construct the validation group
set.seed(154)
vali_proportion <- 0.4
vali_size <- round(nrow(geno) * vali_proportion)
vali_indices <- sample(1:nrow(geno), vali_size, replace = FALSE)
test_indices <- setdiff(1:nrow(geno), vali_indices)
```
**Subsetting of individuals**
```
X_vali <- geno[vali_indices, , drop=FALSE]
y_vali <- phen[vali_indices]
X_test <- geno[test_indices, , drop=FALSE]
y_test <- phen[test_indices]
```
<br>

Step 6: Prepare the regression model input using the CSX-derived AFR and EUR weights
-------------------------------------------------------------------------------------
```
# Read the merged CSX output files
AFR_betas <- fread(file.path("/home/manager/data/Data_Day4/out/csx/combined_AFR_pst_eff.txt"), sep = "\t", header = TRUE)
EUR_betas <- fread(file.path("/home/manager/data/Data_Day4/out/csx/combined_EUR_pst_eff.txt"), sep = "\t", header = TRUE)

# Assuming the beta files have columns: "CHR", "rsid", "pos", "ref", "alt", "beta"
overlap_prs <- merge(AFR_betas, EUR_betas, by = "rsid", suffixes = c("_afr", "_eur"))

# Filter overlap_prs to include only SNPs present in X_vali
overlap_prs <- overlap_prs[rsid %in% colnames(X_vali)]

# Ensure that X_vali and X_test only contain SNPs present in W_afr and W_eur
common_snps <- intersect(colnames(X_vali), overlap_prs$rsid)
X_vali <- X_vali[, common_snps, drop=FALSE]
X_test <- X_test[, common_snps, drop=FALSE]

# Reorder the columns of X_vali and X_test to match the order of SNPs in overlap_prs
X_vali <- X_vali[, match(overlap_prs$rsid, colnames(X_vali)), drop=FALSE]
X_test <- X_test[, match(overlap_prs$rsid, colnames(X_test)), drop=FALSE]

# Extract the overlapping rsid and their corresponding betas
W_afr <- overlap_prs$beta_afr
W_eur <- overlap_prs$beta_eur

# Convert W_afr and W_eur to numeric vectors
W_afr <- as.numeric(W_afr)
W_eur <- as.numeric(W_eur)
```
<br>

Step 7: Prepare the variant weights matrices as vectors
-------------------------------------------------------
```
# Pre-check the alignment between the different objects
if (ncol(X_vali) != length(W_afr) || ncol(X_vali) != length(W_eur)) {
  stop("Dimensions of X_vali and W_afr/W_eur do not match.")
}

# In the validation sample: 

# (i) Compute XWafr_vali
XWafr_vali <- X_vali %*% W_afr
# (ii) Convert XWafr_vali to have zero mean and unit variance
XWafr_vali_z <- scale(XWafr_vali)

# (iii) Compute XWeur_vali
XWeur_vali <- X_vali %*% W_eur
# (iv) Convert XWeur_vali to have zero mean and unit variance
XWeur_vali_z <- scale(XWeur_vali)

# (v) Combine the normalized matrices
XW_vali <- cbind(XWafr_vali_z, XWeur_vali_z)

# Fit the model
model <- lm(scale(y_vali) ~ XWafr_vali_z + XWeur_vali_z - 1) # '- 1' removes the intercept

# Obtain the regression parameters
a_hat <- coef(model)[1]
b_hat <- coef(model)[2]
print(paste("a_hat =", a_hat))
print(paste("b_hat =", b_hat))
```
<br>

Step 8: Predict phenotype on validation and test dataset
--------------------------------------------------------
**Generate a linear combination of AFR and EUR PRSs for each individual**
```
# Each ancestry component is weighted by the regression coicient of that ancestry, in the preceding step
y_hat_vali <- a_hat * XWafr_vali_z + b_hat * XWeur_vali_z

# In the test sample: 

# Compute XWafr_test and XWeur_test
XWafr_test <- X_test %*% W_afr
XWafr_test_z <- scale(XWafr_test)

XWeur_test <- X_test %*% W_eur
XWeur_test_z <- scale(XWeur_test)

# y_hat in the test sample
y_hat <- a_hat * XWafr_test_z + b_hat * XWeur_test_z
```
<br>

Step 9: Plot phenotype distributions of validation and test data:
---------------------------------------------------------------------------
**Check that both distributions are approximately normal**
```
library(ggplot2)

# Create data frames for validation and test sets
vali_data <- data.frame(trait = y_vali, dataset = "Validation")
test_data <- data.frame(trait = y_test, dataset = "Test")

# Combine both data frames
combined_data <- rbind(vali_data, test_data)

# Plot the distributions
ggplot(combined_data, aes(x = trait, fill = dataset)) +
  geom_density(alpha = 0.5) +
  labs(title = "Trait Distributions for Validation and Test Sets", x = "Trait Value", y = "Density") +
  theme_minimal()
```
<br>

Step 10: Plot true values against predicted values
--------------------------------------------------
**The next steps use standard normal phenotype data to reduce scale differences between PRS and trait values**
```
min_true <- min(min(y_vali), min(y_test))
max_true <- max(max(y_vali), max(y_test))
min_pred <- min(min(y_hat_vali), min(y_hat))
max_pred <- max(max(y_hat_vali), max(y_hat))

pdf("true_against_pred.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

plot(scale(y_vali), y_hat_vali, pch = 19, col = rgb(0, 0, 0, 0.5), xlab = 'True Values', ylab = 'Predicted Values', main = 'Validation Dataset')
abline(0, 1, col = 'red', lty = 2)

plot(scale(y_test), y_hat, pch = 19, col = rgb(0, 0, 0, 0.5), xlab = 'True Values', ylab = 'Predicted Values', main = 'Test Dataset')
abline(0, 1, col = 'red', lty = 2)

dev.off()

```
Step 11: Calculate deviance-based _R<sup>2</sup>_
-------------------------------------
```
# Calculate the deviance (SS_res)
deviance <- sum((scale(y_test) - y_hat) ^ 2)
print(paste("deviance =", deviance))

# Calculate the mean of the scaled y_test
y_test_mean <- mean(scale(y_test))

# Calculate the null deviance (SS_tot)
deviance_null <- sum((scale(y_test) - y_test_mean)^ 2)
print(paste("deviance_null =", deviance_null))

# Calculate R2
R2 <- 1 - (deviance / deviance_null)
print(paste("R2 =", R2))
```

## Results

graph:

![](https://github.com/WCSCourses/PRS2024/blob/3f07962be8fc8915e8cde8666a4038f27c48fc83/images/image_afreur_merge.png)

pdf true against pred

![](https://github.com/WCSCourses/PRS2024/blob/3f07962be8fc8915e8cde8666a4038f27c48fc83/images/imagemerge_afrEur_pdf.png)











