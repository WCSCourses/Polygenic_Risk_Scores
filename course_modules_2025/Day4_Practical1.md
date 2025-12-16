 
## Day 4 - Practical 1: Introduction to Admixture analysis

### Module Goals
 The goal of this practical is to provide you with basic understanding of Admixture and the basic elements behind Admixture PRS scores. Upon completion of this practical, you should:
* Gain familiarity with a variety of tools used in the context of admixed population research
* Go through the basic steps of formulating admixture-informed polygenic risk scores


### Part1: Plot Decay of Ancestry LD over time
Here we explore how Admixture LD varies over time and as a function of the genetic distance between loci. This script is entirely generating (simulating) the decay curves to be plotted.In R:
```
library(ggplot2)
library(reshape2)
library(viridis)

# Define a sequence of recombination fractions (r) from 0 to 0.5 in steps of 0.1
r = seq(0, 0.5, 0.1)

# Initialize an empty matrix to store admixture LD values (Dt)
# Rows correspond to different r values (length 6), columns to generations (1–10)
dtmat = matrix(NA, nrow = 6, ncol = 10) # matrix to store values of dt

# Loop over each recombination fraction
for(i in 1:6){
  dt = 0.25
  for(j in 1:10){
	dtmat[i,j] = dt
	dt = dt*(1 -r[i])^j
  }
}

# Melt the matrix into a long data frame with columns: row (r index), column (generation), and Dt
dtmat = reshape2::melt(dtmat)

# Rename the automatically generated columns for clarity
colnames(dtmat) = c("r","g","Dt")
dtmat$r = (dtmat$r - 1)/10

# Build plot
p <- ggplot(dtmat,
            aes(x = g,
                y = Dt,
                group = r,
                colour = as.factor(r))) +
  geom_line(linewidth = 1.2) +
  theme_minimal(base_size = 14) +
  theme(axis.text        = element_text(size = 14),
        axis.title       = element_text(size = 16),
        legend.title     = element_text(size = 14),
        legend.text      = element_text(size = 12),
        plot.title       = element_text(size = 18, face = "bold"),
        plot.subtitle    = element_text(size = 14),
        panel.grid.major = element_line(colour = "gray80")) +
  scale_color_viridis_d() +
  scale_x_continuous(breaks = 1:10) +
  labs(title    = "Decay of Linkage Disequilibrium Over Generations",
       subtitle = "Effect of Recombination Fraction (r) on Admixture LD",
       x        = "Generations since admixture (t)",
       y        = "Admixture LD",
       colour   = "Recombination Fraction (r)")

# Display & save
print(p)

ggsave(filename = "./LD_decay.png",
       plot     = p,
       width    = 8, 
       height   = 6,
       dpi      = 300) 
```
#### **Questions**
##### (i) Describe what happens to admixture LD over time?
##### (ii) How does recombination affect this relationship?
 
### Part 2: Global Ancestry Inference
We will now run an analysis using the software ADMIXTURE to calculate global ancestry proportions across a sample of 28 individuals. Here we perform a supervised analysis. 
 ```
 ./software/admixture ./data/plink/samples_n28_qc_thin.bed 2 --supervised -j4
```
#### **Questions**
##### (i) What do you think the number specified after the input file represents?

The next step is to create a plot the results
```
# Plot Global Ancestry Results
 R
 library(ggplot2)
 library(reshape2)
 
 # Read the data from files
 fam <- read.table("./data/plink/samples_n28_qc_thin.fam", header = FALSE)
 pop <- read.table("./data/plink/samples_n28_qc_thin.pop", header = FALSE)
 Q <- read.table("./samples_n28_qc_thin.2.Q", header = FALSE)
 
 # Merge the data into one data frame
 merged_data <- cbind(fam, pop, Q)
 
 # Extract necessary columns (assuming the structure of the files as in the image)
 colnames(merged_data) <- c("FID", "IID", "PaternalID", "MaternalID", "Sex", "Phenotype", "Population", "AFR", "EUR")
 
 # Order data by decreasing AFR
 merged_data <- merged_data[order(-merged_data$AFR), ]
 
 # Prepare data for plotting
 plot_data <- merged_data[, c("IID", "AFR", "EUR")]
 plot_data$IID <- factor(plot_data$IID, levels = plot_data$IID)  # Ensure order is maintained in plot
 plot_data_melted <- melt(plot_data, id.vars = "IID")
 
 # Create the plot
 p <- ggplot(plot_data_melted, aes(x = IID, y = value, fill = variable)) +
   geom_bar(stat = "identity", position = "stack") +
   scale_fill_manual(values = c("AFR" = "red", "EUR" = "darkgreen")) +
   labs(x = "Individual ID", y = "Ancestry Proportion", fill = "Ancestry", title = "Global Ancestry Proportions") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
         axis.text.y = element_text(face = "bold"),  # Bold y-axis numbers
         axis.title.x = element_text(face = "bold"),
         axis.title.y = element_text(face = "bold"),
         legend.title = element_text(face = "bold"),
         plot.title = element_text(face = "bold", hjust = 0.5))
 
 # Save the plot with specified dimensions
 ggsave("ancestry_plot.png", plot = p, width = 10, height = 6, units = "in", dpi = 300)

# Close R using the command q() and when prompted whether to save workspace image enter "n" for no

# In Bash enter the following command to view the plotted admixture report:
xdg-open ancestry_plot.png
```
#### Questions
##### (i) How do ancestry assignments vary across the 28 individuals shown in the plot? (Use the scale to interpret ancestry proportions for specific individuals)


### Part 3: Local Ancestry Inference 

Next, we would like to annotate each genome with information about the ancestry of an individual chromosomal segment. This procedure is called local ancestry inference and it can be used to incorporate information about the ancestral origin of individual SNPs for more personalised PRS calculation within admixed individuals.
We will use the RFMix software to calculate local ancestry on chromosome 22 for the same set of 28 individuals. The RFMix algorithm uses an unsupervised learning algorithm. In this instance, “unsupervised learning” means that RFMix detects patterns in haplotype structure and groups similar segments together without needing pre-assigned ancestry labels, allowing it to infer each segment’s ancestral origin directly from the genetic data. Here, a haplotype refers to the specific combination of alleles at neighbouring loci along a chromosome that tend to be inherited together as a block.
```
# Run the following UNIX command from the home directory

for i in {22..22}; do
    ./software/rfmix \
        -f ./data/rfmix/chr1-22_phased.bcf.gz \
        -r ./reference/chr22_reference.bcf \
        --analyze-range=26.86-31.80 \
        -m ./data/rfmix/1KG_superpop_vs_ID.txt \
        --chromosome=${i} \
        -g ./reference/1kg_chr1-22.gmap \
        --n-threads=4 \
        -o ./out/rfmix/chr${i}.local_ancestry
done
```
Here’s what each part of that loop and RFMix flag does, and why it matters:

for i in {22..22}; do … done <br/>
What it is: A bash loop that sets i to each integer from 22 to 22 (so just chromosome 22 here).<br/>
Why it’s important: Lets you easily run the exact same command over a range of chromosomes by changing the bounds (e.g. {1..22}), without copy-pasting.<br/>

-f ./data/rfmix/chr1-22_phased.bcf.gz <br/>
What it is: The phased haplotypes file for all samples, in BCF/VCF format.<br/>
Why it’s important: RFMix needs phased data (i.e. which alleles are on each parental chromosome) to model the sequential inheritance patterns that distinguish ancestry tracts.

-r ./reference/chr22_reference.bcf <br/>
What it is: The reference panel BCF/VCF for chromosome 22, containing known‐ancestry haplotypes.<br/>
Why it’s important: These labeled haplotypes train the random‐forest classifier to recognize ancestry‐specific patterns in your target data.<br/>

--analyze-range=26.86-31.80<br/>
What it is: Restricts inference to the genetic map window from 26.86 cM to 31.80 cM on the chromosome.<br/>
Why it’s important: Focusing on a subregion speeds up computation and can improve accuracy if you only care about a particular locus or block.<br/>

-m ./data/rfmix/1KG_superpop_vs_ID.txt<br/>
What it is: The sample‐map file: each line maps a sample (or haplotype) ID to its population label.<br/>
Why it’s important: RFMix uses these labels to teach the Randm Forest classifier which haplotype patterns correspond to which ancestral group.<br/>

--chromosome=${i}<br/>
What it is: Tells RFMix which chromosome number to process—in this case, the value of i (22).<br/>
Why it’s important: Ensures the tool applies the correct genetic map and reference panel for that chromosome.<br/>

-g ./reference/1kg_chr1-22.gmap<br/>
What it is: The genetic map file mapping physical positions (bp) to genetic distances (cM) for chromosomes 1–22.<br/>
Why it’s important: Accurate local ancestry inference depends on knowing recombination rates; the genetic map lets RFMix translate base‐pair positions into recombination distances.<br/>

--n-threads=4<br/>
What it is: Number of CPU threads to use in parallel.<br/>
Why it’s important: Speeds up the random‐forest training and ancestry assignment steps by leveraging multiple cores.<br/>

-o ./out/rfmix/chr${i}.local_ancestry<br/>
What it is: The output prefix for all result files (e.g. .msp.tsv, .fb.tsv).<br/>
Why it’s important: Organizes outputs per chromosome and ensures you can trace results back to the exact input settings.<br/>

#### Questions
##### (i) ...


### Part 4: Plot local admixture on chromosome 22

In this step we will plot the output from the previous RFMix run. Execute the following code from the home directory

```
R
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to read the .msp.tsv file
read_msp_file <- function(msp_file) {
  txt  <- readLines(msp_file)
  hdr  <- txt[grepl("^#", txt)]
  skip <- length(hdr)
  coln <- strsplit(hdr[skip], "\t")[[1]]

  msp  <- read.table(msp_file,
                     header = FALSE,
                     skip   = skip,
                     sep    = "\t",
                     col.names = coln,
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
  return(msp)
}

# Function to read the .Q file
read_Q_file <- function(Q_file) {
  read.table(Q_file,
             header       = TRUE,
             sep          = "\t",
             comment.char = "#",
             col.names    = c("sample", "AFR", "AMR", "EAS", "EUR", "SAS"),
             stringsAsFactors = FALSE)
}

# Specify file paths
msp_file <- './out/rfmix/chr22.local_ancestry.msp.tsv'
Q_file <- './out/rfmix/chr22.local_ancestry.rfmix.Q'

# Read in files
msp_df <- read_msp_file(msp_file)
Q_df <- read_Q_file(Q_file)

# Rename the haplotype columns so they are unique & easy to pivot
haplo_cols_raw <- colnames(msp_df)[7:62]
colnames(msp_df)[7:62] <- paste0("haplo_", haplo_cols_raw)

# 4b. helper to map ancestry codes to labels
determine_ancestry <- function(x) {
  ancestry_map <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  ancestry_map[x + 1]
}

# Long format and metadata columns
plot_data <- msp_df %>%
  pivot_longer(
    cols         = starts_with("haplo_"),
    names_to     = "haplo_col",
    values_to    = "ancestry_code",
    names_repair = "unique"
  ) %>%
  mutate(
    individual = sub("\\..*", "", haplo_col),              # ID prefix
    strand     = ifelse(grepl("\\.0$", haplo_col),
                        "Strand 1", "Strand 2"),
    ancestry   = determine_ancestry(ancestry_code)
  ) %>%
  arrange(individual, strand, sgpos)                       # reproducible order


# Create plot of chromosome 22 across samples
# Spacing logic: give each person a 3-unit band (0 = lower margin, 3 = start of next individual)

plot_data <- plot_data %>%
  mutate(
    ind_index = as.integer(factor(individual,
                                  levels = unique(individual))),
    y_pos = ind_index * 3 +
            ifelse(strand == "Strand 1", 0.4, 1.4)         # within-person gap
  )

# y-axis labels centred between the two haplotypes
label_df <- plot_data %>%
  distinct(individual, ind_index) %>%
  mutate(label_y = ind_index * 3 + 0.9)

p <- ggplot(plot_data,
            aes(x    = sgpos,
                xend = egpos,
                y    = y_pos,
                yend = y_pos,
                colour = ancestry)) +
  geom_segment(linewidth = 4, lineend = "butt") +
  scale_color_manual(values = c(
    "AFR" = "blue",
    "AMR" = "orange",
    "EAS" = "green",
    "EUR" = "red",
    "SAS" = "purple"
  )) +
  scale_y_continuous(
    breaks = label_df$label_y,
    labels = label_df$individual,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x      = "Genetic position (cM)",
    y      = "Individual ID and Haplotype",
    title  = "Local Ancestry Across Chromosome 22",
    colour = "Ancestry"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y  = element_text(size = 7),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9)
  )

# Save plot 
out_dir <- "./out/rfmix"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(
  filename = file.path(out_dir,
                       "local_ancestry_chromosome22_5Mb_subregion.png"),
  plot     = p,
  width    = 10,
  height   = 8,
  dpi      = 300
)

# Navigate to the output folder ./out/rfmix and open the output file using the below command:
xdg-open local_ancestry_chromosome22_5Mb_subregion.png
```
#### Questions
##### (i) How many different continental ancestries do you see represented in this 5 megabase block, across the 58 strands (corresponding to 29 individuals)?


### Part 5: Formatting of admixture files for analysis using PRSice
We now want to integrate the genotype data and local ancestry data for the 29 individuals into the same file. We will do this in two steps. From there we will create a custom set of Plink files and use them to perform an ancestry-informed version of the normal Prsice2 analysis, to evaluate efficacy in our small admixed cohort.

#### Step 1 - Convert phased genotypes to GenomicRange format
```
# Run from the home directory

library(vcfR)
library(memuse)
library(panelr)
library(GenomicRanges)
library(data.table)

clean_memory <- function(vars_to_remove) {
  rm(list = vars_to_remove)
  gc()
}

# Read in phased VCF file
vcfr_inputfile_chr22 <- read.vcfR("./data/rfmix/chr22_phased.vcf.gz", verbose = FALSE)
extracted_haps22 <- extract.haps(vcfr_inputfile_chr22, mask = FALSE, unphased_as_NA = TRUE, verbose = TRUE)
extracted_snp_info22 <- getFIX(vcfr_inputfile_chr22)

# Convert haplotypes to data frame
haps_df <- data.frame(extracted_haps22, check.names = FALSE)
haps_dt <- setDT(haps_df, keep.rownames = "snps")

# Convert SNP info to data frame
snp_info_df <- data.frame(extracted_snp_info22)

# Check for numerical and ordering consistency between haplotypes and SNP info
if (!identical(haps_dt[['snps']], snp_info_df[['ID']])) {
  stop("Numerical and ordering inconsistency between haplotypes and SNP info")
}

# Merge haplotypes and SNP info
merge_chr22 <- cbind(snp_info_df, haps_dt)

# Convert to long format using panelr
chr22_haplo_long <- long_panel(merge_chr22, prefix = "_", begin = 0, end = 1, label_location = "end", as_panel_data = FALSE)

# Insert 'end' column after 'POS' and create 'wave' column
chr22_haplo_long$end <- chr22_haplo_long$POS
chr22_haplo_long$wave <- ifelse(chr22_haplo_long$wave == "0", "+", ifelse(chr22_haplo_long$wave == "1", "-", "Z"))

# Remove row names and drop redundant columns using base R
rownames(chr22_haplo_long) <- NULL
drops <- c("id", "snps", "QUAL", "FILTER")
chr22_haplo_long <- chr22_haplo_long[, !names(chr22_haplo_long) %in% drops]

# Rename columns
names(chr22_haplo_long)[3] <- "start"
names(chr22_haplo_long)[1] <- "strand"
header <- gsub(".*_", "", colnames(chr22_haplo_long)[7:ncol(chr22_haplo_long)])
names(chr22_haplo_long)[7:ncol(chr22_haplo_long)] <- header

# Create GRanges object
gr_obj_chr22 <- makeGRangesFromDataFrame(chr22_haplo_long, ignore.strand = FALSE, keep.extra.columns = TRUE)

# Save the GRanges object
saveRDS(gr_obj_chr22, file = "./out/rfmix/chr22_phased_gr.rds")
```

### Part 5: Step 2 - Merge genotypes from Step 1 with local ancestry calls by RFMix
```
# Read in MSP file
msp <- fread("./out/rfmix/chr22.local_ancestry.msp.tsv")

# Reformat genome-wide local ancestry output as GRanges object
colnames(msp)[1] <- "chm"
msp_gr <- makeGRangesFromDataFrame(msp,
								   seqnames.field = "chm",
								   start.field = "spos",
								   end.field = "epos",
								   keep.extra.columns = TRUE)

# Convert the elements of the GRange object into a dataframe
chr22_rf <- data.frame(seqnames = seqnames(msp_gr),
					   ranges = ranges(msp_gr),
					   strand = strand(msp_gr),
					   mcols(msp_gr), check.names = FALSE)
names(chr22_rf)[1:5] <- c("chr", "start", "end", "width", "strand")
rownames(chr22_rf) <- NULL

# Clean up column headers
header <- gsub(".*_", "", colnames(chr22_rf)[9:ncol(chr22_rf)])
names(chr22_rf)[9:ncol(chr22_rf)] <- header

# Convert local ancestry calls from wide to long format
chr22_rf_long <- long_panel(chr22_rf, prefix = ".", begin = 0, end = 1, label_location = "end", as_panel_data = FALSE)
chr22_rf_long <- as.data.frame(chr22_rf_long, check.names = FALSE)
chr22_rf_long$strand <- ifelse(chr22_rf_long$wave == "0", "+", ifelse(chr22_rf_long$wave == "1", "-", "Z"))

# Drop redundant columns
drops <- c("wave", "id")
chr22_rf_long <- chr22_rf_long[, !(names(chr22_rf_long) %in% drops)]
rownames(chr22_rf_long) <- NULL

# Convert the reconfigured dataframe file into a GRanges object
chr22_msp_gr <- makeGRangesFromDataFrame(chr22_rf_long, ignore.strand = FALSE, keep.extra.columns = TRUE)

# Ensure chr22_haplo_long is a GRanges object before finding overlaps
chr22_haplo_gr <- makeGRangesFromDataFrame(chr22_haplo_long, ignore.strand = FALSE, keep.extra.columns = TRUE)

# Find overlaps and store matching features
# Matching is coordinates based. Base position ("start"/"end")is used to match the two files
matched_regions <- findOverlaps(chr22_haplo_gr, chr22_msp_gr)
chr22_haps_lanc_gr <- chr22_haplo_gr[queryHits(matched_regions)]  #Store matching features in a new dataframe, add metadata from RFmix output.
mcols(chr22_haps_lanc_gr) <- cbind.data.frame(mcols(chr22_haps_lanc_gr), mcols(chr22_msp_gr[subjectHits(matched_regions)]))

# Local ancestry calls are now aligned with genotypic data and positional coordinates
# Convert to dataframe without adding 'X' to numeric column names
genes_df <- as.data.frame(chr22_haps_lanc_gr)
names(genes_df) <- sub('^X', '', names(genes_df))

# Output the dataframe
write.table(genes_df, file = "./out/rfmix/chr22_phased_geno_lanc.txt", quote = FALSE, row.names=F)

# Clean up memory
rm(list = c("chr22_haplo_gr", "chr22_haplo_long", "msp", "msp_gr", "chr22_rf", "header",
			"chr22_rf_long", "chr22_msp_gr", "matched_regions",
			"chr22_haps_lanc_gr"))
gc()
```

#### Part 5: Step 3 - Create separate Plink files
```
# Partition genotype and local ancestry data
chr22_genos <- genes_df[, c(1, 2, 6, 9:36)]
chr22_LA <- genes_df[, c(1, 2, 6, 40:ncol(genes_df))]

# Reintegrate in interleaved format
d <- chr22_genos[, -c(1:3)]  # Retain ID columns only
d2 <- chr22_LA[, -c(1:3)]

# Prepare headers
indx <- rbind(names(d), names(d2))
dmerge <- cbind(d, d2)
dfinal <- dmerge[, indx]

# Convert LAnc data from long to wide format
LA_wide <- lapply(1:ncol(chr22_LA), function(i) as.data.frame(matrix(chr22_LA[, i], ncol = 2, byrow = TRUE)))

# Check length of LA_wide
length(LA_wide)  # Expected number. Each item contains a set of matched strand pairs per individual

LA_final <- as.data.frame(do.call(cbind, LA_wide))
LA_final <- LA_final[, -c(2, 4, 6)]  # Remove duplicates of first 3 cols
colnames(LA_final)[4:ncol(LA_final)] <- colnames(dfinal)  # Apply headers

# Convert genotype calls to horizontal orientation
geno_wide <- lapply(1:ncol(chr22_genos), function(i) as.data.frame(matrix(chr22_genos[, i], ncol = 2, byrow = TRUE)))  # Genotypic part: Get horizontal genotypes
geno_final <- as.data.frame(do.call(cbind, geno_wide))  # Merge horizontal genotypes

# Clean up and finalize geno_final
geno_final <- geno_final[, -c(2, 4, 6)]  # Redundant columns
colnames(geno_final)[4:ncol(geno_final)] <- colnames(dfinal)

# Name columns
colnames(geno_final)[1:3] <- c("CHROM", "BP", "ID")
colnames(LA_final)[1:3] <- c("CHROM", "BP", "ID")

# Write final tables
write.table(geno_final, "./out/rfmix/chr22_geno.txt", row.names = FALSE, quote = FALSE)
write.table(LA_final, "./out/rfmix/chr22_LA.txt", row.names = FALSE, quote = FALSE)

# Close R using the command q() and when prompted whether to save workspace image enter "n" for no
```

#### Part 5: Step 4 - Use custom software (RFTransform) to create Plink files for input into PRSice
```
# This section of the code will generate the custom Plink files. It uses a custom-build software not available anywhere else.
./software/RFTransform/build/RFTransformer ./out/rfmix/chr22_geno.txt ./out/rfmix/chr22_LA.txt ./out/plink/chr22

# Perform general QC ahead of running PRSice on ancestry-deconvolved individuals
# (i) Remove SNPs with low minor allele count (MAC)

# Plink - Remove monomorphic SNPs (minor allele count 0-4).
#AFR
for i in {22..22}; do
		./software/plink2 \
		--bfile ./out/plink/chr${i}-AFR \
		--mac 5 \
		--make-bed \
		--out ./out/plink/chr${i}-AFR.mac
done

#EUR
for i in {22..22}; do
		./software/plink2 \
		--bfile ./out/plink/chr${i}-EUR \
		--mac 5 \
		--make-bed \
		--out ./out/plink//chr${i}-EUR.mac
done

# The chromosome 22 genetic data contained in the above output files has been divided into separate European and African components for each individual in the analysis.
# The next section of code will perform separate plink analyses for each ancestry component using weights from ancestry-matched summary statistics.

# PRSice - Generate ancestry-specific weights
# AFR
Rscript ./software/PRSice.R \
--prsice ./software/PRSice_linux \
--base ./data/plink/AFR-BMI.Phenotype.glm.linear \
--extract ./data/plink/BMI_AFR-base.valid \
--A1 A1 \
--pvalue P \
--stat BETA \
--pheno ./data/plink/pheno.plink \
--beta \
--snp ID \
--score sum \
--target ./out/plink/chr22-AFR.mac \
--binary-target F \
--out ./out/prsice/BMI_AFR-base

# EUR
Rscript ./software/PRSice.R \
--prsice ./software/PRSice_linux \
--base ./data/plink/EUR-BMI.Phenotype.glm.linear \
--extract ./data/plink/BMI_EUR-base.valid \
--A1 A1 \
--pvalue P \
--stat BETA \
--pheno ./data/plink/pheno.plink \
--beta \
--snp ID \
--score sum \
--target ./out/plink/chr22-EUR.mac \
--binary-target F \
--out ./out/prsice/BMI_EUR-base
```

#### Part 5: Step 5 - Evaluate the Admixture-informed PRS
We are now ready to combine the two sets of polygenic risk scores per individual, in order to ascertain and evaluate the full admixture-aware PRS performance for each individual.
```
library(dplyr)

# Read the PRS files into dataframes
file1 <- read.table("./out/prsice/BMI_EUR-base.best", header = TRUE, check.names=F)
file2 <- read.table("./out/prsice/BMI_AFR-base.best", header = TRUE, check.names=F)

# Sum the fourth column ($PRS) of both files
file1$PRS_SUM <- file1$PRS + file2$PRS

# Load the phenotype data
pheno <- read.table("./data/plink/pheno.plink", header = F)
colnames(pheno) <- c("FID", "IID", "phenotype")

# Convert phenotype column to numeric
pheno$phenotype <- as.numeric(pheno$phenotype)

# Merge PRS data with phenotype data
merged_data <- merge(file1[, c("FID", "IID", "PRS", "PRS_SUM")], pheno, by = c("FID", "IID"))
merged_data <- merge(merged_data, file2[, c("FID", "IID", "PRS")], by = c("FID", "IID"))

# Rename columns
names(merged_data)[names(merged_data) == "PRS.x"] <- "PRS1"
names(merged_data)[names(merged_data) == "PRS.y"] <- "PRS2"

# Perform linear regression for each PRS and the combined PRS
model1 <- lm(phenotype ~ PRS1, data = merged_data)
model2 <- lm(phenotype ~ PRS2, data = merged_data)
model_sum <- lm(phenotype ~ PRS_SUM, data = merged_data)

# Extract R-squared values
r_squared1 <- summary(model1)$r.squared
r_squared2 <- summary(model2)$r.squared
r_squared_sum <- summary(model_sum)$r.squared

# Print the R-squared values
cat("R-squared for PRS1: ", r_squared1, "\n")
cat("R-squared for PRS2: ", r_squared2, "\n")
cat("R-squared for PRS_SUM: ", r_squared_sum, "\n")
```

#### **Questions**
##### (i) How does the R-squared of the combined-ancestry PRS perform relative to the 2 partial-genome PRSs?
##### (ii) How do you rationalise the result obtained?
