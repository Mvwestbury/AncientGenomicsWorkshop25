# Introduction
In this exercise we will be covering some analyses of palaeogenomic data, how to simulate ancient DNA data, and evaluate the influences of potential aDNA biases in our results. Specifically we will investigate the impact of reference genome selection, the quality of the data, and the base call approach used in PCA, Neighbour joining tree, and D-statistic analyses.

## Getting started
- If you would like to repeat this on your own server later, data needed for this exercise can be downloaded using the command `wget -O AncientGenomicsWorkshop.tar.gz https://sid.erda.dk/share_redirect/GU3WFFHp8z` and unzipped using `tar -zxvf AncientGenomicsWorkshop.tar.gz`
- Within this you will find:
  - “Results” - Example results from different tasks
  - “Reference_genomes” - reference genomes
  - “Spotted_map_bams” - Files mapped to the spotted hyena
  - “Striped_map_bams” - Files mapped to the striped hyena
- Sample information can be found in "BAM_information.txt" on this Github
- Spotted and Cave hyena data are from https://doi.org/10.1126/sciadv.aay0456
- Striped hyena is from https://doi.org/10.1093/molbev/msy037
- Aardwolf is from https://doi.org/10.1093/molbev/msab055

### **Note**: The paths to the input files for the analyses will need to be checked and adjusted as necessary  



# Task 1: Identification of aDNA damage using mapdamage
The most common approach to infer aDNA damage patterns is to use Mapdamage https://ginolhac.github.io/mapDamage/
- See what the possible parameters are

`mapDamage -h`
### Run mapdamage (optimised for speed)
 * Example individual - **Ccsp015**

`mapDamage -i Spotted_map_bams/Ccsp015.rmdup.sort_RG_Hi1.bam --merge-reference-sequences --no-stats -r Reference_genomes/GWHAZPN00000000.genome_HiC.fasta -d Ccsp015_mapdamage --downsample 1000000`
  - Example individual - **4035**

`mapDamage -i Spotted_map_bams/4035_map_merged_sort_RG_Hi1.bam --merge-reference-sequences --no-stats -r Reference_genomes/GWHAZPN00000000.genome_HiC.fasta -d 4035_mapdamage --downsample 1000000`
- Output directory is defined by -d
  - In this example the results are found in "Ccsp015_mapdamage"
- Look at the output plots of main interest
  - Fragmisincorporation_plot.pdf + Length_plot.pdf

**Task 1 Question:** Which of these individuals is modern and which is ancient? How do you know?



# Task 2: Population genomic analyses
Here we will use some commonly implemented approaches in ancient population genomics that are suitable for low coverage data

**Note:** As these take awhile to run, you can start running it but then cancel it with ctrl +c 
**All results needed can be found in the Results/Task2 directory. Results have been computed using two different reference genomes (Striped and Spotted hyena) and have been split into their own respective directories**

## Run analyses to infer population structure

### PCA (Genotype likelihoods and pseudo haploid base call) - RUN BUT CANCEL (TAKES TOO LONG)
* Make a text file with a list of the bam files you want to use (e.g. ls /home/wpsg/workshop_materials/31_ancient_genomics/Spotted_map_bams/*bam > Bamlist.txt **make sure to remove outgroup bams (Aardwolf_map_merged_sort_RG_Hi1.bam and Striped_hyena_map_merged_sort_RG_Hi1.bam) if not necessary**)
* You can also make one more manually using the BAM_information.txt from this Github
* Make all Bamlist files necessary (Mapped to Spotted hyena, Mapped to Striped hyena, with and without outgroups)
* Perform genotype likelihood (-GL + -Glf) and pseudohaploid (-doIBS) base calls in ANGSD - This example applies filters I commonly use, **if you want to know what all filters mean they are listed in the .arg file output after running the command or visit the website https://www.popgen.dk/angsd/index.php/ANGSD**
  
`angsd -minmapQ 20 -minQ 20 -doCounts 1 -GL 2 -out Spottedmap_minind11 -nThreads 10 -doGlf 2 -doMajorMinor 1 -rmtrans 1 -doMaf 2 -SNP_pval 1e-6 -b Bamlist.txt -rf Reference_genomes/Crocuta_scaffold1.txt -minmaf 0.05 -skiptriallelic 1 -uniqueonly 1 -minind 11 -dohaplocall 2 -doIBS 2 -minminor 2 -docov 1 -makematrix 1 -ref Reference_genomes/GWHAZPN00000000.genome_HiC.fasta -checkbamheaders 0`

You will get a few outputs of interest, they end in `.ibsMat` `.covMat` and `.beagle.gz`

**Note:** If you are changing between references/datasets, pay specific attention to the `-rf` (scaffold list) `-ref` (reference fasta) `-b` (bamlist) `-out` (output prefix) parameters

* Use PCANGSD to compute a covariance matrix from the GL beagle file

`pcangsd -h`

`pcangsd -b Spottedmap_minind11.beagle.gz -t 2 -o Spottedmap_minind11_pcangsd`

* Plot the covariance matrices using R (either ending in .covMat for pseudohaploid or .cov for GL)

  Example plotting code

```R
# Set the working directory
setwd("~/workshop_materials/31_ancient_genomics/Results/")

# Import the covariance matrix (either .covMat for pseudohaploid or .cov for GL)
e=eigen(as.matrix(read.table("Task2/Spotted_map/Spottedmap_minind11.covMat")))

# Extract eigenvalues
eigens=e$values/sum(e$values)*100

# Visualise the eigenvalues from the first ten Principal components
barplot(yaxt="n",names.arg=c(1:10),ylab="Percentage of variation",xlab="Principal component",eigens[1:10])
axis(2,las=2)

# Set the colours (the order of individuals is defined in your bamlist)
colours=c(2,2,2,2,2,2,2,1,1,1,1,1,1,1,2,2,2)

#Plot PCA
plot(e$vectors[,1], e$vectors[,2], lwd=2,
     ylab=sprintf("PC 2 (%.2f%%)", eigens[2]),
     xlab=sprintf("PC 1 (%.2f%%)", eigens[1]),
     col=colours, pch=16, cex=2, cex.lab=1.5, font=2)

#Add label
par(new=T)
labes=c("Cave", "Spotted")
legend("top",labes,cex=1,col=c(1,2),pch=16,bty='n')

```

**Rerun the PCa but using the results from the the different base call methods and reference genomes*

### Pairwise distances/phylogenetic trees (NJ)
Here we will build an unrooted neighbour joining phylogenetic tree from the distance matrix (.ibsMat) output with the above command. This runs quickly so you can run it yourself
You can also rerun the ANGSD command but with an outgroup in the bamfile if you want to construct a rooted tree (output not available in this tutorial)

* Add names (in the same order as in the Bamlist file) to the first column of the ibsMat (Distance matrix file) and add number of individuals to a row at the top

e.g. `cut -f 2 -d "_" Dstats_names.txt |paste - Spottedmap_minind11.ibsMat | cat <(echo "17") - > Spottedmap_minind11.infile` or make your own from the BAM_information.txt

* Convert distance matrix into newick file using FASTME

`fastme -i Spottedmap_minind11.infile -o Spottedmap_minind11.tree`

* The output "Spottedmap_minind11.tree" can then be visualised with your favourite tree visualisation tool (e.g. figtree or online at https://phylo.io/)

**Task 2 Questions:** Is there structure in this dataset? Are there differences between the PCA base call methods?

## Run analysis to infer gene flow (D-statistics)  - RUN BUT CANCEL (TAKES TOO LONG)
**This requires a new bamlist with the outgroup at the bottom (either striped hyena or aardwolf)**

* Compute Dstatistics in 1Mb blocks using a random base call approach in ANGSD

`angsd -minmapQ 20 -minQ 20 -doCounts 1 -out Spottedmap_minind11_stripedH4 -nThreads 5 -doabbababa 1 -rmtrans 1 -b Bamlist_Dstats_striped.txt -rf Reference_genomes/Crocuta_scaffold1.txt -uniqueonly 1 -minind 11 -uselast 1 -blocksize 1000000 -ref Reference_genomes/Crocuta_scaffold1.fasta -checkbamheaders 0`

* Create a text file containing a list of IDs for all individuals apart from the outgroup (Dstats_names.txt) -- For easier filtering I add a population identifier before the name e.g. Cave and Spot
* Perform block jacknifing with the R script as part of the ANGSD toolsuite (Script can be downloaded with `wget -O ANGSD_jackknife.R https://sid.erda.dk/share_redirect/fmXTAv65vF`)
  
`Rscript ANGSD_jackknife.R file=Spottedmap_minind11_stripedH4.abbababa indNames=Dstats_names.txt outfile=Spottedmap_minind11_stripedH4.jack`

* If you have added a population identifier to the indNames file above, you can filter the results using AWK
1. Only consider relevant topologies ((Cave,Cave),Spotted)
2. If the results are negative, flip the H1 and H2 individuals to make it positive for easier comparisons

`awk '$1~/Cave/&&$2~/Cave/&&$3~/Spot/ {print}' Spottedmap_minind11_stripedH4.jack.txt | awk '{if ($6>0) print $1,$2,$3,$6,$9; else print $2,$1,$3,$6*-1,$9*-1;}' | sort -r -k 4 | cat <(echo H1 H2 H3 Dscore Zscore) - | column -t | less`

**Question:** Which individuals have the most gene flow? Which have the least? (remember ABBA-BABA / ABBA+BABA so positive = more ABBA) 

* If comparing between methods, you can also filter and plot using R
```R
# Set the working directory
setwd("~/workshop_materials/31_ancient_genomics/Results/")

# Load the data from each Jackknifed Dstatistics output file (try this with multiple combinations of input files e.g. Different references)
data1 <- read.table("Task2/Spotted_map/Spottedmap_minind11_aardwolfH4.jack.txt", header = TRUE, sep = "\t")
data2 <- read.table("Task2/Spotted_map/Spottedmap_minind11_stripedH4.jack.txt", header = TRUE, sep = "\t")

# Prepend text to the headers of each data frame
names(data1) <- paste("AardwolfH4", names(data1), sep = "_")
names(data2) <- paste("StripedH4", names(data2), sep = "_")

### Further filter to only include the comparisons of interest
## i.e. ((Cave,Cave),Spotted)
filtered_data1 <- data1[grepl("Cave", data1[,1], ignore.case = TRUE) & 
                          grepl("Cave", data1[,2], ignore.case = TRUE) & 
                          grepl("Spot", data1[,3], ignore.case = TRUE), ]
filtered_data2 <- data2[grepl("Cave", data1[,1], ignore.case = TRUE) & 
                          grepl("Cave", data1[,2], ignore.case = TRUE) & 
                          grepl("Spot", data1[,3], ignore.case = TRUE), ]

## OR ((Spotted,Spotted),Cave)
filtered_data1 <- data1[grepl("Spot", data1[,1], ignore.case = TRUE) & 
                          grepl("Spot", data1[,2], ignore.case = TRUE) & 
                          grepl("Cave", data1[,3], ignore.case = TRUE), ]
filtered_data2 <- data2[grepl("Spot", data1[,1], ignore.case = TRUE) & 
                          grepl("Spot", data1[,2], ignore.case = TRUE) & 
                          grepl("Cave", data1[,3], ignore.case = TRUE), ]

## OR you can skip this filtering step and compare all results regardless of topology...

## Combine the data
combined_data <- cbind(filtered_data1, filtered_data2)

## Do a basic plot between the two datasets (Make it look nice if you like)
plot(combined_data$AardwolfH4_Dstat,combined_data$StripedH4_Dstat)

## Add a 1:1 line to represent unbiased results
abline(0,1,col=2)
```
**Question:** Do different outgroups give different D values? What about different mapping references? Any idea what could cause any differences?



# Task 3: Ancient DNA simulation
In this task we will simulate raw sequencing reads from a high quality modern genome with ancient damage using gargammel and map the reads to a reference genome

**Note** This also takes awhile so output can be found in Results/Task3
* Build fasta using consensus base call in ANGSD and unzip it (Takes ~1,000 seconds to run)

`angsd -minq 20 -docounts 1 -minmapq 20 -i Spotted_map_bams/NamCrocuta_map_merged_sort_RG_Hi1.bam -dofasta 2 -setmindepthind 10 -out NamCrocuta -rf Reference_genomes/Crocuta_scaffold1.txt`

`gunzip NamCrocuta.fa.gz`
* Prepare directory you want the output to go for gargammel and then within that directory make three directories “bact” “cont” “endo” 

`mkdir NamCroc`

`cd NamCroc/`

`mkdir bact cont endo`

* Put the fasta you created above into the "endo" directory and index it using SAMtools

`cp NamCrocuta.fa NamCroc/endo`

`samtools faidx NamCrocuta.fa`

* Create a txt file with the proportion of each fragment length based on the mapdamage read lengths from your empirical data (the file lgdistribution.txt in the mapdamage output directory):

`awk '/\+/{sum+=$3; count[$2]+=$3} END{for (i in count) print i"\t"count[i]/sum}' lgdistribution.txt > Fragment_lengths.txt`

* Run gargammel -- for a list of parameters type `gargammel -h` 

`gargammel.pl -c 1 --comp 0,0,1 -f Results/Task1/Ccsp015_mapdamage/Fragment_lengths.txt -mapdamage Results/Task1/Ccsp015_mapdamage/misincorporation.txt single -rl 80 -o NamCroc.damaged NamCroc`

The paired end output fastq of interest will end in _s1.fq.gz _s2.fq.gz

* Map reads using the "Ancient_mapping_PE.sh" script availabe in this github or downloaded with `wget -O Ancient_mapping_PE.sh https://sid.erda.dk/share_redirect/eIdHVdyz5I` **- RUN BUT CANCEL (TAKES TOO LONG)**
```
Ancient_mapping_PE.sh 3 . NamCroc Mapping Reference_genomes/Crocuta_scaffold1.fasta 30 0.01

Parameters are:
## $1 - Threads
## $2 - Raw reads folder
## $3 - Code name for sample
## $4 - Results folder
## $5 - Reference fasta
## $6 - Minimum read length
## $7 - Mismatch parameter
```

* Check the "records" file to see how many reads mapped, average coverage, and the number of bp that mapped

* Check for damage to see if it has worked (mapdamage output)

* **EXTRA If you like** you can downsample and index the final bam file to similar coverage to the lowest coverage ancient sample

`samtools view -s 0.2 -o NamCroc.rmdup.sort_RG_0.2.bam NamCroc.rmdup.sort_RG.bam`

`samtools index NamCroc.rmdup.sort_RG_0.2.bam`


# Task 4: Investigating biases
Repeat the analyses from Task 2 but swap out a single (or multiple) modern individuals with their simulated damaged counterparts in the bamlist. 

* I created three simulated individuals using the mapdamage patterns from Ccsp015 - I used the high coverage ones to ensure the fasta file used in gargammel was highest quality possible.

Again, as the ANGSD commands take awhile to run, all results are found within *Results/Task4* and are split between reference genomes - the file name denotes which individual has been damaged and can open the Bamlist*txt files to see their position in the bamlist for later plotting

* You can simply plot the PCA and NJ tree as before and visually compare them - change the shape or colour of the damaged individuals to make easier visualisation

For the Dstatistics, you can compare results in a similar manner to above but only look at the comparisons comparing the simulated damaged individual to its high quality equivalent

As an example...

```R
# Load the data from each file (Change the data2 for some of the other jack.txt files)
data1 <- read.table("Spottedmap_minind11_stripedH4.jack.txt", header = TRUE, sep = "\t")
data2 <- read.table("Spottedmap_minind11_Namsim_stripedH4.jack.txt", header = TRUE, sep = "\t")

# Prepend text to the headers of each data frame
names(data1) <- paste("Original", names(data1), sep = "_")
names(data2) <- paste("Simulated", names(data2), sep = "_")


## Combine the data
combined_data <- cbind(data1, data2)

## Filter only for rows that contain SIM - I added SIM to the name of the simulated individuals in the Dstats_names.txt file
filtered_data <- combined_data[grepl("SIM", combined_data$Simulated_H1) | 
                          grepl("SIM", combined_data$Simulated_H2) | 
                          grepl("SIM", combined_data$Simulated_H3), ]

## Do a basic plot between the two datasets (Make it look nice if you like)
plot(filtered_data$Original_D,filtered_data$Simulated_D)

## Add a 1:1 line to represent unbiased results
abline(0,1,col=2)
```

**Question:** Do the results change if we replace the high quality modern individual with a low quality ancient equivalent?




# Software
* Mapdamage https://ginolhac.github.io/mapDamage/
* R
* ANGSD https://github.com/ANGSD/angsd and https://www.popgen.dk/angsd/index.php/ANGSD
* PCAngsd https://github.com/Rosemeis/pcangsd
* Gargammel https://github.com/grenaud/gargammel
* SAMtools https://github.com/samtools/samtools
* Fastme http://www.atgc-montpellier.fr/fastme/binaries.php
* Fastp - https://github.com/OpenGene/fastp
* BWA - https://github.com/lh3/bwa

