# Introduction
In this exercise we will be covering some analyses of palaeogenomic data, how to simulate ancient DNA data, and evaluate the influences of potential aDNA biases in our results. Specifically we will investigate the impact of reference genome selection, the quality of the data, and the base call approach used in PCA, Neighbour joining tree, and D-statistic analyses.

## Getting started
- Data needed for this exercise can be downloaded from …
- Within this you will find:
  - “Results” - Example results from different tasks
  - “Reference_genomes” - reference genomes
  - “Spotted_map_bams” - Files mapped to the spotted hyena
  - “Striped_map_bams” - Files mapped to the striped hyena
- Sample information can be found in "BAM_information.txt" on this Github
- Spotted and Cave hyena data are from https://doi.org/10.1126/sciadv.aay0456
- Striped hyena is from https://doi.org/10.1093/molbev/msy037
- Aardwolf is from https://doi.org/10.1093/molbev/msab055

**Note**: Analyses can take some time to run so I recommend starting the subsequent analyses while plotting up results

# Task 1: Identification of aDNA damage using mapdamage
The most common approach to infer aDNA damage patterns is to use Mapdamage https://ginolhac.github.io/mapDamage/
- See what the possible parameters are

`mapDamage -h`
### Run mapdamage (optimised for speed)
 * Example individual - Ccsp015

`mapDamage -i Ccsp015.rmdup.sort_RG_Hi1.bam --merge-reference-sequences --no-stats -r ~/data/References/Crocuta/GWHAZPN00000000.genome_HiC.fasta -d Results/Ccsp015_mapdamage --downsample 1000000`
  - Example individual - 4035

`mapDamage -i 4035_map_merged_sort_RG_Hi1.bam --merge-reference-sequences --no-stats -r ~/data/References/Crocuta/GWHAZPN00000000.genome_HiC.fasta -d Results/4035_mapdamage --downsample 1000000`
- Look at the output plots of main interest - Fragmisincorporation_plot.pdf + Length_plot.pdf

**Question:** Which of these individuals is modern and which is ancient? How do you know?

# Task 2: Population genomic analyses
Here we will use some commonly implemented approaches in ancient population genomics that are suitable for low coverage data
**Note:** It can take awhile to run so while plotting the first outputs make sure the others are running in the background

## Run analyses to infer population structure
### PCA (GL and pseudo haploid base call)
* Make a text file with a list of the bam files you want to use (e.g. Bamlist.txt)
* Perform genotype likelihood (-GL + -Glf) and pseudohaploid (-doIBS) base calls - This example applies filters I commonly use, **if you want to know what all filters mean they are listed in the .arg file output after running the command**
  
`angsd -minmapQ 20 -minQ 20 -doCounts 1 -GL 2 -out Croc_0.1x_mInd13 -nThreads 10 -doGlf 2 -doMajorMinor 1 -rmtrans 1 -doMaf 2 -SNP_pval 1e-6 -b Bamlist.txt -rf ../../../Reference_genomes/Crocuta_scaffold1.txt -minmaf 0.05 -skiptriallelic 1 -uniqueonly 1 -minind 13 -dohaplocall 2 -doIBS 2 -minminor 2 -docov 1 -makematrix 1 -ref References/Crocuta/GWHAZPN00000000.genome_HiC.fasta`

**Note:** If you are changing between references/datasets, pay specific attention to the `-rf` `-ref` `-b` `-out` parameters

* Use PCANGSD to computed a covariance matrix from the GL

`pcangsd -b Spottedmap_minind11.beagle.gz -t 2 -o Spottedmap_minind11_pcangsd`

* Plot the covariance matrices using R

  Example

```R
# Import the covariance matrix (either -covMat for pseudohaploid or .cov for GL)
e=eigen(as.matrix(read.table("Spottedmap_minind11.covMat")))

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

#Print as PDF
dev.copy2pdf(file="Spottedmap_PCA_PH.pdf")
```


### Pairwise distances/phylogenetic trees (NJ)
Here we will build an unrooted neighbour joining phylogenetic tree from the distance matrix output with the above command. You can also rerun the command but with an outgroup in the bamfile if you want to construct a rooted tree
* Add names to the first column of the ibsMat (Distance matrix file) and add number of individuals to a row at the top e.g. `cut -f 2 -d "_" Dstats_names.txt |paste - Spottedmap_minind11.ibsMat | cat <(echo "17") - > Spottedmap_minind11.infile`

* Convert distance matrix into newick file using FASTME

`fastme -i Spottedmap_minind11.infile -o Spottedmap_minind11.tree`

* The output "Spottedmap_minind11.tree" can then be visualised with your favourite tree visualisation tool (e.g. figtree)

**Question:** Is there structure in this dataset? Are there differences between the base call methods?

## Run analysis to infer gene flow (D-statistics)
This requires a new bamlist with the outgroup at the bottom (either striped hyena or aardwolf)

* Compute Dstatistics in 1Mb blocks using a random base call approach in ANGSD

`angsd -minmapQ 20 -minQ 20 -doCounts 1 -out Spottedmap_minind11_stripedH4 -nThreads 5 -doabbababa 1 -rmtrans 1 -b Bamlist_Dstats_striped.txt -rf ../../../Reference_genomes/Crocuta_scaffold1.txt -uniqueonly 1 -minind 11 -uselast 1 -blocksize 1000000 -ref ../../../Reference_genomes/Crocuta_scaffold1.fasta -checkbamheaders 0`

* Perform block jacknifing with the R script as part of the ANGSD toolsuite
  
`Rscript ~/Scripts/ANGSD_jackknife.R file=Spottedmap_minind11_stripedH4.abbababa indNames=Dstats_names.txt outfile=Spottedmap_minind11_stripedH4.jack`

* You can filter the results using AWK
1. Only consider relevant topologies ((Cave,Cave),Spotted)
2. Flip H1 and H2 if negative for easier comparisons

`awk '$1~/Cave/&&$2~/Cave/&&$3~/Spot/ {print}' Spottedmap_minind11_aardwolfH4.jack.txt | awk '{if ($6>0) print $1,$2,$3,$6,$9; else print $2,$1,$3,$6*-1,$9*-1;}' | sort -r -k 5 | less`

**Question:** Which individuals have the most gene flow? Which have the least (remember ABBA-BABA / ABBA+BABA so positive = more ABBA) 

* If comparing between methods, you can also filter and plot using R
```R
# Load the data from each file
data1 <- read.table("Spottedmap_minind11_aardwolfH4.jack.txt", header = TRUE, sep = "\t")
data2 <- read.table("Spottedmap_minind11_stripedH4.jack.txt", header = TRUE, sep = "\t")

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

## Combine the data
combined_data <- cbind(filtered_data1, filtered_data2)

## Do a basic plot between the two datasets (Make it look nice if you like)
plot(combined_data$AardwolfH4_Z,combined_data$StripedH4_Z)

## Add a 1:1 line to represent unbiased results
abline(0,1,col=2)
```
**Question:** Do different outgroups give different D or Z values? What about different mapping references

# Task 3: Ancient DNA simulation
In this task we will simulate raw sequencing reads from a high quality modern genome with ancient damage using gargammel and map the reads to a reference genome
* Build fasta using consensus base call in ANGSD and unzip it

`angsd -minq 20 -docounts 1 -minmapq 20 -i NamCrocuta_map_merged_sort_RG_Hi1.bam -dofasta 2 -setmindepthind 10 -out NamCrocuta -r HiC_scaffold_1`

`gunzip NamCrocuta.fa.gz`
* Prepare directories for gargammel inlcuding three directories “bact” “cont” “endo” 

`mkdir Sequences`

`cd Sequences/`

`mkdir bact cont endo`

* Put the genome into the "endo" directory and index it

`cp NamCrocuta.fa endo`

`samtools faidx NamCrocuta.fa`

* Create a txt file with the proportion of each fragment lengths based on the mapdamage read lengths from your empirical data (file lgdistribution.txt in the mapdamage output directory):

`awk '/\+/{sum+=$3; count[$2]+=$3} END{for (i in count) print i"\t"count[i]/sum}' lgdistribution.txt > Fragment_lengths.txt`

* Run gargammel

`gargammel -h`

`gargammel -c 1 --comp 0,0,1 -f /home/zhc860/data/Cave_hyena/Workshop/Results/Ccsp015_mapdamage/Fragment_lengths.txt -mapdamage /home/zhc860/data/Cave_hyena/Workshop/Results/Ccsp015_mapdamage/misincorporation.txt single -rl 80 -o NamCroc.damaged Sequences`

The paired end output fastq of interest will end in _s1.fq.gz _s2.fq.gz

* Map reads using the script availabe in this github

`Ancient_mapping_PE.sh 3 . NamCroc Mapping ../Crocuta_scaffold1.fasta 30 0.01`	

Check for damage to see if it has worked (mapdamage output)

* Downsample and index final bam file

`samtools view -s 0.2 -o NamCroc.rmdup.sort_RG_0.2.bam NamCroc.rmdup.sort_RG.bam`

`samtools index NamCroc.rmdup.sort_RG.bam`


# Task 4: Investigating biases
Repeat the analyses from Task 2 but swap out a single (or multiple) modern individuals with their simulated damaged counterparts. 

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

