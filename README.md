# Introduction
In this exercise we will be covering some analyses of palaeogenomic data, how to simulate ancient DNA data, and evaluate the influences of potential aDNA biases in our results. Specifically we will investigate the impact of reference genome selection, the quality of the data, and the base call approach used in PCA, Neighbour joining tree, and D-statistic analyses.

## Getting started
- Data needed for this exercise can be downloaded from …
- Within this you will find:
  - “Results” - Example results from different tasks
  - “Reference_genomes” - reference genomes
  - “Spotted_map_bams” - Files mapped to the spotted hyena
  - “Striped_map_bams” - Files mapped to the striped hyena
- Spotted and Cave hyena data are from https://doi.org/10.1126/sciadv.aay0456
- Striped hyena is from https://doi.org/10.1093/molbev/msy037
- Aardwolf is from https://doi.org/10.1093/molbev/msab055


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

### Question: Which of these individuals is modern and which is ancient? How do you know?

# Task 2: Population genomic analyses
Note: It can take awhile to run so while plotting the first outputs make sure the others are running in the background

## Run analyses to infer population structure
### PCA (GL and pseudo haploid base call)
  
* Performed genotype likelihood (-GL + -Glf) and pseudohaploid (-doIBS) base calls `angsd -minmapQ 20 -minQ 20 -doCounts 1 -GL 2 -out Croc_0.1x_mInd13 -nThreads 10 -doGlf 2 -doMajorMinor 1 -rmtrans 1 -doMaf 2 -SNP_pval 1e-6 -b Bamlist.txt -r HiC_scaffold_1 -minmaf 0.05 -skiptriallelic 1 -uniqueonly 1 -minind 13 -dohaplocall 2 -doIBS 2 -minminor 2 -docov 1 -makematrix 1 -ref ~/data/References/Crocuta/GWHAZPN00000000.genome_HiC.fasta`

* Use PCANGSD to computed a covariance matrix from the GL `pcangsd -b Spottedmap_minind11.beagle.gz -t 2 -o Spottedmap_minind11_pcangsd`

* Plot the covariance matrices using R

  Example

```R
e=eigen(as.matrix(read.table("Spottedmap_minind11.covMat")))
eigens=e$values/sum(e$values)*100

barplot(yaxt="n",names.arg=c(1:10),ylab="Percentage of variation",xlab="Principal component",eigens[1:10])
axis(2,las=2)
colours=c(2,2,2,2,2,2,2,1,1,1,1,1,1,1,2,2,2)
plot(e$vectors[,1], e$vectors[,2], lwd=2,
     ylab=sprintf("PC 2 (%.2f%%)", eigens[2]),
     xlab=sprintf("PC 1 (%.2f%%)", eigens[1]),
     col=colours, pch=16, cex=2, cex.lab=1.5, font=2)
par(new=T)
labes=c("Cave", "Spotted")
legend("top",labes,cex=1,col=c(1,2),pch=16,bty='n')

dev.copy2pdf(file="Spottedmap_PCA_PH.pdf")
```


### Pairwise distances/phylogenetic trees (NJ)
* Add names to the first column of the ibsMat (Distance matrix file) and add number of individuals to a row at the top e.g. `cut -f 2 -d "_" Dstats_names.txt |paste - Spottedmap_minind11.ibsMat | cat <(echo "17") - > Spottedmap_minind11.infile`

* Convert distance matrix into newick file using FASTME `fastme -i Spottedmap_minind11.infile -o Spottedmap_minind11.tree`

* The output can then be visualised with your favourite tree visualisation tool (e.g. figtree)

### Question: Is there structure in this dataset? Are there differences between the base call methods?

## Run analysis to infer gene flow (D-statistics)
angsd -minmapQ 20 -minQ 20 -doCounts 1 -out Spottedmap_minind11_stripedH4 -nThreads 5 -doabbababa 1 -rmtrans 1 -b Bamlist_Dstats_striped.txt -rf ../../../Reference_genomes/Crocuta_scaffold1.txt -uniqueonly 1 -minind 11 -uselast 1 -blocksize 1000000 -ref ../../../Reference_genomes/Crocuta_scaffold1.fasta -checkbamheaders 0
Rscript ~/Scripts/ANGSD_jackknife.R file=Spottedmap_minind11_stripedH4.abbababa indNames=Dstats_names.txt outfile=Spottedmap_minind11_stripedH4.jack

