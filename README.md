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
- Run mapdamage (optimised for speed)
 * Example individual - Ccsp015

`mapDamage -i Ccsp015.rmdup.sort_RG_Hi1.bam --merge-reference-sequences --no-stats -r ~/data/References/Crocuta/GWHAZPN00000000.genome_HiC.fasta -d Results/Ccsp015_mapdamage --downsample 1000000`
  - Example individual - 4035

`mapDamage -i 4035_map_merged_sort_RG_Hi1.bam --merge-reference-sequences --no-stats -r ~/data/References/Crocuta/GWHAZPN00000000.genome_HiC.fasta -d Results/4035_mapdamage --downsample 1000000`
- Look at the output plots of main interest - Fragmisincorporation_plot.pdf + Length_plot.pdf

### Question: Which of these individuals is modern and which is ancient? How do you know?
