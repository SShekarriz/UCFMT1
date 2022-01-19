# UCFMT1
Codes and description of methods for UCFMT1 project

## DNA extraction and 16S rRNA gene sequencing
Genomic DNA extraction and PCR amplification of the V3 region of 16S rRNA gene, was conducted using
previously described protocols [9, 51, 52]. Briefly, 0.2 g of fecal matter was mechanically homogenized using
ceramic beads in 800 μL of 200 mM NaPO 4 (pH 8) and 100 μL of guanidine thiocyanate-EDTA-N-lauroyl
sacosine. This was followed by enzymatic lysis of the supernatant using 50 μL of 100 mg/mL lysozyme,
50 μL of 10 U/μL mutanolysin, and 10 μL of 10 mg/mL RNase A for one hour at 37  C. Then, 25
μL of 25% sodium dodecyl sulfate (SDS), 25 μL of 20 mg/mL proteinase K, and 75 μL of 5 M NaCl
was added, and incubated for one hour at 65 C. Supernatants were collected and purified through the
addition of phenol-chloroform-isoamyl alcohol (25:24:1; Sigma, St. Louis, MO, USA). DNA was recovered
using the DNA Clean &amp; Concentrator TM -25 columns, as per manufacturer’s instructions (Zymo,
Irvine, CA, USA) and quantified using the NanoDrop (Thermofisher, Burlington, ON). After genomic
DNA extraction, the V3 region of the 16S rRNA gene was amplified via PCR using these conditions per
reaction well: Total polymerase chain reaction volume of 50 μL (5 μL of 10X buffer, 1.5 μL of 50mM
MgCl 2 , 1 μL of 10 mM dNTPs, 2 μL of 10mg/mL BSA, 5 μL of 1 μM of each primer, 0.25 μL of Taq
polymerase (1.25U/ μL), and 30.25 μL of dH 2 O. Each reaction was divided into triplicate for greater
efficiency. The primers used in this study were developed by Bartram et al.,2011. PCR conditions used
included an initial denaturation at 94 C for 2 minutes, followed by 30 cycles of 94 C for 30s, 50 C for 30s,
72 C for 30s, followed by a final elongation at 72 C for 10 minutes. All samples were sequenced using an
Ilumina MiSeq platform at the McMaster Genome Facility (Hamilton, Ontario, Canada). Samples were
processed in batches, meaning not all samples were extracted and sequenced at the same time.


## 16S rRNA gene sequencing processing pipeline
Cutadapt v. 1.14 [53] was used to filter and trim adapter sequences and PCR primers from the raw reads,
using a quality score cut-off of 30 and a minimum read length of 100 bp. We used DADA2 [24] to resolve
the sequence variants from the trimmed raw reads as follow. DNA sequences were trimmed and filtered
based on the quality of the reads for each Illumina run separately. The Illumina sequencing error rates
were detected, and sequences were denosied to produce ASV count table. The sequence variant tables
from the different Illumina runs were merged to produce a single ASV table. Chimeras were removed
and taxonomy was assigned using the DADA2 implementation of the RDP classifier against the SILVA
database v. 1.3.2 [54], at 50% bootstrap confidence. All downstream analysis was conducted in R v. 4.0.3
[55]. We curated the data and generated plots using phyloseq v. 1.22.3 [56] and the following tidyverse
packages: dplyr v. 0.7.6, tidyr v. 0.8.1, rlang v. 0.2.1, and ggplot2 v. 3.0.0. To visualize sample distances
(beta-diversity), we calculated both Aitchison and Bray–Curtis distances. ASV counts transformed to the centered log-ratio (CLR) using microbiome v.1.12.0 and visualized via Principal Component Analysis
(PCA) for Aitchison distances. We applied PCoA to generate Bray-Curtis distances for ordination plots
and unweighted pair group method with arithmetic mean (UPGMA) for clustering trees using ape package
v. 5.2 [57] and hclust() function in R. Trees were visualized using the stringi package v. 1.2.3 in R and
the Interactive Tree of Life (iTOL) [58]. To measure sample diversity, we calculated Shannon values for
sample or group using phyloseq package and visualize them using ggplot. The variability of microbiota
was tested by PERMANOVA on Bray-curtis distances based on relative-abundance of microbes in each
sample using adonis() function in ape package [57]. The diversity of samples calculated by Shannon
values using phyloseq and the significant changes were measured by Linear Mixed-Effects Models using
lmer package. Those ASVs that were present in   1 sample from donor B selected as donor B’s ASVs.
Then, the donor B’s ASVs were compared in each patient with data from prior and post-FMT, ASVs
with a relative abundance of 0 in a patient before FMT and   0.1% post-FMT were labelled as engrafted.
In order to find the commonly engrafted ASVs, the number of engrafted ASVs was compared across an
increasing number of patients in FMT vs. placebo group. To have an equal number of patients across
these two groups, 20 of 31 patients on placebo treatment were randomly sampled (with 100 re-sampling).

## Library preparation and read-based shotgun metagenomics pipeline
