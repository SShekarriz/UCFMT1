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
We have conducted shotgun metagenomics on 22 samples collected from 11 patients, with 2-time points
each, in this study (4 non-responder, 6 responders patients who received FMT and one patient on placebo).
Genomic DNA was standardized to 5 ng/ μL and sonicated to 500 bp. Using the NEBNext Multiplex
Oligos for Illumina kit (New England Biolabs), DNA ends were blunted, adapter ligated, PCR amplified,
and cleaned as per manufacturers instructions. Library preparations were sent to the McMaster Genome
Facility, and sequenced using the Illumina HiSeq platform. The forward and reverse sequencing runs
were concatanated and trimmed for primer adapters and low quality reads using Trimmomatic [59]. The
taxonomic, and gene-family composition of trimmed shotgun reads identified using Metaphlan2 and
Humann2 pipeline [42]. All downstream analysis was conducted in R v. 4.0.3. The Bray-Curtis distances
calculated based on the relative abundance of known species and gene-families using phyloseq package.
Principal coordinate analysis (PCoA) plots were generated using phyloseq and ggplot2. Unweighted pair
group method with arithmetic mean (UPGMA) trees based on Bray–Curtis distances were generated
using the ape package v. 5.2 [57] and hclust() function in R. Trees were visualized using the stringi
package v. 1.2.3 in R and the Interactive Tree of Life (iTOL). [58]. The diversity of samples measured by
Shannon index using phyloseq and the significant changes were measured by Linear Mixed-Effects Models
using lmer package. For the microbial taxonomy dataset, the engrafted strains were defined as any strains
present in   1 sample from donor B with relative abundance of 0 prior to FMT and 0.1% post-FMT in a
patient. Humann2 uses a detection threshold of 0.01% relative abundance which is equivalent to 0.1x
fold-coverage of a 5 Mbp microbial genome [42]. Given 1000 genes per Mbp, we expect 0.0005% relative
abundance for detection of a gene family. Thus, any gene family with a minimum relative abundance of
0.0005% in donor B samples, 0% before FMT and 0.0005% post- FMT defined as engrafted gene-families.

## Culture-enriched and independent metagenomics on donor B samples
A single fresh, anaerobic fecal sample collected from donor B. The collected sample was cultured using
33 media, and incubation of plates anaerobically and aerobically resulted in 66 culture conditions for
culture-enriched molecular profiling using previously described protocol [17]. The list of media and culture
conditions are described earlier [17]. 16Sr RNA amplicon sequencing was conducted on all the 66 culture
conditions. To determine a representative subset of culture-enriched plates that adequately represent the
sample, the distribution of ASVs in the direct sequencing was compared to the culture-enriched sequencing
per plate pool using the PLCA algorithm [20]. Shotgun metagenomics was conducted on the subset of
plate pools identified by the PLCA algorithm. Genomic DNA was isolated from the thirteen selected plate
pool and shotgun metagenomics conducted as previously described [17, 20]. Direct shotgun metagenomics
conducted on the same fecal sample, which was earlier used for culture-enriched metagenomics as well as
three other fecal samples collected from donor B at different time points (2013, 2017x2).

## Comparison of the culture-enriched metagenomics with direct metagenomics
To build the culture-enriched metagenomic library, the raw shotgun sequences from the selected plate
pools and the original fecal sample collected from donor B co-assembled together as follows. The lowquality
reads and sequencing primers removed using Trimmomatic [59]. The reads decontaminated for
any human DNA utilizing DeconSeq package [60]. The shotgun reads were co-assembled and binned
using metaSPADE [61] and Metabat2 [62] respectively. In addition to CEMG assembly, the microbial
composition of direct metagenomics (DMG) from the fecal sample assembled and binned separately.
These two datasets are labelled as CEMG and DMG in Figure 3. The microbial composition of DMG
and CEMG datasets were then comprehensively evaluated using the following procedure. The single-copy
core genes were identified within each bin using CheckM [63], any bin with minimum 70% completion and
maximum 10% contamination were defined as a metagenome assembled genome (MAG). The shotgun
reads were mapped to the assembled contigs to estimate sequence coverage for all contigs, Bins, MAGs,
and those contigs that were not present in any bin. We used bwa [64] to map reads to assembled contigs
and anvio pipeline [65] to normalize the coverage to the depth of sequencing. The detection values
calculated for each bin using anvio package [65]. The detection value defined as the proportion of a given
MAG that is covered at least 1X; in other words, it estimates the proportion of MAG that recruited reads
to it. We used GTDB-Tk [66] to build a phylogenetic tree, and taxonomic assignment of MAGs. All of
the figures visualized in R v. 4.0.3.

## Microbial engraftment detection in metagenomic data
To detect microbial engraftment, we aimed to construct a comprehensive library of microbial genes and
genomes from donor B. This library contains 4 DMG samples and single CEMG sequencing methods. The
low-quality reads and sequencing primers removed using Trimmomatic [59]. The reads decontaminated
for any human DNA utilizing DeconSeq package [60]. The shotgun reads from both culture-dependent
and independent libraries were co-assembled and binned using metaSPADE [61] and Metabat2 [62]
respectively. The MAGs and MABs were identified using previously described criteria. The single-copy
core genes were identified within each bin using CheckM [63], any bin with minimum 50% completion and
maximum 10% contamination were defined as a metagenome assembled genome (MAG). To include more
number of MAGs in our database, we reduced the completion value of a MAG from 70% to 50%. We used
Prodigal [67] to predict prokaryotic genes and coding DNA sequences (CDS) from the assembled contigs.
The taxonomic labels are assigned to all bins using GTDB-Tk [66]. In total, we were able to assemble
255 metagenomics assembled genomes (minimum 50% completion and maximum 10% contamination)
and 1,130,000 completed prokaryotic genes. After de-novo prediction of genes and MAGs, we mapped
the collected metagenomics samples from before and after FMT (22 samples in total) to the assembled
contigs from donor B. The raw reads from each sample were mapped to the assembled contigs using bwa
mem [64] and the coverage information normalized to the depth of sequencing using anvio package. For
each MAG, the detection and single nucleotide variability measurements calculated using anvio pipeline.
The variability index shows the number of reported single-nucleotide variants per kilobase pair. All the
downstream analyses to detect microbial engraftment at gene and genomic-level were performed in R v.
4.0.3 R.

The assembled MAGs from donor B were classified by comparing the short read mapping coverage
and SNV frequencies from before and after FMT for each patient. Shared category was defined as MAGs
covered above our minimum detection cutoff (  60% proportion of nucleotides in a MAG that has at
least 1X coverage) in a patient both before and after FMT. The MAGs with coverage lower than the
minimum detection cutoff in both time points were classified as Unique to Donor. The MAGs with
0 coverage before FMT and   60% post-FMT were classified as Engrafted and opposite cutoffs were
used for Lost category. To measure variability across detected MAGs, the SNV frequency calculated for
MAGs with minimum 0.6% coverage in both time points using anvio pipeline [65]. The SNV frequency
shows the number of single-nucleotide variants per kilo base pair. The Shared MAGs that showed >=1
SNV per kilobase pair before FMT but their frequencies reduced to <= 0.5 per kbp after FMT were
classified as Strain Replacement. In other words, those MAGs that were present in both time points
but highly similar to donor B’s original MAG only after FMT are defined as strain replacement.

To detect microbial engraftment at the gene-level, we compared the coverage of all the 1,130,000
donor B microbial genes across UC patients before and after FMT. In this model, those genes that their
detection (% of gene covered at least 1X) was 0 before FMT and became at least 0.6 with minimum 5X
coverage after FMT is called engrafted genes (Figure 5 C). To narrow down the number of engrafted
genes, commonly engrafted ones across three patients were labelled as common engraftment. This model
applied to all eleven patients, regardless of their response to FMT. Then, we compared these commonly
engrafted genes against the Uniref90 reference database using Diamond blastp and the identified Uniprot
ids annotated with GO, KEGG, COG, PFAMs, and lineage information.

## Single whole-genome sequencing and comparative genomics
30 Dorea, 1 F.prausnitzii, and 67 Blautia strains were isolated from human gut. The media, culture
conditions for isolation, library preparation, and sequencing protocols as described earlier [68]. In addition,
we have collected 65 Dorea, 98 F.prausnitzii, and 143 Blautia strains available in NCBI RefSeq (May
2020). We have annotated all the genes and CDS using Prokka [69] with default settings. The assembled
genomes were re-classified using GTDB-Tk [66] and phylogenetic trees were constructed within each
genus based on multiple sequence alignment of 120 bacterial marker genes from GTDB database [36]. We
used panaroo [70] with strict mode and mafft aligner to generate core-gene alignment within each species.
We then used approximately-maximum-likelihood model via FastTree [71] to construct phylogenetic trees
for strains within each species. We made a blastn database using all the 402 genomes and tracked the
commonly engrafted genes across these genomes with a minimum  90 pident and qcovhsp  90 cut-offs.
The number of non-redundant positive hits from blastn output were visualized for each genome on the
phylogenetic trees for genus and species collections. All the phylogenetic trees were visualized in v. 1.2.3
R using ggtree, ggtreeExtra, and ape packages.

A single genome with the most number of commonly engrafted genes and fewest contigs were selected
for D. longicatena, F. prausnitzii, and F. saccharivorans as representative strains of commonly engrafted
genes. Then we mapped all the shotgun raw-reads from donor B (5 samples) and UC patients (22
samples) to these three genomes using bwa-mem [64]. The commonly engrafted genes identified for each
representative strain using previously described gene engraftment model. Briefly, genes that were not
present (0% 1X coverage) pre-FMT but present (>0.6% with  5X coverage) post-FMT across  3 patients
were selected as commonly engrafted genes. We used a custom python code to extract all the commonly
engrafted genes and their flanking regions (20,000 bps) from the three representative strains. To find
whether the engrafted genes are the result HGT or strain replacement, we used anvio pipeline [65] for
de-novo characterization and reporting of SNVs for the two selected flanking regions in F.prausnitzii
and F. saccharivorans strains. In short, a table of nucleotide base frequencies for the 80,000 bp gene
clusters, contained commonly engrafted genes, were constructed for F. saccharivorans and F.prausnitzii
representative strains. The consensus nucleotide identified based on anvio’s conservative heuristic model.
We then selected and visualized only those base positions that were identical across all donor B samples
(5 samples). These are donor B’s specific SNVs that we used to see whether the samples collected after
FMT are closer to the donor’s SNV profile (less number of SNV) or contain increased SNVs. The SNV
tables were filtered and visualized in v. 1.2.3 in R using tidyverse packages.

## Species- and strain-specific markers for a metagenomic survey of IBD patients
## and healthy controls

To build species-specific marker, a group of genomes from distinct species were selected from Dorea
sp., Feacalibacterium sp., and Fusicatenibacter sp. collections for pangenome analysis. Species-specific
core-genomes were identified and visualized using Anvio microbial pangenomics workflow. Brierly, gene
calls were annotated with prodigal and MCL algorithm [72] was used to identify gene-cluster across
the pangenome alignment with –mcl-inflation 10 –minbit 0.5 –use-ncbi-blast. To test the accuracy of
species-specific and strain-specific (the commonly engrafted genes) markers, we have used 1200 WGS from
our lab strain collections. These strains are diverse bacterial isolates from all bacterial phyla collected
from the human microbiota. We mapped shotgun reads from 1112 WGS to markers using bwa-mem
[64] with –B 40 –O 60 –E 10 —L 100 parameters to find perfectly aligned reads over their entire length
and samtools [73] to extract coverage information from bam file. We then used the percentage of 1X
coverage to visualize the coverage information for each marker in v. 1.2.3 in R using tidyverse packages.
We used a publicly available metagenomic dataset (PRJNA279196 [39]) to investigate the specificity of D.
longicatena, F. prausnitzii, and F. saccharivorans strains in IBD patients compared to healthy controls.
We downloaded metagenomic samples from the SRA database via Entrez Direct (EDirect) command line.
Metagenomic shotgun reads from all the samples (n=220) were mapped to marker gene-clusters using
bwa-mem [64] with the parameters specified above and samtools [73]. Subsequently, the percentage of 1X
coverage of strain- and species-specific markers were visualized in v. 1.2.3 in R using tidyverse packages.

