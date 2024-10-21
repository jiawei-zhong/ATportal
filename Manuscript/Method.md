## Data collection

To retrieve publications containing adipose transcriptomic datasets, we mined PubMed with species-related, tissue-related, and modalities-related terms using the following search string: `[(human) AND (adipose) AND (transcriptomics OR microarray OR RNA-seq OR proteomics)]`. Next, we manually curated the obtained datasets according to the “data availability” statement of each candidate publication to confirm whether they relate to human white adipose tissue and contain clinical parameters. For the Cell Type and Perturbation module search strings also included the terms `(hMADS OR SGBS)`. Due to the large amount of in vitro data sets, we manually selected the most relevant ones. For single-cell and spatial module, we obtained the integrated human adipose atlas from recently published studies<sup>1,2</sup>. Whenever available, raw data was downloaded. All included data sets are summarized in Table S1.

## Re-annotation and normalization of microarrays

Gencode v43<sup>3</sup> and FANTOM6<sup>5</sup> were obtained and merged to generate a comprehensive annotation aimed to increase coverage of ncRNAs. For overlapping transcripts, only Gencode entries were kept to avoid redundancy. Next, we remapped the probe sets onto hg38 genome to eliminate wrong and incomplete annotations. Briefly, probe sequences were obtained from the manufacture's websites (Table S2) and aligned to hg38 using STAR v2.7.2a<sup>5</sup>. Mapped probe sequences were assigned to corresponding genes by using BEDtools intersect v2.30.0<sup>6</sup>. A probe set was defined as successfully mapped to a gene when more than 50% of the probes intersected with the exons of the corresponding gene. The probe set ID-gene associated information of each array type can be found in Table S3. For Affymetrix arrays, we used robust multiarray averaging normalization from oligo v1.58.0<sup>7</sup> for normalization. For Illumina arrays, we employed the lumiExpresso function from lumi v2.46.0<sup>8</sup>. To facilitate and condense the read-out in the portal, we calculated the mean expression of all probe sets mapping to the same gene.

## Processing of RNA-seq data

For RNA-seq datasets, raw fastq files were mapped to hg38. Mapped reads were then counted by using HTSeq v2.0.3<sup>9</sup> to generate expression matrices. If raw data was not available, matrices were directly downloaded from GEO, and alias gene names were converted to the gene names used in the Gencode and FANTOM6 annotations. Expression values were quantified as reads per million (RPM). Single cell and spatial data in the portal is presented as in our recent publications<sup>1,2</sup>, however, we updated gene alias information as described above. and recalculated marker genes for each cluster using the FindMarkers function in Seurat v 4.3.0<sup>10</sup>.

## Statistical Analysis

Almost all requests will be calculated on the fly based on user inputs, which, among other things, allows to determine settings affecting the results. For continuous variables (e.g., age or BMI), correlations with genes of interest can be calculated using either Pearson or Spearman and results can be freely adjusted for one or multiple traits (Partial Correlation). Categorical variables (e.g., sex or BMI group) allow comparison with either Student’s T-test or a Wilcoxon Rank Sum Test. For comparison across cohorts and categorical variables, users can decide between standard mean difference (SMD, default), or fold change. In the specific case of depot as well as weight loss comparisons, we only included paired samples, and therefore employ paired statistical test (paired T-test or Wilcoxon signed rank test). For heatmaps clustering methods implemented in the hclust package can be chosen under the ´Customization` tab. In the cell type and perturbation modules, for both heatmaps and line charts absolute (allows better comparison between genes) or normalized (for better comparison across conditions) can be plotted. For volcano plots in these modules, fold changes and p-values were calculated using empirical bayes statistics in the limma v3.50.3 package and queried genes are highlighted. For additional cohorts, where only precalculated correlations could be shared due to restrains in the ethical permits, user queries are limited and allow only the selection of the correlation method and adjustments for any combination of age, sex or BMI. For the Gene Finder module, which allows the identification of transcripts fitting user specified criteria, all results were precalculated. For continuous variables, Spearman’s rho was employed, and for categorical variables either standard mean difference (when multiple cohorts were summarized, as for depot enrichment) or log2 fold change (when analyzing only one cohorts, as for FACS enrichment), was calculated. Adipogenesis regulation was assessed using Hotelling’s T-squared distribution. All p-values were corrected for false discovery rate using Benjamini-Hochberg based on all 132,117 transcripts, and not only coding genes. Users can manually define either p or FDR cutoffs in the input sidebar.

## Implementation

The front end of ATportal is a multi-page web application built using in Vue v3.3.11 using HTML, JavaScript, and CSS code with static information. The responsive Clinical, Depots, Cell type, Single-cell, Spatial and Perturbation modules have been developed using shiny v1.7.4 in R and are deployed on a shiny server v1.5.20.1002. The gene finder module was built using Vue.js v3.3.11. JQuery v3.5.1 and DataTables v1.11.3 are used to generate output tables. For this purpose, we use MySQL v8.0.37 for the database, which is accessed using Flask v3.0.0. Nginx v1.18.0 is used as the reverse proxy server. Currently, ATportal was tested on and supports the following browsers: Google Chrome (v128.0.6613.137), Microsoft Edge (v 128.0.2739.67), and Firefox (v102.5 and above). The AT portal is currently not optimized for mobile access.

## Human samples acquisition

Abdominal subcutaneous white adipose tissue biopsies were obtained under local anesthesia from two separate clinical cohorts, NEFA (NCT01727245) and ADIPO-SCAPIS (ethical approval from the Swedish Ethical Review Authority Dnr 2024-01590-02).

## Proteomics

### Sample preparation

Frozen samples were homogenized in 2% SDS buffer (2% SDC, 100 mM Tris-HCl pH=8.5) using a homogenizer and boiled at 95°C, 5min. After determining protein concentration via BCA assay (Thermo, 23225), 25 µg protein of each sample were prepared and frozen for proteomics processing. Protein were later reduced and alkylated using 10 mM TCEP and 40 mM CAA at 40°C in the dark for 10min and digested overnight with 1:50 (protein:enzyme) of trypsin (Sigma, t6567) and LysC (Wako, 129-02541). The next day, samples were acidified by adding 1:1 (v:v) of isopropanol, 2% TFA. After centrifugation for 10 min at 15,000 g, supernatants were loaded onto activated triple layer styrenedivinylbenzene–reversed phase sulfonated STAGE tips (3M Empore). Peptides were washed with 100 µl ethylacetate 1% TFA, 100 µl 30% Methanol 1% TFA and 150 µl 0.2% TFA and eluted with 60 µl elution buffer (80% ACN, 5% NH4OH). Peptides were lyophilized and dissolved in 10 µl MS loading buffer (2% ACN, 0.1% TFA).

### LC-MS/MS

LC-MS/MS analysis 500 ng of peptides was performed on an Orbitrap Exploris 480 (Thermo Fisher Scientific) equipped with a nano-electrospray ion source and FAIMS (CV50) coupled to an EASY-nLC 1200 HPLC (all Thermo Fisher Scientific). Peptides were separated at 60°C on 50 cm columns with an inner diameter of 75 μm packed in-house with ReproSil-Pur C18-AQ 1.9 μm resin (Dr.Maisch GmbH) over 1 h by reversed-phase chromatography using a binary buffer system consisting of buffer A (0.1 formic acid) and buffer B (80% ACN, 0.1% formic acid). Starting with 5% of buffer B, this fraction was increased stepwise to 45% over 45 min followed by a wash-out at 95%, at a constant flow rate of 300 nl/min. Peptides were ionized and transferred from the LC system into to the gas phase using electrospray ionization (ESI). One MS1 scan (300-1650 m/z, max. ion fill time of 45 ms, normalized AGC target= 300%, R= 120.000 at 200 m/z) was followed by 66 MS2 fragment scans of unequally spaced windows (fill time= 22 ms, normalized AGC target = 1000%, normalized HCD collision energy= 30%, R= 15.000). Spectra were acquired in profile mode using positive polarity.

### Data processing

DIA raw files were analyzed using Spectronaut software (v. 18.4.231017.55695, developed by Biognosys AG in Schlieren, Switzerland) with directDIA and searched against the Uniprot human databases: UP000005640_9606 and UP000005640_9606_additional with standard processing parameters (trypsin cleavage with a peptide length ranging from 7 to 52 amino acids, two missed cleavages). Fixed modification settings included carbamidomethylation, variable modifications were methionine oxidation and N-terminal acetylation. The analysis specified a minimum of 3 and a maximum of 6 Best N Fragment ions per peptide. For filtering and quality control, a precursor and protein q-value cutoff of 1% was applied.

### Data analysis

All analysis was conducted in the R statistical computing environment v.4.2.3.R. Before analysis, contaminants and reverse decoy database hits were excluded. Proteins present in less that 25 samples and 5 samples not meeting quality control criteria were excluded from further analysis. To identify blood contamination, single cell cluster markers were exported from the portal (adipose tissue markers) and the list provided by Geyer et al.<sup>11</sup> (blood origin markers) was downloaded. A spearman’s correlation matrix was created for the protein expression of those markers, flowed by ward.D2 hierarchical clustering. Clusters with proteins with high positive correlation with one group and high negative correlation with the other group were used to create a quality control protein list (Table S2). This list was correlated with the proteins present in the data set, and hierarchical clustering was able to pinpoint a distinct cluster of blood contaminants that was removed from further analysis. The remaining proteins were normalized using variance stabilizing normalization. Differential protein expression was examined using two-way Anova with FDR p-value adjustment. Pathway analysis was conducted using the PathfindeR package and data from MSigDB and IntAct databases, while GO enrichment was performed using ClusterProfiler. All plots were generated using ggplot2, with the exception of clustered heatmaps that were visualised using pheatmap.


## Reference
1.	Bäckdahl, J. et al. Spatial mapping reveals human adipocyte subpopulations with distinct sensitivities to insulin. Cell Metab 33, 1869-1882.e6 (2021). 

2.	Massier, L. et al. An integrated single cell and spatial transcriptomic map of human white adipose tissue. Nat Commun 14, 1438 (2023).

3.	Frankish, A. et al. GENCODE reference annotation for the human and mouse genomes. Nucleic Acids Res 47, D766–D773 (2019).

4.	Ramilowski, J. A. et al. Functional annotation of human long noncoding RNAs via molecular phenotyping. Genome Res 30, 1060–1072 (2020).

5.	Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15–21 (2013).

6.	Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841–842 (2010).

7.	Carvalho, B. S. & Irizarry, R. A. A framework for oligonucleotide microarray preprocessing. Bioinformatics 26, 2363–2367 (2010).

8.	Du, P., Kibbe, W. A. & Lin, S. M. lumi: a pipeline for processing Illumina microarray. Bioinformatics (Oxford, England) 24, 1547–1548 (2008).

9.	Anders, S., Pyl, P. T. & Huber, W. HTSeq--a Python framework to work with high-throughput sequencing data. Bioinformatics 31, 166–169 (2015).

10.	Hao, Y. et al. Integrated analysis of multimodal single-cell data. Cell 184, 3573-3587.e29 (2021).

11.	Geyer, P. E. et al. Plasma Proteome Profiling to detect and avoid sample-related biases in biomarker studies. EMBO Mol Med 11, e10427 (2019).
