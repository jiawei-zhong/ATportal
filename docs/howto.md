# How to use and navigate the Adipose Tissue Knowledge Portal

## Gene Summary

The easiest way to use the AT portal, is to enter a gene name into the search bar on the [starting page](adiposetissue.org). This will open the Gene Summary Module (which has the same search bar), and provide a comprehensive overview across all available modules. 

<p align="center">
  <img src="https://raw.githubusercontent.com/jiawei-zhong/ATportal/main/img/SearchBar.PNG">
</p>
<img align="left" src="https://raw.githubusercontent.com/jiawei-zhong/ATportal/main/img/SummarySide.PNG">

The sidebar on the left side displays key information for the queried gene.

| Section | Description |
| --- | --- |
| Gene Symbol | Provided Gene Symbol, Links out to NCBI Gene Summary |
| Gene Name | Official gene name from the HUGO Gene Nomenclature Committee (HGCN) |
| Source | Links to the HGNC |
| Known as | Gene aliases |
| Description | Gene description from NCBI Gene |
| External Resources | Links to <ul><li>Human Protein Atlas</li><li>T2D Knowledge Portal</li><li>LD Knowledge Portal</li></ul> |
| Download | Download All Summay Files as PDF |

<p> All results in the summary section are weighted analysis across all available cohorts for a certain trait, and both proteome and transcriptome data are shown if available. The first 4 sections are from the Clinical module. Under the **Clinical correlation** section, correlation with clinical parameters, separated into 3 groups (antropometric, circulating, tissue specific) is shown, where the Spearman correlation is shown on the x axis, significance level is indicated by asterisks and the color shows the amount of available cohorts. **Sex differences** shows a meta-analysis forestplot for all available cohorts to the left, which indicates both common and random effect size. To the right, a histgram shows the effect size for the queried gene in comparison to all other genes in the AT portal, to allow a better assessment of the effect size. The **Obsese vs. non-obsese** section is build in the manner, but is comparing people living with or without obesity based on the respective BMI cutoffs (<30 kg/m2 vs. >= 30 kg/m2). The **weight loss** section shows barplots for 2 surgery induced and 2 diet induced weight loss cohorts. 
**Depots** allows the comparison of transcription levels in subcutaneous vs omental adipose tissue. Again, a meta-analysis forest plot over the 6 available cohorts is shown.
Next, the **Tissue specificity** section summarizes data from the FANTOM consortium, to compare expression in adipose tissue compared to other human organs and tissues. As the focus here is on adipose tissue, data from whole adipose tissue ("adipose") and isolated mature adipocytes ("adipocytes") is highlighted. The portal also contains a growing variety of **perturbation** data sets. For the summary, fold changes and FDR adjusted p-values across all pertubations (different cell lines and modalities) are shown in one volcano plot to facilitate further explorations in the detailed module. The last section summarizes **Single Cell** data. In the left panels, 2 UMAPs indicate major cell lineages and where the gene of interest is expressed. The information is also summarized in violin charts on the right. The single cell module allows further analysis.  </p>

<p align="center">
  <img src="https://raw.githubusercontent.com/jiawei-zhong/ATportal/main/img/Summary.gif">
</p>
---

## Clinical Relationships

The following image give a general overview about the structure of the modules. 
 <img src="https://raw.githubusercontent.com/jiawei-zhong/ATportal/main/img/Layout_PDF.png" width=600px>

