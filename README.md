# Preprocessing and downstream analysis of scATAC-seq data

Authors: Pia Baronetzky, George Tsitsiridis, Niklas Kemper, Markus Franke
Supervisor: Maria Colomé-Tatché
Advisor: Maria Richter

Author: Pia Baronetzky, George Tsitsiridis, Niklas Kemper, Markus Franke  
Supervisor: Maria Colomé-Tatché  
Advisor: Maria Richter

The progressing development of experimental techniques to investigate chromatin accessibility in single cells demands for suitable computational tools to analyze this data type. Recently, methods have emerged to do so, mainly focusing on quality control and how to best cluster the cells based on their shared chromatin openness. However, in a PCA, it can often be observed that the first component corresponds to library size, which needs to be normalized for. In scRNA-seq data, it has been shown that simple library-size normalization is not sufficient to distinguish technical from biological variation.

To handle this, a plethora of normalization techniques have been developed for scRNA-seq data, improving the ability to differ biological from technical variance. For example, the theory behind one prominent technique is that cells from a cell type will have similar library sizes. In line with that, it has been observed that cells also differ in their chromatin openness, e.g. stem cells showing overall higher chromatin openness than differentiated cells.

Therefore, this project aims to develop a better normalization technique for scATAC-seq data, similar to what is available for scRNA-seq data. We hypothesize that with suitable normalization for scATAC-seq data, we improve the ability to differentiate between cell states and to confidentially call differential openness between cell states. We hope that this can eventually lead to a better understanding of which genes start to be differentially open between cell states before differential gene expression occurs.

**Goals**: Improve library size normalization for single cell ATAC-seq data to reduce technical noise while retaining biologically relevant information

**Data set**: scATAC-seq

**Methods**: Application and adaptation of existing normalization algorithms; potential development of new algorithms

**Literature**:  

Alexandre Gaspar-Maia et al., “Open Chromatin in Pluripotency and Reprogramming,” Nature Reviews Molecular Cell Biology 12, no. 1 (January 2011): 36–47, https://doi.org/10.1038/nrm3036.

Anna Danese et al., “EpiScanpy: Integrated Single-Cell Epigenomic Analysis,” BioRxiv, May 24, 2019, 648097, https://doi.org/10.1101/648097.

Tim Stuart et al., “Multimodal Single-Cell Chromatin Analysis with Signac,” BioRxiv, November 10, 2020, 2020.11.09.373613, https://doi.org/10.1101/2020.11.09.373613

Jeffrey M. Granja et al., “ArchR: An Integrative and Scalable Software Package for Single-Cell Chromatin Accessibility Analysis,” BioRxiv, April 29, 2020, 2020.04.28.066498, https://doi.org/10.1101/2020.04.28.066498

Beate Vieth et al., “A Systematic Evaluation of Single Cell RNA-Seq Analysis Pipelines,” Nature Communications 10, no. 1 (October 11, 2019): 4667, https://doi.org/10.1038/s41467-019-12266-7

Aaron T. L. Lun, Karsten Bach, and John C. Marioni, “Pooling across Cells to Normalize Single-Cell RNA Sequencing Data with Many Zero Counts,” Genome Biology 17, no. 1 (April 27, 2016): 75, https://doi.org/10.1186/s13059-016-0947-7.
