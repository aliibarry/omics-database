---
title: "codedata"
author: "amb"
date: "23/08/2021"
output: html_document
---

Last updated 04/09/2023, AMB

---

**Reports**

> 03.03.2023: added network integration with String DB  
> 03.12.2022: fixed separation by sex, added "clear" button for gene search by file upload.   
> 05.31.2023: added Predictive Pain Score to network integration 

---

**Overview**

User Notes: 
* Searchable by gene symbol
* File uploads for queries available, need to press "clear" between searches if file is uploaded.
* May need to adjust window size manually to ensure they are a useful size on the screen, but downloaded plots will be the correct width.

---

**Data Tables**

Tables are available for download for independent datasets, with group-level data provided for human data. 

---

**Useful Resources**

_Code_  
This app was generated using Shiny ('R') and is hosted through the Interactive Data Network (IDN) at the University of Oxford. 
* A free bookdown tutorial is available at: https://aliibarry.github.io/database-book/  
* Code for this app can be found at: https://github.com/aliibarry/omics-database  
* The code used to generate a previous version of this website is also freely available [here](https://github.com/aliibarry/shiny)  

_Data Interpretation_  
The presence/absence of a gene will depend on various factors as will the likelihood of a gene to reach a significance threshold. When interpreting these data, please read the original manuscript methods, and consider the following points:  
* Sample size (see an example power calculator   [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2191-5)
* Cell type composition  
* Library preparation  
* Sequencing method and depth  
* Count filtering and normalization  
* FDR weighting  (see  [Korthauer et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1716-1))  
* Effect size moderation (eg. [ashr](https://rdrr.io/cran/ashr/) and [apeglm](https://doi.org/10.1093/bioinformatics/bty895))  

Experimental factors may also account for differences, eg:   
* Species  
* Sex  
* Age  
* Model 
* Time 
* Tissue   

For example:
- In Barry _et al_ (2023), differentially expressed genes are calculated 3 days and 4 weeks after SNI in _male and female mice_ (L3-L5). This study uses _8 biological replicates_ per group and FACS purified sensory neuron subtypes. Results were calculated using the Wald test and a *weighted* FDR correction (independent hypothesis weighting, `IHW`). Effect sizes were calculated using Bayesian shrinkage estimators (the `apeglm` method, via `DESeq2`) and are presented as *moderated* log2 fold changes to account for the variability in lowly expressed genes.  
- In Perkins _et al_ (2014), differentially expressed genes are calculated 7 days after SNT in _male rats_ (L5). This study uses _3 biological replicates_ per group on whole DRG. Results were calculated using the Wald test and an *unweighted* FDR correction (Benjamini and Hochberg). Effect sizes were *not moderated* and are presented as log2 fold changes. 

---


**Other Relevant Databases**

As a field, various other databases exist to host relevant -omics data. Please see respective webpages for their associated citations, as these continue to the updated. We to not take responsibility for the data hosted on other sites. 

* [Sensoryomics Web](https://paincenter.utdallas.edu/sensoryomics/)
* [Sensoryomics Shiny App](https://sensoryomics.shinyapps.io/RNA-Data/)
* [Painseq Shiny App](https://painseq.shinyapps.io/tg-painseq/)
* [HPGDB Web](https://humanpaingeneticsdb.ca/hpgdb/)
* [MPI-EM Pain Proteome](http://painproteome.em.mpg.de/)
* [Neuro-immune RNA-seq](https://rna-seq-browser.herokuapp.com/)  
* [Broad Institute](https://singlecell.broadinstitute.org/)

---

**Acknowledgments**

Special thanks to Dr. Lucy McDermott for workshopping url names, and to the rest of the Neural Injury Group for feedback as beta users. This website is hosted by the Interactive Data Network (Oxford).

---

**Funding**

This work was funded in part by the Wellcome Trust (DPhil scholarship to AMB, 215145/Z/18/Z) and a Wellcome Investigator Grant to DB (223149/Z/21/Z), as well as the MRC (MR/T020113/1), and with funding from the MRC and Versus Arthritis to the PAINSTORM consortium as part of the Advanced Pain Discovery Platform (MR/W002388/1). AMB further received a GTC MSDTC Scholarship. GB is funded by Diabetes UK, grant number 19/0005984, MRC and Versus Arthritis through the PAINSTORM consortium as part of the Advanced Pain Discovery Platform (MR/W002388/1) and by the Wellcome Trust (223149/Z/21/Z). 

---