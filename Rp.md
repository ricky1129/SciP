# Title

Ruiqi Shi

---


## 0 Introduction

### 0.1 Background

Coeliac disease is a immune-related disorder of the small intestine, triggered by the ingestion of gluten in specific individuals. The disease mostly affects the duodenum where immune action can lead to significant mucosal damage. The typical cellular identities of the disease include scalloping of the small bowel folds (pictured), fissures, a mosaic pattern to the mucosa, prominence of the submucosa blood vessels, and a nodular pattern to the mucosa[@1]. These cause damage to the normal absorptive functions of the small intestines, and will consequently lead to nutritional deficiencies, which shows as multiple symptoms including abdominal pain, cramping, bloating with abdominal distension, and mouth ulcers. Althought the causes of coeliac diasease is mostly overreacting immune response to gluten, the progression also involves changes in the functions of epithelial cells in the gut. 

The small intestine mucosa is a highly complex and dynamic tissue, composed of several interacting cell groups. As the primary interface between the lumen and the mucosal tissue, the epithelium consists of several specialized cell types including absorptive enterocytes, mucus-producing goblet cells, immune cells (T cells, B cells, Plasma cells, etc.), stem cells and so on. These cells, while performing differenciated functions, all coordinate as a whole within a tightly regulated range to maintain homeostatis; therefore the functional state of the alimentary canal depends deeply on the coordinated responses of different cell types to environmental as well as inflammatory signals. 

### 0.2 Motivations

Although the roles of epithelial and immune systems in the mechanism of coeliac diseases are investigated, it is not well known how stromal compartments contribute in this situation. The intestinal stroma contains a diversed population of cells including fibroblasts, endothelial cells, etc., which contribute to mechanical support as well as epithelial regeneration and communication with immune cells. 

However, the extent to which the changes to stromal cells affect the physcial damages in epithelial cells remain unclear - in particular, whether the disease-associated stromal transcript programs are accompanied by changes in epithelial signalling. This shows an huge limitation in the current understanding of coeliac diseases, as the study about epithelial and immune systems independently might ignore the significance of multi-cell interactions. 

Under this frame, single-cell RNA sequencing (scRNA seq) provides a great opportunity to reveal the transcriptional state of different cell categories and thus investigate their behaviour. This method allows disease-related changes in cellular components be analysed within the scope of a whole. 

### 0.3 Work done

In order to investigate the potential contribution of stromal compartment to epithelial damage, 2 steps are taken. 

#### 0.3.1 Identification of disease-related gene expression programs

The first stage is aimed at identifying specific genes that is differed in expression in coeliac and control groups, so that we can know what genes affect the damaging process of epithelium. 

The analysis focuses on identifying gene expression programs that are associated with intestial barrier integrity and examining their relationships with stromal-derived signalling factors. 

During the process several analysis tools are used: the [NCBI GEO database](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277276) provides data of gene expression in multiple cell groups related to the disease. R package `seurat` is used as a frame to store the data, as well as providing key functions for data analysis such as `logNormalize`, `RunUMAP`, `RunPCA`, etc. By selecting the differently expressed genes among the cells, we can cluster and label cells by their functions, and then we are able to see the difference in expression between infected cells of coeliac disease and healthy ones by using `cnmf` which identifies specific 'programs' in gene expression. 

### 0.3.2 Experimental validation of the genes found

The second stage aimed to test potential effects of the selected genes on the integrity of epithelial barrier. Using a controlled model created with Caco-2 cells (*H. Sapiens* cancer cells of the colon), we investigate the changes in epithelial barrier function. 





[@1]: Calado J, Verdelho Machado M (2022). https://karger.com/pjg/article-pdf/29/2/111/3762388/000514716.pdf. GE - Portuguese Journal of Gastroenterology. 29 (2): 111–124.