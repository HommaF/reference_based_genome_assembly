## De-novo genome assembly and structural annotation using published genomic reads and RNAseq data

graph TD;
	De-novo, reference-guided, whole genome assembly --> Scaffolding --> structural annotation;

# De-novo, reference-guided, whole genome assembly
*Genomic sequencing data*
| Species | SRA ID | Source |
| ------- | ------ | ------ |
| S.peruvianum | ERR418094 | insert paper |

*Workflow based on [Lischer et al. 2017](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1911-6)*
graph TD;
