# HSrats-ADDO-sQTL

This is a pipeline used for the splicing QTL analysis based on the RNA-seq data of HSrats fat and liver tissues.

STEP1. Convert the gene/transcript expression levels to the input format required by ADDO.

STEP2. eQTL mapping based on ADDO using the gene/transcript expression levels of two tissues.

STEP3. Summary the peak eQTLs from the ADDO GWAS results.

STEP4. Summary the corresponding pairs between gene-based eQTLs and transcript-based eQTLs and find the antagonistic/synergistic eQTL pairs.

STEP5. Summary the proportions of different eQTL types and calculate the significance of dominance enrichment among trans-eQTLs.
