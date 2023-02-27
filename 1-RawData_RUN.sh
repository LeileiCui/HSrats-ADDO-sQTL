cp /SAN/mottlab/HSrats/RNASEQ/Pipeline1/Ref7.2/Liver/eqtl_input/transcript_expression.txt ./
mv transcript_expression.txt liver_transcript_expression.txt 
cp /SAN/mottlab/HSrats/RNASEQ/Pipeline1/Ref7.2/Fat/eqtl_input/transcript_expression.txt ./
mv transcript_expression.txt fat_transcript_expression.txt 

cp /SAN/mottlab/HSrats/RNASEQ/Pipeline2/Ref7.2/5_eqtl/Liver/eqtl_input/gene_expression_rint.txt ./
mv gene_expression_rint.txt liver_gene_expression_rint.txt 
cp /SAN/mottlab/HSrats/STACY/Res_Pipeline2/5_eqtl/Fat/eqtl_input/gene_expression_rint.txt ./
mv gene_expression_rint.txt fat_gene_expression_rint.txt 

cp /SAN/mottlab/HSrats/QTLmapping/HAPPY/run_HAPPY/chrall/pruned.hs.chrall.* ./

cp /SAN/mottlab/HSrats/STACY/Res_Pipeline2/5_eqtl/Fat/eqtl_input/covariates.txt ./
mv covariates.txt fat_covariates.txt

cp /SAN/mottlab/HSrats/RNASEQ/Pipeline2/Ref7.2/5_eqtl/Liver/eqtl_input/ori.covariates.un_transposed.txt ./
mv ori.covariates.un_transposed.txt liver_covariates.txt

