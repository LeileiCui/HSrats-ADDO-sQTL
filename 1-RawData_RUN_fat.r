
rm(list=ls());options(stringsAsFactors=FALSE)

phe_gene = read.table("fat_gene_expression_rint.txt", header=T)
phe_trpt = read.table("fat_transcript_expression.txt", header=T)
covs = read.table("fat_covariates.txt", header=T)
geno_id = read.table("pruned.hs.chrall.fam", header=F)

colnames(phe_gene) = gsub("X", "rat_", colnames(phe_gene))
colnames(phe_trpt) = gsub("X", "rat_", colnames(phe_trpt))
colnames(covs) = gsub("X", "rat_", colnames(covs))
geno_id[,1] = paste0("rat_", geno_id[,1]); geno_id[,2] = paste0("rat_", geno_id[,2])

common_id_gene = intersect(intersect(colnames(phe_gene), colnames(covs)), geno_id[,1])
common_id_trpt = intersect(intersect(colnames(phe_trpt), colnames(covs)), geno_id[,1])

rownames(covs) = covs[,1]; covs = covs[,-1]

out_phe_gene = cbind("id" = common_id_gene, t(covs[,common_id_gene]), t(phe_gene[,common_id_gene]))
out_phe_trpt = cbind("id" = common_id_trpt, t(covs[,common_id_trpt]), t(phe_trpt[,common_id_trpt]))
write.table(out_phe_gene, file="HSrats_fat_Genes.phe", row.names=F, col.names=T, quote=F)
write.table(out_phe_trpt, file="HSrats_fat_Trpts.phe", row.names=F, col.names=T, quote=F)

out_geno = as.data.frame(gsub("rat_", "", common_id_gene))
write.table(cbind(out_geno, out_geno), file="tmp_HSrats_fat_Genes.txt", row.names=F, col.names=F, quote=F)
system("plink --bfile pruned.hs.chrall --keep tmp_HSrats_fat_Genes.txt --make-bed --out HSrats_fat_Genes")
write.table(cbind(common_id_gene, common_id_gene, 0, 0, 0, -9), file="HSrats_fat_Genes.fam", row.names=F, col.names=F, quote=F)

out_trpt = as.data.frame(gsub("rat_", "", common_id_trpt))
write.table(cbind(out_trpt, out_trpt), file="tmp_HSrats_fat_Trpts.txt", row.names=F, col.names=F, quote=F)
system("plink --bfile pruned.hs.chrall --keep tmp_HSrats_fat_Trpts.txt --make-bed --out HSrats_fat_Trpts")
write.table(cbind(common_id_trpt, common_id_trpt, 0, 0, 0, -9), file="HSrats_fat_Trpts.fam", row.names=F, col.names=F, quote=F)

out_covs_gene = cbind("measure"=rownames(phe_gene), "covariates"=paste(rownames(covs), collapse=","))
write.table(out_covs_gene, file="HSrats_fat_Genes.covs", row.names=F, col.names=T, quote=F)

out_covs_trpt = cbind("measure"=rownames(phe_trpt), "covariates"=paste(rownames(covs), collapse=","))
write.table(out_covs_trpt, file="HSrats_fat_Trpts.covs", row.names=F, col.names=T, quote=F)

