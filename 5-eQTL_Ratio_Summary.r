
### This script aims to summary the proportions of four eQTL groups across 3 stocks ###

rm(list=ls());options(stringsAsFactors=FALSE)

setwd("/Users/leileicui/Desktop/1_项目资料/UCL联培研究生/Richard学生")

filenames = c("PeakSNP-Fat_Genes.txt", "PeakSNP-Fat_Trpts.txt", "PeakSNP-Liver_Genes.txt", "PeakSNP-Liver_Trpts.txt")

res_all = c(); thre_all = c()
for(i in 1:4){

	data_tmp = read.table(filenames[i], header=T)

	# data_tmp = subset(data_tmp, data_tmp[,"logP_AD"]>data_tmp[,"thrgenome"])
	# sig_thr = 5.5; data_tmp = subset(data_tmp, data_tmp[,"logP_AD"]>sig_thr)
	data_tmp_used = cbind(data_tmp[,"QTL_Type"], abs(as.numeric(data_tmp[,"OD.A"])))

	### Check the ranges of each eQTL groups ###
	# print(range(subset(data_tmp_used, data_tmp_used[,1]=="additive")[,2]))
	# print(range(subset(data_tmp_used, data_tmp_used[,1]=="partial-dominance")[,2]))
	# print(range(subset(data_tmp_used, data_tmp_used[,1]=="complete-dominance")[,2]))
	# print(range(subset(data_tmp_used, data_tmp_used[,1]=="over-dominance")[,2]))

	res11 = nrow(subset(data_tmp_used, data_tmp_used[,1]=="additive"))
	res12 = nrow(subset(data_tmp_used, data_tmp_used[,1]=="additive"))/nrow(data_tmp_used)

	res21 = nrow(subset(data_tmp_used, data_tmp_used[,1]=="partial-dominance"))
	res22 = nrow(subset(data_tmp_used, data_tmp_used[,1]=="partial-dominance"))/nrow(data_tmp_used)

	res31 = nrow(subset(data_tmp_used, data_tmp_used[,1]=="complete-dominance"))
	res32 = nrow(subset(data_tmp_used, data_tmp_used[,1]=="complete-dominance"))/nrow(data_tmp_used)

	res41 = nrow(subset(data_tmp_used, data_tmp_used[,1]=="over-dominance"))
	res42 = nrow(subset(data_tmp_used, data_tmp_used[,1]=="over-dominance"))/nrow(data_tmp_used)

	res_tmp = cbind(c(res11, res21, res31, res41), c(res12, res22, res32, res42))
	res_all = cbind(res_all, res_tmp)

	# thre_tmp = mean(data_tmp[,"thrgenome"])
	thre_tmp = mean(data_tmp[,"thrsuggest"])
	thre_all = c(thre_all, thre_tmp)

}

res_all
thre_all

################################################################################################################################################

### This script aims to summary the trans-acting enrichment across 3 stocks ###

rm(list=ls());options(stringsAsFactors=FALSE)

setwd("/Users/leileicui/Desktop/1_项目资料/UCL联培研究生/Richard学生")

filenames = c("PeakSNP-Fat_Genes.txt", "PeakSNP-Fat_Trpts.txt", "PeakSNP-Liver_Genes.txt", "PeakSNP-Liver_Trpts.txt")

res_all1 = c(); res_all2 = c(); pvalue_all1 = c(); pvalue_all2 = c(); thre_all = c()
for(i in 1:4){

	data_tmp = read.table(filenames[i], header=T)

	# data_tmp = subset(data_tmp, data_tmp[,"logP_AD"]>data_tmp[,"thrgenome"])
	# sig_thr = 5.5; data_tmp = subset(data_tmp, data_tmp[,"logP_AD"]>sig_thr)
	data_tmp_used = data_tmp[, c("QTL_Type", "chr", "pos", "Chromosome", "start_bp", "end_bp")]

	eQTL_data_rm = data_tmp_used[!is.na(data_tmp_used[,"Chromosome"]),]; rownames(eQTL_data_rm) = 1:nrow(eQTL_data_rm)
	eQTL_data_rm = cbind(eQTL_data_rm, "cis_trans"="trans")

	dis_thr = 2000000
	for(j in 1:nrow(eQTL_data_rm)){
		dis_min = min(abs(eQTL_data_rm[j,"start_bp"] - eQTL_data_rm[j,"pos"]), abs(eQTL_data_rm[j,"end_bp"] - eQTL_data_rm[j,"pos"]))

		if(eQTL_data_rm[j,"Chromosome"] == eQTL_data_rm[j,"chr"] & dis_min < dis_thr){
			eQTL_data_rm[j,"cis_trans"] = "cis"
		}
	}

	print(filenames[i])

	eQTL_data_a_cis = subset(eQTL_data_rm, eQTL_data_rm[,"QTL_Type"]=="additive" & eQTL_data_rm[,"cis_trans"]=="cis")
	eQTL_data_a_trans = subset(eQTL_data_rm, eQTL_data_rm[,"QTL_Type"]=="additive" & eQTL_data_rm[,"cis_trans"]=="trans")

	eQTL_data_pd_cis = subset(eQTL_data_rm, eQTL_data_rm[,"QTL_Type"]=="partial-dominance" & eQTL_data_rm[,"cis_trans"]=="cis")
	eQTL_data_pd_trans = subset(eQTL_data_rm, eQTL_data_rm[,"QTL_Type"]=="partial-dominance" & eQTL_data_rm[,"cis_trans"]=="trans")

	eQTL_data_d_cis = subset(eQTL_data_rm, eQTL_data_rm[,"QTL_Type"]=="complete-dominance" & eQTL_data_rm[,"cis_trans"]=="cis")
	eQTL_data_d_trans = subset(eQTL_data_rm, eQTL_data_rm[,"QTL_Type"]=="complete-dominance" & eQTL_data_rm[,"cis_trans"]=="trans")

	eQTL_data_od_cis = subset(eQTL_data_rm, eQTL_data_rm[,"QTL_Type"]=="over-dominance" & eQTL_data_rm[,"cis_trans"]=="cis")
	eQTL_data_od_trans = subset(eQTL_data_rm, eQTL_data_rm[,"QTL_Type"]=="over-dominance" & eQTL_data_rm[,"cis_trans"]=="trans")

	res_tmp1 = cbind(c(nrow(eQTL_data_a_cis), nrow(eQTL_data_pd_cis), nrow(eQTL_data_d_cis), nrow(eQTL_data_od_cis)), c(nrow(eQTL_data_a_trans), nrow(eQTL_data_pd_trans), nrow(eQTL_data_d_trans), nrow(eQTL_data_od_trans)))
	res_tmp2 = cbind(c(nrow(eQTL_data_a_cis)+nrow(eQTL_data_pd_cis), nrow(eQTL_data_d_cis)+nrow(eQTL_data_od_cis)), c(nrow(eQTL_data_a_trans)+nrow(eQTL_data_pd_trans), nrow(eQTL_data_d_trans)+nrow(eQTL_data_od_trans)))

	pvalue_tmp1 = chisq.test(res_tmp1,simulate.p.value=T)$p.value
	# pvalue_tmp2 = chisq.test(res_tmp2,simulate.p.value=T)$p.value

	# pvalue_tmp1 = fisher.test(res_tmp1,simulate.p.value=F)$p.value
	pvalue_tmp2 = fisher.test(res_tmp2,simulate.p.value=F)$p.value

	res_all1 = cbind(res_all1, res_tmp1); pvalue_all1 = c(pvalue_all1, pvalue_tmp1)
	res_all2 = cbind(res_all2, res_tmp2); pvalue_all2 = c(pvalue_all2, pvalue_tmp2)

	# thre_tmp = mean(data_tmp[,"thrgenome"])
	thre_tmp = mean(data_tmp[,"thrsuggest"])
	thre_all = c(thre_all, thre_tmp)

}

res_all1
pvalue_all1

res_all2
pvalue_all2

thre_all

# p-value < 2.2e-16
# p-value < 


### This script aims to summary the significant eQTLs detected by Add model within cis/trans- eQTLs across 3 stocks ###

rm(list=ls());options(stringsAsFactors=FALSE)

setwd("/Users/leileicui/Desktop/1_项目资料/UCL联培研究生/Richard学生")

filenames = c("PeakSNP-Fat_Genes.txt", "PeakSNP-Fat_Trpts.txt", "PeakSNP-Liver_Genes.txt", "PeakSNP-Liver_Trpts.txt")

res_all = c()
for(i in 1:4){

	data_tmp = read.table(filenames[i], header=T)

	# data_tmp = subset(data_tmp, data_tmp[,"logP_AD"]>data_tmp[,"thrgenome"])
	# sig_thr = 5.5; data_tmp = subset(data_tmp, data_tmp[,"logP_AD"]>sig_thr)
	data_tmp_used = data_tmp[, c("QTL_Type", "chr", "pos", "logP_AD", "logP_A", "thrgenome", "thrsuggest", "Chromosome", "start_bp", "end_bp")]

	eQTL_data_rm = data_tmp_used[!is.na(data_tmp_used[,"Chromosome"]),]; rownames(eQTL_data_rm) = 1:nrow(eQTL_data_rm)
	eQTL_data_rm = cbind(eQTL_data_rm, "cis_trans"="trans")

	dis_thr = 2000000
	for(j in 1:nrow(eQTL_data_rm)){
		dis_min = min(abs(eQTL_data_rm[j,"start_bp"] - eQTL_data_rm[j,"pos"]), abs(eQTL_data_rm[j,"end_bp"] - eQTL_data_rm[j,"pos"]))

		if(eQTL_data_rm[j,"Chromosome"] == eQTL_data_rm[j,"chr"] & dis_min < dis_thr){
			eQTL_data_rm[j,"cis_trans"] = "cis"
		}
	}

	print(filenames[i])

	eQTL_data_cis = subset(eQTL_data_rm, eQTL_data_rm[,"cis_trans"]=="cis")
	eQTL_data_trans = subset(eQTL_data_rm, eQTL_data_rm[,"cis_trans"]=="trans")

	res_tmp1 = c(nrow(subset(eQTL_data_cis, eQTL_data_cis[,"logP_A"] > eQTL_data_cis[,"thrsuggest"])), nrow(eQTL_data_cis))
	res_tmp2 = c(nrow(subset(eQTL_data_trans, eQTL_data_trans[,"logP_A"] > eQTL_data_trans[,"thrsuggest"])), nrow(eQTL_data_trans))

	res_all = cbind(res_all, cbind(res_tmp1, res_tmp2))

}

res_all

