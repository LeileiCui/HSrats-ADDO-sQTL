rm(list=ls());options(stringsAsFactors=FALSE)

Liver_2trpts = read.table("Liver_2trpts.txt", header=T)

Liver_res1_1 = subset(Liver_2trpts, is.na(Liver_2trpts[,"G_logP_AD"]) & !is.na(Liver_2trpts[,"T1_logP_AD"]) & !is.na(Liver_2trpts[,"T2_logP_AD"]) & Liver_2trpts[,"T1_chr"]==Liver_2trpts[,"T2_chr"]) # Find examples of antagonistic T-eQTL, which located in the same chromosome #
Liver_res1_2 = subset(Liver_2trpts, is.na(Liver_2trpts[,"G_logP_AD"]) & !is.na(Liver_2trpts[,"T1_logP_AD"]) & !is.na(Liver_2trpts[,"T2_logP_AD"]) & Liver_2trpts[,"T1_chr"]!=Liver_2trpts[,"T2_chr"]) # Find exmapes of antagonistic T-eQTL, which located in different chromosomes #

Liver_res2_1 = subset(Liver_2trpts, !is.na(Liver_2trpts[,"G_logP_AD"]) & !is.na(Liver_2trpts[,"T1_logP_AD"]) & !is.na(Liver_2trpts[,"T2_logP_AD"]) & Liver_2trpts[,"T1_chr"]==Liver_2trpts[,"T2_chr"] & Liver_2trpts[,"G_chr"]==Liver_2trpts[,"T1_chr"]) # Find examples of synergistic T-eQTL, which located in the same chromosome, and also for the G-eQTL #
Liver_res2_2 = subset(Liver_2trpts, !is.na(Liver_2trpts[,"G_logP_AD"]) & !is.na(Liver_2trpts[,"T1_logP_AD"]) & !is.na(Liver_2trpts[,"T2_logP_AD"]) & Liver_2trpts[,"T1_chr"]==Liver_2trpts[,"T2_chr"] & Liver_2trpts[,"G_chr"]!=Liver_2trpts[,"T1_chr"]) # Find examples of synergistic T-eQTL, which located in the same chromosome, but not for the G-eQTL #
Liver_res2_3 = subset(Liver_2trpts, !is.na(Liver_2trpts[,"G_logP_AD"]) & !is.na(Liver_2trpts[,"T1_logP_AD"]) & !is.na(Liver_2trpts[,"T2_logP_AD"]) & Liver_2trpts[,"T1_chr"]!=Liver_2trpts[,"T2_chr"]) # Find examples of synergistic T-eQTL, which located in different chromosomes, and could be the same as the G-eQTL #

Liver_res = rbind(Liver_res1_1[order(Liver_res1_1[,"T1_logP_AD"],decreasing=T),], 
				colnames(Liver_2trpts), 
				Liver_res1_2[order(Liver_res1_2[,"T1_logP_AD"],decreasing=T),], 
				colnames(Liver_2trpts), 
				Liver_res2_1[order(Liver_res2_1[,"G_logP_AD"],decreasing=T),], 
				colnames(Liver_2trpts), 
				Liver_res2_2[order(Liver_res2_2[,"G_logP_AD"],decreasing=T),], 
				colnames(Liver_2trpts), 
				Liver_res2_3[order(Liver_res2_3[,"G_logP_AD"],decreasing=T),])

write.table(Liver_res, file="Liver_examples.txt", row.names=F, col.name=T, quote=F)

