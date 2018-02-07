
# installera paket,
# ladda BSgenome i ny session
# biocLite: refgenomet


library(refgenomet)
library(mutpatterns)

vcf=read_vcfs_as_granges('Data/BTBdata/Annotation/Strelka_P2233_103T_vs_P2233_122N_somatic_snvs.snpeff.ann.canon.vep.ann.pick.PASS.recode.vcf',genome = 'BSgenome.Hsapiens.UCSC.hg38',sample_names = 'test')

type_occurrences=mut_type_occurrences(vcf,'BSgenome.Hsapiens.UCSC.hg38')

plot_spectrum(type_occurrences, CT = TRUE)
mut_mat <- mut_matrix(vcf_list = vcf, ref_genome = 'BSgenome.Hsapiens.UCSC.hg38')

plot_96_profile(mut_mat,condensed=T)

sp_url <- paste("http://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)

# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]


row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])

#Plot mutational profile of the first two COSMIC signatures:
plot_96_profile(cancer_signatures[,1:2], condensed = TRUE, ymax = 0.3)

#Hierarchically cluster the COSMIC signatures based on their similarity with average linkage:
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
plot(hclust_cosmic)

#The cosine similarity between two mutational profiles/signatures can be calculated with :
cos_sim(mut_mat[,1], cancer_signatures[,1])

#Calculate pairwise cosine similarity between mutational profiles and COSMIC signatures:
cos_sim_samples_signatures = rbind(cos_sim_matrix(mut_mat, cancer_signatures),cos_sim_matrix(mut_mat, cancer_signatures))# made a dummy copy
rownames(cos_sim_samples_signatures)[2]='test2'
plot_cosine_heatmap(cos_sim_samples_signatures,col_order = cosmic_order,cluster_rows = TRUE) # needs >1 samples

#Fit mutation matrix to the COSMIC mutational signatures:
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)

# Plot contribution barplot
select <- which(rowSums(fit_res$contribution) > 10)
plot_contribution(fit_res$contribution[select,],cancer_signatures[,select],coord_flip = FALSE,mode = "absolute") # this will likely work once all data has been collected....




