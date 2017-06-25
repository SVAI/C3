library(data.table)
library(deconstructSigs)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ggplot2)
library(topGO)
library(org.Hs.eg.db)
library(gridExtra)

setwd("/media/avantika/My Passport/nf2hack/")

#Read point mutation data from TCGA
files_GBM = grep(".txt", list.files("/media/avantika/My Passport/nf2hack/gdac.broadinstitute.org_GBM.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/", full.names=T), value=T)
files_LGG = grep(".txt", list.files("/media/avantika/My Passport/nf2hack/gdac.broadinstitute.org_LGG.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/", full.names=T), value=T)

files = c(files_GBM, files_LGG)
cols = c("Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Start_position", "End_position", "Strand",
         "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", 
         "GO_Biological_Process",	"GO_Molecular_Function")

#Select relevant fields and concatenate patient files
for(i in 1:length(files)){
  f = files[i]
  patient = strsplit(strsplit(f, split = "TCGA-")[[1]][2], split = ".hg19")[[1]][1]
  type = strsplit(strsplit(f, split = "broadinstitute.org_")[[1]][2], split = ".Mutation_Packager_")[[1]][1]
  patient = paste0(type, "_", patient)
  x = fread(f, header = T)
  x = subset(x, select = cols)
  x = cbind(patient, x)
  write.table(x, file = "muts.txt", append = T, quote = F, col.names = F, row.names = F, sep = "\t")
}

muts = fread("muts.txt")
setnames(muts, c("patient", "Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Start_position", "End_position", "Strand",
                 "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", 
                 "GO_Biological_Process",	"GO_Molecular_Function"))
         
#Generate feature matrix for genes
gene_patient = dcast(unique(muts[,list(patient, Hugo_Symbol)]), Hugo_Symbol~patient, fun.aggregate = length, fill = 0)
write.table(gene_patient, file = "gene_patient.txt", quote = F, sep = "\t")

#Generate feature matrix for GO:biological annotation
go = unique(muts[,list(patient, Hugo_Symbol, Entrez_Gene_Id, GO_Biological_Process,	GO_Molecular_Function)])
biol = cbind(go[, list(patient, Hugo_Symbol, Entrez_Gene_Id)], as.data.table(tstrsplit(go$GO_Biological_Process, split = "\\|")))
biol = na.omit(melt(biol, id.vars = c(1,2,3), value.name = "go"))[,variable := NULL]
biol_patient = dcast(biol[,.N, by = list(go, patient)], go~patient, value.var = "N", fill=0)
write.table(biol_patient, file = "biol_patient.txt", quote=F, sep = "\t")

#Generate feature matrix for GO:molecular annotation
molec = cbind(go[, list(patient, Hugo_Symbol, Entrez_Gene_Id)], as.data.table(tstrsplit(go$GO_Molecular_Function, split = "\\|")))
molec = na.omit(melt(molec, id.vars = c(1,2,3), value.name = "go"))[,variable := NULL]
molec_patient = dcast(molec[,.N, by = list(go, patient)], go~patient, value.var = "N", fill=0)
write.table(molec_patient, file = "molec_patient.txt", quote=F, sep = "\t")

#Generate feature matrix for GO:biological enrichment
biol_go_ngenes = biol[, length(unique(Entrez_Gene_Id)), by=go]
biol_go = dcast(biol[, .N, by = list(patient, go)], go~patient, value.var = "N", fill = 0)
biol_go = merge(biol_go_ngenes, biol_go, by = "go")
setnames(biol_go, "V1", "n_genes_in_go")
biol_go = biol_go[n_genes_in_go>10]

ps = colnames(biol_go)[3:ncol(biol_go)]
n_genes_mut_p = colSums(biol_go[,3:ncol(biol_go)])
n_genes = length(unique(biol$Entrez_Gene_Id))

go_enr = data.table()

for(p in ps){
  n_genes_mut = n_genes_mut_p[p]
  n_genes_not_mut = n_genes - n_genes_mut
  
  x = biol_go[,list(go,n_genes_in_go,get(p))][V3>1]
  
  if(nrow(x)>0){
    fish = c()
    for(j in 1:nrow(x)){
      fish = c(fish, fisher.test(matrix(c(x$V3[j], x$n_genes_in_go[j]-x$V3[j], n_genes_mut, n_genes_not_mut), nrow=2))$p.value)
    }
    x$fish=fish
    x = x[fish<0.05]
    x$patient=p
  }
  
  if(p == ps[1]){
    go_enr = x
  }
  else if(nrow(x)>0){
    go_enr = rbind(go_enr,x)
  }
}

biol_enr_patient = dcast(go_enr, go~patient, value.var = "V3", fill=0)
write.table(biol_enr_patient, file = "biol_enr_patient.txt", quote=F, sep = "\t")

#Feature matrix for molecular GO enrichment
molec_go_ngenes = molec[, length(unique(Entrez_Gene_Id)), by=go]
molec_go = dcast(molec[, .N, by = list(patient, go)], go~patient, value.var = "N", fill = 0)
molec_go = merge(molec_go_ngenes, molec_go, by = "go")
setnames(molec_go, "V1", "n_genes_in_go")
molec_go = molec_go[n_genes_in_go>10]

ps = colnames(molec_go)[3:ncol(molec_go)]
n_genes_mut_p = colSums(molec_go[,3:ncol(molec_go)])
n_genes = length(unique(molec$Entrez_Gene_Id))

go_enr = data.table()

for(p in ps){
  n_genes_mut = n_genes_mut_p[p]
  n_genes_not_mut = n_genes - n_genes_mut
  
  x = molec_go[,list(go,n_genes_in_go,get(p))][V3>1]
  
  if(nrow(x)>0){
    fish = c()
    for(j in 1:nrow(x)){
      fish = c(fish, fisher.test(matrix(c(x$V3[j], x$n_genes_in_go[j]-x$V3[j], n_genes_mut, n_genes_not_mut), nrow=2))$p.value)
    }
    x$fish=fish
    x = x[fish<0.05]
    x$patient=p
  }
  
  if(p == ps[1]){
    go_enr = x
  }
  else if(nrow(x)>0){
    go_enr = rbind(go_enr,x)
  }
}

molec_enr_patient = dcast(go_enr, go~patient, value.var = "V3", fill=0)
write.table(molec_enr_patient, file = "molec_enr_patient.txt", quote=F, sep = "\t")

#Generate feature matrix for sequence based clustering
mut = muts[Variant_Type=="SNP", list(patient, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2)]
setnames(mut, c("sample", "chr", "pos", "ref", "alt"))
sigs_input = mut.to.sigs.input(mut.ref = mut, 
                               sample.id = "sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt", 
                               bsg = BSgenome.Hsapiens.1000genomes.hs37d5)
muttype_patient = data.table(t(sigs_input), keep.rownames = T)
setnames(muttype_patient, "rn", "Variant_Type")
svtype_patient = as.data.table(dcast(muts[Variant_Type %in% c("DEL", "INS"), .N, by = list(patient, Variant_Type)], Variant_Type~patient, value.var = "N", fill=0))
muttype_patient = rbind(muttype_patient, svtype_patient, fill = TRUE)

muttype_patient = melt(muttype_patient, id.vars = 1, variable.name = "patient", value.name = "N")
muttype_patient[is.na(N), N:=0]
muttype_patient[, N.patient := sum(N), by = patient]
muttype_patient[, freq := N/N.patient]

muttype_patient = dcast(muttype_patient, Variant_Type~patient, value.var = "freq")
write.table(muttype_patient, file = "muttype_patient.txt", quote=F, sep="\t")


#Find mutational signatures in TCGA patients
mut2 = as.data.frame(muts[Variant_Type=="SNP", list(patient, Chromosome, Start_position, Reference_Allele, Tumor_Seq_Allele2)])
colnames(mut2) = c("sample", "chr", "pos", "ref", "alt")
sigs_input = mut.to.sigs.input(mut.ref = mut2, 
                               sample.id = "sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt", 
                               bsg = BSgenome.Hsapiens.1000genomes.hs37d5)
sigs = NULL
for(i in 1:nrow(sigs_input)){
  if(i==1){
    sigs = whichSignatures(sigs_input, rownames(sigs_input)[i], signatures.ref = signatures.cosmic, contexts.needed = T)$weights
  }
  else{
    sigs = rbind(sigs, whichSignatures(sigs_input, rownames(sigs_input)[i], signatures.ref = signatures.cosmic, contexts.needed = T)$weights)
    
  }
}

#Generate feature matrix for mutational signatures
sigs_patients = t(sigs)
write.table(sigs_patients, file = "sigs_patients.txt", quote = F, sep = "\t")

#Filter coding mutations for Onno

load("exons.RData")
exons = unique(exons[, list(chrom, start, end)])
mut_o = as.data.table(read.csv("VCF_filtered_GQ_thresh.csv"))
mut_o = mut_o[, c("CHROM", "POS", "REF", "ALT")]
mut_o = mut_o[nchar(as.character(REF))==1 & nchar(as.character(ALT))==1]

mut_o_filt = mut_o[2,]
for(i in 3:nrow(mut_o)){
  if(i%%100==0){print(i)}
  if(nrow(exons[chrom==mut_o[i, CHROM] & start <= mut_o[i, POS] & end >= mut_o[i, POS]])>0){
    mut_o_filt = rbind(mut_o_filt, mut_o[i,])
  }
}

mut_o_filt[, sample := "onno"]
mut_o_filt = as.data.frame(mut_o_filt)

write.table(mut_o_filt, file = "mut_o_filt.txt", quote=F, sep="\t", row.names = F)

#Find signatures for Onno
sigs_input = mut.to.sigs.input(mut.ref = mut_o_filt, 
                               sample.id = "sample", chr = "CHROM", pos = "POS", ref = "REF", alt = "ALT")
sigs = whichSignatures(sigs_input, "onno", signatures.ref = signatures.cosmic, contexts.needed = T)$weights
sigs_onno = t(sigs)
write.table(sigs_onno, file = "sigs_onno.txt", quote = F, sep = "\t")

#Find mutated genes for onno
load("genes.RData")
genes = unique(genes[, list(name2, chrom, cdsStart, cdsEnd)])
mut_o_filt = as.data.table(mut_o_filt)

onno_genes = c()
for(i in 1:nrow(mut_o_filt)){
  x = genes[chrom==mut_o_filt[i, CHROM] & cdsStart <= mut_o_filt[i, POS] & cdsEnd >= mut_o_filt[i, POS]]
  if(nrow(x)>0){
    onno_genes = c(onno_genes, as.character(x$name2))
  }
}
onno_genes = unique(onno_genes)

onno_biol = biol[Hugo_Symbol %in% onno_genes, .N, by = go]
onno_biol = merge(biol_patient, onno_biol, by = "go", all.x=T)[, list(go, N)]
onno_biol[is.na(N), N:=0]
write.table(onno_biol, file = "onno_biol.txt", quote=F, col.names = F, row.names = F, sep="\t")


biol_onno = t(as.data.table(onno_biol))
write.table(biol_onno, file = "biol_onno.txt", quote=F, col.names = F, row.names = F, sep="\t")
