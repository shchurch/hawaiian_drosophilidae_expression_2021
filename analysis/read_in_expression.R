##### READ IN THE DATA #####
# read in dataset from agalma export_expression
export <- fromJSON("analysis/data/June2021_expression_export.json")

# save species names, codes, IDs
species_names <- c('Drosophila_primaeva','Drosophila_sproati','Drosophila_picticornis','Drosophila_macrothrix','Drosophila_mimica','Drosophila_cfdives','Drosophila_nanella','Scaptomyza_varia','Scaptomyza_cyrtandrae','Scaptomyza_varipicta','Drosophila_atroscutellata','Drosophila_tanythrix')
species_id <- c('008D','106A','025A','055A','040C','16_1','002D','CFB','088B','020A','029A','043D');names(species_id) <- species_names
species_codes <- c('Dpri','Dspr','Dpic','Dmac','Dmim','Dcfd','Dnan','Svar','Scyr','Svpt','Datr','Dtan');names(species_codes) <- species_names
species_order <-  c('Drosophila_primaeva','Drosophila_mimica','Drosophila_nanella','Drosophila_atroscutellata','Drosophila_tanythrix','Drosophila_cfdives','Drosophila_sproati','Drosophila_macrothrix','Drosophila_picticornis','Scaptomyza_cyrtandrae','Scaptomyza_varipicta','Scaptomyza_varia')

# save clade information
PNA <- c("Drosophila_sproati","Drosophila_macrothrix","Drosophila_picticornis")
MM <- c("Drosophila_mimica","Drosophila_nanella")
AMC <- c("Drosophila_tanythrix","Drosophila_atroscutellata")
SCAP <- c("Scaptomyza_varipicta","Scaptomyza_cyrtandrae","Scaptomyza_varia")
SCAP2 <- c("Scaptomyza_varipicta","Scaptomyza_cyrtandrae")
PRIM <- c("Drosophila_primaeva")
HAL <- c("Drosophila_cfdives")

# save tissue information (tissue=treatment)
treatment <- c("n","o","b")
names(treatment) <- c("head","ovary","carcass")

# mark which libraries are resequenced for sequencing replicates replicates
resequenced_ids <- c("18_0_14_resequence1","8_0_2_resequence1","40_2_2_resequence1","40_2_2_resequence2","40_2_4_resequence1","18_0_17_resequence1")

##### READ IN THE TREES ##### 
# get species tree, and make ultrametric
species_tree <- phytools::midpoint.root(read.tree(text=export$speciestree))
species_tree$node.label<-NULL
species_tree <- hutan::slide_root_edges( species_tree )
species_tree <- ape::ladderize( species_tree )

species_ultrametric <- chronos( species_tree, lambda=1, model="correlated", quiet=TRUE )
class( species_ultrametric ) <- "phylo"

# reduce species tree to focal species for expression
tips <- species_ultrametric$tip.label
species_tree_focal <- drop.tip(species_ultrametric,tips[which(!(tips %in% species_names))])

# get gene trees, and give them a unique ID hash
genetrees <- lapply(export$genetrees,parse_gene_tree)
genetree_hashes = lapply(genetrees, digest) %>% unlist()

# decompose orthologs with modified agalmar function
new_decompose_orthologs <- function (nhx) {
    duplications = nhx@data$node[(nhx@data$Ev == "D")]
    phy = nhx@phylo
    to_prune = phy$edge[, 2][phy$edge[, 1] %in% duplications]
    subtrees = hutan::decompose_tree(phy, to_prune)
    return(subtrees)
}

orthologs <- lapply(genetrees,new_decompose_orthologs)
ortholog_hashes <- lapply(orthologs,lapply,digest)

# make a dataframe of gene tree hashes to ortholog hashes
gene_ortho_hash <- lapply(seq(1:length(genetrees)),function(x){
		data.frame(genetree = genetree_hashes[x],
					ortholog = unlist(ortholog_hashes[[x]]),stringsAsFactors=F)
	}) %>% bind_rows()

# make a dataframe of ortholog hashes and tip info
ortholog_info <- lapply(orthologs,lapply,get_tip_info) %>% unlist(recursive=F)

seq_ortho_hash <- lapply(seq(1:length(ortholog_info)),function(x){
	data.frame(species = ortholog_info[[x]]$species,
				seq_id = ortholog_info[[x]]$id,
				ortholog = unlist(ortholog_hashes)[x],stringsAsFactors=F)
	}) %>% bind_rows()

# combine sequence, orthology, and genetree datasets
seq_gene_ortho <- left_join(seq_ortho_hash,gene_ortho_hash,by="ortholog")
seq_gene_ortho$species <- seq_gene_ortho$species %>% sub(" ", "_",.)

##### EXTRACT THE EXPRESSION VALUES #####
# wrapper for getting Expression object from export
get_expression <- function(name) {
	id <- species_id[name] # get species ID
	id_export <- export$expression[[id]] # get expression object
	id_expression <- Expression(id_export) # run agalmar Expression function
	return(id_expression) 
	}

# wrapper for calculating TPM10k, using exported TPM value
get_tpm10k <- function(exp) {
	id_numgenes <- summarize_reference(exp)$Genes # get number of genes in reference
	id_tpm <- exp@tpm # get TPM values 
	id_tpm10k <-  log((( id_tpm * id_numgenes ) / 10^4) + 1) # nat log (tpm * length / 1000) + 1)
	rownames(id_tpm10k) <- row.names(exp@tpm) # set names same as TPM
	return(id_tpm10k)
	}

exps <- lapply(species_names,get_expression)
exps10k <- lapply(exps,get_tpm10k)

library_info <- bind_rows(lapply(exps,function(x){data.frame(species = x@species,
											individual = x@individual,
											treatment = x@treatment,
											library_id = x@library_id)}))

# wrapper for calculating differential gene expression, excluding resequenced libraries
wrap_deseq2 <- function(exp) {
	dds <- create_DESeq2(exp,~individual+treatment) # build dds matrix w/ DESeq2
	dds$treatment <- factor(dds$treatment, levels = c("carcass","head","ovary")) # choose treatments
	dds <- dds[rowSums(counts(dds)) >= threshold_min_count,] # filter with threshold
	new_dds <- dds[,which(!colnames(dds) %in% resequenced_ids)]
	dds_res <- DESeq(new_dds) # run DESEq2
	return(dds_res)
}

threshold_min_count <- 0 # arbitrary threshold
significance_threshold <- 0.01
exps_dds <- lapply(exps,wrap_deseq2)

# count number of significantly differentiall expressed genes
n_sig_up <- sapply(exps_dds,function(x){results(x) %>% as.data.frame %>% filter(pvalue < significance_threshold,log2FoldChange > 0) %>% nrow()})
n_sig_down <- sapply(exps_dds,function(x){results(x) %>% as.data.frame %>% filter(pvalue < significance_threshold,log2FoldChange < 0) %>% nrow()})
n_total <- sapply(exps_dds,function(x){results(x) %>% as.data.frame %>% nrow()})

head_n_sig_up <- sapply(exps_dds,function(x){results(x,name="treatment_head_vs_carcass") %>% as.data.frame %>% filter(pvalue < significance_threshold,log2FoldChange > 0) %>% nrow()})
head_n_sig_down <- sapply(exps_dds,function(x){results(x,name="treatment_head_vs_carcass") %>% as.data.frame %>% filter(pvalue < significance_threshold,log2FoldChange < 0) %>% nrow()})
head_n_total <- sapply(exps_dds,function(x){results(x,name="treatment_head_vs_carcass") %>% as.data.frame %>% nrow()})

# Perform a DGE analysis on all libraries, including resequenced ones
wrap_full_deseq2 <- function(exp) {
	dds <- create_DESeq2(exp,~individual+treatment) # build dds matrix w/ DESeq2
	dds$treatment <- factor(dds$treatment, levels = c("carcass","head","ovary")) # choose treatments
	dds <- dds[rowSums(counts(dds)) >= threshold_min_count,] # filter with threshold
	dds_res <- DESeq(dds) # run DESEq2
	return(dds_res)
}

resequenced_species <- which(species_names %in% gsub(" ","_",(library_info %>% filter(library_id %in% resequenced_ids) %>% pull(species))))

full_exps_dds <- lapply(exps[resequenced_species],wrap_full_deseq2)

# plot principal component analyses within species
plot_pca <- function(sp,exps_dds){
	dds <- exps_dds[[sp]]
	vsd <- vst(dds, blind=FALSE) # run variance stabilizing transformation
	pcaData <- plotPCA(vsd, intgroup=c("treatment", "individual"), returnData=TRUE)
	pca_plot <- ggplot(pcaData,aes(x=PC1, y=PC2, color=treatment, shape=individual)) + 
			geom_point(size = 5, alpha=0.75) + 
			scale_color_manual(values = c("#4C505B","#5179BD","#D05B5B")) + 
			guides(color = guide_legend(override.aes=list(shape=15))) + 
			xlab(paste("PC1 (",round(attr(pcaData,'percentVar')[1] * 100,2),"%)",sep="")) + 
			ylab(paste("PC2 (",round(attr(pcaData,'percentVar')[2] * 100,2),"%)",sep="")) + 
			ggtitle(species_codes[[sp]]) +
			theme(plot.title = element_text(hjust = 0.5)) +
			theme(text = element_text(size=10)) +
			theme(legend.position = "none")

	return(pca_plot)
}

pcas <- lapply(order(match(names(species_id),species_order)),plot_pca,exps_dds=exps_dds)

pdf(file = "figures_and_panels/panel_pca_exemplar.pdf",width=3,height=2,useDingbats=F)
print(plot_pca(1,exps_dds) + theme(legend.position = "right"))
dev.off()

pdf(file = "figures_and_panels/panel_pca_grid.pdf",width=9,height=11,useDingbats=F)
grid.arrange(grobs = pcas,ncol=3)
dev.off()

# plot principal component analyses within species
plot_pca <- function(sp,exps_dds){
	dds <- exps_dds[[sp]]
	vsd <- vst(dds, blind=FALSE) # run variance stabilizing transformation
	pcaData <- plotPCA(vsd, intgroup=c("treatment", "individual"), returnData=TRUE)
	pca_plot <- ggplot(pcaData,aes(x=PC1, y=PC2, fill=treatment, shape=individual)) + 
			geom_point(size = 5, alpha=0.5, color="black") + 
			scale_fill_manual(values = c("#4C505B","#5179BD","#D05B5B")) + 
			scale_shape_manual(values = c(21,22,24)) + 
			guides(color = guide_legend(override.aes=list(shape=15))) + 
			xlab(paste("PC1 (",round(attr(pcaData,'percentVar')[1] * 100,2),"%)",sep="")) + 
			ylab(paste("PC2 (",round(attr(pcaData,'percentVar')[2] * 100,2),"%)",sep="")) + 
			ggtitle(species_codes[[sp]]) +
			theme(plot.title = element_text(hjust = 0.5)) +
			theme(text = element_text(size=10)) +
			theme(legend.position = "none")

	return(pca_plot)
}

res_pcas <- lapply(seq(1:length(resequenced_species)),plot_pca,exps_dds=full_exps_dds)

pdf(file = "figures_and_panels/panel_resequenced_pca_grid.pdf",width=9,height=3,useDingbats=F)
grid.arrange(grobs = res_pcas,ncol=3)
dev.off()

### READING IN HOMOLOGY INFERENCE AND ANNOTATING THE DATA

# commands for getting homology info out of agalma

#sqlite3 agalma/data/backup.aglama.sqlite
#.output run_seqid_header.csv
#.mode csv
#select run_id ,model_id, header from agalma_sequences where type = 'n';

#.output homology_54_table.csv
#.mode csv
#select * from agalma_homology_models where homology_id in (select id from agalma_homology where run_id == 54);

hom_table <- read.delim("analysis/data/homology_54_table.csv",sep=",",header=F,stringsAsFactors=F,col.names=c("homolog","seq_id"))
full_run_seq_header <- read.delim("analysis/data/run_seqid_header.csv",sep=",",header=F,stringsAsFactors=F,col.names=c("run_sp","seq_id","header")) 

dmel_runid <- 31 # specific to agalma run
exp_runid <- c(13,14,15,16,17,18,20,21,22,23,24,25) # specific to agalma run

# get information on agalma run, sequence ID, and trinity header
run_seq_header <- full_run_seq_header %>% filter(run_sp %in% exp_runid) %>% separate(header, c("header",NA), sep = "\ ",extra="drop")

# get Drosophila melanogaster annotations from the header of the Dmel transcriptome
dmel_seqs <- full_run_seq_header %>% 
	filter(run_sp == dmel_runid) %>%
	extract(header,into=c("transcript","name","annotation","parent_gene"),"^(.*?)\ .*name=(.*?)-.*;.*FlyBase_Annotation_IDs:(.*?)[,;].*parent=(.*?);.*$") %>% 
	select(seq_id,name,transcript,annotation,parent_gene) %>% 
	mutate(annotation = paste("Dmel_",gsub("-.*","",annotation),sep=""))

get_dmel_seq_id <- function(gt_num) {
	gt_seqs <- gsub("Drosophila_melanogaster@","",genetrees[[gt_num]]@phylo$tip.label[grepl("melanogaster",genetrees[[gt_num]]@phylo$tip.label)])
	return(data.frame(genetree=rep(genetree_hashes[[gt_num]],length(gt_seqs)),seq_id=as.integer(gt_seqs)))
}

# annotate genetrees with D melanogaster header and name
genetree_dmel_seqs <- pblapply(seq(1:length(genetrees)),get_dmel_seq_id)
genetree_dmel_id <- left_join(bind_rows(genetree_dmel_seqs),dmel_seqs,by="seq_id")

# annotate genetrees with combined headers and names from D melanogaster
genetree_all_dmel_id <- genetree_dmel_id %>% group_by(genetree) %>% arrange(name) %>% summarize(name = paste(name,collapse=";"),annotation = paste(annotation,collapse=";"),parent_gene = paste(parent_gene,collapse=";"))

# use the results of HI species -> Dmel BLAST search to annotate HI transcriptomes
blast_match <- function(sp_id){
	blast_dmel_transcript <- read.delim(paste("analysis/BLAST/BLAST_results/",sp_id,"_dmel-transcript.tsv",sep=""),header=F,stringsAsFactors=F,col.names=c("qseqid","qlen","sseqid","slen","frames","pident","nident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","stitle"))
	blast_dmel_translation <- read.delim(paste("analysis/BLAST/BLAST_results/",sp_id,"_dmel-translation.tsv",sep=""),header=F,stringsAsFactors=F,col.names=c("qseqid","qlen","sseqid","slen","frames","pident","nident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","stitle"))
	blast_dros_transcript  <- read.delim(paste("analysis/BLAST/BLAST_results/",sp_id,"_dros-transcript.tsv",sep=""),header=F,stringsAsFactors=F,col.names=c("qseqid","qlen","sseqid","slen","frames","pident","nident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","stitle"))
	blast_dros_translation <- read.delim(paste("analysis/BLAST/BLAST_results/",sp_id,"_dros-translation.tsv",sep=""),header=F,stringsAsFactors=F,col.names=c("qseqid","qlen","sseqid","slen","frames","pident","nident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","stitle"))

	dmel_blastx <- blast_dmel_translation %>% 
		rename(qseqid = "header",stitle = "title") %>% 
		extract(title,into=c("prot","name","parent_gene","transcript","sp_match"),"^(.*?)\ .*name=(.*?)-.*;.*parent=(.*?),(.*?);.*species=(.*?);.*$") %>% 
		select(header,name,prot,transcript,parent_gene,sp_match) %>% 
		distinct(header,.keep_all=T) 

	dros_blastx <- blast_dros_translation %>% 
		rename(qseqid = "header",stitle = "title") %>% 
		extract(title,into=c("prot","name","parent_gene","transcript","sp_match"),"^(.*?)\ .*name=(.*?)-.*;.*parent=(.*?),(.*?);.*species=(.*?);.*$") %>% 
		select(header,name,prot,transcript,parent_gene,sp_match) %>% 
		distinct(header,.keep_all=T)

	dmel_blastn <- blast_dmel_transcript %>% 
		rename(qseqid = "header",stitle = "title") %>% 
		extract(title,into=c("transcript","name","parent_gene","sp_match"),"^(.*?)\ .*name=(.*?)-.*;.*parent=(.*?);.*species=(.*?);.*$") %>% 
		select(header,name,transcript,parent_gene,sp_match) %>% 
		distinct(header,.keep_all=T)

	dros_blastn <- blast_dros_transcript %>% 
		rename(qseqid = "header",stitle = "title") %>% 
		extract(title,into=c("transcript","name","parent_gene","sp_match"),"^(.*?)\ .*name=(.*?)-.*;.*parent=(.*?);.*species=(.*?);.*$") %>% 
		select(header,name,transcript,parent_gene,sp_match) %>% 
		distinct(header,.keep_all=T)

	blast_list <- list(dmel_blastx,dros_blastx,dmel_blastn,dros_blastn)

	return(blast_list)
}

blast_list <- pblapply(species_id,blast_match)
