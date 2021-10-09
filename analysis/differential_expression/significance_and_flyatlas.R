#source("read_in_expression.R")

### READ IN FLYATLAS ###

flyatlas <- read.delim("analysis/data/120721_FlyAtlas2_sqldump.csv",sep=",",stringsAsFactors=F,header=T) %>% 
				filter(type == "protein coding")
flyatlas_f <- flyatlas %>% filter(stage == "Adult",sex == "Female") %>% # female only values
				select(-stage,-sex)
flyatlas_ref <- flyatlas_f %>% filter(tissue_reference == "Yes") %>% # from FlyAtlas instructions, enrichment is calculated of reference tissue
				mutate(refFPKM = ifelse(FPKM>2,FPKM,2)) %>% 
				select(parent_gene,refFPKM)
flyatlas_en <- flyatlas_f %>% left_join(.,flyatlas_ref,by="parent_gene") %>% # calculate tissue enrichment values
				filter(tissue_reference == "No") %>% 
				mutate(enrichment = FPKM / refFPKM) %>% 
				select(-tissue_reference)
flyatlas_ov <- flyatlas_en %>% filter(tissue == "Ovary") # select only ovary values
flyatlas_max <- flyatlas_en %>% # check which tissue is the maximum
				filter(!(tissue %in% c("Head","Brain","Eye"))) %>% # exclude these values because in our data they were sequenced separetely
				group_by(parent_gene) %>% 
				arrange(desc(enrichment)) %>% 
				filter(row_number()==1)

### ANNOTATE GENES BY PARENT GENE ###

annotate_genes <- function(edds,sp){
	dds <- edds # grab DEseq2 results
	res <- results(dds,name="treatment_ovary_vs_carcass") # select contrast of interest: ovary vs carcass
	res_df <- res %>% 
		as.data.frame %>% 
		tibble::rownames_to_column() %>% 
		mutate(seq_id = as.numeric(rowname)) %>% 
		left_join(.,run_seq_header,by="seq_id") %>% 
		left_join(.,blast_list[[sp]][[1]],by="header") %>% 
		mutate(significant = ifelse(padj > significance_threshold,"NO","YES"),
			direction = ifelse(log2FoldChange > 0,"UP","DOWN"),
			blast = ifelse(is.na(parent_gene),"NO","YES"),
			in_flyatlas = ifelse(parent_gene %in% flyatlas_en$parent_gene,"YES","NO"),
			flyatlas_en = ifelse(parent_gene %in% (flyatlas_ov %>% filter(enrichment>1) %>% pull(parent_gene)),"YES","NO"),
			flyatlas_max = ifelse(parent_gene %in% (flyatlas_max  %>% filter(tissue == "Ovary") %>% pull(parent_gene)),"YES","NO"),
			species = sp,
			species_code = species_codes[species]) %>%
		left_join(.,seq_gene_ortho,by=c("species","seq_id"))
	return(res_df)
}

annotated_genes <- pblapply(seq(1:length(species_names)),function(x){annotate_genes(exps_dds[[x]],species_names[[x]])})

### IDENTIFY CORE GENES BY PARENT GENE ###

core_genes_threshold <- 10

core_parent_genes <- annotated_genes %>% bind_rows() %>% 
	filter(significant == "YES",direction == "UP") %>% 
	filter(!(is.na(parent_gene))) %>% distinct(species,parent_gene) %>% 
	group_by(parent_gene) %>% summarize(n = n()) %>% arrange(desc(n)) %>% 
	filter(n > core_genes_threshold) %>% ungroup() %>% left_join(.,annotated_genes %>% bind_rows(),by="parent_gene") 

annotated_genes <- pblapply(annotated_genes,function(x){x %>% mutate(core_parent_gene = ifelse(parent_gene %in% core_parent_genes$parent_gene,"YES","NO"))}) 

frac_sig_core <- annotated_genes %>% bind_rows() %>% filter(significant=="YES",direction=="UP") %>%
	group_by(species,core_parent_gene) %>% tally() %>% pivot_wider(names_from=core_parent_gene,values_from=n) %>%
	group_by(species) %>% summarize(val=YES/NO)

### PLOT DGE VOLCANOS SHOWING SIGNIFICANCE

# generate a volcano plot comparing the 
plot_volcanos <- function(sp,ann_genes){
	res_df <- ann_genes[[sp]] %>% as.data.frame()
	volcano_plot <- ggplot(res_df, aes (x = log2FoldChange, y = -log10(padj))) +
				geom_blank() + 
				geom_point(data = res_df %>% filter(significant == "NO"),shape=16,size = 0.1,color = "gray") + 
				geom_point(data = res_df %>% filter(significant == "YES"),shape=16,size = 0.1,color = "black") + 
				theme(legend.position="none") +
				scale_x_continuous(limits=c(-30,30)) +
				ggtitle(res_df$species_code[1]) +
				theme(text = element_text(size=6)) +
				ylab("-log10 (adjusted p-value)") +
				xlab("log2 fold change") +
				geom_vline(xintercept=0,linetype="dashed",size=0.1)

	return(volcano_plot)
}

volcanos <- lapply(seq(1:length(annotated_genes)),plot_volcanos,ann_genes=annotated_genes)

	pdf(file = "figures_and_panels/panel_volcano_exemplar.pdf",width=3,height=1.1,useDingbats=F)
	print(plot_volcanos(1,annotated_genes) + theme(legend.position = "right"))
	dev.off()

plot_grid_volcanos <- function(sp,ann_genes){
	res_df <- ann_genes[[sp]] %>% as.data.frame()
	volcano_plot <- ggplot(res_df, aes (x = log2FoldChange, y = -log10(padj))) +
				geom_blank() + 
				geom_point(data = res_df %>% filter(significant == "NO"),shape=16,size = 0.25,color = "gray") + 
				geom_point(data = res_df %>% filter(significant == "YES"),shape=16,size = 0.25,color = "black") + 
				theme(legend.position="none") +
				ggtitle(res_df$species_code[1]) +
				theme(plot.title = element_text(hjust = 0.5)) +
				theme(text = element_text(size=10)) +
				ylab("-log10 (adjusted p-value)") +
				xlab("log2 fold change") +
				geom_vline(xintercept=0,linetype="dashed",size=0.1)

	return(volcano_plot)
}

grid_volcanos <- lapply(seq(1:length(annotated_genes)),plot_grid_volcanos,ann_genes=annotated_genes)

	pdf(file = "figures_and_panels/panel_volcano_grid.pdf",width=9,height=11,useDingbats=F)
	grid.arrange(grobs = grid_volcanos[match(species_order,species_names)],ncol=3)
	dev.off()

### PLOT REFERENCE TREE

	pdf(file="figures_and_panels/panel_species_ultrametric.pdf",width=6,height=5,useDingbats=F)
	panel_species_ultrametric <- species_ultrametric
	panel_species_ultrametric$tip.label <- gsub("Drosophila_","D. ",panel_species_ultrametric$tip.label)
	panel_species_ultrametric$tip.label <- gsub("Scaptomyza_","S. ",panel_species_ultrametric$tip.label)
	plot(panel_species_ultrametric)
	dev.off()

### PLOT JITTER PLOT SHOWING SIGNIFICANCE ###

r <- annotated_genes %>% 
	bind_rows()
significance_jitter <- ggplot(r,
	aes(x=log2FoldChange,y=factor(species_code,levels=rev(species_codes[species_order])))) + 
	geom_blank() + 
	geom_jitter(data = r %>% filter(significant == "NO"),shape=16,size = 0.1,color = "gray") + 
	geom_jitter(data = r %>% filter(significant == "YES"),shape=16,size = 0.1,color = "black") + 
	scale_x_continuous(limits=c(-30,30)) +
	theme(legend.position="none") +
	theme(text = element_text(size=6)) +
	theme(axis.title.y=element_blank()) +
	xlab("log2 fold change") +
	geom_vline(xintercept=0,linetype="dashed",size=0.1)


	pdf(file = "figures_and_panels/panel_significance_jitter.pdf",width=3,height=2,useDingbats=F)
	print(significance_jitter)
	dev.off()

	
### PLOT JITTER PLOT WITH ATLAS OVERLAY ###

atlas_jitter <- ggplot(r,
	aes(x=log2FoldChange,y=factor(species_code,levels=rev(species_codes[species_order])))) + 
	geom_blank() + 
	geom_jitter(data = r %>% filter(significant == "NO"),shape=16,size = 0.1,color = "gray") + 
	geom_jitter(data = r %>% filter(significant == "YES",blast == "NO"),shape=16,size = 0.1,color = "black") + 
	geom_jitter(data = r %>% filter(significant == "YES",blast == "YES",flyatlas_max == "NO"),shape=16,size = 0.1,color = "dark cyan") + 
	geom_jitter(data = r %>% filter(significant == "YES",blast == "YES",flyatlas_max == "YES"),shape=16,size = 0.1,color = "dark orange") + 
	theme(legend.position="none") +
	theme(text = element_text(size=6)) +
	theme(axis.title.y=element_blank()) +
	scale_x_continuous(limits=c(-30,30)) +
	geom_vline(xintercept=0,linetype="dashed",size=0.1)
	xlab("log2 fold change")

	pdf(file = "figures_and_panels/panel_atlas_jitter.pdf",width=3,height=2,useDingbats=F)
	print(atlas_jitter)
	dev.off()

### PLOT JITTER PLOT WITH CORE PARENT GENE OVERLAY

core_parent_gene_jitter <- ggplot(r,
	aes(x=log2FoldChange,y=factor(species_code,levels=rev(species_codes[species_order])))) + 
	geom_blank() + 
	geom_jitter(data = r %>% filter(core_parent_gene == "NO"),shape=16,size = 0.1,color = "gray") + 
	geom_jitter(data = r %>% filter(core_parent_gene == "YES"),shape=16,size = 0.1,color = "#C12756") + 
	theme(legend.position="none") +
	theme(text = element_text(size=6)) +
	theme(axis.title.y=element_blank()) +
	scale_x_continuous(limits=c(-30,30)) +
	geom_vline(xintercept=0,linetype="dashed",size=0.1)
	xlab("log2 fold change")

	pdf(file = "figures_and_panels/panel_atlas_core_jitter.pdf",width=3,height=2,useDingbats=F)
	print(core_parent_gene_jitter)
	dev.off()


### PLOT SCATTER PLOT WITH ENRICHMENT VS MEAN LOG FOLD CHANGE ###

mean_fold_change_enrichment <- r %>% 
 	group_by(parent_gene,species) %>% 
 	summarize(meanl2fc = mean(log2FoldChange)) %>% 
 	group_by(parent_gene)  %>% 
 	summarize(mean_ml2fc = mean(meanl2fc)) %>% 
 	left_join(.,flyatlas_ov,by="parent_gene") %>% 
 	left_join(.,r %>% distinct(parent_gene,core_parent_gene),by="parent_gene")

mean_fold_change_enrichment_scatter <- ggplot(mean_fold_change_enrichment,aes(x = mean_ml2fc,y= enrichment)) + 
	geom_blank() + 
	geom_point(data=mean_fold_change_enrichment%>%filter(core_parent_gene=="NO"),shape=16,size = 0.25,color="gray") +
	geom_point(data=mean_fold_change_enrichment%>%filter(core_parent_gene=="YES"),shape=16,size = 0.25,color="#C12756") +
	scale_color_manual(values = c("gray","red")) +
	theme(legend.position="none") +
	theme(text = element_text(size=6)) +
	geom_vline(xintercept=0,linetype="dashed",size=0.1) +
	ylab("Dmel ovary expression enrichment") +
	xlab("mean log2 fold change across Hawaiian species and transcripts")

	pdf(file = "figures_and_panels/panel_mean_fold_change_enrichment.pdf",width=3,height=2,useDingbats=F)
	print(mean_fold_change_enrichment_scatter)
	dev.off()

### PLOT LABELED CORE PARENT GENES ###

core_parent_gene_names_scatter <- ggplot(mean_fold_change_enrichment %>% filter(core_parent_gene == "YES"),aes(x = mean_ml2fc,y= enrichment,label = name)) +
	geom_point(color="#C12756",size=0.1) +
	geom_text_repel(max.overlaps=100,size=1,segment.size=0.1,box.padding=0.05) + 
	theme(legend.position="none") +
	theme(text = element_text(size=6)) +
	ylab("Dmel ovary expression enrichment") +
	xlab("mean log2 fold change across Hawaiian species and transcripts")

	pdf(file = "figures_and_panels/panel_core_parent_gene_names.pdf",width=3,height=3,useDingbats=F)
	print(core_parent_gene_names_scatter)
	dev.off()


