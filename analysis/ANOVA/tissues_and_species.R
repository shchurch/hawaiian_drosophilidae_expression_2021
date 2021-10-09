#source("read_in_expression.R")

##### BUILD FULL EXPRESSION COMPARISON DATA FRAME #####
build_expression10k_dataframe <- function(n,exps10k,exps) {
	dat <- exps10k[[n]] # grab one tpm10k set at a time
	info <- data.frame(species = sub(" ", "_",exps[[n]]@species),
		treatment = as.character(exps[[n]]@treatment),
		library_id = as.character(exps[[n]]@library_id),
		individual = as.character(exps[[n]]@individual),stringsAsFactors=F) # build data frame of sample information
	info_dat <- dat %>% as.data.frame() %>% 
		mutate(seq_id = as.numeric(rownames(dat))) %>% 
		gather(library_id,expression,info$library_id) %>% 
		left_join(.,info,by="library_id") # combine expression and species info
	return(info_dat)
}

# build expression tables for each species and combine them
full_expression <- bind_rows(lapply(seq(1:length(species_names)),build_expression10k_dataframe,exps10k=exps10k,exps=exps))
expression_genetree <- full_join(seq_gene_ortho,full_expression %>% select(-species),by="seq_id") 

# calculate the average epxression across homology group and biological replicate
ave_ave_expression <- expression_genetree %>% 
	filter(!(library_id %in% resequenced_ids)) %>% # choose unique sequence libraries
	mutate(id = paste(species,treatment,sep="_"),key=genetree) %>% # choose ID variable
	filter(!is.na(expression),!is.na(genetree)) %>% # remove missing data
	group_by(key,library_id,id,species,treatment) %>% summarize(mean_val = mean(expression)) %>% ungroup() %>%
	group_by(key,id,species,treatment) %>% summarize(val = mean(mean_val)) %>% ungroup()

# build a plot that compares PCA, expects 'key','id', and 'val'
build_pca_plot <- function(dat) {
	pca_df <- dat %>% mutate(species = species_codes[species]) %>% 
		tidyr::spread(.,key,val) %>%
		select_if(~ !any(is.na(.))) 
	pca_plot <- autoplot(prcomp(pca_df %>% select(-id,-species,-treatment)),data=pca_df,colour="species",shape="treatment",size=2.5) + 
		scale_color_brewer(palette="Paired") +
		theme(text = element_text(size=7)) +
		theme(legend.position = "none") 
		#coord_fixed()
	return(pca_plot)
}

homolog_pca <- build_pca_plot(ave_ave_expression)

	pdf(file = "figures_and_panels/panel_homolog_pca.pdf",width=3,height=3,useDingbats=F)
	print(homolog_pca)
	dev.off()

panel_species_tree_focal <- species_tree_focal
panel_species_tree_focal$tip.label <- species_codes[panel_species_tree_focal$tip.label]

	pdf(file="figures_and_panels/panel_species_focal.pdf",width=3,height=5,useDingbats=F)
	plot(panel_species_tree_focal)
	dev.off()

### COUNT REPRESENTATION

count_homology <- left_join(expression_genetree %>% distinct(genetree) %>% filter(!is.na(genetree)),
	(expression_genetree %>% 
		select(species,genetree,expression) %>% 
		na.omit %>% 
		distinct(species,genetree)  %>% 
	group_by(genetree) %>% tally()),by="genetree") %>%
	mutate(n = replace_na(n,0))

count_orthology <- left_join(expression_genetree %>% distinct(ortholog) %>% filter(!is.na(ortholog)),
	(expression_genetree %>% 
		select(species,ortholog,expression) %>% 
		na.omit %>% 
		distinct(species,ortholog)  %>% 
	group_by(ortholog) %>% tally()),by="ortholog") %>%
	mutate(n = replace_na(n,0))

	pdf(file="figures_and_panels/panel_homolog_representation.pdf",width=5,height=1.5,useDingbats=F)
	print(ggplot(count_homology,aes(x=n)) + 
		geom_bar(fill="black") + 
		scale_x_continuous(breaks = c(0,seq(1:12))) + 
		xlab("species with expression data per homology group") + 
		ylab("count"))
	dev.off()

	pdf(file="figures_and_panels/panel_ortholog_representation.pdf",width=5,height=1.5,useDingbats=F)
	print(ggplot(count_orthology,aes(x=n)) + 
		geom_bar(fill="black") + 
		scale_x_continuous(breaks = c(0,seq(1:12))) + 
		xlab("species with expression data per orthology group") + 
		ylab("count"))
	dev.off()

### BUILD ANOVA DATASETS FOR LINEAR MODELING

anova_ov_table <- expression_genetree %>% 
	filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
	filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
	mutate(key = paste("sample_",library_id,sep=""),id=genetree) %>% # use lib ID as future column name, genetree as row name
	filter(!is.na(expression)) %>% # filter out seqIDs that didn't have expression measured
	filter(!is.na(genetree)) %>% # filter out seqIDs that didn't match a genetree
	group_by(key,id) %>% # unique lib ID and genetree
	summarize(val = mean(expression)) %>% # take the average of all transcripts in that group
	ungroup() %>%
	tidyr::spread(.,key,val) # pivot to long form data

### FULL ANOVA, WITH MISSING VALUES WHICH WILL BE REPLACED WITH 0

anova_ov_tab <- anova_ov_table %>% select(-id) %>% as.data.frame() 
rownames(anova_ov_tab) <- anova_ov_table %>% pull(id)

	write.table(anova_ov_tab,file="analysis/ANOVA/ANOVA_results/anova_ovary_table.tsv",row.names=T,col.names=T,sep="\t")

### REDUCED ANOVA, NO MISSING VALUES

anova_ov_red_tab <- na.omit(anova_ov_tab)
	write.table(anova_ov_red_tab,file="analysis/ANOVA/ANOVA_results/anova_ovary_reduced_table.tsv",row.names=T,col.names=T,sep="\t")

### PRINT ANOVA METADATA 

anova_metadata <- library_info %>% filter(!library_id %in% resequenced_ids) %>%
	mutate(id = paste("sample_",library_id,sep=""),
	species = species_codes[gsub(" ","_",species)]) %>% 
	select(species,treatment,id) 
	write.table(anova_metadata,file="analysis/ANOVA/ANOVA_results/anova_metadata.tsv",sep="\t",row.names=F,col.names=T,quote=F)

### RUN FULL ANOVA

anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ovary_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ovary_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species -r TRUE"
	system(anova_command)
anova_ov_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ovary_table_results.tsv",header=T,stringsAsFactors=F)

### RUN REDUCED ANOVA

anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ovary_reduced_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ovary_reduced_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species"
	system(anova_command)
anova_ov_red_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ovary_reduced_table_results.tsv",header=T,stringsAsFactors=F)

### TMP
### BUILD AND RUN REDUCED ANOVA ON ALL SPECIES -  SVAR
anova_ov_red2_table <- expression_genetree %>% 
	filter(species != "Scaptomyza_varia") %>% 
	filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
	filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
	mutate(key = paste("sample_",library_id,sep=""),id=genetree) %>% # use lib ID as future column name, genetree as row name
	filter(!is.na(expression)) %>% # filter out seqIDs that didn't have expression measured
	filter(!is.na(genetree)) %>% # filter out seqIDs that didn't match a genetree
	group_by(key,id) %>% # unique lib ID and genetree
	summarize(val = mean(expression)) %>% # add up all transcripts in that group
	ungroup() %>%
	tidyr::spread(.,key,val) # pivot to long form data

anova_ov_red2_tab <- anova_ov_red2_table %>% select(-id) %>% as.data.frame() 
rownames(anova_ov_red2_tab) <- anova_ov_red2_table %>% pull(id)
anova_ov_red2_tab <- na.omit(anova_ov_red2_tab)
	write.table(anova_ov_red2_tab,file="analysis/ANOVA/ANOVA_results/anova_ovary_reduced2_table.tsv",row.names=T,col.names=T,sep="\t")

anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ovary_reduced2_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ovary_reduced2_table.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species"
	system(anova_command)
anova_ov_red2_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ovary_reduced2_table.tsv",header=T,stringsAsFactors=F)

### BUILD AND RUN REDUCED ANOVA ON HI DROSOPHILA CLADE
anova_ov_HIDROS_table <- expression_genetree %>% 
	filter(species %in% c(PNA,PRIM,HAL,AMC,MM)) %>% 
	filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
	filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
	mutate(key = paste("sample_",library_id,sep=""),id=genetree) %>% # use lib ID as future column name, genetree as row name
	filter(!is.na(expression)) %>% # filter out seqIDs that didn't have expression measured
	filter(!is.na(genetree)) %>% # filter out seqIDs that didn't match a genetree
	group_by(key,id) %>% # unique lib ID and genetree
	summarize(val = mean(expression)) %>% # add up all transcripts in that group
	ungroup() %>%
	tidyr::spread(.,key,val) # pivot to long form data

anova_ov_HIDROS_tab <- anova_ov_HIDROS_table %>% select(-id) %>% as.data.frame() 
rownames(anova_ov_HIDROS_tab) <- anova_ov_HIDROS_table %>% pull(id)
anova_ov_HIDROS_red_tab <- na.omit(anova_ov_HIDROS_tab)
	write.table(anova_ov_HIDROS_red_tab,file="analysis/ANOVA/ANOVA_results/anova_ovary_HIDrosophila_reduced_table.tsv",row.names=T,col.names=T,sep="\t")

anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ovary_HIDrosophila_reduced_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ovary_HIDrosophila_reduced_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species"
	system(anova_command)
anova_ov_HIDROS_red_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ovary_HIDrosophila_reduced_table_results.tsv",header=T,stringsAsFactors=F)

### BUILD AND RUN REDUCED ANOVA ON PNA CLADE

anova_ov_PNA_table <- expression_genetree %>% 
	filter(species %in% PNA) %>% 
	filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
	filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
	mutate(key = paste("sample_",library_id,sep=""),id=genetree) %>% # use lib ID as future column name, genetree as row name
	filter(!is.na(expression)) %>% # filter out seqIDs that didn't have expression measured
	filter(!is.na(genetree)) %>% # filter out seqIDs that didn't match a genetree
	group_by(key,id) %>% # unique lib ID and genetree
	summarize(val = mean(expression)) %>% # add up all transcripts in that group
	ungroup() %>%
	tidyr::spread(.,key,val) # pivot to long form data

anova_ov_PNA_tab <- anova_ov_PNA_table %>% select(-id) %>% as.data.frame() 
rownames(anova_ov_PNA_tab) <- anova_ov_PNA_table %>% pull(id)
anova_ov_PNA_red_tab <- na.omit(anova_ov_PNA_tab)
	write.table(anova_ov_PNA_red_tab,file="analysis/ANOVA/ANOVA_results/anova_ovary_PNA_reduced_table.tsv",row.names=T,col.names=T,sep="\t")

anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ovary_PNA_reduced_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ovary_PNA_reduced_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species"
	system(anova_command)
anova_ov_PNA_red_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ovary_PNA_reduced_table_results.tsv",header=T,stringsAsFactors=F)

### BUILD AND RUN REDUCED ANOVA ON PW CLADE

anova_ov_PW_table <- expression_genetree %>% 
	filter(species %in% c("Drosophila_sproati","Drosophila_macrothrix")) %>% 
	filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
	filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
	mutate(key = paste("sample_",library_id,sep=""),id=genetree) %>% # use lib ID as future column name, genetree as row name
	filter(!is.na(expression)) %>% # filter out seqIDs that didn't have expression measured
	filter(!is.na(genetree)) %>% # filter out seqIDs that didn't match a genetree
	group_by(key,id) %>% # unique lib ID and genetree
	summarize(val = mean(expression)) %>% # add up all transcripts in that group
	ungroup() %>%
	tidyr::spread(.,key,val) # pivot to long form data

anova_ov_PW_tab <- anova_ov_PW_table %>% select(-id) %>% as.data.frame() 
rownames(anova_ov_PW_tab) <- anova_ov_PW_table %>% pull(id)
anova_ov_PW_red_tab <- na.omit(anova_ov_PW_tab)
	write.table(anova_ov_PW_red_tab,file="analysis/ANOVA/ANOVA_results/anova_ovary_PW_reduced_table.tsv",row.names=T,col.names=T,sep="\t")

anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ovary_PW_reduced_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ovary_PW_reduced_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species"
	system(anova_command)
anova_ov_PW_red_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ovary_PW_reduced_table_results.tsv",header=T,stringsAsFactors=F)


### SUMMARIZE AND PLOT ANOVA RESULTS
summarize_anova_results <- function(anova_results){
	anova_res_sum <- anova_results %>% 
	mutate(genetree = rownames(anova_results)) %>% 
	rowwise() %>% 
	mutate(total_SS = sum(treatment_SS,species_SS,Residuals_SS),
		prop_species = species_SS/total_SS,
		prop_treatment = treatment_SS /total_SS,
		sum_prop = sum(prop_treatment,prop_species),
		VG = ifelse(sum_prop > 0.75,
				ifelse(prop_species > (2*prop_treatment),"SVG",
					ifelse(prop_treatment > (2*prop_species),"TVG",
					"NO")),"NO")) %>% as.data.frame()
	return(anova_res_sum)
}

anova_phyl_dist <- lapply(list(anova_ov_red_results,anova_ov_HIDROS_red_results,anova_ov_PNA_red_results,anova_ov_PW_red_results),summarize_anova_results)

mean_prop_treatment <- sapply(anova_phyl_dist,function(x){x %>% pull(prop_treatment) %>% mean()})
mean_prop_species <- sapply(anova_phyl_dist,function(x){x %>% pull(prop_species) %>% mean()})
n_homologs <- sapply(anova_phyl_dist,nrow)
all_node_depths <- node.depth.edgelength(species_tree_focal)
scaled_all_node_depths <- (1-all_node_depths * (1 / max(all_node_depths)))
node_depths <- round(scaled_all_node_depths[c(getMRCA(species_tree_focal,c(PNA,SCAP)),
				getMRCA(species_tree_focal,c(PNA,AMC)),
				getMRCA(species_tree_focal,c(PNA)),
				getMRCA(species_tree_focal,c("Drosophila_sproati","Drosophila_macrothrix")))],2)

anova_phyl_dist <- lapply(seq(1:length(node_depths)),function(x){anova_phyl_dist[[x]] %>% mutate(node_depth = node_depths[x])})

prop_phyl_dist <- ggplot(data = data.frame(prop = mean_prop_treatment,depth = node_depths),aes(y = prop, x = depth)) + 
	geom_point(color="dark orange") + 
	geom_line(color="dark orange") + 
	geom_point(data = data.frame(prop = mean_prop_species,depth = node_depths),color="dark cyan") +
	geom_line(data = data.frame(prop = mean_prop_species,depth = node_depths),color="dark cyan") + 
	theme(legend.position = "none") + 
	theme(text = element_text(size=6)) +
	xlab("scaled evolutionary distance") + 
	ylab("mean proportion\nof variance") +
	scale_x_continuous(limits=c(0,1.04))

	pdf(file = "figures_and_panels/panel_prop_variation_by_phylogenetic_distance.pdf",width=3,height=1.75,useDingbats=F)
	print(prop_phyl_dist)
	dev.off()

anova_phyl_dist2 <- lapply(list(anova_ov_red2_results,anova_ov_HIDROS_red_results,anova_ov_PNA_red_results,anova_ov_PW_red_results),summarize_anova_results)

mean_prop_treatment2 <- sapply(anova_phyl_dist2,function(x){x %>% pull(prop_treatment) %>% mean()})
mean_prop_species2 <- sapply(anova_phyl_dist2,function(x){x %>% pull(prop_species) %>% mean()})
prop_phyl_dist2 <- ggplot(data = data.frame(prop = mean_prop_treatment2,depth = node_depths),aes(y = prop, x = depth)) + 
	geom_point(color="dark orange") + 
	geom_line(color="dark orange") + 
	geom_point(data = data.frame(prop = mean_prop_species2,depth = node_depths),color="dark cyan") +
	geom_line(data = data.frame(prop = mean_prop_species2,depth = node_depths),color="dark cyan") + 
	theme(legend.position = "none") + 
	theme(text = element_text(size=6)) +
	xlab("scaled evolutionary distance") + 
	ylab("mean proportion\nof variance") +
	scale_x_continuous(limits=c(0,1.04))

	pdf(file = "figures_and_panels/panel_prop_variation_by_phylogenetic_distance2.pdf",width=3,height=1.75,useDingbats=F)
	print(prop_phyl_dist2)
	dev.off()


n_homologs_phyl_dist <- ggplot(bind_rows(anova_phyl_dist),aes(x=node_depth,fill=VG)) + 
	geom_bar(width=0.05) +
	scale_fill_manual(values = c("dark gray","dark cyan","dark orange")) +
	theme(legend.position = "none") + 
	theme(text = element_text(size=6)) +
	xlab("scaled evolutionary distance") + 
	ylab("homolog groups") +
	scale_x_continuous(limits=c(0,1.04))

	pdf(file = "figures_and_panels/panel_homologs_by_phylogenetic_distance.pdf",width=3,height=1.6,useDingbats=F)
	print(n_homologs_phyl_dist)
	dev.off()

plot_anova_ov_PW_red <- ggplot(summarize_anova_results(anova_ov_PW_red_results),aes(y=prop_treatment,x=prop_species,color=VG)) + 
	geom_point(size=0.25) + 
	scale_color_manual(values = c("dark gray","dark cyan","dark orange")) +
	coord_fixed() +
	theme(legend.position="none") +
	theme(text = element_text(size=6)) +
	ylab("proportion of variance across tissue") + 
	xlab("proportion of variance across species")

	pdf(file = "figures_and_panels/panel_anova_plot_ovary_PictureWing_reduced.pdf",width=3,height=3,useDingbats=F)
	print(plot_anova_ov_PW_red)
	dev.off()

plot_anova_ov_red <- ggplot(summarize_anova_results(anova_ov_red_results),aes(y=prop_treatment,x=prop_species,color=VG)) + 
	geom_point(size=0.25) + 
	scale_color_manual(values = c("dark gray","dark cyan","dark orange"))  +
	scale_y_continuous(limits=c(0,1))  +
	theme(text = element_text(size=6)) +
	theme(legend.position="none") +
	ylab("proportion of variance across tissue") + 
	xlab("proportion of variance across species")

	pdf(file = "figures_and_panels/panel_anova_plot_ovary_reduced.pdf",width=3,height=3,useDingbats=F)
	print(plot_anova_ov_red)
	dev.off()

anova_ov_red_meanfoldchange <- summarize_anova_results(anova_ov_red_results) %>% select(VG,prop_treatment,prop_species,genetree) %>% 
	left_join(.,(r %>% filter(species %in% PNA) %>%
		filter(!is.na(genetree)) %>% 
		group_by(genetree) %>% 
		summarize(meanl2fc = mean(log2FoldChange)) %>% 
		ungroup()),
	by="genetree") 

plot_anova_meanfoldchange <- ggplot(anova_ov_red_meanfoldchange,aes(y = prop_treatment, x = meanl2fc,color=VG)) + 
	geom_vline(xintercept=0,linetype="dashed") +
	geom_point(size=0.25) + 
	scale_color_manual(values = c("dark gray","dark cyan","dark orange")) +
	theme(text = element_text(size=6)) +
	theme(legend.position="none") +
	ylab("proportion of variance across tissue") + 
	xlab("mean log2 fold change")

	pdf(file = "figures_and_panels/panel_anova_vs_meanlog2foldchange.pdf",width=3,height=3,useDingbats=F)
	print(plot_anova_meanfoldchange)
	dev.off()

plot_anova_meanfoldchange_TVGlabel <- ggplot(data=anova_ov_red_meanfoldchange %>% left_join(.,genetree_all_dmel_id,by="genetree") %>% filter(VG == "TVG",meanl2fc > -2),aes(y=prop_treatment,x=meanl2fc,label=name)) + 
	geom_point(size=0.25,color="dark orange") +
	geom_text_repel(max.overlaps=25,size=0.85,segment.size=0.1,box.padding=0.05) +
	theme(text = element_text(size=6)) +
	ylab("proportion of variance across tissue") + 
	xlab("mean log2 fold change") + 
	scale_x_continuous(limits=c(-1,8))

	pdf(file = "figures_and_panels/panel_anova_vs_meanlog2foldchange_TVGlabel.pdf",width=3,height=3,useDingbats=F)
	print(plot_anova_meanfoldchange_TVGlabel)
	dev.off()
