### BUILD AND RUN ANOVA DATASETS FOR LINEAR MODELING WITH ORTHOLOGS
anova_ov_table <- expression_genetree %>% 
	filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
	filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
	mutate(key = paste("sample_",library_id,sep=""),id=ortholog) %>% # use lib ID as future column name, genetree as row name
	filter(!is.na(expression),!is.na(ortholog)) %>% # filter out seqIDs that didn't match an ortholog group
	mutate(val = expression) %>% select(key,id,val) %>% 
	tidyr::spread(.,key,val) # pivot to long form data

anova_ov_tab <- anova_ov_table %>% select(-id) %>% as.data.frame() 
rownames(anova_ov_tab) <- anova_ov_table %>% pull(id)
anova_ov_red_tab <- na.omit(anova_ov_tab)
	
# RUN AFTER REMOVING NAs
	write.table(anova_ov_red_tab,file="analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_reduced_table.tsv",row.names=T,col.names=T,sep="\t")
anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_reduced_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_reduced_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species"
	system(anova_command)
anova_ov_red_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_reduced_table_results.tsv",header=T,stringsAsFactors=F)

# PLOT REDUCED ORTHOLOG DATASET
plot_anova_ov <- ggplot(summarize_anova_results(anova_ov_red_results),aes(y=prop_treatment,x=prop_species,color=VG)) + 
	geom_point(size=0.25) + 
	scale_color_manual(values = c("dark gray","dark cyan","dark orange"))  +
	scale_y_continuous(limits=c(0,1))  +
	theme(text = element_text(size=6)) +
	theme(legend.position="none") +
	ylab("proportion of variance across tissue") + 
	xlab("proportion of variance across species")

	pdf(file = "figures_and_panels/panel_ortholog_anova_plot_ovary.pdf",width=3,height=3,useDingbats=F)
	print(plot_anova_ov)
	dev.off()

### BUILD AND RUN REDUCED ANOVA ON HI DROSOPHILA CLADE
anova_ov_HIDROS_table <- expression_genetree %>% 
	filter(species %in% c(PNA,PRIM,HAL,AMC,MM)) %>% 
	filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
	filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
	mutate(key = paste("sample_",library_id,sep=""),id=ortholog) %>% # use lib ID as future column name, genetree as row name
	filter(!is.na(expression),!is.na(ortholog)) %>% # filter out seqIDs that didn't match a genetree
	mutate(val = expression) %>% select(key,id,val) %>% 
	tidyr::spread(.,key,val) # pivot to long form data

anova_ov_HIDROS_tab <- anova_ov_HIDROS_table %>% select(-id) %>% as.data.frame() 
rownames(anova_ov_HIDROS_tab) <- anova_ov_HIDROS_table %>% pull(id)
anova_ov_HIDROS_red_tab <- na.omit(anova_ov_HIDROS_tab)
	write.table(anova_ov_HIDROS_red_tab,file="analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_HIDROS_reduced_table.tsv",row.names=T,col.names=T,sep="\t")

anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_HIDROS_reduced_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_HIDROS_reduced_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species"
	system(anova_command)
anova_ov_HIDROS_red_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_HIDROS_reduced_table_results.tsv",header=T,stringsAsFactors=F)

### BUILD AND RUN REDUCED ANOVA ON PNA CLADE
anova_ov_PNA_table <- expression_genetree %>% 
	filter(species %in% PNA) %>% 
	filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
	filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
	mutate(key = paste("sample_",library_id,sep=""),id=ortholog) %>% # use lib ID as future column name, genetree as row name
	filter(!is.na(expression),!is.na(ortholog)) %>% # filter out seqIDs that didn't match a genetree
	mutate(val = expression) %>% select(key,id,val) %>% 
	tidyr::spread(.,key,val) # pivot to long form data

anova_ov_PNA_tab <- anova_ov_PNA_table %>% select(-id) %>% as.data.frame() 
rownames(anova_ov_PNA_tab) <- anova_ov_PNA_table %>% pull(id)
anova_ov_PNA_red_tab <- na.omit(anova_ov_PNA_tab)
	write.table(anova_ov_PNA_red_tab,file="analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_PNA_reduced_table.tsv",row.names=T,col.names=T,sep="\t")

anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_PNA_reduced_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_PNA_reduced_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species"
	system(anova_command)
anova_ov_PNA_red_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_PNA_reduced_table_results.tsv",header=T,stringsAsFactors=F)

### BUILD AND RUN REDUCED ANOVA ON PW CLADE
anova_ov_PW_table <- expression_genetree %>% 
	filter(species %in% c("Drosophila_sproati","Drosophila_macrothrix")) %>% 
	filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
	filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
	mutate(key = paste("sample_",library_id,sep=""),id=ortholog) %>% # use lib ID as future column name, genetree as row name
	filter(!is.na(expression),!is.na(ortholog)) %>% # filter out seqIDs that didn't match a genetree
	mutate(val = expression) %>% select(key,id,val) %>% 
	tidyr::spread(.,key,val) # pivot to long form data

anova_ov_PW_tab <- anova_ov_PW_table %>% select(-id) %>% as.data.frame() 
rownames(anova_ov_PW_tab) <- anova_ov_PW_table %>% pull(id)
anova_ov_PW_red_tab <- na.omit(anova_ov_PW_tab)
	write.table(anova_ov_PW_red_tab,file="analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_PW_reduced_table.tsv",row.names=T,col.names=T,sep="\t")

anova_command <- "Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_PW_reduced_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_PW_reduced_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species"
	system(anova_command)
anova_ov_PW_red_results <- read.table("analysis/ANOVA/ANOVA_results/anova_ortholog_ovary_PW_reduced_table_results.tsv",header=T,stringsAsFactors=F)

# PLOT ANOVA RESULTS WITH ORTHOLOGS IN PW DATASET
plot_anova_ov_PW_red <- ggplot(summarize_anova_results(anova_ov_PW_red_results),aes(y=prop_treatment,x=prop_species,color=VG)) + 
	geom_point(size=0.25) + 
	scale_color_manual(values = c("dark gray","dark cyan","dark orange"))  +
	scale_y_continuous(limits=c(0,1))  +
	theme(text = element_text(size=6)) +
	theme(legend.position="none") +
	ylab("proportion of variance across tissue") + 
	xlab("proportion of variance across species")

	pdf(file = "figures_and_panels/panel_ortholog_anova_plot_ovary_PictureWing_reduced.pdf",width=3,height=3,useDingbats=F)
	print(plot_anova_ov_PW_red)
	dev.off()

### SUMMARIZE AND PLOT ANOVA RESULTS
anova_phyl_dist <- lapply(list(anova_ov_red_results,anova_ov_HIDROS_red_results,anova_ov_PNA_red_results,anova_ov_PW_red_results),summarize_anova_results)

mean_prop_treatment <- sapply(anova_phyl_dist,function(x){x %>% pull(prop_treatment) %>% mean()})
mean_prop_species <- sapply(anova_phyl_dist,function(x){x %>% pull(prop_species) %>% mean()})
n_homologs <- sapply(anova_phyl_dist,nrow)

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

	pdf(file = "figures_and_panels/panel_ortholog_prop_variation_by_phylogenetic_distance.pdf",width=3,height=1.75,useDingbats=F)
	print(prop_phyl_dist)
	dev.off()

# calculate the average epxression across orthology group and biological replicate
ave_ortholog_expression <- expression_genetree %>% 
	filter(!(library_id %in% resequenced_ids)) %>% # choose unique sequence libraries
	mutate(id = paste(species,treatment,sep="_"),key=ortholog) %>% # choose ID variable
	filter(!is.na(expression),!is.na(ortholog)) %>% # remove missing data
	group_by(key,id,species,treatment) %>% summarize(val = mean(expression)) %>% ungroup()

# calculate the average ratio of expression across orthology group and biological replicate
ave_ortholog_ratio <- expression_genetree %>% 
	filter(!(library_id %in% resequenced_ids)) %>% select(-library_id) %>% 
	filter(!is.na(expression)) %>% # remove missing data
	mutate(expression = expression + pseudocount) %>% 
	pivot_wider(names_from = treatment,values_from = expression) %>% 
	mutate(ov_ratio = log(ovary / carcass), hd_ratio = log(head / carcass)) %>% 
	pivot_longer(c(ov_ratio,hd_ratio),names_to="treatment",values_to="ratio") %>% 
	filter(!is.na(ratio),!is.na(genetree)) %>%
	mutate(id = paste(species,treatment,sep="_"),key=ortholog) %>% # choose ID variable
	group_by(key,id,species,treatment) %>% summarize(val = mean(ratio)) %>% ungroup()

ortholog_ratio_pca <- build_pca_plot(ave_ortholog_ratio)

pdf(file = "figures_and_panels/panel_ortholog_ratio_pca.pdf",width=3,height=3,useDingbats=F)
print(ortholog_ratio_pca)
dev.off()

species_ultrametric_focal <- drop.tip(species_ultrametric,species_ultrametric$tip.label[!(species_ultrametric$tip.label %in% species_names)])

# label nodes on tree
species_ultrametric_labeled <- species_ultrametric_focal
species_ultrametric_labeled$node.labels <- seq(1:species_ultrametric_labeled$Nnode) + length(species_ultrametric_labeled$tip.label)
species_tip_labels <- data.frame(species = species_ultrametric_labeled$tip.label, node = as.character(seq(1:length(species_ultrametric_labeled$tip.label))),stringsAsFactors=F)

# find trees with minimum number of tips, above threshold
filtered_genetree_hashes <- ave_ortholog_expression %>% 
	distinct(key,species,.keep_all=T) %>% 
	group_by(key) %>% tally() %>% ungroup() %>% 
	filter(n >= min_tips_for_asr) %>% 
	pull(key)

ortho_list <- unlist(orthologs,recursive=F)
get_dmel_seq_id_ortho <- function(o_num) {
	ortho_seqs <- gsub("Drosophila_melanogaster@","",ortho_list[[o_num]]$tip.label[grepl("melanogaster",ortho_list[[o_num]]$tip.label)])
	return(data.frame(ortholog=rep(unlist(ortholog_hashes,recursive=F)[[o_num]],length(ortho_seqs)),seq_id=as.integer(ortho_seqs)))
}
# annotate genetrees with D melanogaster header and name
ortholog_dmel_seqs <- pblapply(seq(1:length(unlist(orthologs,recursive=F))),get_dmel_seq_id_ortho)
ortholog_dmel_id <- left_join(bind_rows(ortholog_dmel_seqs),dmel_seqs,by="seq_id")

# annotate genetrees with combined headers and names from DMEL
ortholog_all_dmel_id <- ortholog_dmel_id %>% group_by(ortholog) %>% arrange(name) %>% summarize(name = paste(name,collapse=";"),annotation = paste(annotation,collapse=";"),parent_gene = paste(parent_gene,collapse=";"))

# perform each contrast calculation for absolute expression values for the ovary, and for the ovary-carcass ratio 
ovary_abs_changes <- get_evolutionary_changes(ave_ortholog_expression %>% 
	filter(treatment == "ovary") %>% 
	select(-treatment,-id),filtered_genetree_hashes,species_ultrametric_labeled) %>%
	left_join(.,ortholog_all_dmel_id %>% rename(ortholog = "key"),by="key")

ovary_rel_changes <- get_evolutionary_changes(ave_ortholog_ratio %>% 
	filter(treatment == "ov_ratio") %>% 
	select(-id,-treatment),filtered_genetree_hashes,species_ultrametric_labeled) %>%
	left_join(.,ortholog_all_dmel_id %>% rename(ortholog = "key"),by="key")

# for downstream analyses, use the ratio based values
ovary_changes <- ovary_rel_changes

node_order <- c("13","22","14","23","20","21","15","16","17","19","18","3","2","1","6","5","4","12","11","10","9","8","7")
dir_changes <- ovary_changes %>% mutate(change_type = ifelse(child_val < 0 & parent_val > 0,"FROM_OVARY",ifelse(child_val > 0 & parent_val < 0,"TO_OVARY","SAME")))

changes_by_nodes <- ggplot(dir_changes,aes(y=scaled_change,x=factor(child_node,levels=node_order),color=change_type)) + 
	geom_blank() +
	geom_jitter(size=0.25,data=dir_changes %>% filter(change_type == "SAME"),color="gray") +
	geom_jitter(size=0.25,data=dir_changes %>% filter(change_type == "FROM_OVARY"),color="blue") +
	geom_jitter(size=0.25,data=dir_changes %>% filter(change_type == "TO_OVARY"),color="red") +
	scale_y_continuous(limits=c(-35,35)) +
	theme(text = element_text(size=6)) +
	xlab("branch in phylogeny") + 
	ylab("scaled evolutionary change") + 
	theme(legend.position="none")

pdf(file = "figures_and_panels/panel_ortholog_changes_by_nodes.pdf",width=5.2,height=2.5,useDingbats=F)
print(changes_by_nodes)
dev.off()

