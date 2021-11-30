# calculate the average ratio of expression across homology group and biological replicate
pseudocount <- 0.01
ave_ave_ratio <- expression_genetree %>% 
	filter(!(library_id %in% resequenced_ids)) %>% select(-library_id) %>% 
	filter(!is.na(expression)) %>% # remove missing data
	mutate(expression = expression + pseudocount) %>% 
	pivot_wider(names_from = treatment,values_from = expression) %>% 
	mutate(ov_ratio = log(ovary / carcass), hd_ratio = log(head / carcass)) %>% 
	pivot_longer(c(ov_ratio,hd_ratio),names_to="treatment",values_to="ratio") %>% 
	filter(!is.na(ratio),!is.na(genetree)) %>%
	mutate(id = paste(species,treatment,sep="_"),key=genetree) %>% # choose ID variable
	group_by(key,individual,id,species,treatment) %>% summarize(mean_val = mean(ratio)) %>% ungroup() %>% 
	group_by(key,id,species,treatment) %>% summarize(val = mean(mean_val)) %>% ungroup()

ave_ratio <- expression_genetree %>% 
	filter(!(library_id %in% resequenced_ids)) %>% select(-library_id) %>% 
	filter(!is.na(expression)) %>% # remove missing data
	mutate(expression = expression + pseudocount) %>% 
	pivot_wider(names_from = treatment,values_from = expression) %>% 
	mutate(ov_ratio = log(ovary / carcass), hd_ratio = log(head / carcass)) %>% 
	pivot_longer(c(ov_ratio,hd_ratio),names_to="treatment",values_to="ratio") %>% 
	filter(!is.na(ratio),!is.na(genetree)) %>%
	mutate(id = paste(species,treatment,sep="_"),key=seq_id) %>% # choose ID variable
	group_by(key,id,species,treatment,genetree) %>% summarize(val = mean(ratio)) %>% ungroup()

build_pca_plot <- function(dat) {
	pca_df <- dat %>% mutate(species = species_codes[species]) %>% 
		tidyr::spread(.,key,val) %>%
		select_if(~ !any(is.na(.))) 
	pca_plot <- autoplot(prcomp(pca_df %>% select(-id,-species,-treatment)),data=pca_df,colour="species",shape="treatment",size=2.5) + 
		scale_color_brewer(palette="Paired") +
		#scale_shape_manual(values=c(17,15)) +
		theme(text = element_text(size=7)) +
		theme(legend.position = "none") 
		#coord_fixed()
	return(pca_plot)
}

homolog_ratio_pca <- build_pca_plot(ave_ave_ratio)

	pdf(file = "figures_and_panels/panel_homolog_ratio_pca.pdf",width=3,height=3,useDingbats=F)
	print(homolog_ratio_pca)
	dev.off()

species_ultrametric_focal <- drop.tip(species_ultrametric,species_ultrametric$tip.label[!(species_ultrametric$tip.label %in% species_names)])

# label nodes on tree
species_ultrametric_labeled <- species_ultrametric_focal
species_ultrametric_labeled$node.labels <- seq(1:species_ultrametric_labeled$Nnode) + length(species_ultrametric_labeled$tip.label)
species_tip_labels <- data.frame(species = species_ultrametric_labeled$tip.label, node = as.character(seq(1:length(species_ultrametric_labeled$tip.label))),stringsAsFactors=F)

# function for reconstructing expression contrasts on the species tree
# takes a hash for a given gene cluster (genetree_hash) used to pick homolog. genes
species_tree_asr <- function(genetree_hash,data,tree) {
	dat <- data %>% filter(key == genetree_hash) # select homologous genes in  cluster
	Rphylopars_dat <- dat %>% pull(val) # pull value for expression
	names(Rphylopars_dat) <- dat$species # build phylopars dataset
	drop_tree <- drop.tip(tree,tree$tip.label[!(tree$tip.label %in% dat$species)]) # select only species with data
	asr <- anc.recon(Rphylopars_dat,drop_tree,CI=TRUE) # perform ASR with confidence interval
	rownames(asr$Yhat) <- drop_tree$node.label
	rownames(asr$lowerCI) <- drop_tree$node.label
	rownames(asr$upperCI) <- drop_tree$node.label
	return(asr)
}


get_evolutionary_changes <- function(data,filt_hash,tree) {
	species_tip_labels <- data.frame(species = tree$tip.label, node = as.character(seq(1:length(tree$tip.label))),stringsAsFactors=F)

	# perform ASR on each gene cluster
	tree_cluster_asr <- lapply(filt_hash,species_tree_asr,
		data = data,
		tree = tree)

	# combine ASR results
	asr_node_vals <- lapply(seq(1:length(filt_hash)),function(x){data.frame(key = filt_hash[x], 
			val = tree_cluster_asr[[x]]$Yhat[,"V1"],
			node = rownames(tree_cluster_asr[[x]]$Yhat),stringsAsFactors=F)}) %>% 
		bind_rows()

	# get node and tip values
	tip_vals <- left_join(data %>% filter(key %in% filt_hash),species_tip_labels,by="species")
	all_node_vals <- bind_rows(asr_node_vals,tip_vals %>% select(-species))

	# organize values by 'child' and 'parent'
	# tree$edge lists nodes by parent[,1] to child[,2]
	# nodes can appear in both columns of $edge
	child_vals <- all_node_vals %>% 
		filter(node %in% tree$edge[,2]) %>% 
		rename(node = "child_node",val = "child_val")

	parent_vals <- all_node_vals %>% 
		filter(node %in% tree$edge[,1]) %>% 
		rename(node = "parent_node",val = "parent_val")

	# use the tree structure to create a set of child-parent pairs
	child_vals$parent_node <- plyr::mapvalues(child_vals$child_node,
		tree$edge[,2],
		tree$edge[,1])

	# use the tree structure to pull edge lengths by child node values
	branch_lengths <- data.frame(child_node = as.character(tree$edge[,2]), branch_length = tree$edge.length)
	# calculate contrasts as parent - child
	changes <- left_join(child_vals,parent_vals,by=c("key","parent_node")) %>% 
		left_join(.,branch_lengths,by="child_node") %>% 
		mutate(change = child_val - parent_val,
			scaled_change = change / branch_length) 

	return(changes)
}

# find trees with minimum number of tips, above threshold
min_tips_for_asr <- 3 
filtered_genetree_hashes <- ave_ave_expression %>% 
	distinct(key,species,.keep_all=T) %>% 
	group_by(key) %>% tally() %>% ungroup() %>% 
	filter(n >= min_tips_for_asr) %>% 
	pull(key)

# perform each contrast calculation for absolute expression values for the ovary, and for the ovary-carcass ratio 
ovary_abs_changes <- get_evolutionary_changes(ave_ave_expression %>% 
	filter(treatment == "ovary") %>% 
	select(-treatment,-id),filtered_genetree_hashes,species_ultrametric_labeled) %>%
	left_join(.,genetree_all_dmel_id %>% rename(genetree = "key"),by="key")

ovary_rel_changes <- get_evolutionary_changes(ave_ave_ratio %>% 
	filter(treatment == "ov_ratio") %>% 
	select(-id,-treatment),filtered_genetree_hashes,species_ultrametric_labeled) %>%
	left_join(.,genetree_all_dmel_id %>% rename(genetree = "key"),by="key") %>% 
	mutate(cd = as.numeric(child_node)) %>%
	mutate(branch = as.character(ifelse(cd < 13,cd, cd - 1))) %>%
	select(-cd)

# for downstream analyses, use the ratio based values
ovary_changes <- ovary_rel_changes 
full_keys <- ovary_changes %>% 
	filter(as.integer(child_node) < 13) %>% 
	group_by(key) %>% tally() %>% filter(n == 12) %>% 
	pull(key)

full_key_samp <- sample(full_keys,100)

panel_dist_changes <- ggplot(ovary_changes %>% filter(key %in% full_key_samp),aes(y = scaled_change,x=key)) + 
	geom_hline(yintercept = 0,color= "dark gray") + 
	geom_point(size=0.25) + 
	geom_point(data=ovary_changes %>% filter(child_node == 21,key %in% full_key_samp),color="red",size=0.25) + 
	theme(text = element_text(size=6)) +
	xlab("genes") + 
	ylab("scaled evolutionary change") + 
	scale_y_continuous(limits=c(-35,35)) +
	theme(axis.text.x = element_blank()) +
	theme(axis.ticks.x = element_blank()) 

	pdf(file = "figures_and_panels/panel_dist_changes.pdf",width=4,height=2,useDingbats=F)
	print(panel_dist_changes)
	dev.off()

panel_hist_changes <- ggplot(ovary_changes,aes(y=scaled_change)) + geom_histogram(binwidth=0.1,fill="black") + scale_x_log10() + geom_histogram(data=ovary_changes %>% filter(child_node == 21),fill="red",binwidth=0.1) +
	scale_y_continuous(limits=c(-35,35)) +
	theme(text = element_text(size=6)) +
	xlab("count, log10 transformed") + 
	ylab("scaled evolutionary change") 

	pdf(file = "figures_and_panels/panel_hist_changes.pdf",width=1.5,height=2.1,useDingbats=F)
	print(panel_hist_changes)
	dev.off()

node_order <- c("21","13","22","19","20","14","15","16","18","17","3","2","1","6","5","4","12","11","10","9","8","7")
dir_changes <- ovary_changes %>% mutate(change_type = ifelse(child_val < 0 & parent_val > 0,"FROM_OVARY",ifelse(child_val > 0 & parent_val < 0,"TO_OVARY","SAME"))) %>% 
	mutate(cd = as.numeric(child_node)) %>%
	mutate(branch = as.character(ifelse(cd < 13,cd, cd - 1))) %>%
	select(-cd)

changes_by_nodes <- ggplot(dir_changes,aes(y=scaled_change,x=factor(branch,levels=node_order),color=change_type)) + 
	geom_blank() +
	geom_jitter(size=0.25,data=dir_changes %>% filter(change_type == "SAME"),color="gray") +
	geom_jitter(size=0.25,data=dir_changes %>% filter(change_type == "FROM_OVARY"),color="blue") +
	geom_jitter(size=0.25,data=dir_changes %>% filter(change_type == "TO_OVARY"),color="red") +
	scale_y_continuous(limits=c(-35,35)) +
	theme(text = element_text(size=6)) +
	xlab("branch in phylogeny") + 
	ylab("scaled evolutionary change") + 
	theme(legend.position="none")

	pdf(file = "figures_and_panels/panel_changes_by_nodes.pdf",width=5.2,height=2.5,useDingbats=F)
	print(changes_by_nodes)
	dev.off()

change1 <- dir_changes %>% filter(change_type == "TO_OVARY") %>% arrange(desc(change)) %>% filter(row_number() == 1)
change2 <- dir_changes %>% filter(change_type == "TO_OVARY") %>% arrange(desc(change)) %>% filter(row_number() == 2)
change3 <- dir_changes %>% filter(change_type == "FROM_OVARY") %>% arrange(change) %>% filter(row_number() == 1)
change4 <- dir_changes %>% filter(change_type == "FROM_OVARY") %>% arrange(change) %>% filter(row_number() == 2)

ch4_name <- seq_gene_ortho %>% filter(genetree == change4$key) %>% left_join(.,run_seq_header,by="seq_id") %>% filter(species %in% species_names)
ch4_name <- left_join(data.frame(sp = seq(1:length(species_names)),species = species_names),ch4_name,by="species")
ch4_name <- lapply(seq(1:nrow(ch4_name)),function(x){blast_list[[ch4_name[x,]$sp]][[1]] %>% filter(header == ch4_name[x,]$header)}) %>% bind_rows %>% pull(name) %>% unique()

parent_child_vals <- ggplot(dir_changes,aes(y=child_val,x=parent_val,color=change_type)) + 
	geom_point(size=0.25) + 
	geom_text_repel(data=change1,label="a",box.padding=0.3,segment.size=0.15,min.segment.length = 0) + 
	geom_text_repel(data=change2,label="b",box.padding=0.5,segment.size=0.15,min.segment.length = 0) + 
	geom_text_repel(data=change3,label="c",box.padding=0.25,segment.size=0.15,min.segment.length = 0) + 
	geom_text_repel(data=change4,label="d",box.padding=0.3,segment.size=0.15,min.segment.length = 0) + 
	scale_color_manual(values = c("blue","gray","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	geom_hline(yintercept=0,linetype="dashed",size=0.25) + 
	theme(legend.position="none") + 
	theme(text = element_text(size=6)) +
	xlab("ancestral bias in tissue expression") + 
	ylab("descendant bias in tissue expression") 

	pdf(file = "figures_and_panels/panel_parent_child_vals.pdf",width=3,height=3,useDingbats=F)
	print(parent_child_vals)
	dev.off()

ave_ratio_summary <- ave_ratio %>% filter(treatment == "ov_ratio") %>% 
	mutate(bias = ifelse(val > 0,"ovary","carcass"),
	sp_code = species_codes[species]) %>% 
	select(genetree,species,sp_code,val,bias) 

top_bias_change1 <- ggplot(ave_ratio_summary %>% filter(genetree == change1$key),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + 
	geom_point(size=1) +
	scale_color_manual(values = c("blue","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	scale_x_continuous(limits = c(-5.5,5.5)) +
	ggtitle(change1$name) +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(text = element_text(size=6)) +
	theme(axis.title.y = element_blank()) +
	theme(axis.title.x = element_blank()) +
	theme(legend.position="none") 
top_bias_change2 <- ggplot(ave_ratio_summary %>% filter(genetree == change2$key),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + 
	geom_point(size=1) +	
	scale_color_manual(values = c("blue","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	scale_x_continuous(limits = c(-5.5,5.5)) +
	ggtitle(change2$name) +
	theme(plot.title = element_text(hjust = 0.5)) +
		theme(text = element_text(size=6)) +
	theme(axis.title.y = element_blank()) +
	theme(axis.title.x = element_blank()) +
	theme(legend.position="none") 
top_bias_change3 <- ggplot(ave_ratio_summary %>% filter(genetree == change3$key),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + 
	geom_point(size=1) +	
	scale_color_manual(values = c("blue","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	scale_x_continuous(limits = c(-5.5,5.5)) +
	ggtitle(change3$name) +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(text = element_text(size=6)) +
	theme(axis.title.y = element_blank()) +
	xlab("bias in tissue expression") +
	theme(legend.position="none") 
top_bias_change4 <- ggplot(ave_ratio_summary %>% filter(genetree == change4$key),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + 
	geom_point(size=1) +
	scale_color_manual(values = c("blue","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	scale_x_continuous(limits = c(-5.5,5.5)) +
	ggtitle(paste(ch4_name,"*",sep="")) +
	theme(text = element_text(size=6)) +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(axis.title.y = element_blank()) +
	xlab("bias in tissue expression") +
	theme(legend.position="none") 

	pdf(file = "figures_and_panels/panel_parent_examples.pdf",width=3,height=3,useDingbats=F)
	grid.arrange(top_bias_change1,top_bias_change2,top_bias_change3,top_bias_change4,ncol=2)
	dev.off()

change_summary <- ovary_changes %>% group_by(key) %>% filter(!is.na(scaled_change)) %>% 
	summarize(change_sd = sd(scaled_change),
		change_max = max(scaled_change),
		change_min = min(scaled_change),
		change_diff = change_max - change_min,
		change_var = var(scaled_change),
		change_n = n()) %>%
	left_join(.,genetree_all_dmel_id %>% rename(genetree = "key"),by="key")

full_cormat <- ovary_changes %>% 
	filter(key %in% full_keys) %>% 
	select(key,child_node,scaled_change) %>% 
	spread(.,key,scaled_change)  %>% 
	select(-child_node) %>% 
	cor(.,method="pearson",use="everything")


mitab <- read.delim("analysis/data/physical_interactions_mitab_fb_2021_04.tsv",header=F,stringsAsFactors=F,comment.char="#") %>% 
	mutate(parent_gene1 = gsub("flybase:(.*)","\\1",V1),
	parent_gene2 = gsub("flybase:(.*)","\\1",V2),
	interaction_type1 = gsub("psi-mi:MI:[0-9]+\\((.*)\\)","\\1",V21),
	interaction_type2 = gsub("psi-mi:MI:[0-9]+\\((.*)\\)","\\1",V22)) %>% 
	distinct(parent_gene1,parent_gene2,interaction_type1,interaction_type2)
physical_interactions <- mitab %>% left_join(.,(genetree_dmel_id %>% 
		mutate(genetree1 = genetree,parent_gene1 = parent_gene) %>% 
		select(genetree1,parent_gene1)),by="parent_gene1") %>% 
	left_join(.,(genetree_dmel_id %>% 
		mutate(genetree2 = genetree,parent_gene2 = parent_gene) %>% 
		select(genetree2,parent_gene2)),by="parent_gene2") %>%
	filter(genetree1 %in% rownames(full_cormat),genetree2 %in% rownames(full_cormat)) %>%
	filter(genetree1 != genetree2) %>% 
	na.omit 

build_sample_unknown_cormat <- function(fc,interaction_dat){
	sample_cormat1 <- sample(seq(1:nrow(fc)),10000,replace=T)
	sample_cormat2 <- sample(seq(1:nrow(fc)),10000,replace=T)
	sample_cormat <- data.frame(genetree1 = rownames(fc)[sample_cormat1],genetree2 = colnames(fc)[sample_cormat2]) %>% 
		filter(genetree1 != genetree2) %>% 
		distinct(genetree1,genetree2) %>%
		mutate(key = paste(genetree1,genetree2,sep="_")) %>%
		filter(!key %in% (interaction_dat %>% 
			mutate(key = paste(genetree1,genetree2,sep="_")) %>%
			pull(key))) %>%
		select(-key)
	unknown_cormat <- sample(sapply(seq(1:nrow(sample_cormat)),function(x){fc[sample_cormat$genetree1[x],sample_cormat$genetree2[x]]}),5000)
	return(unknown_cormat)
}

physical_cor <- sapply(seq(1:nrow(physical_interactions)),function(x){full_cormat[physical_interactions$genetree1[x],physical_interactions$genetree2[x]]})
unknown_physical_cor <- build_sample_unknown_cormat(full_cormat,physical_interactions)

ggi <- read.delim("analysis/data/gene_genetic_interactions_fb_2021_03.tsv",header=F,stringsAsFactors=F,comment.char="#")
names(ggi) <- c("name1","parent_gene1","name2","parent_gene2","interaction","rf")
gene_genetic_interactions <- ggi %>% separate_rows(.,c(name1,parent_gene1),sep="\\|") %>% separate_rows(.,c(name2,parent_gene2),sep="\\|")

genetic_interactions <- gene_genetic_interactions %>% left_join(.,(genetree_dmel_id %>% 
		mutate(genetree1 = genetree,parent_gene1 = parent_gene) %>% 
		select(genetree1,parent_gene1)),by="parent_gene1") %>% 
	left_join(.,(genetree_dmel_id %>% 
		mutate(genetree2 = genetree,parent_gene2 = parent_gene) %>% 
		select(genetree2,parent_gene2)),by="parent_gene2")

enhance <- genetic_interactions %>% filter(interaction == "enhanceable") %>% 
	filter(genetree1 %in% rownames(full_cormat),genetree2 %in% rownames(full_cormat)) %>% 
	filter(genetree1 != genetree2) %>% 
	na.omit
suppress <- genetic_interactions %>% filter(interaction == "suppressible") %>% 
	filter(genetree1 %in% rownames(full_cormat),genetree2 %in% rownames(full_cormat)) %>% 
	filter(genetree1 != genetree2) %>% 
	na.omit


enhance_cor <- sapply(seq(1:nrow(enhance)),function(x){full_cormat[enhance$genetree1[x],enhance$genetree2[x]]})
suppress_cor <- sapply(seq(1:nrow(suppress)),function(x){full_cormat[suppress$genetree1[x],suppress$genetree2[x]]})

unknown_genetic_cormat <- build_sample_unknown_cormat(full_cormat,bind_rows(enhance,suppress))

physical_interaction_comparison <- ggplot(bind_rows(data.frame(V1 = unknown_physical_cor,V2 = "unknown"),
	data.frame(V1=physical_cor,V2="physical")),aes(y=V1,x=factor(V2,levels=c("unknown","physical")))) + 
	geom_jitter(size=0.2,color="dark gray",width=0.2) +
	geom_boxplot(outlier.color = NA,color="black",fill=NA,lwd=0.2) +
	scale_y_continuous(limits=c(-1,1)) + 
	xlab("protein interactions") +
	ylab("correlation coefficient") +
	theme(text = element_text(size=6)) +
	theme(legend.position="none")

genetic_interaction_comparison <- ggplot(bind_rows(data.frame(V1 = unknown_genetic_cormat,V2 = "unknown"),
	data.frame(V1 = enhance_cor,V2 = "enhance"),
	data.frame(V1=suppress_cor,V2="suppress")),aes(y=V1,x=factor(V2,levels=c("unknown","suppress","enhance")),color=V2)) + 
	geom_jitter(size=0.2,color="dark gray",width=0.2) +
	geom_boxplot(outlier.color = NA,color="black",fill=NA,lwd=0.2) +
	scale_y_continuous(limits=c(-1,1)) + 
	xlab("genetic interactions") +
	ylab("correlation coefficient") +
	theme(text = element_text(size=6)) +
	theme(legend.position="none")

	pdf(file = "figures_and_panels/panel_physical_interaction_comparison.pdf",width=1.2,height=2,useDingbats=F)
	print(physical_interaction_comparison)
	dev.off()

	pdf(file = "figures_and_panels/panel_genetic_interaction_comparison.pdf",width=1.5,height=2,useDingbats=F)
	print(genetic_interaction_comparison)
	dev.off()

physical_fcs <- replicate(100,build_sample_unknown_cormat(fc=full_cormat,interaction_dat=physical_interactions),simplify=F)
physical_t_test <- sapply(physical_fcs,function(x){t.test(physical_cor,x)$p.value})

genetic_fcs <- replicate(100,build_sample_unknown_cormat(fc=full_cormat,interaction_dat=bind_rows(enhance,suppress)),simplify=F)
enhance_t_test <- sapply(genetic_fcs,function(x){t.test(enhance_cor,x)$p.value})
suppress_t_test <- sapply(genetic_fcs,function(x){t.test(suppress_cor,x)$p.value})
enhance_suppress_t_test <- t.test(enhance_cor,suppress_cor)

# yolk protein gene family == 352dd19ce3c67e0ebb9b12e5d9189bb5
target <- genetree_dmel_id %>% filter(genetree == "352dd19ce3c67e0ebb9b12e5d9189bb5") %>% pull(parent_gene)
target_path <- mitab %>% filter(parent_gene1 %in% target | parent_gene2 %in% target) 
target_genetree <- genetree_dmel_id %>% filter(parent_gene %in% c(target_path$parent_gene1,target_path$parent_gene2)) %>%
	mutate(name = gsub(";.*","",name))
target_changes <- ovary_changes %>% filter(key %in% target_genetree$genetree)
target_cormat <- target_changes %>% select(name,child_node,scaled_change) %>% 
	spread(.,name,scaled_change)  %>% 
	select(-child_node) %>% 
	cor(.,method="pearson",use="pairwise.complete.obs")

target_name <- "Yp1;Yp2;Yp3"
mat <- target_cormat
mat[mat < 100] <- 0
mat[target_name,] <- target_cormat[target_name,]
mat[,target_name] <- target_cormat[,target_name]
mat[target_name,target_name] <- 0

cor_color_option <- "B"
cor_color_range <- viridis::viridis(option=cor_color_option,n=4)

	pdf(file = "figures_and_panels/panel_real_target_interaction_correlations.pdf",width=3,height=3,useDingbats=F)
	qgraph(mar=c(5,5,5,5),borders=F,vTrans=0,mat,maximum=1,layout="circle",posCol=cor_color_range[3],negCol=cor_color_range[2],labels=colnames(mat),label.cex=0.4,label.font=2,label.scale=F)
	dev.off()

target_of_int <- target_genetree %>% filter(parent_gene %in% target) %>% pull(genetree) %>% unique
ovary_changes_cor <- left_join(ovary_changes,data.frame(key = names(full_cormat[target_of_int,]),correlation = full_cormat[target_of_int,]),by="key") %>% 
	na.omit %>% 
	mutate(cd = as.numeric(child_node)) %>%
	mutate(branch = as.character(ifelse(cd < 13,cd, cd - 1))) %>%
	select(-cd)

target_correlated_values <- ggplot(ovary_changes_cor,aes(x=factor(branch,levels=node_order),y=scaled_change,group=name,color=correlation)) +
	viridis::scale_color_viridis(option = "B") + 
	geom_jitter(size=0.25,pch=16,aes(alpha=abs(correlation))) + 
	geom_point(data=ovary_changes %>% filter(name == target_name),color="white",fill="black",size=1.25,stroke=0.25,pch=21) + 
	theme(legend.position="none") + 
	theme(text = element_text(size=6)) +
	ylab("scaled evolutionary change") + 
	xlab("phylogenetic branch")

	pdf(file = "figures_and_panels/panel_target_correlated_values.pdf",width=3.25,height=2,useDingbats=F)
	print(target_correlated_values)
	dev.off()

correlation_threshold <- 0.825
target_int_changes <- ovary_changes_cor %>% filter(abs(correlation) > correlation_threshold)
target_int_cormat <- target_int_changes %>% select(name,child_node,scaled_change) %>% 
	spread(.,name,scaled_change)  %>% 
	select(-child_node) %>% 
	cor(.,method="pearson",use="pairwise.complete.obs")

tarmat <- target_int_cormat
tarmat[tarmat < 10] <- NA 
tarmat[target_name,] <- target_int_cormat[target_name,]
tarmat[,target_name] <- target_int_cormat[,target_name]
tarmat[target_name,target_name] <- NA

	pdf(file = "figures_and_panels/panel_observed_target_interaction_correlations.pdf",width=2,height=2,useDingbats=F)
	qgraph(mar=c(5,5,5,5),borders=F,vTrans=0,tarmat,maximum=1,layout="circle",posCol=cor_color_range[3],negCol=cor_color_range[2],labels=colnames(tarmat),label.cex=0.4,label.font=2,label.scale=F)
	dev.off()

