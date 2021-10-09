# perform each contrast calculation for absolute expression values for the head, and for the head-carcass ratio 
head_abs_changes <- get_evolutionary_changes(ave_ave_expression %>% 
	filter(treatment == "head") %>% 
	select(-treatment,-id),filtered_genetree_hashes,species_ultrametric_labeled) %>%
	left_join(.,genetree_all_dmel_id %>% rename(genetree = "key"),by="key")

head_rel_changes <- get_evolutionary_changes(ave_ave_ratio %>% 
	filter(treatment == "hd_ratio") %>% 
	select(-id,-treatment),filtered_genetree_hashes,species_ultrametric_labeled) %>%
	left_join(.,genetree_all_dmel_id %>% rename(genetree = "key"),by="key")

# for downstream analyses, use the ratio based values
head_changes <- head_rel_changes
full_keys <- head_changes %>% 
	filter(as.integer(child_node) < 13) %>% 
	group_by(key) %>% tally() %>% filter(n == 12) %>% 
	pull(key)

full_key_samp <- sample(full_keys,100)

panel_dist_changes <- ggplot(head_changes %>% filter(key %in% full_key_samp),aes(y = scaled_change,x=key)) + 
	geom_hline(yintercept = 0,color= "dark gray") + 
	geom_point(size=0.25) + 
	geom_point(data=head_changes %>% filter(child_node == 21,key %in% full_key_samp),color="red",size=0.25) + 
	theme(text = element_text(size=6)) +
	xlab("genes") + 
	ylab("scaled evolutionary change") + 
	scale_y_continuous(limits=c(-35,35)) +
	theme(axis.text.x = element_blank()) +
	theme(axis.ticks.x = element_blank()) 

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_dist_changes.pdf",width=4,height=2,useDingbats=F)
print(panel_dist_changes)
dev.off()

panel_hist_changes <- ggplot(head_changes,aes(y=scaled_change)) + geom_histogram(binwidth=0.1,fill="black") + scale_x_log10() + geom_histogram(data=head_changes %>% filter(child_node == 21),fill="red",binwidth=0.1) +
	scale_y_continuous(limits=c(-35,35)) +
	theme(text = element_text(size=6)) +
	xlab("count, log10 transformed") + 
	ylab("scaled evolutionary change") 

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_hist_changes.pdf",width=1.5,height=2.1,useDingbats=F)
print(panel_hist_changes)
dev.off()

node_order <- c("13","22","14","23","20","21","15","16","17","19","18","3","2","1","6","5","4","12","11","10","9","8","7")
dir_changes <- head_changes %>% mutate(change_type = ifelse(child_val < 0 & parent_val > 0,"FROM_HEAD",ifelse(child_val > 0 & parent_val < 0,"TO_HEAD","SAME")))

changes_by_nodes <- ggplot(dir_changes,aes(y=scaled_change,x=factor(child_node,levels=node_order),color=change_type)) + 
	geom_blank() +
	geom_jitter(size=0.25,data=dir_changes %>% filter(change_type == "SAME"),color="gray") +
	geom_jitter(size=0.25,data=dir_changes %>% filter(change_type == "FROM_HEAD"),color="blue") +
	geom_jitter(size=0.25,data=dir_changes %>% filter(change_type == "TO_HEAD"),color="red") +
	scale_y_continuous(limits=c(-35,35)) +
	theme(text = element_text(size=6)) +
	xlab("branch in phylogeny") + 
	ylab("scaled evolutionary change") + 
	theme(legend.position="none")

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_changes_by_nodes.pdf",width=5.2,height=2.5,useDingbats=F)
print(changes_by_nodes)
dev.off()

change1 <- dir_changes %>% filter(change_type == "TO_HEAD") %>% arrange(desc(change)) %>% filter(row_number() == 1)
change2 <- dir_changes %>% filter(change_type == "TO_HEAD") %>% arrange(desc(change)) %>% filter(row_number() == 2)
change3 <- dir_changes %>% filter(change_type == "FROM_HEAD") %>% arrange(change) %>% filter(row_number() == 1)
change4 <- dir_changes %>% filter(change_type == "FROM_HEAD") %>% arrange(change) %>% filter(row_number() == 2)

ch4_name <- seq_gene_ortho %>% filter(genetree == change4$key) %>% left_join(.,run_seq_header,by="seq_id") %>% filter(species %in% species_names)
ch4_name <- left_join(data.frame(sp = seq(1:length(species_names)),species = species_names),ch4_name,by="species")
ch4_name <- lapply(seq(1:nrow(ch4_name)),function(x){blast_list[[ch4_name[x,]$sp]][[1]] %>% filter(header == ch4_name[x,]$header)}) %>% bind_rows %>% pull(name) %>% unique()

parent_child_vals <- ggplot(dir_changes,aes(y=child_val,x=parent_val,color=change_type)) + 
	geom_point(size=0.25) + 
	geom_text_repel(data=change1,label="a",box.padding=0.25,segment.size=0.1,min.segment.length = 0) + 
	geom_text_repel(data=change2,label="b",box.padding=0.25,segment.size=0.1,min.segment.length = 0) + 
	geom_text_repel(data=change3,label="c",box.padding=0.25,segment.size=0.1,min.segment.length = 0) + 
	geom_text_repel(data=change4,label="d",box.padding=0.25,segment.size=0.1,min.segment.length = 0) + 
	scale_color_manual(values = c("blue","gray","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	geom_hline(yintercept=0,linetype="dashed",size=0.25) + 
	theme(legend.position="none") + 
	theme(text = element_text(size=6)) +
	xlab("ancestral bias in tissue expression") + 
	ylab("descendant bias in tissue expression") 

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_parent_child_vals.pdf",width=3,height=3,useDingbats=F)
print(parent_child_vals)
dev.off()

a <- ave_ave_ratio %>% filter(treatment == "hd_ratio") %>% 
	mutate(bias = ifelse(val > 0,"head","carcass"),
	sp_code = species_codes[species]) %>% 
	select(key,sp_code,val,bias) 

top_bias_change1 <- ggplot(a %>% filter(key == change1$key),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + geom_point() +
	scale_color_manual(values = c("blue","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	scale_x_continuous(limits = c(-5,5)) +
	ggtitle(change1$name) +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(text = element_text(size=6)) +
	theme(axis.title.y = element_blank()) +
	theme(axis.title.x = element_blank()) +
	theme(legend.position="none") 
top_bias_change2 <- ggplot(a %>% filter(key == change2$key),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + geom_point() +
	scale_color_manual(values = c("blue","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	scale_x_continuous(limits = c(-5,5)) +
	ggtitle(change2$name) +
	theme(plot.title = element_text(hjust = 0.5)) +
		theme(text = element_text(size=6)) +
	theme(axis.title.y = element_blank()) +
	theme(axis.title.x = element_blank()) +
	theme(legend.position="none") 
top_bias_change3 <- ggplot(a %>% filter(key == change3$key),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + geom_point() +
	scale_color_manual(values = c("blue","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	scale_x_continuous(limits = c(-5,5)) +
	ggtitle(change3$name) +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(text = element_text(size=6)) +
	theme(axis.title.y = element_blank()) +
	xlab("bias in tissue expression") +
	theme(legend.position="none") 
top_bias_change4 <- ggplot(a %>% filter(key == change4$key),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + geom_point() +
	scale_color_manual(values = c("blue","red")) + 
	geom_vline(xintercept=0,linetype="dashed",size=0.25) + 
	scale_x_continuous(limits = c(-5,5)) +
	ggtitle(paste(ch4_name,"*",sep="")) +
	theme(text = element_text(size=6)) +
	theme(plot.title = element_text(hjust = 0.5)) +
	theme(axis.title.y = element_blank()) +
	xlab("bias in tissue expression") +
	theme(legend.position="none") 

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_parent_examples.pdf",width=3,height=3,useDingbats=F)
grid.arrange(top_bias_change1,top_bias_change2,top_bias_change3,top_bias_change4,ncol=2)
dev.off()

change_summary <- head_changes %>% group_by(key) %>% filter(!is.na(scaled_change)) %>% 
	summarize(change_sd = sd(scaled_change),
		change_max = max(scaled_change),
		change_min = min(scaled_change),
		change_diff = change_max - change_min,
		change_var = var(scaled_change),
		change_n = n()) %>%
	left_join(.,genetree_all_dmel_id %>% rename(genetree = "key"),by="key")

full_cormat <- head_changes %>% 
	filter(key %in% full_keys) %>% 
	select(key,child_node,scaled_change) %>% 
	spread(.,key,scaled_change)  %>% 
	select(-child_node) %>% 
	cor(.,method="pearson",use="everything")

physical_cor <- sapply(seq(1:nrow(physical_interactions)),function(x){full_cormat[physical_interactions$genetree1[x],physical_interactions$genetree2[x]]})
unknown_physical_cor <- build_sample_unknown_cormat(full_cormat,physical_interactions)

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

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_physical_interaction_comparison.pdf",width=1.2,height=2,useDingbats=F)
print(physical_interaction_comparison)
dev.off()

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_genetic_interaction_comparison.pdf",width=1.5,height=2,useDingbats=F)
print(genetic_interaction_comparison)
dev.off()

physical_fcs <- replicate(100,build_sample_unknown_cormat(fc=full_cormat,interaction_dat=physical_interactions),simplify=F)
physical_t_test <- sapply(physical_fcs,function(x){t.test(physical_cor,x)$p.value})

genetic_fcs <- replicate(100,build_sample_unknown_cormat(fc=full_cormat,interaction_dat=bind_rows(enhance,suppress)),simplify=F)
enhance_t_test <- sapply(genetic_fcs,function(x){t.test(enhance_cor,x)$p.value})
suppress_t_test <- sapply(genetic_fcs,function(x){t.test(suppress_cor,x)$p.value})
enhance_suppress_t_test <- t.test(enhance_cor,suppress_cor)


head_phys_t_test <- physical_t_test
head_enh_t_test <- enhance_t_test
head_sup_t_test <- suppress_t_test
head_enh_sup_t_test <- enhance_suppress_t_test
save(head_phys_t_test,file="figures_and_panels/head_figures_and_panels/physical_t_test.RData")
save(head_enh_t_test,file="figures_and_panels/head_figures_and_panels/enhance_t_test.RData")
save(head_sup_t_test,file="figures_and_panels/head_figures_and_panels/suppress_t_test.RData")
save(head_enh_sup_t_test,file="figures_and_panels/head_figures_and_panels/enhance_suppress_t_test.RData")

# yolk protein gene family == 352dd19ce3c67e0ebb9b12e5d9189bb5
target <- genetree_dmel_id %>% filter(genetree == "352dd19ce3c67e0ebb9b12e5d9189bb5") %>% pull(parent_gene)
target_path <- mitab %>% filter(parent_gene1 %in% target | parent_gene2 %in% target) 
target_genetree <- genetree_dmel_id %>% filter(parent_gene %in% c(target_path$parent_gene1,target_path$parent_gene2)) %>%
	mutate(name = gsub(";.*","",name))
target_changes <- head_changes %>% filter(key %in% target_genetree$genetree)
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

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_real_target_interaction_correlations.pdf",width=3,height=3,useDingbats=F)
qgraph(mar=c(5,5,5,5),borders=F,vTrans=0,mat,maximum=1,layout="circle",posCol=cor_color_range[3],negCol=cor_color_range[2],labels=colnames(mat),label.cex=0.4,label.font=2,label.scale=F)
dev.off()

target_of_int <- target_genetree %>% filter(parent_gene %in% target) %>% pull(genetree) %>% unique
head_changes_cor <- left_join(head_changes,data.frame(key = names(full_cormat[target_of_int,]),correlation = full_cormat[target_of_int,]),by="key") %>% 
	na.omit

target_correlated_values <- ggplot(head_changes_cor,aes(x=factor(child_node,levels=node_order),y=scaled_change,group=name,color=correlation)) +
	viridis::scale_color_viridis(option = "B") + 
	geom_jitter(size=0.25,pch=16,aes(alpha=abs(correlation))) + 
	geom_point(data=head_changes %>% filter(name == target_name),color="white",fill="black",size=1.25,stroke=0.25,pch=21) + 
	theme(legend.position="none") + 
	theme(text = element_text(size=6)) +
	ylab("scaled evolutionary change") + 
	xlab("phylogenetic branch")

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_target_correlated_values.pdf",width=3.25,height=2,useDingbats=F)
print(target_correlated_values)
dev.off()

target_int_changes <- head_changes_cor %>% filter(abs(correlation) > correlation_threshold)
target_int_cormat <- target_int_changes %>% select(name,child_node,scaled_change) %>% 
	spread(.,name,scaled_change)  %>% 
	select(-child_node) %>% 
	cor(.,method="pearson",use="pairwise.complete.obs")

tarmat <- target_int_cormat
tarmat[tarmat < 10] <- NA 
tarmat[target_name,] <- target_int_cormat[target_name,]
tarmat[,target_name] <- target_int_cormat[,target_name]
tarmat[target_name,target_name] <- NA

pdf(file = "figures_and_panels/head_figures_and_panels/panel_head_observed_target_interaction_correlations.pdf",width=2,height=2,useDingbats=F)
qgraph(mar=c(5,5,5,5),borders=F,vTrans=0,tarmat,maximum=1,layout="circle",posCol=cor_color_range[3],negCol=cor_color_range[2],labels=colnames(tarmat),label.cex=0.4,label.font=2,label.scale=F)
dev.off()

