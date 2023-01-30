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
	filter(genetree1 %in% rownames(cormat_no_outliers),genetree2 %in% rownames(cormat_no_outliers)) %>%
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

physical_cor <- sapply(seq(1:nrow(physical_interactions)),function(x){cormat_no_outliers[physical_interactions$genetree1[x],physical_interactions$genetree2[x]]})
unknown_physical_cor <- build_sample_unknown_cormat(cormat_no_outliers,physical_interactions)

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
	filter(genetree1 %in% rownames(cormat_no_outliers),genetree2 %in% rownames(cormat_no_outliers)) %>% 
	filter(genetree1 != genetree2) %>% 
	na.omit
suppress <- genetic_interactions %>% filter(interaction == "suppressible") %>% 
	filter(genetree1 %in% rownames(cormat_no_outliers),genetree2 %in% rownames(cormat_no_outliers)) %>% 
	filter(genetree1 != genetree2) %>% 
	na.omit

enhance_cor <- sapply(seq(1:nrow(enhance)),function(x){cormat_no_outliers[enhance$genetree1[x],enhance$genetree2[x]]})
suppress_cor <- sapply(seq(1:nrow(suppress)),function(x){cormat_no_outliers[suppress$genetree1[x],suppress$genetree2[x]]})

unknown_genetic_cormat <- build_sample_unknown_cormat(cormat_no_outliers,bind_rows(enhance,suppress))

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

	pdf(file = "figures_and_panels/panel_physical_interaction_comparison_no_outliers.pdf",width=1.2,height=2,useDingbats=F)
	print(physical_interaction_comparison)
	dev.off()

	pdf(file = "figures_and_panels/panel_genetic_interaction_comparison_no_outliers.pdf",width=1.5,height=2,useDingbats=F)
	print(genetic_interaction_comparison)
	dev.off()

physical_fcs <- replicate(100,build_sample_unknown_cormat(fc=cormat_no_outliers,interaction_dat=physical_interactions),simplify=F)
physical_t_test <- sapply(physical_fcs,function(x){t.test(physical_cor,x)$p.value})

genetic_fcs <- replicate(100,build_sample_unknown_cormat(fc=cormat_no_outliers,interaction_dat=bind_rows(enhance,suppress)),simplify=F)
enhance_t_test <- sapply(genetic_fcs,function(x){t.test(enhance_cor,x)$p.value})
suppress_t_test <- sapply(genetic_fcs,function(x){t.test(suppress_cor,x)$p.value})
enhance_suppress_t_test <- t.test(enhance_cor,suppress_cor)

# yolk protein gene family == 5daa431b904a92eb1d01639824a469b6
target <- genetree_dmel_id %>% filter(genetree == "5daa431b904a92eb1d01639824a469b6") %>% pull(parent_gene)
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

	pdf(file = "figures_and_panels/panel_real_target_interaction_correlations_no_outliers.pdf",width=3,height=3,useDingbats=F)
	qgraph(mar=c(5,5,5,5),borders=F,vTrans=0,mat,maximum=1,layout="circle",posCol=cor_color_range[3],negCol=cor_color_range[2],labels=colnames(mat),label.cex=0.4,label.font=2,label.scale=F)
	dev.off()

target_of_int <- target_genetree %>% filter(parent_gene %in% target) %>% pull(genetree) %>% unique
ovary_changes_cor <- left_join(ovary_changes,data.frame(key = names(cormat_no_outliers[target_of_int,]),correlation = cormat_no_outliers[target_of_int,]),by="key") %>% 
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

	pdf(file = "figures_and_panels/panel_target_correlated_values_no_outliers.pdf",width=3.25,height=2,useDingbats=F)
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

	pdf(file = "figures_and_panels/panel_observed_target_interaction_correlations_no_outliers.pdf",width=2,height=2,useDingbats=F)
	qgraph(mar=c(5,5,5,5),borders=F,vTrans=0,tarmat,maximum=1,layout="circle",posCol=cor_color_range[3],negCol=cor_color_range[2],labels=colnames(tarmat),label.cex=0.4,label.font=2,label.scale=F)
	dev.off()
