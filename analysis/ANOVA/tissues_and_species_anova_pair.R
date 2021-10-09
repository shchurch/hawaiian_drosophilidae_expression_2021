# LINEAR MODELING WITH TWO SPECIES ONLY 

# get species pair combinations
species_pairs <- combn(species_names,2,simplify=F)
pair_names <- lapply(combn(species_codes,2,simplify=F),paste,collapse="_")

# create and run anova
run_pair_anova <- function(pair,name){
	anova_ov_pair_table <- expression_genetree %>% 
		filter(species %in% pair) %>% 
		filter(treatment %in% c("ovary","carcass")) %>% # focus on comparison of ovary and body
		filter(!library_id %in% resequenced_ids) %>% # remove resequenced libraries
		mutate(key = paste("sample_",library_id,sep=""),id=genetree) %>% # use lib ID as future column name, genetree as row name
		filter(!is.na(expression)) %>% # filter out seqIDs that didn't have expression measured
		filter(!is.na(genetree)) %>% # filter out seqIDs that didn't match a genetree
		group_by(key,id) %>% # unique lib ID and genetree
		summarize(val = mean(expression)) %>% # add up all transcripts in that group
		ungroup() %>%
		tidyr::spread(.,key,val) # pivot to long form data

	anova_ov_pair_tab <- anova_ov_pair_table %>% select(-id) %>% as.data.frame() 
	rownames(anova_ov_pair_tab) <- anova_ov_pair_table %>% pull(id)
	anova_ov_pair_red_tab <- na.omit(anova_ov_pair_tab)
		write.table(anova_ov_pair_red_tab,file=paste("analysis/ANOVA/ANOVA_results/anova_ovary_",name,"_reduced_table.tsv",sep=""),row.names=T,col.names=T,sep="\t")

	anova_command <- paste("Rscript analysis/ANOVA/Breschi_anova.R -i analysis/ANOVA/ANOVA_results/anova_ovary_",name,"_reduced_table.tsv -o analysis/ANOVA/ANOVA_results/anova_ovary_",name,"_reduced_table_results.tsv -p 0.01 -m analysis/ANOVA/ANOVA_results/anova_metadata.tsv --merge_mdata_on id -F treatment+species",sep="")
		system(anova_command)
	anova_ov_pair_red_results <- read.table(paste("analysis/ANOVA/ANOVA_results/anova_ovary_",name,"_reduced_table_results.tsv",sep=""),header=T,stringsAsFactors=F)
	return(anova_ov_pair_red_results)
}

anova_pair_results <- lapply(seq(1:length(species_pairs)),function(x){run_pair_anova(species_pairs[[x]],pair_names[[x]])})

anova_phyl_dist <- lapply(anova_pair_results,summarize_anova_results)

mean_prop_treatment <- sapply(anova_phyl_dist,function(x){x %>% pull(prop_treatment) %>% mean()})
mean_prop_species <- sapply(anova_phyl_dist,function(x){x %>% pull(prop_species) %>% mean()})
n_homologs <- sapply(anova_phyl_dist,nrow)
all_node_depths <- node.depth.edgelength(species_tree_focal)
scaled_all_node_depths <- (1-all_node_depths * (1 / max(all_node_depths)))

get_node_depths <- function(pair){
	return(round(scaled_all_node_depths[getMRCA(species_tree_focal,pair)],2))
}

pair_depths <- lapply(species_pairs,get_node_depths)

mean_props <- bind_rows(data.frame(depth=simplify(pair_depths),type="species",mean_prop=mean_prop_species),
	data.frame(depth=simplify(pair_depths),type="treatment",mean_prop=mean_prop_treatment))

ggplot(data = mean_props,aes(y = mean_prop, x = depth, color = type)) + 
	geom_point()  +
	geom_smooth(method='lm')+ 
	theme(legend.position = "none") + 
	theme(text = element_text(size=6)) +
	xlab("scaled evolutionary distance") + 
	ylab("mean proportion\nof variance") +
	scale_x_continuous(limits=c(0,1.04)) + 
	scale_color_manual(values=c("dark cyan","dark orange"))

anova_phyl_dist <- lapply(seq(1:length(pair_depths)),function(x){anova_phyl_dist[[x]] %>% mutate(node_depth = pair_depths[x])})

prop_phyl_dist <- ggplot(data = mean_props,aes(y = mean_prop, x = depth, color = type)) + 
	geom_point()  +
	geom_smooth(method='lm')+ 
	theme(legend.position = "none") + 
	theme(text = element_text(size=6)) +
	xlab("scaled evolutionary distance") + 
	ylab("mean proportion\nof variance") +
	scale_x_continuous(limits=c(0,1.04)) + 
	scale_color_manual(values=c("dark cyan","dark orange"))

pdf(file = "figures_and_panels/panel_pairwise_prop_variation_by_phylogenetic_distance.pdf",width=3,height=1.75,useDingbats=F)
print(prop_phyl_dist)
dev.off()



