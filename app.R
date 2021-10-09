library(shiny)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ape)
library(qgraph)
theme_set(theme_classic())
set.seed(84095)
load(file="analysis/results_Robjects/genetree_all_dmel_id.RData")
load(file="analysis/results_Robjects/ovary_changes.RData")
load(file="analysis/results_Robjects/full_correlation_matrix.RData")
load(file="analysis/results_Robjects/ovary_average_average_ratio.RData")
source("analysis/phylogenetic_expression/qgraph_as_ggraph.R")

species_names <- c('Drosophila_primaeva','Drosophila_sproati','Drosophila_picticornis','Drosophila_macrothrix','Drosophila_mimica','Drosophila_cfdives','Drosophila_nanella','Scaptomyza_varia','Scaptomyza_cyrtandrae','Scaptomyza_varipicta','Drosophila_atroscutellata','Drosophila_tanythrix')
species_id <- c('008D','106A','025A','055A','040C','16_1','002D','CFB','088B','020A','029A','043D');names(species_id) <- species_names
species_codes <- c('Dpri','Dspr','Dpic','Dmac','Dmim','Dcfd','Dnan','Svar','Scyr','Svpt','Datr','Dtan');names(species_codes) <- species_names
species_order <-  c('Drosophila_primaeva','Drosophila_mimica','Drosophila_nanella','Drosophila_atroscutellata','Drosophila_tanythrix','Drosophila_cfdives','Drosophila_sproati','Drosophila_macrothrix','Drosophila_picticornis','Scaptomyza_cyrtandrae','Scaptomyza_varipicta','Scaptomyza_varia')


selected_gene_families <- c("Yp1;Yp2;Yp3","Octbeta1R;Octbeta2R","Doa","Haspin","nos")
gene_family_names <- genetree_all_dmel_id %>% filter(genetree %in% rownames(full_cormat)) %>% filter(!name %in% selected_gene_families) %>% pull(name) %>% unique %>% sort

gene_family_names <- c(selected_gene_families,gene_family_names)


# Define UI for miles per gallon app ----
ui <- pageWithSidebar(

  # App title ----
  headerPanel("Evolutionary changes in ovary-biased expression across 12 species of Hawaiian Drosophilidae flies"),

  # Sidebar panel for inputs ----
  sidebarPanel(

  		helpText("This R::Shiny app allows you to explore the correlation of evolutionary changes in the ovary-biased expression using RNA sequence data from twelve species of Hawaiian Drosophilidae flies.  You can select any gene family for which to display expression data."),
  		selectizeInput(inputId="variable",choices=gene_family_names,label="Select a family of genes",options = list(create = TRUE)),
  		helpText("Note: not all gene families will be represented. For the purposes of this app we have restricted to only gene families represented in all the datasets of twelve species."),
  		helpText("Gene families have been defined as groups of homologous genes, as inferred with the software agalma (https://bitbucket.org/caseywdunn/agalma), which uses sequence similarity to cluster genes."),

  ),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Plot of the requested variable against mpg ----
      plotOutput("qplot")

    )
)

# Data pre-processing ----
# Tweak the "am" variable to have nicer factor labels -- since this
# doesn't rely on any user inputs, we can do this once at startup
# and then use the value throughout the lifetime of the app


cor_color_option <- "B"
cor_color_range <- viridis::viridis(option=cor_color_option,n=4)


# Define server logic to plot various variables against mpg ----
server <- function(input, output) {

  # Compute the formula text ----
  # This is in a reactive expression since it is shared by the
  # output$caption and output$mpgPlot functions
  target_genetree <- reactive({
    #target <- input$variable
	#target_all_dmel <- dmel_seqs %>% filter(grepl(target,name))
	#target_parent_gene <- dmel_seqs %>% filter(name == target) %>% pull(parent_gene) %>% .[1] 
	#target_genetree <- genetree_dmel_id %>% filter(parent_gene == target_parent_gene) %>% pull(genetree) %>% unique
	#target_name <- genetree_all_dmel_id %>% filter(genetree == target_genetree) %>% pull(name)
	
  })

  # Generate a plot of the requested variable against mpg ----
  # and only exclude outliers if requested
  output$qplot <- renderPlot({

	target <- input$variable
  target_genetree <- genetree_all_dmel_id %>% filter(name == target) %>% pull(genetree) %>% unique
  ovary_changes_cor <- left_join(ovary_changes,data.frame(key = names(full_cormat[target_genetree,]),
  															correlation = full_cormat[target_genetree,]),by="key") %>% 
  							na.omit
  bias_colors <- setNames(c("red","blue"),c("ovary","carcass"))

	plot_exp <- ggplot(a %>% filter(key == target_genetree),aes(y=factor(sp_code,levels=species_codes[species_order]),x=val,color = bias)) + geom_point(size=3) +
		scale_color_manual(values = bias_colors) + 
		geom_vline(xintercept=0,linetype="dashed",size=0.5) + 
		scale_x_continuous(limits = c(-5,5)) +
		ggtitle("ovary:carcass bias") +
		theme(plot.title = element_text(hjust = 0.5)) +
		theme(text = element_text(size=15)) +
		theme(axis.title.y = element_blank()) +
		theme(axis.title.x = element_blank()) +
		theme(legend.position="none") 


	correlation_threshold <- 0.825
	target_changes <- ovary_changes_cor %>% filter(abs(correlation) > correlation_threshold)
	target_cormat <- target_changes %>% select(name,child_node,scaled_change) %>% 
		spread(.,name,scaled_change)  %>% 
		select(-child_node) %>% 
		cor(.,method="pearson",use="pairwise.complete.obs")

	node_order <- c("13","22","14","23","20","21","15","16","17","19","18","3","2","1","6","5","4","12","11","10","9","8","7")

	plot_tree <- grid::rasterGrob(png::readPNG("figures_and_panels/labeled_ultrametric_tree-01.png"))
	plot_changes <- ggplot(ovary_changes_cor,aes(x=factor(child_node,levels=node_order),y=scaled_change,group=name,color=correlation)) +
		viridis::scale_color_viridis(option = "B") + 
		geom_jitter(size=1.5,pch=16,aes(alpha=abs(correlation))) + 
		geom_point(data=ovary_changes %>% filter(name == target),color="white",fill="black",size=3,stroke=1,pch=21) + 
		theme(legend.position="none") + 
		theme(text = element_text(size=10)) +
		ylab("scaled evolutionary\nchange") + 
		xlab("phylogenetic branch")

	tarmat <- target_cormat
	tarmat[tarmat < 10] <- NA 
	tarmat[target,] <- target_cormat[target,]
	tarmat[,target] <- target_cormat[,target]
	tarmat[target,target] <- NA

	plot_network <- as.ggraph(qgraph(mar=c(2,2,2,2),border.width=0.1,tarmat,maximum=1,layout="circle",esize=10,posCol=cor_color_range[3],
		negCol=cor_color_range[2],labels=colnames(tarmat),color="white",
		label.font=2,label.scale=F))

	layout <- rbind(c(1,3),c(2,3),c(4,4),c(5,5),c(6,6),c(6,6))

	t1 <- grid::textGrob("Above is the phylogeny of the 12 species\nstudied here, labeled by phylogenetic branch.\n\nTo the right shows the average expression\nratio between the ovary and\ncarcass for each species\nRed=more ovary expression\nBlue=more carcass epxression.",just="left",x=0,y=0.6)
	t2 <- grid::textGrob("Above shows the distribution of evolutionary changes along each branch in the phylogeny.\nThe black point indicates the selected gene family. Colored points are changes\nin other gene families.\nBelow is the network of genes with a strong evolutionary correlation of expression to the selected gene family.\nFor both panels, orange=strong positive correlation, purple=strong negative correlation.",just="left",x=0,y=0.5)

	grid.arrange(grobs = list(plot_tree,t1,plot_exp,plot_changes,t2,plot_network), layout_matrix = layout)
  },height=1000,width=500)

}

shinyApp(ui, server)