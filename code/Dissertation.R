##### ORGANISATION #####
#### Data entry ####
library(readxl) # for reading xlsx
library(vegan) # for ecological analyses
library(tidyr) # for tidying
library(plyr) # for organising data
library(dplyr) # for organising data
library(wrapr) # for merging datasets
library(FSA) # for Dunn test
library(car) # for Levene test
library(bbmle) # for AIC test
library(factoextra)
library(Factoshiny)
library(MASS) # for StepAIC
# library(wCorr/weights/WGCNA) # for weighted correlation
library(multcomp)
library(goeveg)
library(ordinal) # for ordinal regression
library(ape)
library(phytools)
library(viridis)


# For nice plots
# library(BiodiversityR)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(ggforce)
library(ggpubr)
# library(ggord)
library(ggbiplot)
# library(Cairo)
library(gridExtra)
library(grid)
library(ggfortify)
library(ggtree)

# windowsFonts(LMR12=windowsFont("LMRoman12-Regular"))

# Read Data
speciesAbundances<-read_xlsx("../data/Data Analysis.xlsx", sheet = 1)
speciesBinomial<-read_xlsx("../data/Data Analysis.xlsx", sheet = 2)
specGen<-read_xlsx("../data/Data Analysis.xlsx", sheet = 3)
habitatVariables<-read_xlsx("../data/Data Analysis.xlsx", sheet = 4)
specChar.final<-read_xlsx("../data/Data Analysis.xlsx", sheet=7)
specTraits<-read_xlsx("../data/Data Analysis.xlsx", sheet=8)
newSpecTraits<-read_xlsx("../data/Traits.xlsx", sheet=3)

# attach(habitatVariables)
# attach(speciesAbundances)
# attach(speciesBinomial)
# attach(specGen)
# attach(specTraits)
# attach(specChar.final)
# attach(newSpecTraits)

# habitatVariables$TreeShadeArea<-asin(sqrt(habitatVariables$TreeShadeArea/100)) # Arcsine transformation as values are percentages
# habitatVariables$RelativeDryness<-asin(sqrt(habitatVariables$RelativeDryness/100)) # Arcsine transformation as values are percentages
# habitatVariables$RelativeDisturbance<-asin(sqrt(habitatVariables$RelativeDisturbance/100))
# Helper: 0–100 %  -> logit, with small clamp for 0/1
logit01 <- function(x, eps = 1e-6) {
  p <- x / 100
  p <- pmin(pmax(p, eps), 1 - eps)
  qlogis(p)
}

cols <- c("TreeShadeArea", "RelativeDryness", "RelativeDisturbance")
habitatVariables[cols] <- lapply(habitatVariables[cols], logit01)
habitatVariables$Group<-factor(habitatVariables$Group)

#### Organisation ####
speciesAbundances$Habitats<-NULL
rownames(speciesAbundances)<-c("Evergreen",	"Deciduous",	"Olive-Grove",	"Maquis",	"Phrygana",	"Successional-Meadow",	"Agricultural",	"Dunes",	
                               "Marshlands",	"Lowland-River",	"Highland-River",	"Lagoon",	"Salt-Flats",	"Large-Town",	"Small-Town",	"SoEBSP")

speciesBinomial$Habitats<-NULL
rownames(speciesBinomial)<-c("Evergreen",	"Deciduous",	"Olive-Grove",	"Maquis",	"Phrygana",	"Successional-Meadow",	"Agricultural",	"Dunes",	
                             "Marshlands",	"Lowland-River",	"Highland-River",	"Lagoon",	"Salt-Flats",	"Large-Town",	"Small-Town",	"SoEBSP")

specGen$Altitude<-NULL # Duplicate
colnames(habitatVariables)<-c("Habitats","Group", "Alt","DTW","TVH","SVH", "FTP", "SFPP","OP","RS","TSA","RDr", "RDi")

specGen<-left_join(specGen,habitatVariables, by="Habitats")
habitatVariables$Habitats<-NULL
rownames(habitatVariables)<-c("Evergreen",	"Deciduous",	"Olive-Grove",	"Maquis",	"Phrygana",	"Successional-Meadow",	"Agricultural",	"Dunes",	
                              "Marshlands",	"Lowland-River",	"Highland-River",	"Lagoon",	"Salt-Flats",	"Large-Town",	"Small-Town",	"SoEBSP")

# Get mean of altitude and number per species with variables
speciesMean <- ddply(specChar.final, c("Species"), summarise,
                     mean_Altitude  = mean(Altitude),
                     mean_Number = mean(Number))
species.Mean.Traits<-left_join(speciesMean,specTraits, by="Species")

species.Mean.Traits$Gregariousness<-as.factor(species.Mean.Traits$Gregariousness)
species.Mean.Traits$Myrmecophily<-as.factor(species.Mean.Traits$Myrmecophily)
species.Mean.Traits$Generation_numbers<-as.factor(species.Mean.Traits$Generation_numbers)
species.Mean.Traits$Host_plant_apparency<-as.factor(species.Mean.Traits$Host_plant_apparency)
species.Mean.Traits$Migration<-as.factor(species.Mean.Traits$Migration)
species.Mean.Traits$Overwintering_stage<-as.factor(species.Mean.Traits$Overwintering_stage)
species.Mean.Traits$Range_size<-as.factor(species.Mean.Traits$Range_size)
species.Mean.Traits$Range_type<-as.factor(species.Mean.Traits$Range_type)
species.Mean.Traits$Voltinism<-as.factor(species.Mean.Traits$Voltinism)
species.Mean.Traits$PuL_shrub_layer_Sl<-as.factor(species.Mean.Traits$PuL_shrub_layer_Sl)
species.Mean.Traits$PuL_canopy_layer_Cl<-as.factor(species.Mean.Traits$PuL_canopy_layer_Cl)
species.Mean.Traits$PuL_buried_Bu<-as.factor(species.Mean.Traits$PuL_buried_Bu)
species.Mean.Traits$PuL_attended_At<-as.factor(species.Mean.Traits$PuL_attended_At)
species.Mean.Traits$PuL_fieldlayer_Gl<-as.factor(species.Mean.Traits$PuL_fieldlayer_Gl)
species.Mean.Traits$PuL_groundlayer_Gl<-as.factor(species.Mean.Traits$PuL_groundlayer_Gl)




#### ORDINATION ANALYSES ####
#### Log-transformed Abundance Analysis #### 
speciesAbundances_log<-log(speciesAbundances+1) #[c(1,3,5,6,7,9,10,11,12,14,16,18,20,22,23,25,26,32,35,39,40,44,46),]
decorana(speciesAbundances_log)

rda.1<- rda(speciesAbundances_log ~ Alt + DTW + TVH + SVH + FTP + SFPP + OP + RS + TSA + RDr + RDi, data=habitatVariables) # Maximal Model
rda.1<- ordistep(rda.1,direction="both", Pout = 0.05, permutations = how(nperm = 9000)) # Perform simplification

summary(rda.1)
anova.cca(rda.1, step=999, by = "axis") # 66.77% Explained variance
round(100*(summary(rda.1)$cont$importance[2, 1:2]), 2) # 37.68 & 12.07 = 49.75% 1st & 2nd axes.

# Plotting
# par(family="LMRoman10-Regular")
# ordipointlabel(rda.1, choices = c(1,2), col = c("black", "#FC4E07"), pch = 1:2, font = c(1,3), scaling = 1, yaxt = "n")
# axis(side = 2, las = 2)
# bp <- scores(rda.1, display = 'bp', scaling=1)
# mul <- ordiArrowMul(bp, fill = 0.75)
# arrows(0, 0, mul * bp[,1], mul * bp[,2], length = 0.05, col = "#07B5FC")
# text(ordiArrowTextXY(mul * bp, rownames(bp)), rownames(bp), col = "#07B5FC")


# Extract scores for sites, species, and biplot (bp)
site_scores <- vegan::scores(rda.1, display = "sites", scaling = 1)
species_scores <- vegan::scores(rda.1, display = "species", scaling = 1)
bp_scores <- vegan::scores(rda.1, display = "bp", scaling = 1)

site_df <- as.data.frame(site_scores)
site_df$Habitat <- rownames(site_df)  # Assuming habitat names are rownames

species_df <- as.data.frame(species_scores)
species_df$Label <- rownames(species_df)

bp_df <- as.data.frame(bp_scores)
bp_df$Label <- rownames(bp_df)

# Determine the arrow multiplication factor
mul <- ordiArrowMul(bp_scores, fill = 1)
bp_df$RDA1 <- mul * bp_df$RDA1
bp_df$RDA2 <- mul * bp_df$RDA2

# Create the plot
plot_rda_habitats <- ggplot() +
  # Plot site scores
  geom_point(data = site_df, aes(x = RDA1, y = RDA2), color = "#5f136e", shape = 16, size = 1.5, alpha = 1) +
  geom_text_repel(data = site_df, aes(x = RDA1, y = RDA2, label = Habitat), 
                  family = "arial", fontface = "italic", size = 3, color = "#5f136e") +
  # Plot species scores
  geom_point(data = species_df, aes(x = RDA1, y = RDA2), color = "#e65d2f", shape = 17, size = 1.5, alpha = 0.3) +
  # geom_text_repel(data = species_df, aes(x = RDA1, y = RDA2, label = Label), 
  #                 family = "LMRoman10-Regular", fontface = "italic", size = 3.5, color = "#FC4E07", alpha = 0.3) +
  # Plot biplot arrows and labels
  geom_segment(data = bp_df, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = bp_df, aes(x = RDA1, y = RDA2, label = Label), 
                  color = "black", size = 3.5, family = "arial") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  # Customize plot appearance
  theme_classic() +
  theme(text = element_text(family = "arial", size = 12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(title = NULL, x = "RDA1 (37.7%)", y = "RDA2 (12.1%)")
# Create the plot
plot_rda_species <- ggplot() +
  # Plot site scores
  geom_point(data = site_df, aes(x = RDA1, y = RDA2), color = "#5f136e", shape = 16, size = 1.5, alpha = 0.3) +
  # geom_text_repel(data = site_df, aes(x = RDA1, y = RDA2, label = Habitat), 
  #                 family = "LMRoman10-Regular", size = 3, color = "#00AFBB") +
  # Plot species scores
  geom_point(data = species_df, aes(x = RDA1, y = RDA2), color = "#e65d2f", shape = 17, size = 1.5, alpha = 1) +
  geom_text(data = species_df, aes(x = RDA1, y = RDA2 + 0.07, label = Label),
                  family = "arial", fontface = "italic", size = 3.5, color = "#e65d2f", alpha = 1) +
  # Plot biplot arrows and labels
  geom_segment(data = bp_df, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = bp_df, aes(x = RDA1, y = RDA2, label = Label), 
                  color = "black", size = 3.5, family = "arial") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  # Customize plot appearance
  theme_classic() +
  theme(text = element_text(family = "arial", size = 12), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(title = NULL, x = "RDA1 (37.7%)", y = "RDA2 (12.1%)")

plot_rda = ggarrange(plot_rda_habitats+ rremove("xlab"), plot_rda_species + rremove("xylab"), align = "hv", common.legend = TRUE, legend = 'right')
annotate_figure(plot_rda, bottom = text_grob("RDA1 (37.7%)", family = "arial"))
  

  

# #### Trait Ordination ####
# #Organisation
# newSpecTraits<-data.frame(newSpecTraits)
# rownames(newSpecTraits) <- newSpecTraits$Species
# newSpecTraits$Species<-NULL
# options(ggrepel.max.overlaps = Inf)


# #PCA

# #K-Means Cluster
# res.PCA<-PCA(newSpecTraits,ncp=Inf,scale.unit=FALSE,graph=FALSE)
# res.HCPC<-HCPC(res.PCA,graph=FALSE)

# plot2<-fviz_cluster(res.HCPC,font.family="LMRoman10-Regular")+
#   theme_classic()+ 
#   theme(text=element_text(family="LMRoman10-Regular", size=12),legend.position = c(0.92,0.85))+
#   scale_color_manual('Cluster', values=c("#00AFBB","#8E44AD", "#FC4E07"),labels=c("Intermediate","Specialist","Generalist"))+
#   scale_fill_manual('Cluster', values=c("#00AFBB","#8E44AD", "#FC4E07"),labels=c("Intermediate","Specialist","Generalist"))+
#   scale_shape_manual('Cluster', values=c(19,19,19),labels=c("Intermediate","Specialist","Generalist"))+
#   rremove("xlab") + rremove("ylab")



# rda.3<-prcomp(newSpecTraits,center = TRUE)

# species.colours = join(plot[["data"]], plot2[["data"]], by = "name")
# species.colours <- species.colours %>%
#   mutate(colours = case_when(
#     cluster == 1 ~ "#00AFBB",
#     cluster == 2 ~ "#8E44AD",
#     cluster == 3 ~ "#FC4E07"))

# biplot<-fviz_pca_biplot(rda.3,repel=TRUE, col.var = "black", col.ind = species.colours$cluster,font.family="LMRoman10-Regular", select.var=list(contrib = 10), label = "var", alpha.ind = 0.4)+
#   theme_classic()+
#   scale_color_manual('Cluster', values=c("#00AFBB", "#8E44AD", "#FC4E07"),labels=c("Specialist","Intermediate","Generalist"))+
#   scale_fill_manual('Cluster', values=c("#00AFBB","#8E44AD", "#FC4E07"),labels=c("Specialist","Intermediate","Generalist"))+
#   scale_shape_manual('Cluster', values=c(16,15,17),labels=c("Specialist", "Intermediate","Generalist"))+
#   theme(text=element_text(family="LMRoman10-Regular", size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.92,0.85))+
#   labs(title = NULL, x = "PC1 (26.1%)", y = "PC2 (12.4%)")+
#   geom_text_repel(label = species.colours$name, alpha = .5, color = alpha(alpha(species.colours$colours), .8), family = "LMRoman10-Regular", fontface = "italic")
# biplot 



# Levins Niche Breadth
levins<-read_xlsx("../data/Traits.xlsx", sheet=4)
levins<-data.frame(levins)
library(MicroNiche)

# Calculate Levins Breadth based on individual habitats
sampleInfo <- colnames(levins[,-1])
levins_output = levins.Bn(levins,16, sampleInfo)

# Calculate Levins Breadth based on grouped abundances (forest, grassland, aquatic, urban)
levins_categorised = levins
levins_categorised$forest = rowSums(levins_categorised[,c(2:4)])
levins_categorised$grassland = rowSums(levins_categorised[,c(5:8)])
levins_categorised$aquatic = rowSums(levins_categorised[,c(9:14)])
levins_categorised$urban = rowSums(levins_categorised[,c(15:16)])
levins_categorised = levins_categorised[,c(1,18:21)]
sampleInfo = colnames(levins_categorised[,-1])
levins_categorised_output = levins.Bn(levins,4, sampleInfo, q = )



# Phylogenetics
tree = read.tree('../data/EUROPEAN_BUTTERFLIES_FULLMCC_DROPTIPED.nwk')

tree$tip.label = sub("_", " ", tree$tip.label)
species_names = c("Iphiclides podalirius", "Papilio alexanor", "Papilio machaon", "Carcharodus alceae", "Ochlodes sylvanus", "Pyrgus armoricanus", 
                  "Pyrgus malvae", "Spialia orbifer", "Thymelicus acteon", "Thymelicus sylvestris", "Colias crocea", "Euchloe ausonia", "Gonepteryx cleopatra", 
                  "Gonepteryx rhamni", "Leptidea sinapis", "Pieris brassicae", "Pieris mannii", "Pieris rapae", "Pontia edusa", "Aricia agestis", 
                  "Cacyreus marshalli", "Celastrina argiolus", "Favonius quercus", "Lampides boeticus", "Leptotes pirithous", "Lycaena ottomana", 
                  "Lycaena phlaeas", "Polyommatus icarus", "Satyrium acaciae", "Satyrium ilicis", "Argynnis pandora", "Argynnis paphia", "Brintesia circe", 
                  "Charaxes jasius", "Coenonympha pamphilus", "Hipparchia syriaca", "Hipparchia volgensis", "Kirinia roxelana", "Lasiommata maera", 
                  "Lasiommata megera", "Libythea celtis", "Limenitis reducta", "Maniola jurtina", "Melanargia larissa", "Melitaea cinxia", "Melitaea didyma", 
                  "Nymphalis polychloros", "Pararge aegeria", "Polygonia egea", "Vanessa atalanta", "Vanessa cardui")
species_abbreviations = c("M.jur",	"G.cleo",	"M.lar",	"P.bra",	"P.ica",	"T.act",	"C.croc",	"B.cir",	"O.sylv",	"L.sin",	"C.arg",	
                                      "K.rox",	"T.sylv",	"H.volg",	"L.phle",	"A.ages", 	"M.did",	"A.paph",	"P.rap",	"L.red",	"V.car",	"L.cel",	
                                      "S.aca",	"G.rham",	"L.piri",	"H.syr",	"I.pod",	"L.meg",	"S.ili",	"M.cin",	"V.ata",	"P.aeg",	"N.pol",	
                                      "P.egea",	"C.jas",	"P.alex",	"C.alc",	"P.edu",	"C.mars",	"C.pam",	"L.otto",	"L.mae",	"F.quer",	"L.boe",
                                      "P.malv","P.man",	"A.pand",	"P.mach",	"P.armo",	"S.orb",	"E.aus")

tree_subset = keep.tip(tree, tree$tip.label[tree$tip.label %in% unlist(species_names)])
plot(tree_subset)

# Step 1: Generate short_species_names from species_names
short_species_names <- sapply(
  species_names,
  function(x) {
    parts <- strsplit(x, " ")[[1]]
    paste0(substr(parts[1], 1, 1), ".", parts[2])
  }
)

# Step 2: Build a correctly ordered data frame
sort_order <- order(short_species_names)
mapping_df <- data.frame(
  species_name = species_names[sort_order],
  short_name = short_species_names[sort_order],
  abbreviation = sort(species_abbreviations),
  stringsAsFactors = FALSE
)

# Step 3: Replace tree tip labels by matching species_name
tree_subset$tip.label <- trimws(tree_subset$tip.label)
label_map <- setNames(mapping_df$abbreviation, mapping_df$species_name)
tree_subset$tip.label <- label_map[tree_subset$tip.label]

plot(tree_subset)

# Step 1: Count non-zero values for each species
non_zero_counts <- colSums(speciesAbundances > 0)

# Step 2: Create a data frame mapping tree labels to counts
tip_data <- data.frame(
  label = names(non_zero_counts),
  total_abundance = non_zero_counts
)

# Step 3: Color the full tree by non-zero abundance counts
p <- ggtree(tree_subset) %<+% tip_data +
  geom_tree(aes(color = as.integer(total_abundance)), size = 0.8) +
  geom_tiplab(size = 2.5) +
  scale_color_gradient(name = "Habitats", low = "#fcfdbf", high = "#e24d66") +
  theme_tree2() +
  # ggtitle("Tree Colored by Total Abundance per Species") +
  theme(legend.position = "right") +
  geom_text2(aes(subset = !isTip), label = NA)  # <- hides internal node labels


p + theme_void()


trait_vector <- colSums(speciesAbundances > 0)[tree_subset$tip.label]
lambda_result <- phylosig(tree_subset, trait_vector, method = "lambda", test = TRUE)
print(lambda_result)





plotBranchbyTrait(tree_subset,non_zero_counts,"tips",palette=colorRampPalette(c("#fcfdbf","#e24d66")),show.tip.label=TRUE)


library(ggtree)
library(phytools)
library(dplyr)

tip_trait <- non_zero_counts[tree_subset$tip.label]
anc       <- fastAnc(tree_subset, tip_trait)                 # internal nodes

node_df <- data.frame(
  node  = 1:(Ntip(tree_subset) + tree_subset$Nnode),
  trait = c(tip_trait, anc)
)


ggtree(tree_subset, layout = "fan") %<+% node_df +
  geom_tree(aes(colour = trait), size = 0.8) +
  geom_tiplab(size = 2.5, colour = "darkgray") +           # solid label colour
  scale_colour_viridis_c(
    option    = "inferno",
    direction = -1,              # ← reverses the palette
    name      = "Habitats",
    limits    = range(non_zero_counts)   # keeps legend at 1–16
  ) +
  theme_tree2() + theme_void() +
  theme(legend.position = "right")



































# #### TRAIT ANALYSES - NUMERIC#### 
#   #Reload speciesAbundances & habitatVariables before running
# speciesTraits.long<-left_join(speciesAbundances,habitatVariables,"Habitats")
# speciesTraits.long<-reshape(data = speciesTraits.long,
#                             idvar = "Habitats",
#                             varying = c(2:52),
#                             v.name = c("Number"),
#                             times = c("M.jur",	"G.cleo",	"M.lar",	"P.bra",	"P.ica",	"T.act",	"C.croc",	"B.cir",	"O.sylv",	"L.sin",	"C.arg",	
#                                       "K.rox",	"T.sylv",	"H.volg",	"L.phle",	"A.ages", 	"M.did",	"A.paph",	"P.rap",	"L.red",	"V.car",	"L.cel",	
#                                       "S.aca",	"G.rham",	"L.piri",	"H.syr",	"I.pod",	"L.meg",	"S.ili",	"M.cin",	"V.ata",	"P.aeg",	"N.pol",	
#                                       "P.egea",	"C.jas",	"P.alex",	"C.alc",	"P.edu",	"C.mars",	"C.pam",	"L.otto",	"L.mae",	"F.quer",	"L.boe",
#                                       "P.malv","P.man",	"A.pand",	"P.mach",	"P.armo",	"S.orb",	"E.aus"),
#                             new.row.names = 1:2000,
#                             direction = "long")
# colnames(speciesTraits.long)[14]<-"Species"
# speciesTraits.long<-left_join(speciesTraits.long, specTraits, "Species")
# speciesTraits.long$Group<-NULL
# speciesTraits.long$AL1<-NULL   
# speciesTraits.long$AL2<-NULL 
# speciesTraits.long$AL3<-NULL 
# speciesTraits.long$FL1<-NULL 
# speciesTraits.long$FL2<-NULL 
# speciesTraits.long$FL3<-NULL 
# speciesTraits.long$FL4<-NULL 
# speciesTraits.long$FL5<-NULL 
# speciesTraits.long$`LFM_–_flower`<-NULL 
# speciesTraits.long$`LFM_–_leaf`<-NULL 
# speciesTraits.long<-filter(speciesTraits.long, Number>0)



# #### Species Wingspan ####
# speciesTraits.long1<-na.omit(speciesTraits.long)
# speciesTraits.long1<-ddply(speciesTraits.long1, c("Habitats"            ,         "Altitude"              ,      
#                                                   "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                   "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                   "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                   "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                   "RelativeDryness"             , "RelativeDisturbance" ), summarise,
#                            Wingspan = mean(Wingspan))
# ggboxplot(data=speciesTraits.long1, x= "OrnamentalPresence", y = "Wingspan")

# shapiro.test(speciesTraits.long1$Wingspan) # Independent variables all p<0.05 (except RelativeDryness), Wingspan p>0.05

# m1<-glm(speciesTraits.long1$Wingspan ~ speciesTraits.long1$Altitude)
# m2<-glm(speciesTraits.long1$Wingspan ~ 1)
# anova(m1,m2,test="Chi")

# p1<-ggplot(speciesTraits.long1,aes(Altitude,Wingspan))+
#   geom_smooth(method="lm",se=FALSE)+
#   geom_point()+
#   theme_bw()+
#   theme(text=element_text(family="LMR12", size=12, face="bold"))



# #### Species Feeding Index ####
# speciesTraits.long1<-na.omit(speciesTraits.long)
# speciesTraits.long1<-ddply(speciesTraits.long1, c("Habitats"            ,         "Altitude"              ,      
#                                                   "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                   "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                   "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                   "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                   "RelativeDryness"             , "RelativeDisturbance" ), summarise,
#                            Feeding_index = mean(Feeding_index))


# shapiro.test(speciesTraits.long1$Feeding_index) # Independent variables all p<0.05 (except RelativeDryness), Wingspan p>0.05

# m1<-glm(speciesTraits.long1$Feeding_index ~ speciesTraits.long1$StrongFloweringPlantPresence)
# m2<-glm(speciesTraits.long1$Feeding_index ~ 1)
# anova(m1,m2,test="Chi")

# ggboxplot(data=speciesTraits.long1, x= "RockStructures", y = "Feeding_index")
# p2<-ggplot(speciesTraits.long1,aes(Altitude,Feeding_index))+
#   geom_smooth(method="lm",se=FALSE)+
#   geom_point()+
#   theme_bw()+
#   theme(text=element_text(family="LMR12", size=12, face="bold"))
# p2<-p2+labs(y="Feeding Index")
# p2

# #### TRAIT ANALYSES - ORDINAL ####
# #### Voltinism ####
# speciesTraits.long1<-ddply(speciesTraits.long, c("Habitats"            ,         "Altitude"              ,      
#                                                  "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                  "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                  "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                  "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                  "RelativeDryness"             , "RelativeDisturbance" , "Voltinism"), summarise,
#                            NumberOfVoltinism=sum(Voltinism),
#                            SumVoltinism = sum(Number))

# speciesTraits.long1$NumberOfVoltinism<-speciesTraits.long1$NumberOfVoltinism/speciesTraits.long1$Voltinism
# speciesTraits.long1$Voltinism<-as.factor(speciesTraits.long1$Voltinism)
# speciesTraits.long1$Voltinism<-speciesTraits.long1$Voltinism

# p3<-ggplot(speciesTraits.long1,aes(Altitude,NumberOfVoltinism, group=Voltinism))+
#   geom_smooth(method="glm", method.args = list(family = 'poisson'),aes(colour=Voltinism),se=FALSE)+
#   geom_point()+
#   theme_bw()+
#   theme(text=element_text(family="LMR12", size=12, face="bold"),legend.position = c(0.5,0.82), legend.key.size = unit(0.25,'cm'))
# p3<-p3+labs(color="Voltinism",y= "Number of Species")
# p3
# ggboxplot(data=speciesTraits.long1, x= "Voltinism", y = "NumberOfVoltinism")

# m1<-glm(NumberOfVoltinism~as.factor(Voltinism) * StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# m2<-glm(NumberOfVoltinism~as.factor(Voltinism) + StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# anova(m1,m2, test="Chi")


# #### Hostplant Apparency ####
# speciesTraits.long1<-ddply(speciesTraits.long, c("Habitats"            ,         "Altitude"              ,      
#                                                  "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                  "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                  "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                  "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                  "RelativeDryness"             , "RelativeDisturbance" , "Host_plant_apparency"), summarise,
#                            NumberOfHostplantApparency=sum(Host_plant_apparency),
#                            SumHostplantApparency = sum(Number))

# speciesTraits.long1$NumberOfHostplantApparency<-speciesTraits.long1$NumberOfHostplantApparency/speciesTraits.long1$Host_plant_apparency
# speciesTraits.long1$Host_plant_apparency<-as.factor(speciesTraits.long1$Host_plant_apparency)

# ggplot(speciesTraits.long1,aes(Altitude,NumberOfHostplantApparency, group=Host_plant_apparency))+
#   geom_smooth(method="glm", method.args = list(family = 'poisson'),aes(colour=Host_plant_apparency),se=FALSE)+
#   geom_point()+
#   theme_bw()
# ggboxplot(data=speciesTraits.long1, x= "Host_plant_apparency", y = "NumberOfHostplantApparency")

# m1<-glm(NumberOfHostplantApparency~as.factor(Host_plant_apparency) * RockStructures, data=speciesTraits.long1, family="poisson")
# m2<-glm(NumberOfHostplantApparency~as.factor(Host_plant_apparency) + RockStructures, data=speciesTraits.long1, family="poisson")
# anova(m1,m2, test="Chi")

# #### Range Size ####
# speciesTraits.long1<-ddply(speciesTraits.long, c("Habitats"            ,         "Altitude"              ,      
#                                                  "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                  "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                  "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                  "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                  "RelativeDryness"             , "RelativeDisturbance" , "Range_size"), summarise,
#                            NumberOfRangeSize=sum(Range_size),
#                            SumRangeSize = sum(Number))

# speciesTraits.long1$NumberOfRangeSize<-speciesTraits.long1$NumberOfRangeSize/speciesTraits.long1$Range_size
# speciesTraits.long1$Range_size<-as.factor(speciesTraits.long1$Range_size)

# ggplot(speciesTraits.long1,aes(Altitude,NumberOfRangeSize, group=Range_size))+
#   geom_smooth(method="glm", method.args = list(family = 'poisson'),aes(colour=Range_size),se=FALSE)+
#   geom_point()+
#   theme_bw()
# ggboxplot(data=speciesTraits.long1, x= "Range_size", y = "NumberOfRangeSize")

# m1<-glm(NumberOfRangeSize~as.factor(Range_size) * StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# m2<-glm(NumberOfRangeSize~as.factor(Range_size) + StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# anova(m1,m2, test="Chi")

# #### Overwintering Stage ####
# speciesTraits.long1<-ddply(speciesTraits.long, c("Habitats"            ,         "Altitude"              ,      
#                                                  "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                  "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                  "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                  "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                  "RelativeDryness"             , "RelativeDisturbance" , "Overwintering_stage"), summarise,
#                            NumberOfOverwintering=sum(Overwintering_stage),
#                            SumOverwintering = sum(Number))

# speciesTraits.long1$NumberOfOverwintering<-speciesTraits.long1$NumberOfOverwintering/speciesTraits.long1$Overwintering_stage
# speciesTraits.long1$Overwintering_stage<-as.factor(speciesTraits.long1$Overwintering_stage)

# ggboxplot(data=speciesTraits.long1, x= "Overwintering_stage", y = "NumberOfOverwintering")

# m1<-glm(NumberOfOverwintering~as.factor(Overwintering_stage) * StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# m2<-glm(NumberOfOverwintering~as.factor(Overwintering_stage) + StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# anova(m1,m2, test="Chi")

# #### Ovum Placement ####
# speciesTraits.long1<-ddply(speciesTraits.long, c("Habitats"            ,         "Altitude"              ,      
#                                                  "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                  "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                  "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                  "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                  "RelativeDryness"             , "RelativeDisturbance" , "Ovum_placement"), summarise,
#                            NumberOfOvumPlacement=sum(Ovum_placement),
#                            SumOvumPlacement = sum(Number))

# speciesTraits.long1$NumberOfOvumPlacement<-speciesTraits.long1$NumberOfOvumPlacement/speciesTraits.long1$Ovum_placement
# speciesTraits.long1$Ovum_placement<-as.factor(speciesTraits.long1$Ovum_placement)

# ggboxplot(data=speciesTraits.long1, x= "Ovum_placement", y = "NumberOfOvumPlacement")

# m1<-glm(NumberOfOvumPlacement~as.factor(Ovum_placement) * as.factor(StrongFloweringPlantPresence), data=speciesTraits.long1, family="poisson")
# m2<-glm(NumberOfOvumPlacement~as.factor(Ovum_placement) + as.factor(StrongFloweringPlantPresence), data=speciesTraits.long1, family="poisson")
# anova(m1,m2, test="Chi")



# #### TRAIT ANALYSES - BINOMIAL ####
# #### Gregariousness #### 
# speciesTraits.long1<-na.omit(speciesTraits.long)
# speciesTraits.long1<-ddply(speciesTraits.long, c("Habitats"            ,         "Altitude"              ,      
#                                                  "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                  "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                  "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                  "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                  "RelativeDryness"             , "RelativeDisturbance" , "Gregariousness"), summarise,
#                            NumberGregariousness=sum(Gregariousness))

# speciesTraits.long1$NumberGregariousness<-speciesTraits.long1$NumberGregariousness/speciesTraits.long1$Gregariousness
# speciesTraits.long1$Gregariousness<-as.factor(speciesTraits.long1$Gregariousness)

# ggboxplot(data=speciesTraits.long1, x= "Gregariousness", y = "NumberGregariousness")

# m1<-glm(NumberGregariousness~as.factor(Gregariousness) * StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# m2<-glm(NumberGregariousness~as.factor(Gregariousness) + StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# anova(m1,m2, test="Chi")

# ggplot(speciesTraits.long1,aes(Altitude,NumberGregariousness, group=Gregariousness))+
#   geom_smooth(method="glm", method.args = list(family = 'poisson'),aes(colour=Gregariousness),se=FALSE)+
#   geom_point()+
#   theme_bw()

# #### Myrmecophily #### 
# speciesTraits.long1<-ddply(speciesTraits.long, c("Habitats"            ,         "Altitude"              ,      
#                                                  "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                  "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                  "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                  "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                  "RelativeDryness"             , "RelativeDisturbance" , "Myrmecophily"), summarise,
#                            NumberMyrmecophily=sum(Myrmecophily))

# speciesTraits.long1$NumberMyrmecophily<-speciesTraits.long1$NumberMyrmecophily/speciesTraits.long1$Myrmecophily

# speciesTraits.long1$Myrmecophily<-as.factor(speciesTraits.long1$Myrmecophily)

# ggboxplot(data=speciesTraits.long1, x= "StrongFloweringPlantPresence", y = "NumberMyrmecophily")

# m1<-glm(NumberMyrmecophily~as.factor(Myrmecophily) * StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# m2<-glm(NumberMyrmecophily~as.factor(Myrmecophily) + StrongFloweringPlantPresence, data=speciesTraits.long1, family="poisson")
# anova(m1,m2, test="Chi")

# ggplot(speciesTraits.long1,aes(ShortVegHeight,NumberMyrmecophily, group=Myrmecophily))+
#   geom_smooth(method="glm", method.args = list(family = 'poisson'),aes(colour=Myrmecophily),se=FALSE)+
#   geom_point()+
#   theme_bw()

# #### Migratory #### 
# speciesTraits.long1<-ddply(speciesTraits.long, c("Habitats"            ,         "Altitude"              ,      
#                                                  "DistanceToWater"              ,"TallVegHeight"        ,        
#                                                  "ShortVegHeight"      ,         "FruitTreePresence"     ,       
#                                                  "StrongFloweringPlantPresence", "OrnamentalPresence"  ,          
#                                                  "RockStructures"      ,         "TreeShadeArea"         ,       
#                                                  "RelativeDryness"             , "RelativeDisturbance" , "Migration"), summarise,
#                            NumberMigration=sum(Migration))

# speciesTraits.long1$NumberMigration<-speciesTraits.long1$NumberMigration/speciesTraits.long1$Migration
# speciesTraits.long1$Migration<-as.factor(speciesTraits.long1$Migration)

# ggboxplot(data=speciesTraits.long1, x= "Migration", y = "Altitude")
# ggplot(speciesTraits.long1,aes(Altitude,NumberMigration, group=Migration))+
#   geom_smooth(method="glm", method.args = list(family = 'poisson'),aes(colour=Migration),se=FALSE)+
#   geom_point()+
#   theme_bw()

# m1<-glm(NumberMigration[Migration==1]~Altitude[Migration==1], data=speciesTraits.long1, family="poisson")
# m2<-glm(NumberMigration[Migration==1]~1, data=speciesTraits.long1, family="poisson")
# anova(m1,m2, test="Chi")

# #Presentation
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(2, 2)))
# define_region <- function(row, col){
#   viewport(layout.pos.row = row, layout.pos.col = col)} 
# print(p1, vp=define_region(2, 1))
# print(p2, vp = define_region(1, 1))
# print(p3, vp = define_region(1, 2))



