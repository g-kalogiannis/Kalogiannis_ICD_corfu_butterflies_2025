library(readxl)
library(vegan)
library(tidyr)
library(plyr)
library(dplyr)
library(wrapr)
library(FSA)
library(car)
library(bbmle)
library(factoextra)
library(Factoshiny)
library(MASS)
library(multcomp)
library(goeveg)
library(ordinal)
library(ape)
library(phytools)
library(viridis)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(ggforce)
library(ggpubr)
library(ggbiplot)
library(gridExtra)
library(grid)
library(ggfortify)
library(ggtree)
library(tibble)
library(betapart)

# Read Data
speciesAbundances<-read.csv("../data/species_abundances.csv")
habitatVariables<-read.csv("../data/habitat_variables.csv")

rownames(speciesAbundances) = speciesAbundances$Habitats
rownames(habitatVariables) = habitatVariables$Habitats
colnames(habitatVariables)<-c("Habitats", "Group", "Alt","DTW","TVH","SVH", "FTP", "SFPP","OP","RS","TSA","RDr", "RDi")

# Habitat Summaries
sp_cols <- setdiff(names(speciesAbundances), "Habitats")
sp_mat <- as.matrix(speciesAbundances[sp_cols])
abund    <- rowSums(sp_mat, na.rm = TRUE)
rich     <- rowSums(sp_mat > 0, na.rm = TRUE)
shannonH <- vegan::diversity(sp_mat, index = "shannon", MARGIN = 1)  # ln-base
evenJ    <- ifelse(rich > 1, shannonH / log(rich), NA_real_)         # J' = H'/ln(S)

summary_by_habitat <- tibble::tibble(
  Habitats   = speciesAbundances$Habitats,
  Abundance  = abund,
  Richness   = rich,
  Diversity  = shannonH,
  Evenness  = evenJ
) %>%
  left_join(habitatVariables, by = "Habitats")
write.csv(summary_by_habitat[,c(1:5)], '../results/habitat_summaries.csv', quote=FALSE, row.names = FALSE)

speciesAbundances$Habitats<-NULL
habitatVariables$Habitats<-NULL

# Logit transform percentages
logit01 <- function(x, eps = 1e-6) {
  p <- x / 100
  p <- pmin(pmax(p, eps), 1 - eps)
  qlogis(p)
}
cols <- c("TSA", "RDr", "RDi")
habitatVariables[cols] <- lapply(habitatVariables[cols], logit01)
habitatVariables$Group<-factor(habitatVariables$Group)

#### Abundance-based Bray-Curtis Beta-Diversity ####
betaAbundResults <- round(beta.pair.abund(speciesAbundances, index.family = "bray")$beta.bray, 2)
betaAbundResults <- as.matrix(betaAbundResults)
write.csv(betaAbundResults, "../results/abundance_bray_beta_diversity.csv", quote=FALSE)

#### Log-transformed Abundance Analysis #### 
speciesAbundancesLog<-log(speciesAbundances+1)
decorana(speciesAbundancesLog)

rda.1<- rda(speciesAbundancesLog ~ Alt + TVH + SVH + FTP + SFPP + OP + RS + TSA, data=habitatVariables) # Maximal Model
rda.1<- ordistep(rda.1,direction="backward", Pout = 0.05, permutations = how(nperm = 9000)) # Perform simplification

summary(rda.1)
anova.cca(rda.1, step=999, by = "axis")
round(100*(summary(rda.1)$cont$importance[2, 1:2]), 2)

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
  geom_point(data = site_df, aes(x = RDA1, y = RDA2), color = "#5f136e", shape = 16, size = 1.5, alpha = 1) +
  geom_text_repel(data = site_df, aes(x = RDA1, y = RDA2, label = Habitat), fontface = "italic", size = 3, color = "#5f136e") +
  geom_point(data = species_df, aes(x = RDA1, y = RDA2), color = "#e65d2f", shape = 17, size = 1.5, alpha = 0.3) +
  geom_segment(data = bp_df, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = bp_df, aes(x = RDA1, y = RDA2, label = Label), color = "black", size = 3.5, family = "arial") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = NULL, x = paste0("RDA1 (", round(100*(summary(rda.1)$cont$importance[2, 1:2]), 2)[[1]], ")"), y = paste0("RDA2 (", round(100*(summary(rda.1)$cont$importance[2, 1:2]), 2)[[2]], ")"))
plot_rda_species <- ggplot() +
  geom_point(data = site_df, aes(x = RDA1, y = RDA2), color = "#5f136e", shape = 16, size = 1.5, alpha = 0.3) +
  geom_point(data = species_df, aes(x = RDA1, y = RDA2), color = "#e65d2f", shape = 17, size = 1.5, alpha = 1) +
  geom_text(data = species_df, aes(x = RDA1, y = RDA2 + 0.07, label = Label), family = "arial", fontface = "italic", size = 3.5, color = "#e65d2f", alpha = 1) +
  geom_segment(data = bp_df, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text_repel(data = bp_df, aes(x = RDA1, y = RDA2, label = Label), color = "black", size = 3.5, family = "arial") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  theme(text = element_text(family = "arial", size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = NULL, x = paste0("RDA1 (", round(100*(summary(rda.1)$cont$importance[2, 1:2]), 2)[[1]], ")"), y = paste0("RDA2 (", round(100*(summary(rda.1)$cont$importance[2, 1:2]), 2)[[2]], ")"))

plot_rda = ggarrange(plot_rda_habitats+ rremove("xlab"), plot_rda_species + rremove("xylab"), align = "hv", common.legend = TRUE, legend = 'right')
annotate_figure(plot_rda, bottom = text_grob(paste0("RDA1 (", round(100*(summary(rda.1)$cont$importance[2, 1:2]), 2)[[1]], ")")))


#### Significant Variable Analysis ####
vars     <- c("Richness", "Abundance", "Shannon_H", "Shannon_J")
col_pts  <- "grey30"
col_line <- "#5f136e"
box_cols <- c("#d6e4ff", "#ffd6d6")   # SFPP=0,1

plot_lm_alt <- function(df, yvar, ylim = NULL, show_x_axis = TRUE, show_y_axis = TRUE) {
  ok <- is.finite(df$Alt) & is.finite(df[[yvar]])
  d  <- df[ok, , drop = FALSE]
  x  <- d$Alt; y <- d[[yvar]]
  if (is.null(ylim)) ylim <- range(y, finite = TRUE)
  plot(x, y, pch = 16, col = col_pts, xlab = if (show_x_axis) "Altitude" else "", ylab = if (show_y_axis) yvar else "", xaxt = if (show_x_axis) "s" else "n", yaxt = if (show_y_axis) "s" else "n", ylim = ylim)
  m <- lm(y ~ x)
  newx <- data.frame(x = seq(min(x), max(x), length.out = 100))
  pr   <- predict(m, newdata = newx, interval = "confidence")
  lines(newx$x, pr[, "fit"], col = col_line, lwd = 2)
  lines(newx$x, pr[, "lwr"], col = col_line, lty = 2)
  lines(newx$x, pr[, "upr"], col = col_line, lty = 2)
  s <- summary(m)
  legend("topleft", bty = "n", legend = paste0("R² = ", round(s$r.squared, 3), "\n", "p = ", round(signif(coef(s)[2, 4], 3), 3)))
}
plot_wilcox_sfpp <- function(df, yvar, ylim = NULL, show_x_axis = TRUE, show_y_axis = TRUE) {
  ok <- is.finite(df$SFPP) & is.finite(df[[yvar]])
  d  <- df[ok, , drop = FALSE]
  stopifnot(all(sort(unique(d$SFPP)) %in% c(0,1)))
  d$SFPP <- factor(d$SFPP, levels = c(0,1), labels = c("0","1"))
  y <- d[[yvar]]
  if (is.null(ylim)) ylim <- range(y, finite = TRUE)
  boxplot(y ~ SFPP, data = d, col = box_cols, outline = FALSE, axes = FALSE, xlab = if (show_x_axis) "SFPP" else "", ylab = if (show_y_axis) yvar else "", ylim = ylim)
  stripchart(y ~ SFPP, data = d, vertical = TRUE, method = "jitter",pch = 16, col = col_pts, add = TRUE)
  if (show_x_axis) axis(1, at = 1:2, labels = levels(d$SFPP))
  if (show_y_axis) axis(2)
  box(bty = "o")
  wt  <- wilcox.test(y ~ SFPP, data = d, exact = FALSE, conf.int = TRUE)
  eff <- if (!is.null(wt$estimate)) paste0(", Δ̃ = ", round(unname(wt$estimate), 2)) else ""
  legend("topleft", bty = "n", legend = c(paste0("Wilcoxon: W = ", wt$statistic, eff), paste0("p = ", round(signif(wt$p.value, 3), 3))))
}
vars   <- c("Richness", "Abundance", "Shannon_H", "Shannon_J")
gap_w <- -0.05   # was 0.06
lay <- matrix(c( 1, 0,  2,
                 3, 0,  4,
                 5, 0,  6,
                 7, 0,  8), nrow = 4, byrow = TRUE)
op <- par(no.readonly = TRUE)
layout(lay, widths = c(1, gap_w, 1), heights = rep(1, 4))
left_margin   <- 3.0   # was 4
right_margin  <- 0.6   # was 1
top_margin    <- 0.01   # was 1.5
bottom_margin <- 2.5   # was 3 (still enough for bottom x labels)
par(oma = c(0,0,0,0), mgp = c(1.7, 0.45, 0), tcl = -0.2)

for (i in seq_along(vars)) {
  v <- vars[i]
  ylim_row  <- range(summary_by_habitat[[v]], finite = TRUE)
  is_bottom <- (i == length(vars))
  par(mar = c(bottom_margin, left_margin, top_margin, right_margin))
  plot_lm_alt(summary_by_habitat, v, ylim = ylim_row,
              show_x_axis = is_bottom, show_y_axis = TRUE)
  par(mar = c(bottom_margin, left_margin, top_margin, right_margin))
  plot_wilcox_sfpp(summary_by_habitat, v, ylim = ylim_row,
                   show_x_axis = is_bottom, show_y_axis = FALSE)
}
par(mfrow=c(1,1))


#### Phylogenetics ####
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

# Pagel's Lambda
trait_vector <- colSums(speciesAbundances > 0)[tree_subset$tip.label]
lambda_result <- phylosig(tree_subset, trait_vector, method = "lambda", test = TRUE)
print(lambda_result)

# FastAnc Tree
non_zero_counts <- colSums(speciesAbundances > 0)
tip_data <- data.frame(
  label = names(non_zero_counts),
  total_abundance = non_zero_counts
)
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
