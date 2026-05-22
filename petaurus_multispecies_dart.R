# Script for Abhishek by Eilish 
# three species:  P. breviceps and P. notatus
# maybe some hybrids
# 

# Map of the samples locations



# Relationships between the three species
# pca
# tsne
# splitstree

# Hybrids 

# Gene flow 

# Diversity stats

# devtools::install_github("rhysnewell/fast-nnt", subdir = "fastnntr")
# remotes::install_github("https://github.com/eilishmcmaster/fastDiversity")
# BiocManager::install("ggtree")
# remotes::install_github("KlausVigo/tanggle")
# remotes::install_github("jokergoo/ComplexHeatmap")
# remotes::install_github('https://github.com/jasongbragg/RRtools')

library(fastnntr) # devtools::install_github("rhysnewell/fast-nnt", subdir = "fastnntr")
library(fastDiversity) #remotes::install_github("https://github.com/eilishmcmaster/fastDiversity")
library(tanggle) # git install
library(adegenet)
library(circlize)
library(ComplexHeatmap) #git install? 
library(data.table)
library(dplyr)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggtree)
library(igraph)
library(openxlsx)
library(RColorBrewer)
library(RRtools) # remotes::install_github('https://github.com/jasongbragg/RRtools')
library(scales)
library(SNPRelate)
library(stringr)
library(tidyr)
library(vegan)

source('https://github.com/eilishmcmaster/SoS_functions/blob/33bd7065baa91b7e2a1800b6481d63573fb38d28/dart2svdquartets.r?raw=TRUE')
devtools::source_url("https://github.com/eilishmcmaster/SoS_functions/blob/main/sos_functions.R?raw=TRUE")

topskip   <- 6
nmetavar  <- 18
basedir <- "" #main directory
species <- "petaurus" #species name
dataset <- "DPetas26-11610" #dart order
missingness <- 0.3

setwd('/Users/abhishek/honours')

d1        <- new.read.dart.xls.onerow(basedir,species,dataset,topskip, nmetavar, euchits=FALSE)
hist(colMeans(is.na(d1$gt)))  ### histogram of missingness of each locus before filtering

d1$sample_names[which(!d1$sample_names %in% meta$sample_names)] # samples in dart not in meta

meta      <- read.meta.data.full.analyses.df(d1, basedir, species, dataset)

d2        <- dart.meta.data.merge(d1, meta) 

# # example make dms with only one species
# d2 <- remove.by.list(d2, d2$sample_names[which(d2$meta$analyses[,'sp']=='breviceps')])
# length(d2$sample_names)

# up to 80% can be missing
d3        <- remove.poor.quality.snps(d2, min_repro=0.96, max_missing=0.8) %>% remove.fixed.snps()

d4        <- sample.one.snp.per.locus.random(d3, seed=12345) 


hist(colMeans(is.na(d1$gt)))  ### histogram of missingness of each locus before filtering

d5 <- remove.by.missingness(d4, 0.8) # remove samples with high missingness
length(d4$sample_names)
length(d5$sample_names)
# odd_samples %in% d5$sample_names

dms <- d5

length(dms$locus_names) # how many polymorphic loci remain
length(dms$sample_names) # how many samples remain

dms_maf5 <- remove.by.maf(dms, 0.05) # make minor allele frequency filtered version 

m2 <- dms$meta$analyses %>% as.data.frame() #get metadata 
m2$lat <- as.numeric(m2$lat)
m2$long <- as.numeric(m2$long)

# Define colors and shapes for species (instead of sites)
species_levels <- factor(m2$sp)  # Use species column 
sp_colours <- named_list_maker(levels(species_levels), palette = 'Set3',4)
sp_shapes <- setNames(1:length(levels(species_levels)), levels(species_levels))  # Assign shapes

# Define colors and shapes for sites)
site_levels <- factor(m2$site)  # Use site column 
site_colours <- named_list_maker(levels(site_levels), palette = 'Spectral',11)
site_shapes <- setNames(1:length(levels(site_levels)), levels(site_levels))  # Assign shapes
 
#### RAW GENOTYPE HEATMAP ####
submat <- dms$gt[,sample(1:ncol(dms$gt), 1000, replace = FALSE)] %>% as.matrix
submat[is.na(submat)] <- -9

row_ann <- rowAnnotation(Group = m2$sp[match(rownames(submat), m2$sample)])

# Define a set of discrete color values
discrete_colors <- c("-9" = "white", "0" = "lightblue", "1" = "coral4", "2" = "orange")

# Apply the discrete color scheme to the heatmap
gt <- Heatmap(submat, 
              col = discrete_colors,
              use_raster = FALSE,
              right_annotation = row_ann,
              cluster_columns = TRUE, cluster_rows = TRUE,
              column_names_gp = gpar(fontsize = 0),
              row_names_gp = gpar(fontsize = 6))

# Set up the PNG device
pdf("petaurus/outputs/petaurus_gt.pdf", width = 25, height = 20)

# Draw the plot
draw(gt)

# Turn off the PNG device
dev.off()

#### KINSHIP ####
kin <- individual_kinship_by_pop(dms, basedir, species, dataset, dms$meta$analyses[,"sp"],
                                 maf=0.05, mis=0.1, as_bigmat=TRUE)

adj <- kin > 1/2^(3/2) # kinship for clones 
diag(adj) <- FALSE
g <- graph_from_adjacency_matrix(adj, mode = "undirected")
components <- components(g)
clones <- as.data.frame(cbind(genet = as.vector(components$membership), sample = as.vector(names(components$membership))))

##### GET CLONES ####

clones_out <- merge(clones, m2[,c("sample","lat","long","site", "sp")], by="sample", all.x = TRUE) #add some metadata

clones_out <- clones_out %>%
  arrange(site, genet)
clones_out$genet_original <- clones_out$genet
old_genet_unique <- unique(clones_out$genet)
new_genet_unique <- 1:length(old_genet_unique)
clones_out$genet <- new_genet_unique[match(clones_out$genet,old_genet_unique)]

write.xlsx(clones_out, "petaurus/outputs/clones_out.xlsx", rowNames=FALSE, colNames=TRUE)

# Remove duplicate clones from working data
#close clones with highest quality data (least missingess)
missingness_gt <- as.data.frame(rowSums(is.na(dms[["gt"]]))/ncol(dms[["gt"]])) #get missingness of all samples
colnames(missingness_gt) <- "missingness"
clone_missing <- merge(clones_out, missingness_gt, by.x="sample", by.y=0) # merge the clone data
clone_missing <- clone_missing[order(clone_missing$missingness),] #order by missingness low to high

unique_genets <- distinct(clone_missing, genet, .keep_all=TRUE) #keeps the top result for each genet (lowest missingness)

#make a list removing duplicate clones
non_clones <- dms$sample_names[which(!dms$sample_names %in% clones_out$sample)]
approved_clones <- unique_genets$sample

clones_dereplicated <- c(non_clones, approved_clones) # list of sample to keep

# remove clones from dms
dms_no_clones <- remove.by.list(dms, clones_dereplicated)

##### KINSHIP HEATMAP ####

col_fun <- colorRamp2(c(0,0.45,0.5), c("grey95", "red",'black'))

kinship_df <- kin %>% as.data.frame()
kinship_df$sample <- rownames(kinship_df)

hm_sites2 <- merge(kinship_df, m2[,c("sample","site", "lat", "long","sp")],
                   by="sample", all.x=TRUE, all.y=FALSE)
hm_sites2 <- hm_sites2[match(rownames(kinship_df),hm_sites2$sample),]
rownames(hm_sites2) <- hm_sites2[,"sample"]

hm_sites <- hm_sites2
hm_sites2[,"sample"] <- NULL

# site_ann <- HeatmapAnnotation(Site = hm_sites2$site,
#                               col=list(Site=site_colours),
#                               annotation_name_gp = gpar(fontsize = 10),
#                               annotation_name_side="right")

sp_ann <- HeatmapAnnotation(Sp = hm_sites2$sp, col = list(Sp = sp_colours),
                            annotation_name_gp = gpar(fontsize = 10),
                            annotation_name_side="right")

# site_ann2 <- rowAnnotation(Site = hm_sites2$site,
#                            col=list(Site=site_colours),
#                            annotation_name_gp = gpar(fontsize = 0))


sp_ann2 <- rowAnnotation(Sp = hm_sites2$sp, col = list(Sp = sp_colours),
                         annotation_name_gp = gpar(fontsize = 0))


hma <- Heatmap( as.matrix(hm_sites2[ , c(1:(nrow(hm_sites2)))]),
                col=col_fun,
                bottom_annotation = sp_ann, #c(site_ann, sp_ann),
                right_annotation = sp_ann2, #c(site_ann2, sp_ann2),
                name = "Kin", #title of legend
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_names_max_width = unit(15, "cm"),
                border_gp = gpar(col = "black", lty = 1)
                # column_order=order(hm_sites2$sp),
                # row_order=order(hm_sites2$sp)
)

draw(hma, merge_legend = TRUE)

# Set up the PNG device
png('petaurus/outputs/petaurus_kin.png', width = 30, height = 30, units = 'cm', res = 300)

# Draw the plot
draw(hma, merge_legend = TRUE)

# Turn off the PNG device
dev.off()


#### PCA ####
##### PCA OF ALL SAMPLES ####

# dms_maf1 <- remove.by.maf(dms_no_clones, 0.01) # if there are a lot of loci, minor allele filter is a good way to reduce them while keeping information esp. in PCA

dms_pca <- remove.poor.quality.snps(dms_no_clones, min_repro=0.96, max_missing=0.3) %>% 
  remove.by.maf(., 0.05) %>%
  remove.loci.randomly(.,2000)

# if there are a lot of loci, minor allele filter is a good way to reduce them while keeping information esp. in PCA
# duplicate clones should be excluded from PCA
length(dms_pca$locus_names)

# Convert filtered genotypes to genlight object
gen_d5 <- new("genlight", dms_pca$gt)

# Perform PCA
gen_pca <- glPca(gen_d5, parallel = TRUE, nf = 7) #can take a while


# Extract PCA scores
g_pca_df <- gen_pca[["scores"]]

# Merge with metadata
g_pca_df2 <- merge(g_pca_df, d2$meta$analyses, 
                   by.x = "row.names", by.y = "sample", 
                   all.x = FALSE, all.y = FALSE)

# Create axis labels with % variance explained
pcnames <- paste0(colnames(g_pca_df), " (", 
                  round(gen_pca[["eig"]][1:6] / sum(gen_pca[["eig"]]) * 100, 2), 
                  "%)")


# Common PCA plot elements
base_pca_plot <- function(data, x, y, xlabel, ylabel) {
  ggplot(data, aes_string(x = x, y = y, colour = "sp", shape = "sp")) +
    geom_vline(xintercept = 0, alpha = 0.2) +
    geom_hline(yintercept = 0, alpha = 0.2) +
    geom_point(size = 1.5) +
    xlab(xlabel) + ylab(ylabel) +
    theme_few() +
    theme(
      legend.key.size = unit(0, 'lines'),
      legend.position = "right",
      legend.text = element_text(face = "italic"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    ) +
    scale_colour_manual(values = sp_colours) +
    scale_shape_manual(values = sp_shapes)
}

# PCA Plot 1: PC1 vs PC2 (with convex hull)
pca_plot1 <- base_pca_plot(g_pca_df2, "PC1", "PC2", pcnames[1], pcnames[2]) +
  guides(
    colour = guide_legend(title.position = "left", ncol = 3),
    shape = guide_legend(title.position = "left", ncol = 3)
  ) 
# PCA Plot 2: PC3 vs PC4
pca_plot2 <- base_pca_plot(g_pca_df2, "PC3", "PC4", pcnames[3], pcnames[4]) +
  guides(colour = guide_legend(title.position = "top"))

# PCA Plot 3: PC5 vs PC6
pca_plot3 <- base_pca_plot(g_pca_df2, "PC5", "PC6", pcnames[5], pcnames[6]) +
  guides(colour = guide_legend(title.position = "top"))

# Arrange all three PCA plots
all3_pca_plots <- ggarrange(pca_plot1, pca_plot2, pca_plot3,
                            labels = "AUTO",
                            font.label = list(face = "plain"),
                            common.legend = TRUE, 
                            legend = "bottom", ncol = 3)

# Display and save combined plot
print(all3_pca_plots)

ggsave("petaurus/outputs/petaurus_pca.png",
       all3_pca_plots, width = 22, height = 10,
       units = "cm", dpi = 300, bg = 'white')


#### DAPC ####
# Find optimal clusters for up to k = 20
grp <- find.clusters(gen_d5, max.n.clust = 20) # retained 90 PCs
                                               # chose k=3 clusters 

# table(pop(gen_d5), grp$grp)

# Run the DAPC 
dapc1 <- dapc(gen_d5, grp$grp) #retained 150 PCs, 2 discriminant functions 
dapc2 <- dapc(gen_d5, grp$grp) #retained 100 PCs, 2 discriminant functions
dapc3 <- dapc(gen_d5, grp$grp) #retained 200 PCs, 2 discriminant functions 
dapc4 <- dapc(gen_d5, grp$grp) #retained 250 PCs, 2 discriminant functions

dapc1
dapc2

# Scatterplot
scatter(dapc1, posi.da = "bottomleft")

# Trying to label by metadata species 
testingmeta <- read_xlsx('/Users/abhishek/honours/petaurus/meta/petaurus_DPetas26-11610_meta.xlsx')

# Reorder metadata to match genlight individuals
testing2 <- testingmeta[match(indNames(gen_d5), testingmeta$sample), ]

# Verify alignment
all(indNames(gen_d5) == testing2$sample)

# Create grouping factor
pop_data <- factor(testing2$sp)

# Colours
custom_colours <- adegenet::seasun(nlevels(pop_data))

# Plot
scatter(dapc1, grp = pop_data, posi.da = "topright")

#### AMOVA ####
library(poppr)

# Reorder metadata to match genlight
testingmeta <- testingmeta[match(indNames(gen_d5), testingmeta$sample), ]
all(indNames(gen_d5) == testingmeta$sample)

# Assign populations to genlight object for hierarchy 
strata(gen_d5) <- data.frame(
  popspecies = as.factor(testingmeta$sp), 
  popbioregion = as.factor(testingmeta$bioregion)
) 

# Run AMOVA
amova_result <- poppr.amova(gen_d5, ~popspecies/popbioregion)
amova_result

# Test significance
# amova_sig <- randtest(amova_result, nrepet = 999) #taking rly long 
amova_sig_quick <- randtest(amova_result, nrepet = 99)
amova_sig_quick

#### DIVERSITY STATS #####

stats_0.3 <- fastDiversity::faststats(gt = dms_no_clones$gt, # diversity stats
                                      genetic_group_variable = dms_no_clones$meta$analyses[,'sp'],
                                      site_variable = dms_no_clones$meta$analyses[,'site'], 
                                      # run_HWE_test = TRUE, 
                                      # get_CI_loci_resampling = TRUE,
                                      # loci_resample_n = 1000, 
                                      # boots = 1000, 
                                      max_missingness = 0.3, maf=0.05)

write.xlsx(stats_0.3, paste0(species,'/outputs/petaurus_stats.xlsx'))

# the individual observed heterozygosity
m2$ind_ho <- fastDiversity::individual_Ho(dms$gt, 
                                          genetic_group_variable = dms$meta$analyses[,'sp'],
                                          maf=0.05, max_missingness = 0.3)

ind_ho <- fastDiversity::individual_Ho(dms$gt, 
                                          genetic_group_variable = dms$meta$analyses[,'sp'],
                                          maf=0.05, max_missingness = 0.3)
ind_ho[odd_samples]

ind_miss <- rowMeans(is.na(dms$gt))  # the missing proportion of individuals
m2$miss <- ind_miss[m2$sample]

library(ggplot2)

ggplot(m2, aes(x=ind_ho, fill=sp))+
  geom_histogram(bins=50)+
  facet_wrap(~sp)+
  theme_few()


ggplot(m2, aes(x=miss, y=ind_ho, color=sp, label=sample))+
  geom_point()+
  theme_few()+
  ggrepel::geom_text_repel()


 #### SPLITSTREE ####
gt_dist <- dist(dms$gt, method='euclidean') %>% as.matrix #distance matrix

Nnet <- fastnntr::run_neighbornet_networkx(gt_dist) #phylogenetic network

# Prepare plotting data
x <- data.frame(x = Nnet$.plot$vertices[,1], 
                y = Nnet$.plot$vertices[,2], 
                sample = rep(NA, nrow(Nnet$.plot$vertices)))

x[Nnet$translate$node, "sample"] <- Nnet$translate$label

# Merge with metadata
x <- merge(x, m2, by.x = "sample", by.y = "sample", all.x = TRUE, all.y = FALSE)
x_tips <- x[!is.na(x$sample) & !is.na(x$site), ]

# Adjust axis limits
net_x_axis <- max(x_tips$x) - min(x_tips$x)
net_y_axis <- max(x_tips$y) - min(x_tips$y)

# Update tip labels to show site only
# Nnet$translate$label <- gsub("\\\\n", "\n", x_tips[match(Nnet$tip.label, x_tips$sample), "pop_label"])

hull <- x %>% group_by(sp) %>%
  slice(chull(x, y)) %>%
  filter(!is.na(sp))




# Create the plot with species-based colors and shapes
splitstree_plot <- ggplot(Nnet, aes(x = x, y = y)) +
  geom_shape(data = hull, alpha = 0.2, expand = 0.01, radius = 0.01,
             aes(fill = sp), color = 'transparent') +
  geom_splitnet(layout = "slanted", size = 0.2) +
  geom_point(data = x_tips, aes(x = x, y = y, colour = sp_original, shape = sp_original), size = 1) +
  scale_fill_manual(values = sp_colours, na.translate = FALSE, guide = "none") +    
  scale_colour_manual(values = sp_colours, na.translate = FALSE) +    
  scale_shape_manual(values = sp_shapes) +            
  geom_tiplab2(size = 2, hjust = -0.15, fontface = "italic", lineheight = 0.8) + #angle = 0, 
  theme_void() +
  expand_limits(x = c(min(x_tips$x) - 0.2 * net_x_axis, max(x_tips$x) + 0.3 * net_x_axis),
                y = c(min(x_tips$y) - 0.1 * net_y_axis, max(x_tips$y) + 0.1 * net_y_axis)) +
  theme(
    legend.text = element_text(face = "plain"),
    legend.key.size = unit(0, 'lines')
  ) +
  coord_fixed()

splitstree_plot

# Save the plot
ggsave(paste0(species,'/outputs/petaurus_splitstree_original.pdf'), 
       splitstree_plot, width=30, height=30, units='cm', dpi=300)

# Display the plot
pca_splits_plot <- ggarrange(all3_pca_plots, splitstree_plot, nrow=2, labels=c("","D"), font.label=list(face="plain"))

# Save the plot
ggsave(paste0(species,'/outputs/petaurus_pca_splitstree.png'), 
       pca_splits_plot, width=22, height=18, units='cm', dpi=300, bg='white')



#### splitstree missing ####


x_tips$miss <- ind_miss[x_tips$sample]

# 
odd_samples <- c('WGM244_203', 'WGM244_120', 'WGM244_121', 'WGM244_186')

ind_miss[odd_samples]
summary(ind_miss)
hist(ind_miss)

# Create the plot with species-based colors and shapes
splitstree_plot_missing <- ggplot(Nnet, aes(x = x, y = y)) +
  geom_shape(data = hull, alpha = 0.2, expand = 0.01, radius = 0.01,
             aes(fill = sp), color = 'transparent') +
  geom_splitnet(layout = "slanted", size = 0.2) +
  geom_point(data = x_tips, aes(x = x, y = y, colour = miss), size = 1) +
  # scale_fill_manual(values = sp_colours, na.translate = FALSE, guide = "none") +    
  # scale_colour_manual(values = sp_colours, na.translate = FALSE) +    
  # scale_shape_manual(values = sp_shapes) +            
  geom_tiplab2(size = 2, hjust = -0.15, fontface = "italic", lineheight = 0.8) + #angle = 0, 
  theme_void() +
  expand_limits(x = c(min(x_tips$x) - 0.2 * net_x_axis, max(x_tips$x) + 0.3 * net_x_axis),
                y = c(min(x_tips$y) - 0.1 * net_y_axis, max(x_tips$y) + 0.1 * net_y_axis)) +
  theme(
    legend.text = element_text(face = "plain"),
    legend.key.size = unit(2, 'lines')
  ) +
  coord_fixed()


splitstree_plot_missing

# Save the plot
ggsave(paste0(species,'/outputs/petaurus_splitstree_missing.pdf'), 
       splitstree_plot_missing, width=30, height=30, units='cm', dpi=300)




#### FST ANALYSIS ####
# Remove samples with NA in site column (not species, since all are Grevillea petaurus)
# Filter sites with n < 2

dms_no_clones <- dms #TEMPORARY!!! until you do kinship

site_df <- dms_no_clones$meta$analyses[,'site'] %>% as.data.frame()
site_df[site_df==0] <- NA

sitepop_freq <- as.data.frame(table(site_df))
not_n1_sites <- as.vector(sitepop_freq[sitepop_freq$Freq < 3, 1])
keep_samples <- dms_no_clones$sample_names[which(dms_no_clones$meta$analyses[,'site'] %in% sitepop_freq[sitepop_freq$Freq >= 3, 1])]

fst_dms <- remove.by.list(dms_no_clones, keep_samples)

length(fst_dms$sample_names)
table(fst_dms$meta$site)


# Check filtered data
print("After filtering sites with n < 2:")
print(length(fst_dms$sample_names))
print(table(fst_dms$meta$analyses[,'site']))

# Calculate FST and geographic distance
gds_file <- dart2gds(fst_dms, basedir, species, dataset)
pFst <- population.pw.Fst(fst_dms, fst_dms$meta$analyses[,'site'], basedir, species, dataset,
                          maf_val=0.05, miss_val=0.3)

pS <- population.pw.spatial.dist(fst_dms, fst_dms$meta$analyses[,'site'])

# Make self comparisons NA
diag(pFst$Fst) <- NA
diag(pS$S) <- NA

##### FST MANTEL ####

meta_agg_df <- data.frame(lat = fst_dms$meta$lat, 
                          long = fst_dms$meta$long,
                          site = fst_dms$meta$analyses[,"site"],
                          sp = fst_dms$meta$analyses[,"sp"])
meta_agg <- meta_agg_df %>%
  group_by(sp, site) %>%
  summarise(lat = mean(lat, na.rm = TRUE),
            long = mean(long, na.rm = TRUE),
            .groups = 'drop')


# Mantel test
man <- mantel(xdis = pS$S, ydis = pFst$Fst, permutations = 10000, na.rm = TRUE)
print(man)

# Mantel plot preparation
Fst_sig <- as.data.frame(as.table(pFst$Fst))
colnames(Fst_sig) <- c("sp1", "sp2", "Fst")
dist_long <- as.data.frame(as.table(pS$S))
Fst_sig$Geo_dist <- dist_long$Freq
Fst_sig$Geo_dist2 <- Fst_sig$Geo_dist / 1000  # Convert to km


# colour by pairwise groupings
# Mantel plot
fstp1 <- ggplot(Fst_sig, aes(x = Geo_dist2, y = Fst)) +
  geom_point(size = 1, alpha = 0.3) +
  labs(x = "Distance (km)", y = "FST", colour = "Comparison") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom") +
  ggtitle("Mantel Test: Fst vs Geographic Distance",
          subtitle = paste("r =", round(man$statistic, 3), "p =", sprintf("%.3f", man$signif)))

# Display and save plot
print(fstp1)

ggsave("petaurus/outputs/petaurus_fst_mantel.png", 
       fstp1, width = 15, height = 10, units = "cm", dpi = 600)

##### FST HEATMAP ####
# Print Mantel result
print(paste("Mantel statistic r is", round(man$statistic, 3), ", P =", man$signif))

# Make heatmaps
# Geo distance
geo_d <- pS$S
mat_dist <- geo_d / 1000  # Convert to km

# FST
mat_fst <- pFst$Fst

# Order heatmaps consistently
order_hm <- Heatmap(mat_fst, cluster_rows = TRUE, cluster_columns = TRUE)
od <- colnames(mat_fst)[column_order(order_hm)]
mat_fst <- mat_fst[od, od]
mat_dist <- mat_dist[od, od]

mat_dist <- merge(mat_dist, meta_agg[,c('sp','site')], by.x=0, by.y="site", all.y=FALSE) #add aggregated df to mat2 (fst)
rownames(mat_dist) <- mat_dist$Row.names
mat_dist$Row.names <- NULL
mat_dist <- mat_dist[match(colnames(mat_dist)[1:nrow(mat_dist)],rownames(mat_dist)),]

# Specify heatmap colors
library(circlize)
gene_col <- colorRamp2(c(0, 0.5, 1), c("#8DD3C7", "white", "#FB8072"))
palette <- colorRamp2(c(0, max(geo_d/1000, na.rm = TRUE)), c("white", "#80B1D3"))


# Define species colors (replace named_list_maker if needed)
row_ann <- rowAnnotation(Group = mat_dist$sp,
                         col=list(Group=sp_colours),  show_legend = FALSE,
                         na_col="white",
                         annotation_legend_param = list(labels_gp=gpar(fontface="italic",fontsize=8),
                                                        title_gp=gpar(fontsize=10)),
                         annotation_name_gp = gpar(fontsize = 0),
                         annotation_name_side="top")

col_ann <- HeatmapAnnotation(Group = mat_dist$sp, 
                             col = list(Group = sp_colours),  show_legend = FALSE,
                             annotation_name_gp = gpar(fontsize = 0),
                             annotation_legend_param = list(labels_gp=gpar(fontface="italic", fontsize=8),
                                                            title_gp=gpar(fontsize=10)),
                             annotation_name_side="left",
                             na_col = "white")

# Geo heatmap
geo <- Heatmap(as.matrix(mat_dist[,1:nrow(mat_dist)]), 
               rect_gp = gpar(type = "none"),
               bottom_annotation = col_ann,
               width = nrow(mat_dist) * unit(6, "mm"),
               height = nrow(mat_dist) * unit(6, "mm"),
               col = palette,
               na_col = "white",
               column_names_gp = gpar(fontsize = 8),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               name = "Distance (km)",
               heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                           labels_gp = gpar(fontsize = 6)),
               cell_fun = function(j, i, x, y, w, h, fill) {
                 if (i > j) {
                   grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                   grid.text(sprintf("%.0f", mat_dist[i, j]), x, y, gp = gpar(fontsize = 6))
                 }
               })

geo

# FST heatmap
gene <- Heatmap(as.matrix(mat_fst), 
                right_annotation = c(row_ann),
                rect_gp = gpar(type = "none"),
                width = nrow(mat_fst) * unit(6, "mm"),
                height = nrow(mat_fst) * unit(6, "mm"),
                col = gene_col,
                na_col = "grey",
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 0),
                border_gp = gpar(col = "black", lty = 1),
                name = "FST",
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                            labels_gp = gpar(fontsize = 5)),
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if (i <= j) {
                    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                    grid.text(sprintf("%.2f", mat_fst[i, j]), x, y, gp = gpar(fontsize = 6))
                  }
                })

gene
# Draw combined heatmap
gene_width <- nrow(mat_fst) * unit(6, "mm")
draw(geo + gene, ht_gap = -gene_width, merge_legend=TRUE)

# Set the file name and parameters
filename <- "petaurus/outputs/petaurus_fst.png"
width <- 24
height <- 20
dpi <- 600
units <- "cm"

# Set up the PNG device
png(filename, width = width, height = height, units = units, res = dpi)

# Draw the plot
draw(geo + gene, ht_gap = -gene_width)

# Turn off the PNG device
dev.off()

#### MAKE ADMIXTURE INPUT FILES ####
igds_file <- dart2gds(dms, '', species, dataset)

igds <- snpgdsOpen(igds_file)

snpgdsGDS2BED(igds, bed.fn=paste0('',species, "_", dataset,"_bed"))


fam <- data.frame(fam = 1:length(dms$sample_names), ind = dms$sample_names, 
                  fat = rep(0, length(dms$sample_names)), mot = rep(0, length(dms$sample_names)), 
                  sex = rep(0, length(dms$sample_names)),
                  pheno = rep(-9, length(dms$sample_names)), stringsAsFactors = FALSE)

write.table(fam, 
            paste0('',species, "_", dataset,"_bed.fam"), sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

snpgdsClose(igds)

# #### LINUX COMMAND ADMIXTURE ####
# for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; \
# do admixture/dist/admixture_macosx-1.3.0/admixture --cv=20 petaurus_DPetas26-11610_bed.bed $K -j20 | tee ADMIXTURE_log${K}.out; done
# grep -h CV ADMIXTURE_log*.out > ADMIXTURE_CV_results.tsv

# Plot cross validation error for most optimal k value
library(data.table)
library(ggplot2)
library(stringr)
library(ggthemes)


cv <- fread("/Users/abhishek/honours/petaurus/outputs/ADMIXTURE_CV_results.tsv")

cv$K <- str_extract(cv$V3, "\\d+") |> as.numeric()

cv$CV_error <- cv$V4

head(cv)

cv_plot <- ggplot(cv, aes(x = as.numeric(K), y = CV_error)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = seq(min(cv$K), max(cv$K), by = 1)
  ) +
  labs(
    x = "Number of clusters (K)",
    y = "Cross-validation error",
    title = "ADMIXTURE CV error"
  ) +
  theme_bw()
cv_plot

# ancestry prediction for all individuals 

q_file <- "/Users/abhishek/honours/petaurus/outputs/petaurus_DPetas26-11610_bed.3.Q"
fam_file <- "/Users/abhishek/honours/petaurus/outputs/petaurus_DPetas26-11610_bed.fam"

opt_k <-3   # set K value here

admix_P <- read.table(q_file)
admix_fam <- read.table(fam_file, sep = "\t")

admix_P$sample <- admix_fam$V2

# merge with metadata (meta must exist already)
admix_P <- merge(admix_P, m2[,c('sample','sp')], by = "sample")

# convert to long format for plotting
admix_P2 <- admix_P %>%
  pivot_longer(
    cols = -c("sample","sp"),
    names_to = "population",
    values_to = "Q"
  ) %>% filter(Q>0.00001)
admix_P2$population <- gsub("V","", admix_P2$population) %>% as.numeric


admix_colours_df <- admix_P2 %>%
  group_by(sp, population) %>%
  summarise(mean_Q = mean(Q), .groups = "drop") %>%
  arrange(desc(mean_Q)) %>%
  group_by(sp) %>%
  slice_max(mean_Q, n = 1) %>% filter(mean_Q>0.90)

# admix_colours <- sp_colours[admix_colours_df$sp]
# names(admix_colours) <- admix_colours_df$population

admix_P2$population <- factor(
  admix_P2$population,
  levels = unique(admix_P2$population[order(nchar(as.character(admix_P2$population)))])
)

plot_admix <- ggplot(admix_P2, aes(x = sample, y = Q, fill = factor(population))) +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  theme_few() +
  labs(
    y = "ADMIXTURE Ancestry",
    fill = "Source\npopulation"
  ) +
  facet_grid(
    ~ factor(sp),#, levels = site_order
    scales = "free_x",
    space = "free_x"
  ) +
  theme(
    axis.text.x = element_text(size=5, angle=90, vjust=0.5),
    axis.title.x = element_blank(),
    strip.text.x = element_text(angle=90, face='italic', hjust=0.01, size=10),
    strip.background = element_blank(),
    panel.spacing = unit(0.07, "lines"),
    legend.position = "none",
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.5, "cm")
  ) +
  scale_y_continuous(limits = c(0, 1.001), expand = c(0, 0))
  # scale_fill_manual(values = admix_colours)

plot_admix
ggsave(paste0(species,'/outputs/petaurus_plot_admix.pdf'), 
       plot_admix, width=40, height=30, units='cm', dpi=300)

#### RAXML ################################################################################
dt <- as.data.table(dms_no_clones$gt) #replaced dms with dms_no_clones

# remove "invariant" sites (0,1 or 1,2)
cols_idx <- which(
  dt[, lapply(.SD, function(x) {
    ux <- unique(na.omit(x))
    length(ux) > 0 && (all(ux %in% c(0, 1)) || all(ux %in% c(1, 2)))
  })] |> unlist()
)

dms_raxml <- remove.snps.from.dart.data(dms_no_clones, snps_to_remove = cols_idx, input_as_names = FALSE)

dart2svdquartets(dms_raxml, basedir, species, paste0(dataset,'_raxml'), add_pop=TRUE, pop=dms_raxml$sample_names)
genotype_matrix <- read.nexus.data('petaurus/popgen/raw_SNPFilt_1SNPperClone/svdq/petaurus_DPetas26-11610_raxml.nex')

dnabin <- ape::nexus2DNAbin(genotype_matrix)
ape::write.FASTA(dnabin, file='petaurus/popgen/raw_SNPFilt_1SNPperClone/svdq/petaurus_DPetas26-11610_raxml.fasta')

# cipres
# galaxy

## recreate dart plate ## 
library(purrr)
# read dart csv
df_dart <- read.csv("petaurus/dart_raw/Report_DPetas26-11610_SNP_mapping_2.csv", header = FALSE)
View(df_dart)

# Keep 21st column onward 
df_u <- df_dart[, 21:ncol(df_dart)]

# Extract relevant rows
plate_row   <- as.character(unlist(df_u[3, ]))   # row 3 = plate
letter_row  <- as.character(unlist(df_u[4, ]))   # row 4 = letter
number_row  <- as.numeric(unlist(df_u[5, ]))     # row 5 = number
sample_row <- as.character(unlist(df_u[7, ]))   # row 7 = sample ID

# Build tidy dataframe
dart_tidy <- data.frame(plate = plate_row, letter = letter_row, number = number_row, sample = sample_row)

# Reconstruct plates
plate_table <- dart_tidy |>
  group_split(plate) |>
  setNames(unique(dart_tidy$plate)) |>
  map(~ .x |>
        arrange(letter, number) |>
        pivot_wider(id_cols = letter, names_from = number, values_from = sample))

plate_table["1"]
plate_table["2"]
plate_table["3"]





