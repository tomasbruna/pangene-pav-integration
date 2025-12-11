################################################################################
# 0. Package installation
################################################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
# List of required packages
packages <- c(
  "data.table",
  "reshape2",
  "ggplot2",
  "cowplot",
  "vegan",
  "RColorBrewer",
  "parallel"
)

# Install any missing packages
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

sapply(packages, require, character.only = TRUE)

################################################################################
################################################################################
# 1. Analyze PAV distances wrt pankmer distances in Soybean
################################################################################
################################################################################

################################################################################
# 1.0 ad hoc functions
read_pankmerDist <- function(path){
  tmp <- fread(path)
  
  mat <- data.matrix(tmp[,-1,with = F])
  ns <- colnames(mat)

  rownames(mat) <- tmp$V1
  
  mat <- mat[colnames(mat),]
  return(mat)
}

convert_cb2pav <- function(genome, orthogroupID){
  # remove duplicates (so no CNV)
  x <- data.table(genome = genome, orthogroupID = orthogroupID)
  x <- subset(x, !duplicated(x))
  
  # convert to a matrix
  y <- dcast(x, orthogroupID ~ genome, fun.aggregate = length)
  mat <- as.matrix(y[,-1])
  dis <- dist(t(mat), method = "binary")
  mat <- as.matrix(dis)

  return(mat)
}

get_hclustOrd <- function(mat, method = "average"){
  hck <- hclust(as.dist(mat), method = method)
  return(hck$labels[hck$order])
}

melt_matrix <- function(x, genomeOrd = NULL){
  dt <- data.table(reshape2::melt(x))
  setnames(dt, 1:2, c("genome1", "genome2"))
  dt[,`:=`(genome1 = as.character(genome1), 
           genome2 = as.character(genome2),
           scaled = scale_between(value, 0, 1))]
  if(!is.null(genomeOrd)){
    dt[,`:=`(genomeFac1 = factor(genome1, levels = genomeOrd),
             genomeFac2 = factor(genome2, levels = genomeOrd),
             x = as.integer(factor(genome1, levels = genomeOrd)),
             y = as.integer(factor(genome2, levels = genomeOrd)))]
  }
  return(dt)
}

scale_between <- function(x, min, max, scale1toMean = TRUE){
  if(length(unique(x)) > 1){
    return((x - min(x)) / (max(x) - min(x)) * (max - min) + min)
  }else{
    if(scale1toMean){
      return(mean(c(min, max)))
    }else{
      return(max)
    }
  }
}

################################################################################
# 1.1 parameters and paths
dir.create("R_figures", showWarnings = FALSE)
out_dir <- file.path(getwd(), "R_figures")
# -- working directory
wd <- file.path(getwd())
print(file.path(wd, "Soybean/pankmer/pangrowth_growth.txt"))
# -- soybean
rawSoyFile <- file.path(
  wd, "Soybean/genespace/external/gs/results/combBed.txt")
mdSoyFile <- file.path(
  wd, "Soybean/soybean_md.txt")
pankmerSoyFile <- file.path(
  wd, "Soybean/pankmer/pankmer_dist_matrix.tsv")
igcSoyFile <- file.path(
  wd, "Soybean/genespace/jgi/gs/results/combBed.txt")
kmerSoy <- as.numeric(unlist(fread(
  file.path(wd, "Soybean/pankmer/pangrowth_growth.txt"))))

# -- cotton
rawCotFile <- file.path(
  wd, "Cotton/genespace/external/gs/results/combBed.txt")
pankmerCotFile <- file.path(
  wd, "Cotton/pankmer/pankmer_dist_matrix.tsv")
igcCotFile <- file.path(
  wd, "Cotton/genespace/jgi/gs/results/combBed.txt")
kmerCot <- as.numeric(unlist(fread(
  file.path(wd, "Cotton/pankmer/pangrowth_growth.txt"))))


# -- parameters
minChrSize <- 0
dropFromNames <- ".gnm.*|.genome$"

################################################################################
# 1.2 Load and process soybean data
# -- read metadata
mdSoy <- fread(mdSoyFile)

# -- convert columns so they match other data structures
mdSoy[,id := gsub(".gnm.*|glyma.", "", basename(assFile))]
mdSoy[,id := gsub("-", "_", id, fixed = T)]

# -- get grouping variables
mdSoy[,grp := ifelse(grepl("Liu", name), "Liu2020",
                     ifelse(grepl("Chu,", name), "Chu2021",
                            ifelse(id %in% c("Wm82", "FiskebyIII"), 
                                   "Phytozome", "ref")))]

# -- strip off other columns
mdSoy <- mdSoy[,c("id", "grp")]
mdSoy[,refLab := ifelse(grepl("wm82", tolower(id)), "Wm82",
                        ifelse(grepl("zh13", tolower(id)), "Zh13", "other"))]

# -- get vector of groupings
soyGrp <- mdSoy$grp; names(soyGrp) <- mdSoy$id
refGrp <- mdSoy$refLab; names(refGrp) <- mdSoy$id

# -- read in pankmer matrix and rename to id
pkSoy <- read_pankmerDist(pankmerSoyFile)

# -- read in the combBed files and rename genomes
rawSoyIn <- fread(rawSoyFile)
rawSoyIn[,genome := gsub(".gnm.*", "", genome)]

igcSoyIn <- fread(igcSoyFile)
igcSoyIn[,genome := gsub(".gnm.*", "", genome)]

# -- remove any chromosomes with too few genes (scaffolds, contigs, etc.)
# -- minChrSize is set to 0 by default, adjust for this to take effect
rawSoyIn[,ngenes := .N, by = c("genome", "chr")]
rawSoyIn <- subset(rawSoyIn, ngenes > minChrSize)

igcSoyIn[,ngenes := .N, by = c("genome", "chr")]
igcSoyIn <- subset(igcSoyIn, ngenes > minChrSize)

################################################################################
# 1.3 convert soybean files to PAV distance matrices

# -- igc globHOGs
igcSoy <- with(igcSoyIn, convert_cb2pav(genome = genome, orthogroupID = globHOG))

# -- raw globHOGs
rawSoy <- with(rawSoyIn, convert_cb2pav(genome = genome, orthogroupID = globHOG))

# -- pankmer (needs to be renamed)
colnames(pkSoy) <- rownames(pkSoy) <- 
  gsub("-","_",gsub(".genome$", "", rownames(pkSoy)))

# -- make heatmaps
heatmap(pkSoy, scale = "none")
heatmap(rawSoy, scale = "none")
heatmap(igcSoy, scale = "none")

# -- get hclust vectors
pkSoyOrd <- get_hclustOrd(pkSoy)
rawSoyOrd <- get_hclustOrd(rawSoy)
igcSoyOrd <- get_hclustOrd(igcSoy)

# -- reshape and add metadata
pkSoyTp <- melt_matrix(pkSoy)
rawSoyTp <- melt_matrix(rawSoy)
igcSoyTp <- melt_matrix(igcSoy)

################################################################################
# 1.4 reformat and plot heatmaps

# -- melt matrices, reorder genomes and combine
tpSoy <- rbind(
  data.table(melt_matrix(pkSoy, genomeOrd = pkSoyOrd), type = "pankmer"),
  data.table(melt_matrix(rawSoy, genomeOrd = rawSoyOrd), type = "hogsRaw"),
  data.table(melt_matrix(igcSoy, genomeOrd = igcSoyOrd), type = "hogsIGC"))


# -- add the rug and line colors
tpSoy[,`:=`(refx = refGrp[genome1], refy = refGrp[genome2],
            rugColx = soyGrp[genome1], rugColy = soyGrp[genome2])]

# -- add plot ordering
tpSoy[,type := factor(type, levels = c("pankmer", "hogsRaw", "hogsIGC"))]

# -- make the heatmap plots
tpRugSoy <- subset(tpSoy, genome1 == genome2)
pltSoyHeatmaps <- ggplot()+
  geom_tile(
    data = tpSoy, 
    aes(x = x, 
        y = y, 
        fill = value))+
  facet_grid(.~type)+
  scale_fill_viridis_c(option = "A")+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.key.size = unit(.4, 'cm'),
    legend.title = element_blank(), 
    legend.text=element_text(size=6),
    title = element_text(size = 7),
    strip.text = element_text(size = 6))+
  geom_vline(data = subset(tpSoy, refx == "Wm82"), aes(xintercept = x),
             alpha = .5, linewidth = .25, col = "black")+
  geom_vline(data = subset(tpSoy, refx == "Zh13"), aes(xintercept = x),
             alpha = .5, linewidth = .25, col = "white")+
  geom_hline(data = subset(tpSoy, refy == "Zh13"), aes(yintercept = y),
             alpha = .5, linewidth = .25, col = "white")+
  geom_hline(data = subset(tpSoy, refy == "Wm82"), aes(yintercept = y),
             alpha = .5, linewidth = .25, col = "black")+
  geom_rug(data = tpRugSoy, aes(x = x, col = rugColx), linewidth = 1, outside = T, sides = "bl")+
  coord_fixed(clip = "off")+
  scale_color_manual(values = c("green3", "dodgerblue2", "cyan", "black"), guide = "none")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))

################################################################################
# 1.5 Make the xy plots and do correlations

# -- calculate mantel correlation coefs
vegan::mantel(pkSoy, igcSoy)
vegan::mantel(pkSoy, rawSoy)

# -- reorder genomes so no duplicate
tpSoy[,`:=`(g1 = ifelse(x > y, genome1, genome2),
            g2 = ifelse(x > y, genome2, genome1))]

# -- pull out pankmer data
x <- with(subset(tpSoy, type == "pankmer"), data.table(
  genome1 = g1, 
  genome2 = g2, 
  pkdist = value, 
  pkgrp = ifelse(rugColx == rugColy, rugColx, "diff")))

# -- pull out pav data
y <- with(subset(tpSoy, type != "pankmer"), data.table(
  genome1 = g1, 
  genome2 = g2, 
  pavdist = value, 
  type = type, 
  pavgrp = ifelse(rugColx == rugColy, rugColx, "diff")))

# -- merge
# First merge
corRawSoy1 <- merge(
  subset(x, !duplicated(x)), 
  subset(y, !duplicated(y)), 
  by = c("genome1", "genome2")
)
# Second merge (swapping the columns)
corRawSoy2 <- merge(
  subset(x, !duplicated(x)),
  subset(y, !duplicated(y)),
  by.x = c("genome2", "genome1"),
  by.y = c("genome1", "genome2")
)
# Combine the results
corRawSoy <- rbind(corRawSoy1, corRawSoy2)
# Filter out self-pairs if needed
corRawSoy <- subset(corRawSoy, genome1 != genome2)

# -- reorder for coloring
corRawSoy[,typeGrp := ifelse(pavgrp == pkgrp, as.character(pavgrp), "diff")]
corRawSoy[,typeGrp := factor(typeGrp, levels = c("Chu2021", "Liu2020", "Phytozome", "ref", "diff"))]

# -- plot
pltCotCorrs <- ggplot(subset(corRawSoy, pkdist > 0 & pavdist > 0), aes(x = pkdist, y = pavdist, col = typeGrp, group = 1))+
  geom_point(size = .25)+
  facet_grid(type ~ ., switch = "y", scale = "free", space = "free")+

  # coord_fixed()+
  scale_color_manual(values = c("green3", "dodgerblue2", "cyan", "black", "lightgrey"), guide = "none")+
  theme(axis.text = element_text(size = 6,color = "black"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "black", linewidth = .25, linetype = 2),
        strip.text = element_text(size = 7, color = "black"),
        strip.placement = "outside",
        strip.clip = "off",
        strip.text.y.left = element_text(angle=0, vjust=1),
        strip.text.y = element_text(margin = margin(t=-10, r=-30), size = 8),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0.01,0.01))+
  scale_x_continuous(expand = c(0.01,0.01))


################################################################################
################################################################################
# 1.6 Combine plots and save

pdf(file.path(out_dir, "fig2_v0.6_heatmapsCorSoyRaw.pdf"),
    height = 2.2, width = 9)
cowplot::plot_grid(pltSoyHeatmaps, pltCotCorrs, rel_widths = c(3,1))
dev.off()
################################################################################
################################################################################

################################################################################
################################################################################
# 2. Analyze PAV distances wrt pankmer distances in Cotton
################################################################################
################################################################################

################################################################################
# 2.1 parameters and paths


################################################################################
# 2.2 Cotton metadata
mdCot <- data.table(id = unique(fread(igcCotFile, select = "genome")[[1]]))
mdCot[,grp := ifelse(id %in% c("FM958", "UA48", "U1HAP1", "TM1"), "Phytozome",
                     ifelse(id %in% c("B713", "Bar32"), "homology", "integrative"))]

# -- get vector of groupings
cotGrp <- mdCot$grp; names(cotGrp) <- mdCot$id

# -- read in pankmer matrix and rename to id
pkCot <- read_pankmerDist(pankmerCotFile)
rownames(pkCot) <- colnames(pkCot) <- gsub(".genome$", "", colnames(pkCot))

# -- read in the combBed files
rawCotIn <- fread(rawCotFile)
igcCotIn <- fread(igcCotFile)

# -- subset to only the big chrs
igcCotIn[,ngenes := .N, by = c("genome", "chr")]
igcCotIn <- subset(igcCotIn, ngenes > minChrSize)

rawCotIn[,ngenes := .N, by = c("genome", "chr")]
rawCotIn <- subset(rawCotIn, ngenes > minChrSize)

################################################################################
# 2.3 convert files to PAV distance matrices

# -- igc globHOGs
igcCot <- with(igcCotIn, convert_cb2pav(genome = genome, orthogroupID = globHOG))

# -- raw globHOGs
rawCot <- with(rawCotIn, convert_cb2pav(genome = genome, orthogroupID = globHOG))

# -- pankmer 
pkCot <- read_pankmerDist(pankmerCotFile)
colnames(pkCot) <- rownames(pkCot) <- 
  gsub("-","_",gsub(".genome$", "", rownames(pkCot)))

# -- make heatmaps
heatmap(pkCot, scale = "none")
heatmap(rawCot, scale = "none")
heatmap(igcCot, scale = "none")

# -- get hclust vectors
pkCotOrd <- get_hclustOrd(pkCot)
rawCotOrd <- get_hclustOrd(rawCot)
igcCotOrd <- get_hclustOrd(igcCot)

# -- reshape and add metadata
pkCotTp <- melt_matrix(pkCot)
rawCotTp <- melt_matrix(rawCot)
igcCotTp <- melt_matrix(igcCot)


################################################################################
# 2.4 reformat and plot heatmaps

# -- melt matrices, reorder genomes and combine
tpCot <- rbind(
  data.table(melt_matrix(pkCot, genomeOrd = pkCotOrd), type = "pankmer"),
  data.table(melt_matrix(rawCot, genomeOrd = rawCotOrd), type = "hogsRaw"),
  data.table(melt_matrix(igcCot, genomeOrd = igcCotOrd), type = "hogsIGC"))

# -- add the line colors
tpCot[,`:=`(rugColx = cotGrp[genome1], rugColy = cotGrp[genome2])]

# -- add plot ordering
tpCot[,type := factor(type, levels = c("pankmer", "hogsRaw", "hogsIGC"))]

# -- make the heatmap plots
tpRugCot <- subset(tpCot, genome1 == genome2)
pltCotHeatmaps <- ggplot()+
  geom_tile(
    data = tpCot, 
    aes(x = x, y = y, fill = value))+
  facet_grid(.~type)+
  scale_fill_viridis_c(option = "A")+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.key.size = unit(.4, 'cm'),
    legend.title = element_blank(), 
    legend.text=element_text(size=6),
    title = element_text(size = 7),
    strip.text = element_text(size = 6))+
  geom_rug(data = tpRugCot, aes(x = x, col = rugColx), linewidth = 2, outside = T, sides = "bl")+
  coord_fixed(clip = "off")+
  scale_color_manual(values = c("green3", "dodgerblue2", "cyan", "black"), guide = "none")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0))

################################################################################
# 2.5 Make the xy plots and do correlations

# -- calculate mantel correlation coefs
vegan::mantel(pkCot, rawCot)
vegan::mantel(pkCot, igcCot)

# -- reorder genomes so no duplicate
tpCot[,`:=`(g1 = ifelse(x > y, genome1, genome2),
            g2 = ifelse(x > y, genome2, genome1))]

# -- pull out pankmer data
x <- with(subset(tpCot, type == "pankmer"), data.table(
  genome1 = g1, 
  genome2 = g2, 
  pkdist = value, 
  pkgrp = ifelse(rugColx == rugColy, rugColx, "diff")))

# -- pull out pav data
y <- with(subset(tpCot, type != "pankmer"), data.table(
  genome1 = g1, 
  genome2 = g2, 
  pavdist = value, 
  type = type, 
  pavgrp = ifelse(rugColx == rugColy, rugColx, "diff")))

# -- merge
# First merge
corRawCot1 <- merge(
  subset(x, !duplicated(x)),
  subset(y, !duplicated(y)),
  by = c("genome1", "genome2")
)
# Second merge (swapping the columns)
corRawCot2 <- merge(
  subset(x, !duplicated(x)),
  subset(y, !duplicated(y)),
  by.x = c("genome2", "genome1"),
  by.y = c("genome1", "genome2")
)
# Combine the results
corRawCot <- rbind(corRawCot1, corRawCot2)
# Filter out self-pairs if needed
corRawCot <- subset(corRawCot, genome1 != genome2)

# -- reorder for coloring
corRawCot[,typeGrp := ifelse(pavgrp == pkgrp, as.character(pavgrp), "diff")]
corRawCot[,typeGrp := factor(typeGrp, levels = c("homology", "integrative", "diff"))]

# -- plot
pltCotCorrs <- ggplot(subset(corRawCot, pkdist > 0 & pavdist > 0), aes(x = pkdist, y = pavdist, col = typeGrp, group = 1))+
  # stat_smooth(
  #   method = "lm", 
  #   se = F, linetype = 1, linewidth = 2)+
  geom_point(size = 1)+
  facet_grid(type ~ ., switch = "y", scale = "free", space = "free")+
  
  # coord_fixed()+
  scale_color_manual(values = c("green3", "dodgerblue2", "lightgrey"), guide = "none")+
  theme(axis.text = element_text(size = 6,color = "black"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "black", linewidth = .25, linetype = 2),
        strip.text = element_text(size = 7, color = "black"),
        strip.placement = "outside",
        strip.clip = "off",
        strip.text.y.left = element_text(angle=0, vjust=1),
        strip.text.y = element_text(margin = margin(t=-10, r=-30), size = 8),
        strip.background = element_blank())+
  scale_y_continuous(expand = c(0.05,0.05))+
  scale_x_continuous(expand = c(0.02,0.02))

################################################################################
################################################################################
# 2.6 Combine multipanel plots and save

pdf(file.path(out_dir, "figS1_v0.6_heatmapsCorCotRaw.pdf"),
    height = 2.2, width = 9)
cowplot::plot_grid(pltCotHeatmaps, pltCotCorrs, rel_widths = c(3,1))
dev.off()
################################################################################
################################################################################

################################################################################
################################################################################
# 3. Build pangenome expansion curves
################################################################################
################################################################################

################################################################################
# 3.0 ad hoc functions
calc_expCurv <- function(cb,             # GENESPACE combBed.txt read in as a data.table
                         maxNSamp = 1000,# maximum number of samples to simulate
                         nCores = 4){    # number of parallel processes
  ig <- unique(cb$genome)
  ng <- length(ig)
  nclist <- lapply(1:ng, function(i){
    cat(sprintf("%s / %s: ", i, ng))
    gc <- lapply(1:(maxNSamp*10), function(j)
      sample(ig, i, replace = F))
    isd <- sapply(gc, function(x) paste(x[order(x)], collapse = ""))
    gc <- gc[!duplicated(isd)]
    if(length(gc) > maxNSamp) gc <- gc[1:maxNSamp]
    nc <- unlist(mclapply(gc, mc.cores = nCores, function(j)
      with(cb, uniqueN(globHOG[genome %in% j]))))
    cat(round(mean(nc)), "\n")
    return(nc)
  })
  return(nclist)
}

calc_crvStats <- function(x, 
                          grpID, 
                          quants = c(.05, .5, .95)){
  out <- rbindlist(lapply(1:length(x), function(i){
    y <- x[[i]]
    qs <- quantile(y, quants)
    return(data.table(
      grp = grpID,
      nGenomes = i,
      propGenomes = i/length(x),
      quantile = quants,
      ord = ifelse(quants == 0.05, -i, i),
      nOGs = as.numeric(qs)))
  }))
  setkey(out, quantile, ord)
  return(out)
}

################################################################################
# 3.1 Read in the kmer expansion curves
nsamp <- 500

pggrow <- rbind(
  data.table(species = "cotton", nGenomes = 1:length(kmerCot), nBp = kmerCot),
  data.table(species = "soybean", nGenomes = 1:length(kmerSoy), nBp = kmerSoy))

################################################################################
# 3.2 split out Liu and Chu for within method curves
chuIds <- mdSoy$id[mdSoy$grp == "Chu2021"]
liuIds <- mdSoy$id[mdSoy$grp == "Liu2020"]

chuSoyIn <- subset(rawSoyIn, genome %in% chuIds)
liuSoyIn <- subset(rawSoyIn, genome %in% liuIds)
chuSoyIgcIn <- subset(igcSoyIn, genome %in%chuIds)
liuSoyIgcIn <- subset(igcSoyIn, genome %in% liuIds)

################################################################################
# 3.3 combine all together in a list
cbList <- list(
  cottonRaw = subset(rawCotIn, !duplicated(paste(genome, globHOG))), 
  cottonIGC = subset(igcCotIn, !duplicated(paste(genome, globHOG))), 
  soybeanRaw = subset(rawSoyIn,  !duplicated(paste(genome, globHOG))), 
  soybeanIGC = subset(igcSoyIn, !duplicated(paste(genome, globHOG))), 
  soybeanLiuRaw = subset(liuSoyIn,  !duplicated(paste(genome, globHOG))), 
  soybeanLiuIGC = subset(liuSoyIgcIn, !duplicated(paste(genome, globHOG))), 
  soybeanChuRaw = subset(chuSoyIn,  !duplicated(paste(genome, globHOG))), 
  soybeanChuIGC = subset(chuSoyIgcIn,  !duplicated(paste(genome, globHOG))))

################################################################################
# 3.4 Pangenome expansion curves
expList <- sapply(cbList, simplify = F, USE.NAMES = T, function(x)
  calc_expCurv(cb = x, maxNSamp = nsamp, nCores = 4))

# convert to a data.table
exps <- rbindlist(lapply(names(expList), function(i)
  calc_crvStats(expList[[i]], grpID = i)))

# add metadata for plotting
exps[,`:=`(species = ifelse(grepl("cotton", grp), "cotton", "soybean"),
           consort = ifelse(grepl("Chu", grp), "Chu2021", 
                            ifelse(grepl("Liu", grp), "Liu2020", "all")),
           method = ifelse(grepl("Raw", grp), "Original", "IGC"))]

################################################################################
# 3.5 Plot panel A-B: cotton and soybean all

# -- scale pangenome expansion curves so that the maximum = mean(max) OG curves
exps[,n1kOGs := nOGs / 1000]
mxs <- with(exps, tapply(n1kOGs, paste(species, method), function(x) max(x)))
mxs <- c(cotton = mean(mxs[grep("cotton", names(mxs))]),
         soybean = mean(mxs[grep("soybean", names(mxs))]))
pggrow[,mx := mxs[species]]
pggrow[,sclbp := scale_between(x = c(0, nBp), 0, mx[1])[-1], by = "species"]

# --  get colors etc. together
colsRaw <- c("#472974","#ce4071")

################################################################################
# -- 3.6 Plot full expansion curves (Fig. 2A-B)
tpa <- subset(exps, consort == "all")
p2ab <- ggplot(tpa)+
  geom_bar(
    data = pggrow, 
    aes(x = nGenomes, y = sclbp),
    fill = "black", alpha = 0.2, stat = "identity")+
  facet_grid(. ~ species, scale = "free")+
  geom_polygon(
    data = subset(tpa, quantile != 0.5),
    aes(x = nGenomes, y = n1kOGs, fill = method),
    alpha = .5, col = NA)+
  geom_line(
    data = subset(tpa, quantile == 0.5),
    aes(x = nGenomes, y = n1kOGs, col = method),
    # linetype = 2,
    linewidth = 1)+
  scale_y_continuous(
    expand = c(0,0))+
  scale_color_manual(values = colsRaw, guide = "none")+
  scale_fill_manual(values = colsRaw, guide = "none")+
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#b5b3b3", linetype = 2, linewidth = .1),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 8, color = "black"),
    strip.background = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7))

################################################################################
# -- 3.7 Plot within-method soybean expansion curves (Fig. 2A-B)
tpc <- subset(exps, consort != "all")
p2c <- ggplot(tpc)+
  facet_grid(consort ~ ., scale = "free", space = "free", switch = "y")+
  geom_polygon(
    data = subset(tpc, quantile != 0.5),
    aes(x = nGenomes, y = n1kOGs, fill = method),
    alpha = .5, col = NA)+
  geom_line(
    data = subset(tpc, quantile == 0.5),
    aes(x = nGenomes, y = n1kOGs, col = method),
    # linetype = 2,
    linewidth = 1)+
  scale_y_continuous(
    expand = c(0,0), limits = c(0, max(tpc$n1kOGs)))+
  scale_color_manual(values = colsRaw, guide = "none")+
  scale_fill_manual(values = colsRaw, guide = "none")+
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#b5b3b3", linetype = 2, linewidth = .1),
    axis.ticks = element_blank(),
    strip.placement = "outside",
    strip.clip = "off",
    strip.text.y.left = element_text(angle=0, vjust=1),
    strip.text.y = element_text(margin = margin(t=0, r=-40), size = 8, color = "black"),
    strip.background = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7))
p2abc <- cowplot::plot_grid(p2ab, p2c, rel_widths = c(3, 1), nrow = 1)


################################################################################
################################################################################
# 4. Build pangenome categorizations
################################################################################
################################################################################

################################################################################
# 4.0 ad hoc functions
calc_pgCats <- function(genome, og){
  x <- data.table(genome = genome, og = og)
  x <- subset(x, !duplicated(x))
  x[,nGenomes := uniqueN(genome), by = "og"]
  x[,propGenome := nGenomes / uniqueN(genome)]
  x[,cls := ifelse(propGenome == 1, "core",
                   ifelse(propGenome >= .9, "nearcore",
                          ifelse(propGenome >= .5, "shell",
                                 ifelse(nGenomes == 1, "private", "cloud"))))]
  return(x)
}

################################################################################
# 4.1 calculate categories by gene
catList <- sapply(cbList, simplify = F, USE.NAMES = T, function(x)
  calc_pgCats(genome = x$genome, og = x$globHOG))

################################################################################
# 4.2 prepare data for barplots
cntList <- rbindlist(lapply(names(catList), function(i){
  x <- catList[[i]]
  y <- x[,list(n1kOGs = uniqueN(og)/1000), by = "nGenomes"]
  y[,grp := i]
  return(y)
}))

# add metadata for plotting
cntList[,`:=`(
  species = ifelse(grepl("cotton", grp), "cotton", "soybean"),
  consort = ifelse(grepl("Chu", grp), "Chu2021", 
                   ifelse(grepl("Liu", grp), "Liu2020", "all")),
  method = ifelse(grepl("Raw", grp), "Original", "IGC"))]

################################################################################
# 4.3 make full species barplots
tpc <- subset(cntList, consort == "all")
p2de <- ggplot(tpc)+
  facet_grid(. ~ species, scale = "free")+
  geom_bar(
    data = tpc,
    aes(x = nGenomes, y = n1kOGs, fill = method),
    position = "dodge", stat = "identity", col = NA)+
  scale_fill_manual(values = colsRaw, guide = "none")+
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#b5b3b3", linetype = 2, linewidth = .1),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 8, color = "black"),
    strip.background = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7))

################################################################################
# 4.4 make within method barplots
tpf <- subset(cntList, consort != "all")
p2f <- ggplot(tpf)+
  facet_grid(consort ~ ., scale = "free", space = "free", switch = "y")+
  geom_bar(
    data = tpf,
    aes(x = nGenomes, y = n1kOGs, fill = method),
    position = "dodge", stat = "identity", col = NA)+
  scale_fill_manual(values = colsRaw, guide = "none")+
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#b5b3b3", linetype = 2, linewidth = .1),
    axis.ticks = element_blank(),
    strip.placement = "outside",
    strip.clip = "off",
    strip.text.y.left = element_text(angle=0, vjust=1),
    strip.text.y = element_text(margin = margin(t=0, r=-40), size = 8, color = "black"),
    strip.background = element_blank(),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7))
p2def <- cowplot::plot_grid(p2de, p2f, rel_widths = c(3, 1), nrow = 1)

################################################################################
# 4.5 make pies
# -- calcuate n genes
pieTp <- rbindlist(lapply(names(catList), function(i){
  x <- catList[[i]]
  y <- x[,list(nHOGs = .N), by = "cls"]
  y[,grp := i]
  return(y)
}))

# -- reclass and color categories
pieTp[,cls := factor(cls, levels = c("private", "cloud", "shell", "nearcore", "core"))]
pieCols <- RColorBrewer::brewer.pal(n = 5, name = "RdYlBu")

# -- get metadata in order
pieTp[,`:=`(
  species = ifelse(grepl("cotton", grp), "cotton", "soybean"),
  consort = ifelse(grepl("Chu", grp), "Chu2021", 
                   ifelse(grepl("Liu", grp), "Liu2020", "all")),
  method = ifelse(grepl("Raw", grp), "Original", "IGC"),
  speciesCons = gsub("IGC|raw|Raw","", grp))]
pieTp[,percHOGs := nHOGs/sum(nHOGs), by = "grp"]
pieTp[,grp := factor(grp, levels = unique(pieTp$grp))]
p2defPies <- ggplot(pieTp, aes(x="", y=percHOGs, fill=cls, label = paste0(round(percHOGs * 100), "%"))) +
  geom_bar(stat="identity", width=1) +
  geom_text(position = position_stack(vjust = 0.5), size = 6/.pt, color = "black")+
  coord_polar("y", start=0)+
  facet_grid(. ~ grp, scale = "free")+
  scale_fill_manual(values = pieCols, guide = "none")+
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 8, color = "black"),
    strip.background = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank())

################################################################################
################################################################################
# 4.7 save the multipanel plot

pdf(file.path(out_dir, "fig2_v0.6raw.pdf"), height = 5, width = 9)
cowplot::plot_grid(p2abc, p2def, p2defPies, nrow = 3, rel_heights = c(2, 2, 1.2))
dev.off()
################################################################################
################################################################################


################################################################################
################################################################################
# 5. Build fixed order category curves (Fig. S2)
################################################################################ 
################################################################################

################################################################################
# 5.0 ad hoc functions
countPavTypes_splitCB <- function(cb, genomes){
  out <- rbindlist(lapply(1:length(genomes), function(i){
    j <- genomes[i]
    gs <- genomes[1:i]
    cbi <- subset(cb, genome %in% gs)
    cntByOg <- cbi[,list(inNGenomes = uniqueN(genome)), by = "globHOG"]
    cntByCnt <- cntByOg[,list(n = .N), by = "inNGenomes"]
    return(data.table(nGenomes = i, genome = j, cntByCnt))
  }))
  return(out)
}

################################################################################
# 5.1 get counts of genes for x axis ordering
nr <- with(rawCotIn, table(genome))
mdCot[,nRawGenes := nr[id]]
ni <- with(igcCotIn, table(genome))
mdCot[,nIGCGenes := ni[id]]

nr <- with(rawSoyIn, table(genome))
mdSoy[,nRawGenes := nr[id]]
ni <- with(igcSoyIn, table(genome))
mdSoy[,nIGCGenes := ni[id]]

################################################################################
# 5.2 get ordering vectors
mdCot[,grpFac := factor(grp, levels = c("homology", "integrative"))]
setorder(mdCot, grpFac, -nRawGenes)
rawCotOrd <- mdCot$id
mdCot[,barOrd := 1:.N]

mdSoy[,grpFac := factor(grp, levels = c("ref", "Phytozome", "Chu2021", "Liu2020"))]
setorder(mdSoy, grpFac, -nRawGenes)
rawSoyOrd <- mdSoy$id
mdSoy[,barOrd := 1:.N]

################################################################################
# 5.3 split combBed files and calculate PAV categories in each grouping
splListByn <- list(
  cottonRaw = countPavTypes_splitCB(cb = rawCotIn, genomes = rawCotOrd), 
  cottonIGC = countPavTypes_splitCB(cb = igcCotIn, genomes = rawCotOrd), 
  soybeanRaw = countPavTypes_splitCB(cb = rawSoyIn, genomes = rawSoyOrd), 
  soybeanIGC = countPavTypes_splitCB(cb = igcSoyIn, genomes = rawSoyOrd))

################################################################################
# 5.4 make plotting datasets for each species
frdCot <- rbindlist(lapply(names(splListByn)[1:2], function(i)
  data.table(grp = i, splListByn[[i]])))
frdSoy <- rbindlist(lapply(names(splListByn)[3:4], function(i)
  data.table(grp = i, splListByn[[i]])))

frdCot[,pg := inNGenomes / nGenomes]
frdSoy[,pg := inNGenomes / nGenomes]
################################################################################
# 5.5 get metadata in order
# frdCot[,cls := factor(cls, levels = c("private", "cloud", "shell", "nearcore", "core"))]
# frdSoy[,cls := factor(cls, levels = c("private", "cloud", "shell", "nearcore", "core"))]

# -- add in consortium
cotGrp <- mdCot$grp; names(cotGrp) <- mdCot$id
soyGrp <- mdSoy$grp; names(soyGrp) <- mdSoy$id
frdCot[,`:=`(
  consort = factor(cotGrp[as.character(genome)], levels = c("homology", "integrative")),
  method = ifelse(grepl("Raw", grp), "Original", "IGC"))]
frdSoy[,`:=`(
  consort = factor(soyGrp[as.character(genome)], levels = c("ref", "Phytozome", "Chu2021", "Liu2020")),
  method = ifelse(grepl("Raw", grp), "Original", "IGC"))]

# -- xaxis labeling 
cntxCot <- subset(frdCot, !duplicated(paste(genome, grp)))
cntxSoy <- subset(frdSoy, !duplicated(paste(genome, grp)))

# -- coloring
pieCols <- RColorBrewer::brewer.pal(n = 5, name = "RdYlBu")
conColsCot <- c("green3", "dodgerblue2", "lightgrey")
conColsSoy <- c("black","cyan", "green3", "dodgerblue2")
mxCot <- with(frdCot, max(tapply(n, paste(genome, method), sum)))
mxSoy <- with(frdSoy, max(tapply(n, paste(genome, method), sum)))


################################################################################
# 5.6 plot
# -- cotton
plt4barCot <- ggplot()+
  geom_bar(data = frdCot, aes(x = nGenomes, y = n/1000, fill = pg, group = pg), 
           stat = "identity", position = "stack")+
  geom_rug(data = cntxCot, aes(x = nGenomes, col = consort), linewidth = 8, outside = F, sides = "b")+
  facet_grid(.~grp, scale = "free", space = "free")+
  scale_x_continuous(breaks = cntxCot$nGenomes, labels = cntxCot$genome)+
  scale_fill_gradientn(colors = pieCols, guide = "none")+
  scale_color_manual(values = conColsCot, guide = "none")+
  scale_y_continuous(expand = c(0, 0), limits = c(-5, 5+mxCot/1000))+
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "black", linetype = 2, linewidth = .25),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 7))+
  labs(x = "genome ID ordered by consortium and decreasing number of genes",
       y = "n. phylogenetic orthogroups (x1000)")

plt4barSoy <- ggplot()+
  geom_bar(data = frdSoy, aes(x = nGenomes, y = n/1000, fill = pg, group = pg), 
           stat = "identity", position = "stack")+
  geom_rug(data = cntxSoy, aes(x = nGenomes, col = consort), linewidth = 3, outside = F, sides = "b")+
  facet_grid(.~grp, scale = "free", space = "free")+
  scale_x_continuous(breaks = cntxSoy$nGenomes, labels = cntxSoy$genome)+
  scale_fill_gradientn(colors = pieCols, guide = "none")+
  scale_color_manual(values = conColsSoy, guide = "none")+
  scale_y_continuous(expand = c(0, 0), limits = c(-5, 5+mxSoy/1000))+
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "black", linetype = 2, linewidth = .25),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 7))+
  labs(x = "genome ID ordered by consortium and decreasing number of genes",
       y = "n. phylogenetic orthogroups (x1000)")

pdf(file.path(out_dir, "Splfig2_v0.6raw.pdf"), height = 5, width = 6)
cowplot::plot_grid(plt4barCot, plt4barSoy, nrow = 2)
dev.off()
