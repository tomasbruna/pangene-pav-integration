library(data.table)
library(ggplot2)

################################################################################
# ad hoc function to pull orthologs and categories
compile_orthoCounts <- function(orthoFile,
                                queryGenome,
                                targetGenome,
                                combBedList){
  bq <- combBedList[[queryGenome]][,c("genome", "id")]
  bt <- combBedList[[targetGenome]][,c("genome", "id")]
  ort <- parse_orthologues(orthoFile)
  ort[,n1 := .N, by = "id1"]
  ort[,n2 := .N, by = "id2"]
  ort[,cat := ifelse(n1 == n2 & n2 == 1, "1x",
                     ifelse(n1 == n2, "multi", "CNV"))]
  allid <- rbind(bq, bt)
  orthcat <- with(ort, data.table(
    genome = c(gen1, gen2),
    id = c(id1, id2),
    cat = c(cat, cat)))
  orthcat <- subset(orthcat, !duplicated(orthcat))

  qout <- merge(allid, orthcat, by = c("genome", "id"), all = T)
  qout$cat[is.na(qout$cat)] <- "PAV"
  out <- qout[,list(n = .N), by = c("genome", "cat")]
  out[,`:=`(query = queryGenome, target = targetGenome)]
  out[,tot := sum(n), by = "genome"]
  out[,prop := n/tot]
  return(out)
}

parse_divtime <- function(disFile, parseNames = FALSE){
  dis <- fread(disFile)
  if(parseNames){
    din <- paste0(substr(dis[[1]], 1,1), gsub(".* ", "", dis[[1]]))
    din[din == "CPawnee"] <- "CillinoinensisPawnee"
    dis[[2]] <- din
    setnames(dis,2, "GenomeID")
  }
  ns <- dis$GenomeID
  matstart <- which(colnames(dis) == "RefSeq_year") + 1
  setnames(dis, matstart:ncol(dis), ns)
  dis <- dis[,c("GenomeID", as.character(ns)), with = F]
  dim <- melt(dis, id.vars = "GenomeID")
  setnames(dim, c("query", "target", "mya"))
  dim[,`:=`(mya = as.numeric(mya))]
  return(dim)
}

plot_meanCat <- function(areaDat, linMod, grpName){
  plt <- ggplot()+
    geom_area(data = areaDat, aes(x = mya, y = meanPerc, fill = cat))+
    stat_smooth(
      data = subset(areaDat, cat == "PAV"),
      aes(x = mya, y = meanPerc),
      method = "lm", col = "black", linewidth = .5, linetype = 2, formula=y~x)+
    scale_fill_manual(values = c("skyblue","dodgerblue", "dodgerblue4", "black"))+
    labs(title = sprintf("%s: PAV ~ %s%% + %s%% (per 100mya)",
                         grpName, round(coef(linMod)[1], 1), round(coef(linMod)[2]*100, 1)),
         x = "Divergence time (M ybp)",
         y = "Mean % of genes (per phylogentic node)")+
    theme(
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      # axis.ticks = element_blank(),
      axis.text = element_text(family = "Helvetica", size = 6, color = "black"),
      axis.title = element_text(family = "Helvetica", size = 7, color = "black"),
      title = element_text(family = "Helvetica", size = 8, color = "black"))+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100))
  return(plt)
}

parse_orthologues <- function(filepath){
  orthID <- id1 <- id2 <- NULL
  x <- fread(filepath, showProgress = F, header = TRUE, check.names = FALSE)
  refID <- colnames(x)[2]
  altID <- colnames(x)[3]
  setnames(x, c("og", "id1", "id2"))
  x1 <- subset(x, !grepl(",", paste(id1, id2)))
  x2 <- subset(x, grepl(",", paste(id1, id2)))
  x2[,orthID := 1:.N]
  x2r <- x2[,list(id1 = unique(strsplit(id1, ",")[[1]])),
            by = "orthID"]
  x2a <- x2[,list(id2 = unique(strsplit(id2, ",")[[1]])),
            by = "orthID"]
  x2 <- merge(x2r, x2a, by = "orthID", all = T, allow.cartesian = T)
  x1[,orthID := (1:.N)+max(x2$orthID)]
  x <- rbind(x1[,colnames(x2), with = F], x2)
  x[,`:=`(gen1 = refID, gen2 = altID,
          id1 = gsub(" ", "", id1), id2 = gsub(" ", "", id2))]
  return(x)
}

################################################################################
# Download data
options(timeout = 3600)
download.file(
  url = "https://zenodo.org/records/17902589/files/pav_cnv_plants_amniotes_data.tgz",
  destfile = "pav_cnv_plants_amniotes_data.tgz",
  mode = "wb"
)
untar("pav_cnv_plants_amniotes_data.tgz", exdir = ".")

################################################################################
# File paths
# -- working directory
wd <- getwd()

# -- divergence time estimates from timetree
plantDivTimeMatrix <- file.path(
  wd, "pav_cnv_plants_amniotes_data/PlantDivergenceTimes.csv")
animalDivTimeMatrix <-file.path(
  wd, "pav_cnv_plants_amniotes_data/AmnioteDivergenceTimes.csv")

# -- combined bed files from genespace
plantCombBed <- file.path(
  wd, "pav_cnv_plants_amniotes_data/plants/combBed.txt")
animalCombBed <- file.path(
  wd, "pav_cnv_plants_amniotes_data/amniotes/combBed.txt")

# -- pairwise ortholog files from orthofinder
animalOrthos <- list.files(
  path = file.path(wd, "pav_cnv_plants_amniotes_data/amniotes"),
  pattern = "__v__", full.names = T)
plantOrthos <- list.files(
  path = file.path(wd, "pav_cnv_plants_amniotes_data/plants"),
  pattern = "__v__", full.names = T)

################################################################################
# Animals
# -- parse the ortholog files to get genome IDs
fmd <- data.table(f = animalOrthos)
fmd[,c("query", "target") := tstrsplit(gsub(".tsv", "", basename(f)), "__v__")]

# -- read in the combined bed file and split by genome (for faster access)
cb <- fread(animalCombBed)
cbs <- split(cb, by = "genome")

# -- compile orthologs into pav/cnv categories
pavcAnimal <- rbindlist(lapply(1:nrow(fmd), function(i){
  qu <- fmd$query[i]
  ta <- fmd$target[i]
  of <- fmd$f[i]
  if(qu != ta){
    outCnts <- compile_orthoCounts(
      orthoFile = of,
      queryGenome = qu,
      targetGenome = ta,
      combBedList = cbs)
    return(outCnts)
  }
}))

# -- parse divergence time matrix
dim <- parse_divtime(animalDivTimeMatrix, parseNames = FALSE)

# -- merge with pav counts and get mean counts by divergence time
dimAnimal <- merge(dim, pavcAnimal, by = c("query", "target"))
dilAnimal <- dimAnimal[,list(meanPerc = mean(prop) * 100), by = c("mya", "cat")]

# -- run linear model and correlations
ctAnimal <- with(subset(dilAnimal, cat == "PAV"), cor.test(x = mya, y = meanPerc))
modAnimal <- with(subset(dilAnimal, cat == "PAV"), lm(meanPerc ~ mya))

# -- make the plot
pltAnimals <- plot_meanCat(areaDat = dilAnimal, linMod = modAnimal, grpName = "Amniotes")
print(pltAnimals)

################################################################################
# Plants
# -- parse the ortholog files to get genome IDs
fmd <- data.table(f = plantOrthos)
fmd[,c("query", "target") := tstrsplit(gsub(".tsv", "", basename(f)), "__v__")]

# -- read in the combined bed file and split by genome (for faster access)
cb <- fread(plantCombBed)
cbs <- split(cb, by = "genome")

# -- compile orthologs into pav/cnv categories
pavcPlant <- rbindlist(lapply(1:nrow(fmd), function(i){
  qu <- fmd$query[i]
  ta <- fmd$target[i]
  of <- fmd$f[i]
  if(qu != ta){
    outCnts <- compile_orthoCounts(
      orthoFile = of,
      queryGenome = qu,
      targetGenome = ta,
      combBedList = cbs)
    return(outCnts)
  }
}))

# -- parse divergence time matrix (need to reformat names, see ad hoc above)
dim <- parse_divtime(plantDivTimeMatrix, parseNames = TRUE)

# -- merge with pav counts and get mean counts by divergence time
dimPlant <- merge(dim, pavcPlant, by = c("query", "target"))
dilPlant <- dimPlant[,list(meanPerc = mean(prop) * 100), by = c("mya", "cat")]

# -- run linear model and correlations
ctPlant <- with(subset(dilPlant, cat == "PAV"), cor.test(x = mya, y = meanPerc))
modPlant <- with(subset(dilPlant, cat == "PAV"), lm(meanPerc ~ mya))

# -- make the plot
pltPlants <- plot_meanCat(areaDat = dilPlant, linMod = modPlant, grpName = "Plants")
print(pltPlants)

################################################################################
# Pie charts
# -- make data for pie chart
pietp <- rbind(
  data.table(subset(pavcPlant, query == "Osativa" & target == "Sbicolor"), grp = "plant", rw = 1),
  data.table(subset(pavcPlant, query == "Osativa" & target == "Athaliana"),grp = "plant", rw = 2),
  data.table(subset(pavcPlant, query == "Osativa" & target == "Cjaponica"),grp = "plant", rw = 3),
  data.table(subset(pavcAnimal, query == "Hsapiens" & target == "Ggorilla"),grp = "animal", rw = 1),
  data.table(subset(pavcAnimal, query == "Hsapiens" & target == "Mmusculus"),grp = "animal", rw = 2),
  data.table(subset(pavcAnimal, query == "Hsapiens" & target == "Ggallus"), grp = "animal", rw = 3))
pietp[,lab := paste(query, target)]

# -- make the plot
piePlts <- ggplot(pietp, aes(x = "", y = prop, fill = cat))+
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("skyblue","dodgerblue", "dodgerblue4", "grey"), guide = "none")+
  facet_grid(rw ~ grp)+
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    title = element_blank())
print(piePlts)

################################################################################
# Output

dir.create("R_figures", showWarnings = FALSE)
out_dir <- file.path(getwd(), "R_figures")

# -- combine area maps into a single plot
pdf(file.path(out_dir, "Fig1.pdf"),
    width = 9, height = 3)
cowplot::plot_grid(pltAnimals, pltPlants, nrow = 1)
dev.off()

# -- write pie charts
pdf(file.path(out_dir, "Fig1_pies.pdf"),
    width = 2, height = 3)
print(piePlts)
dev.off()

# -- combine plant and animal pav cats
outc <- rbind(
  data.table(dimPlant, group = "Plants"),
  data.table(dimAnimal, group = "Animals"))
fwrite(outc, file = file.path(getwd(), "SIData1_v1.1_plantAnimalOGCats.csv"))
