## This script can be used to test out the functions defined and is mainly for application development purposes

# Source to the plottingFunctionsZack.R script

source("~/GenomicProfilingApp/R/plottingFunctionsZack.R")

# Read in the region .bed file as a GRanges object

bed_gr <- import.bed("~/GenomicProfilingApp/data/top3000.bed")

# Get flanks for the genes

genes <- getFeature(bed_gr)
gr_flank <- getFeature(bed_gr, start_flank = 20, end_flank = 20)

# Import the bigwig files

fbw <- list.files("~/Research Project/wb3Feb/Test Data", pattern = ".f.bw$", full.names = TRUE)
names(fbw) <- fbw %>% str_remove("/home/s2274585/Research Project/wb3Feb/Test Data/") %>% str_remove(".f.bw")
bwf <- importBWlist(fbw, names(fbw), selection = gr_flank)
rbw <- list.files("~/Research Project/wb3Feb/Test Data", pattern = ".r.bw$", full.names = TRUE)
names(rbw) <- rbw %>% str_remove("/home/s2274585/Research Project/wb3Feb/Test Data/") %>% str_remove(".r.bw")
bwr <- importBWlist(rbw, names(rbw), selection = gr_flank)

# Make the grl from the genes

grl <- list("genes" = genes)

# Make the matrix list

matl <- matList(bwf = bwf, bwr = bwr, grl = grl, names = names(fbw), extend = 20, w = 1, strand = "for")

# Defining the windows and colour palettes for the heatmaps

colmap = c(wes_palette("Darjeeling2")[2],wes_palette("Darjeeling1")[4],wes_palette("Zissou1")[1],wes_palette("Darjeeling1")[3])
colmap = rep(colmap,3)

wins <- c("Upstream" = 20, "Gene" = 40, "Downstream" = 20)

# Creation of the enrichedHeatmap heatmap list

hml <- hmList(matl = matl, wins = wins, show_row_names = FALSE, col_fun = "red0", axis_labels = c("-20b", "TSS", "TES", "+20b"))

