# Function used to select the subset of arguments found in dots (ellipsis) that are
# specific for the function \code{fn}
get_args_from_dots = function(dots, fn){
  mask = names(dots) %in% methods::formalArgs(fn)
  return(dots[mask])
}


### Graphic: ggplot2 Theme Setting ------------------------------------------------------
# nowaklab theme for ggplot2

### Graphic: Color Choices ------------------------------------------------------
# Colors: set final colors for organs
sample_col <- c("PRL" = "#65A7F3", "PRR" = "#65A7F3", # Prostate Left, Prostate Right
                   "LVL" = "#BA361E", "LVR" = "#BA361E", "LVM" = "#BA361E",  "LVC" = "#BA361E", #  Liver Left Lobe, Liver Right Lobe, Liver Median Lobe, Liver Caudate Lobe   
                   "LGL" = "#4D0D29", "LGR" = "#4D0D29", # Lung Left, Lung Right
                   "ADL" = "#4387B7", "ADR" = "#4387B7", # Adrenal Left, Adrenal Right
                   "KDL" = "#405484", "KDR" = "#405484", # Kidney Left, Kidney Right
                   "BDR" = "A970DB", # Bladder
                   "BRN" = "#F1D351", # Brain
                   "BLD" = "#DB8C9B", # Blood
                   "SPN" = "#1A3F13", # Spleen
                   "THM" = "#6CD24B", # Thymus
                   # LN
                   "LNCA" = "#669F6C", # Lymph Node Caudal
                   # LN
                   "SCL" = "#E7943B", "SCR" = "#E7943B", # Scapula Left, Scapula Right
                   "HML" = "#E7943B", "HMR" = "#E7943B", # Humerus Left, Humerus Right
                   "TBL" = "#6A4925", "TBR" = "#6A4925", # Tibia Left, Tibia, Right
                   "PVL" = "#E8B985" ,  "PVR" =  "#E8B985", # Pelvis Left/Pelvis Right,
                   "RBL" = "#E7B985", "RBR" = "#E7B985") # Rib Left, Rib Right

# Color Option 1 
# create multiple colors
colony_25xcols <- c("#b2df8a", "#1F78C8", "#ff0000", "#ff7f00", "#36648B", "#FFD700", "#FB6496", "#a6cee3", "#33a02c", "#CAB2D6", 
                    "#FDBF6F", "#999999", "#EEE685", "#C8308C", "#FF83FA", "#C814FA", "#0000FF", "#6A33C2", "#00E2E5", "#00FF00", 
                    "#778B00", "#BEBE00", "#8B3B00", "#A52A3C", "#D9D9D9FF")
# set order of colors
colony_25xcols[-which(colony_25xcols %in% "#D9D9D9FF")]
# add greseqtab_df for cut-off
colony_col_max <- c("#D9D9D9FF", colony_25xcols)

# Theme for ggplot2 [tgood one]
barplot_nowaklab_theme <- function(axis.title.font = "Helvetica", axis.title.col = "black",
                                   axis.text.font = "Helvetica", axis.text.col = "black",  
                                   legend.text.font = "Helvetica")
{
  # General setings
  theme(
    
    # Axis lines
    #axis.line = element_line(colour = "black"),
    axis.line = element_line(), # for 
    
    # Tick axis x and y axes
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(1.5, "mm"),
    
    # X axis text
    # X text straight
    axis.text.x  = element_text(size = 10, color = axis.text.col, angle = 0, vjust = 0, hjust= 0.5, family = axis.text.font),
    axis.title.x = element_text(size = 12, color = axis.title.col, face = "bold", margin = unit(c(3, 3, 3, 3), "mm"), family = axis.title.font),
    # X text angled
    #axis.text.x  = element_text(color="black", angle=55, vjust=1, hjust=1, size=14),
    #axis.title.x = element_text(size=14, color="black", face = "plain"),
    
    # Y axis text  
    axis.text.y = element_text(size = 10, color = axis.text.col, angle = 0, vjust = 0.5, hjust= 1, family = axis.text.font),
    axis.title.y = element_text(size = 12, color = axis.title.col , face = "bold", margin = unit(c(3, 3, 3, 3), "mm"), family = axis.title.font),    
    
    # Legend   
    legend.title = element_text(size = 12, face = "bold"),
    legend.title.align = c(0),
    legend.text = element_text(size = 12, color="black", face = "plain", hjust= 1, margin = margin(l=0.1, r=0.2, unit="cm"), family = legend.text.font),
    legend.text.align = c(0),
    legend.position = "right", 
    legend.box = "vertical",
    legend.background = element_rect(colour = NA),
    legend.key = element_rect(colour=NA, fill="white"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    #legend.key.size = unit(0.5, "cm"),
    #legend.spacing = unit(0.1, "cm"),
    legend.spacing.x = unit(0.1, "cm"), # space between key and text in the legend
    legend.spacing.y = unit(0.2,"cm"), # space between different legends
    #legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    #legend.box.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
    
    # Guides

    # Plot Margins  
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    
    # Border
    #panel.border = element_rect(colour = "gray25", size=0.8, fill=NA)
    panel.border = element_blank(),
    
    # Strip
    #strip.background=element_blank(),
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text.x = element_text(size = 12, color="black", face = "bold"),
    strip.switch.pad.grid=unit(0, "cm"),
    strip.switch.pad.wrap=unit(0, "cm"),
    
    # Panels  
    # Background of panel  
    panel.background=element_rect(fill = "white", colour = "white", size = 1, linetype = "solid"),
    #panel.grid = element_blank(),
    #panel.grid.major = element_line(size = 0.5, linetype = "dashed", colour = "white"), 
    #panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "white"),    
    # X and Y Major
    panel.grid.major.y = element_blank(),
    # Alternative: panel.grid.major.y = element_line(size = 0.5, linetype = "dotted", colour = "#999999"),
    panel.grid.major.x = element_blank(),
    # X and Y Ninor
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
    #panel.spacing.y = unit(-0.5, "lines"),
    #panel.spacing.x = unit(-0.5, "lines"))
  )
  
}



### Identify Mutations ------------------------------------------------------ 
identifyMutations <- function(seq, ref) {
  
  indicesToKeep <- as.vector(as.matrix(ref)!="-")
  seqSubset <- seq[indicesToKeep]
  refSubset <- ref[indicesToKeep]
  # vectors for translating the coordinates
  barcodpToBarcode <- rep(NA, length(seq))
  barcodpToBarcode[which(indicesToKeep)] <- 1:length(which(indicesToKeep))  
  barcodpToBarcode <- nafill(barcodpToBarcode, type = "locf")
  
  barcodeToBarcodep <- rep(NA, length(indicesToKeep))  
  barcodeToBarcodep <-  which(!is.na(barcodpToBarcode))
  
  # barcode + scale
  # Create a dataframe by first creating a vector of T/F in each position of the sequence
  # which tells whether there is mutation in each position wrt the reference
  # then, with the function rle obtain the count of consecutive values (e.g., if the sequence is 1 1 0 0 then, rle 
  # returns length=2 2 values = 1 0, which means that it found 2 ones followed by two zeros)
  # In this way, this code treats each subsequent mutatio as a unique mutation
  # This is ok if that's composed only by deletions, because deletions are acuired together
  # However, if we have multiple substitutions we want to treat them separate. Also,
  # if we find first a substitutions and then a deletion (or the opposite, or a substitution between two deletions)
  # we want to separate the two kinds of mutations.
  # For substitutions, we want to assign an ID to them based on the insertion
  arle.inf.df <- data.frame(unclass(rle(as.numeric(seq == ref)))) 
  # find start and end coordinates of each segment
  arle.inf.df$end.barcodep <- cumsum(arle.inf.df$lengths)
  if (nrow(arle.inf.df)==1) {
    arle.inf.df$start.barcodep <- 1
  } else {
    arle.inf.df$start.barcodep <- c(0, arle.inf.df$end.barcodep[1:(nrow(arle.inf.df)-1)])+1
  }
  # mismatches mean mutations - subset to mutations
  arle.inf.df <- subset(arle.inf.df, values==0)
  
  # if any mutations found, add their sequence and type
  if (nrow(arle.inf.df)>0) {
    # annotate each mutations with alt and ref alleles, type and sequence ID
    for (mi in 1:nrow(arle.inf.df)) {
      # First extract reference and alternative sequences
      arle.inf.df$ref[mi]<-paste(ref[ arle.inf.df$start.barcodep[mi]:arle.inf.df$end.barcodep[mi]], collapse='')
      arle.inf.df$alt[mi]<-paste(seq[ arle.inf.df$start.barcodep[mi]:arle.inf.df$end.barcodep[mi]], collapse='')
    }
    arle.inf.df$type <- "sub"
    arle.inf.df[arle.inf.df$lengths>1, "type"] <- "complex"
    # if alt sequence is made entirely of hyphens, it is a deletions
    arle.inf.df[gsub("-", "", arle.inf.df$alt)=="", "type"] <- "del"
    arle.inf.df[gsub("-", "", arle.inf.df$ref)=="", "type"] <- "ins"
    # add the barcode scale
    arle.inf.df <- subset(arle.inf.df, type=="ins")
    arle.inf.df$start.barcode <- barcodpToBarcode[arle.inf.df$start.barcodep]
    arle.inf.df$end.barcode <- arle.inf.df$start.barcode + 1
    # from the barcodep scale, only select insertions
    
  }    
  
  
  ### barcode scale with a scale ------------------------------------------------------ 
  arle.df <- data.frame(unclass(rle(as.numeric(as.vector(seqSubset) == as.vector(refSubset)))))      
  # find start and end coordinates of each segment
  arle.df$end.barcode <- cumsum(arle.df$lengths)
  if (nrow(arle.df)==1) {
    arle.df$start.barcode <- 1
  } else {
    arle.df$start.barcode <- c(0, arle.df$end.barcode[1:(nrow(arle.df)-1)])+1
  }   
  # mismatches mean mutations - subset to mutations
  arle.df <- subset(arle.df, values==0)
  
  # if any mutations found, add their sequence and type
  if (nrow(arle.df)>0) {
    # annotate each mutations with alt and ref alleles, type and sequence ID
    for (mi in 1:nrow(arle.df)) {
      arle.df$ref[mi]<-paste(refSubset[ arle.df$start.barcode[mi]:arle.df$end.barcode[mi]], collapse='')
      arle.df$alt[mi]<-paste(seqSubset[ arle.df$start.barcode[mi]:arle.df$end.barcode[mi]], collapse='')
    }
    arle.df$type <- "sub"
    arle.df[arle.df$lengths>1, "type"] <- "complex"
    # if alt sequence is made entirely of hyphens, it is a deletions
    arle.df[gsub("-", "", arle.df$alt)=="", "type"] <- "del"
    # there should not be any insertions here: they are dealt 
  }    
  # add the barcodep scale
  arle.df$start.barcodep <- barcodeToBarcodep[arle.df$start.barcode]
  arle.df$end.barcodep <- barcodeToBarcodep[arle.df$end.barcode]
  
  if ((nrow(arle.inf.df)>0) &&(nrow(arle.df)>0)) {
    all.muts <- rbind(arle.inf.df[,colnames(arle.df)], arle.df)  
  } else if ((nrow(arle.inf.df)>0)) {
    all.muts <- arle.inf.df
  }  else {
    all.muts <- arle.df
  }
  
  return(all.muts)
  
}

