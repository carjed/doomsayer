samples_1kg <- read_tsv("/mnt/norbert/data/1kg/20130606_sample_info.txt")
comp_1kg <- read_tsv("/mnt/norbert/home/jedidiah/projects/doomsayer/demo/1kg_pca/W_components.txt")
outliers_1kg <- read_tsv("/mnt/norbert/home/jedidiah/projects/doomsayer/demo/1kg_pca/doomsayer_drop.txt", col_names=F)

# Get table of family/population IDs for outlier samples
outliers_1kg_desc <- merge(outliers_1kg, samples_1kg, by.x="X1", by.y="Sample") %>%
  dplyr::select(X1, `Family ID`, Population, `DNA Source from Coriell`)

samples_1kg_2 <- samples_1kg %>% 
  dplyr::select(ID=Sample, Population, `EBV Coverage`, `DNA Source from Coriell`)

comp_1kg_merge <- merge(comp_1kg, samples_1kg_2, by="ID")

# Plot PC2 vs PC3, with points colored by DNA source
comp_1kg_merge %>% 
  mutate(outlier=ifelse(ID %in% outliers_1kg$X1, T, F)) %>%
  mutate(`outlier cluster`=ifelse(ID %in% c("HG00186", "NA18498", "HG00373", "HG01149"), "2|3 (C>G-enriched)", ifelse(outlier==TRUE, "1 (T>G-enriched)", "none"))) %>%
  mutate(`DNA Source from Coriell`=ifelse(is.na(`DNA Source from Coriell`), "Unknown", `DNA Source from Coriell`)) %>%
  # dplyr::filter(!is.na(`DNA Source from Coriell`)) %>%
  dplyr::filter(S2<15) %>%
  ggplot(aes(x=S2, y=S3, colour=`DNA Source from Coriell`, shape=`outlier cluster`, size=`outlier cluster`))+
    geom_point(alpha=0.5)+
    scale_shape_manual(values=c(15,17,1))+
    scale_colour_manual(values=c("red", "blue", "grey60"))+
    scale_size_manual(values=c(4,4,2))+
    xlab("PC2")+
    ylab("PC3")+
    theme_bw()+
    theme(axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16),
          legend.text=element_text(size=12),
          legend.title=element_text(size=16),
          legend.position="bottom",
          legend.box = "vertical")

ggsave("ERV_qc/figs/1kg_dna_source.png", width=8, height=8)

c1kg2<-comp_1kg_merge[,1:4]


