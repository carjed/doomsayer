#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------

# Read signatures output by Doomsayer run
fhs_sigs <- read_tsv("/mnt/norbert/home/jedidiah/projects/doomsayer/demo/ramachandran_pca_filter/W_components.txt")
fhs_spectra <- read_tsv("/mnt/norbert/home/jedidiah/projects/doomsayer/demo/ramachandran_pca_filter/subtype_count_matrix.txt")

# List of 7 FHS outliers
fhs_outliers <- c("NWD219552", "NWD431332", "NWD730596", "NWD968087", "NWD952803", "NWD686119", "NWD270841")



fhs_source <- read_csv("/mnt/norbert/data/topmed/fram_dna_source.csv")

fhs_merge <- merge(fhs_sigs, fhs_source, by.x="ID", by.y="NWD_ID")

topmed_fams <- read_tsv("/mnt/norbert/data/topmed/freeze5.64960.4606.ped", col_names=F)

# get ids of offspring 
fhs_probands <- topmed_fams %>% 
  dplyr::filter(X3 %in% fhs_outliers | X4 %in% fhs_outliers) %>% 
  dplyr::select(ID=X2)

# get ids of non-outlier parents
fhs_inliers <- topmed_fams %>% 
  dplyr::filter(X3 %in% fhs_outliers | X4 %in% fhs_outliers) %>% 
  dplyr::select(X3, X4) %>% 
  gather(X3:X4) %>% 
  dplyr::filter(!(value %in% fhs_outliers)) %>% 
  dplyr::select(ID=value)

# 1-mer spectra
fhs_outlier_1mer_spectra <- fhs_spectra %>%
  dplyr::filter(ID %in% fhs_outliers) %>%
  gather(subtype, count, 2:97) %>%
  separate(subtype, c("category", "motif"), sep = "[.]") %>%
  mutate(subtype=paste0(substr(motif,1,1), "[", category, "]", substr(motif,3,3))) %>%
  mutate(subtype=gsub("_", ">", subtype)) %>%
  mutate(category=gsub("_", ">", category)) %>%
  # mutate(category=factor(category, levels=orderedcatsnc_short)) %>%
  arrange(category) %>%
  group_by(category) %>%
  summarise(n=sum(count)) %>%
  mutate(prop=n/sum(n))

fhs_inlier_1mer_spectra <- fhs_spectra %>%
  dplyr::filter(ID %in% fhs_inliers$ID) %>%
  gather(subtype, count, 2:97) %>%
  separate(subtype, c("category", "motif"), sep = "[.]") %>%
  mutate(subtype=paste0(substr(motif,1,1), "[", category, "]", substr(motif,3,3))) %>%
  mutate(subtype=gsub("_", ">", subtype)) %>%
  mutate(category=gsub("_", ">", category)) %>%
  # mutate(category=factor(category, levels=orderedcatsnc_short)) %>%
  arrange(category) %>%
  group_by(category) %>%
  summarise(n=sum(count)) %>%
  mutate(prop=n/sum(n))

# 3-mer spectra
fhs_outlier_3mer_spectra <- fhs_spectra_outliers %>%
  gather(subtype, count, 2:97) %>%
  separate(subtype, c("category", "motif"), sep = "[.]") %>%
  mutate(subtype=paste0(substr(motif,1,1), "[", category, "]", substr(motif,3,3))) %>%
  mutate(subtype=gsub("_", ">", subtype)) %>%
  mutate(category=gsub("_", ">", category)) %>%
  # mutate(category=factor(category, levels=orderedcatsnc_short)) %>%
  arrange(category) %>%
  group_by(category, motif) %>%
  summarise(n=sum(count)) %>%
  ungroup() %>%
  mutate(prop=n/sum(n))

fhs_inlier_3mer_spectra <- fhs_spectra %>%
  dplyr::filter(ID %in% fhs_inliers$ID) %>%
  gather(subtype, count, 2:97) %>%
  separate(subtype, c("category", "motif"), sep = "[.]") %>%
  mutate(subtype=paste0(substr(motif,1,1), "[", category, "]", substr(motif,3,3))) %>%
  mutate(subtype=gsub("_", ">", subtype)) %>%
  mutate(category=gsub("_", ">", category)) %>%
  # mutate(category=factor(category, levels=orderedcatsnc_short)) %>%
  arrange(category) %>%
  group_by(category, motif) %>%
  summarise(n=sum(count)) %>%
  ungroup() %>%
  mutate(prop=n/sum(n))


transmitted <- read_tsv("/mnt/norbert/data/topmed/fhs_outliers_subtype_count_matrix.txt")

transmitted_1mer_spectra <- transmitted %>%
  dplyr::filter(!(ID %in% fhs_outliers)) %>%
  gather(subtype, count, 2:97) %>%
  separate(subtype, c("category", "motif"), sep = "[.]") %>%
  mutate(subtype=paste0(substr(motif,1,1), "[", category, "]", substr(motif,3,3))) %>%
  mutate(subtype=gsub("_", ">", subtype)) %>%
  mutate(category=gsub("_", ">", category)) %>%
  # mutate(category=factor(category, levels=orderedcatsnc_short)) %>%
  arrange(category) %>%
  group_by(category) %>%
  summarise(n=sum(count)) %>%
  mutate(prop=n/sum(n))

transmitted_3mer_spectra <- transmitted %>%
  dplyr::filter(!(ID %in% fhs_outliers)) %>%
  gather(subtype, count, 2:97) %>%
  separate(subtype, c("category", "motif"), sep = "[.]") %>%
  mutate(subtype=paste0(substr(motif,1,1), "[", category, "]", substr(motif,3,3))) %>%
  mutate(subtype=gsub("_", ">", subtype)) %>%
  mutate(category=gsub("_", ">", category)) %>%
  # mutate(category=factor(category, levels=orderedcatsnc_short)) %>%
  arrange(category) %>%
  group_by(category, motif) %>%
  summarise(n=sum(count)) %>%
  ungroup() %>%
  mutate(prop=n/sum(n))

merge(fhs_outlier_3mer_spectra, transmitted_3mer_spectra, by=c("category", "motif")) %>%
  mutate(d1=paste0(n.x, " (", round(prop.x*100, 2), "%)")) %>%
  mutate(d2=paste0(n.y, " (", round(prop.y*100, 2), "%)")) %>%
  dplyr::select(category, motif, d1, d2)


#-----------------------------------------------------------------------------
# DEPRECATED
# compare spectra of outlier/inlier parents with that of de novos
#-----------------------------------------------------------------------------

# combined parent spectra
fhs_parent_spectra <- bind_rows(fhs_inlier_1mer_spectra, fhs_outlier_1mer_spectra) %>% 
  group_by(category) %>%
  summarise(n=sum(n)) %>%
  mutate(prop=n/sum(n))

data_path <- "/mnt/norbert/data/topmed/dnms/new_sorted/"   # path to the data
dnm_files <- dir(data_path, pattern = "*.txt") # get file names

getChrNum <- function(x){
  chr <- gsub("\\..*", "", x)
  chr <- as.numeric(gsub("chr", "", chr))
}

dnm_files <- data.frame(name = dnm_files) %>%
  mutate(chr=getChrNum(name)) %>%
  arrange(chr)

dnms <- dnm_files$name %>%
  # read in all the files, appending the path before the filename
  map(~ read_tsv(file.path(data_path, .), col_names=TRUE)) %>%
  reduce(rbind)

dnms <- dnms %>%
  dplyr::filter(ID %in% fhs_probands$ID)

dnms <- dnms %>% 
  dplyr::select(TYPE, MOTIF, ID) %>%
  mutate(MOTIF=substr(MOTIF,3,5))

dnms %>% 
  group_by(TYPE) %>% 
  count() %>%
  ungroup() %>%
  mutate(prop=n/sum(n))

dnm_spectra <- dnms %>% 
  group_by(TYPE) %>% 
  count() %>%
  ungroup() %>%
  mutate(prop=n/sum(n)) %>%
  mutate(TYPE=gsub("_", ">", TYPE)) %>%
  mutate(TYPE=fct_recode(TYPE, "T>G" = "A>C", "T>C" = "A>G", "T>A" = "A>T")) %>%
  mutate(TYPE=as.character(TYPE)) %>%
  arrange(TYPE)

dnm_3mer_spectra <- dnms %>% 
  group_by(TYPE, MOTIF) %>% 
  count() %>%
  ungroup() %>%
  mutate(prop=n/sum(n)) %>%
  mutate(TYPE=gsub("_", ">", TYPE)) %>%
  mutate(TYPE=fct_recode(TYPE, "T>G" = "A>C", "T>C" = "A>G", "T>A" = "A>T")) %>%
  mutate(TYPE=as.character(TYPE)) %>%
  arrange(TYPE)
