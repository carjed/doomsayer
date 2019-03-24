bridges_pca <- read_tsv("/mnt/norbert/home/jedidiah/projects/doomsayer/demo/bridges_pca_3_05/doomsayer_drop.txt", col_names=F)
bridges_pca4 <- read_tsv("/mnt/norbert/home/jedidiah/projects/doomsayer/demo/bridges_pca_4_05/doomsayer_drop.txt", col_names=F)
bridges_pca5 <- read_tsv("/mnt/norbert/home/jedidiah/projects/doomsayer/demo/bridges_pca_5_05/doomsayer_drop.txt", col_names=F)
bridges_nmf <- read_tsv("/mnt/norbert/home/jedidiah/projects/doomsayer/demo/bridges_nmf_3_05/doomsayer_drop.txt", col_names=F)

bridges_ped <- read_tsv("/mnt/norbert/data/bridges/bridges.ped")
cluster_mems <- read_tsv("/mnt/norbert/home/jedidiah/projects/doomsayer/demo/bridges_pca_3_05/R/cluster_mems.txt")
bridges_ped_ol <- merge(cluster_mems, bridges_ped, by="ID")

cnv_outliers <- read_tsv("/mnt/norbert/data/bridges/bridges_cnv_LQfilter.eigen.outliers", col_names=F)

bridges_pca %>% dplyr::filter(!(X1 %in% bridges_pca4$X1))
