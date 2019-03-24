#-----------------------------------------------------------------------------
# The FHS spectra was originally generated using a version of Doomsayer that 
# used a different subtype naming convention, so columns of the spectra matrix
# were in a different order than required for the current version. 
#
# This script updates the subtype names and rearranges the columns of this matrix
#-----------------------------------------------------------------------------

revcomp <- function(DNAstr) {
  step1 <- chartr("ACGT","TGCA",DNAstr)
  step2 <- unlist(strsplit(step1, split=""))
  step3 <- rev(step2)
  step4 <- paste(step3, collapse="")
  return(step4)
}

# Subtypes in format 
m_in <- read_tsv("/mnt/norbert/data/topmed/subtype_counts/ramachandran_NMF_M_spectra_no_ellinor.txt")

m_sort <- m_in %>% dplyr::select(c("ID", 
                                   names(m_in)[50:97], names(m_in)[49:2]))

# custom.sort <- function(x){x[order(as.numeric(substring(x, 7)))]}
m_sort2 <- m_sort %>% dplyr::select(c("ID",
                                      names(m_sort[2:49]),
                                      names(m_sort)[50:97][1:16][c(rep(0:3, each=4) + rep(seq(1,13,4), 4))],
                                      names(m_sort)[50:97][17:32][c(rep(0:3, each=4) + rep(seq(1,13,4), 4))],
                                      names(m_sort)[50:97][33:48][c(rep(0:3, each=4) + rep(seq(1,13,4), 4))]))

names(m_sort2)[50:97] <- paste0(rep(c("T_A", "T_C", "T_G"), each=16), ".", unlist(lapply(gsub(".*[.]", "", names(m_sort2)[50:97]), revcomp)))

names(m_sort2)[-1] == sort(names(m_sort2)[-1])


write_tsv(m_sort2, "/mnt/norbert/data/topmed/subtype_counts/ramachandran_NMF_M_spectra_no_ellinor_reform.txt")
