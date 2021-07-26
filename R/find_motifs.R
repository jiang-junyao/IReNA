#' Find motifs in the footprints through fimo software.
#'
#' @param motif motif file, you can choose our bulit-in motif database of 'mus musculus', 'homo sapiens', 'zebrafish' and 'chicken' by 'motif = Tranfac201803_Mm_MotifTFsF', 'motif = Tranfac201803_Hs_MotifTFsF', 'motif = Tranfac201803_Zf_MotifTFsF', 'motif = Tranfac201803_Ch_MotifTFsF' respectively, or you can upload your own motif data base, but the formata use be the same as our built-in motif database.
#' @param step numeric, indicating the numbers of motif in each group. Because there are so many motifs need to be calculated, so we divided every 20 motifs into groups and then calculated each group using FIMO at the same time by nohup
#' @param Dir characater, indicating the directory of output files
#' @param sequence_dir sequence file directory
#'
#' @return return scripts to run Fimo in command line
#' @export
#'
#' @examples
find_motifs <- function(motif, step = 20, Dir, sequence_dir) {
  con1 <- motif
  no1 <- 1
  out2 <- c()
  for (i in seq(1, nrow(con1), by = step)) {
    out1 <- c()
    if (nrow(con1) - i > step) {
      for (j in i:(i + step - 1)) {
        col1 <- paste0(
          "fimo --parse-genomic-coord --max-stored-scores 2000000  --text >",
          Dir, "/mmATACPhxW_Fimo/", rownames(con1[j, ]), ".txt ",
          sequence_dir
        )
        out1 <- c(out1, col1)
      }
      out1 <- as.data.frame(out1)
      write.table(out1, paste0(Dir, "/Fimo", no1, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    } else {
      for (k in i:nrow(con1)) {
        col1 <- paste0(
          "fimo --parse-genomic-coord --max-stored-scores 2000000  --text >",
          Dir, "/mmATACPhxW_Fimo/", rownames(con1[k, ]), ".txt ",
          sequence_dir
        )
        out1 <- c(out1, col1)
      }
      out1 <- as.data.frame(out1)
      write.table(out1, paste0(Dir, "/Fimo", no1, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    # write.table(out1,paste0(Dir1,'/Programs/Fimo',no1,'.txt'),row.names=F,col.names = F,quote = F)
    no1 <- no1 + 1
    col2 <- paste0("nohup sh ", Dir, "/Fimo", no1, ".txt")
    out2 <- c(out2, col2)
  }
  out2 <- as.data.frame(out2)
  write.table(out2, paste0(Dir, "/Fimo", "_All", ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}
