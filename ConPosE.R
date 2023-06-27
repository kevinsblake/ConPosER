
conposer_id <- function(seqs.file, msa=c("ClustalW", "ClustalOmega", "Muscle"), gap.lim=0.30){
  
  # Import fasta file  
  seqs.fasta <- readAAStringSet(seqs.file)
  
  # Generate MSA
  seqs.aln <- msa(seqs.fasta, msa) # can also use "Muscle" etc
  
  # Create empty dataframe template
  rows = c("-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  df <- data.frame(matrix(nrow=length(rows), ncol=1))
  rownames(df) = rows
  
  # Make consensus matrix + sum num AAs at each position (except gap)
  cons.mat <- consensusMatrix(seqs.aln)
  cons.mat_v2 <- merge(df, cons.mat, by="row.names", all=TRUE)
  rownames(cons.mat_v2) <- cons.mat_v2[,1]
  cons.mat_v2[,1:2] <- NULL
  cons.mat_v2[is.na(cons.mat_v2)] <- 0
  cons.mat2 <- as.data.frame(t(rbind(cons.mat_v2, colSums(cons.mat_v2 != 0)))) 
  rownames(cons.mat2) <- NULL
  
  cons.mat2 <- cons.mat2 %>% rename(gap = "-") %>%
    tibble::rownames_to_column("pos")
  
  # Filter output
  cons.mat.fil <- cons.mat2 %>% filter(gap < length(seqs.fasta)) %>%
    select(-pos) %>%
    tibble::rownames_to_column("pos")
  cons.fil.cons <- cons.mat.fil %>% filter(cons.mat.fil[,ncol(cons.mat.fil)] == 1 & gap == 0) %>% select(-gap) 
  
  # Print AA name at conserved position
  cons.fil.cons$AA <- names(cons.fil.cons)[-1][max.col(cons.fil.cons[-1] !=0, 'first')]
  cons.find <- cons.fil.cons[, c('pos', 'AA')]
  
  return(cons.find)
}


conposer_plot <- function(seqs.file, filename="geneplot.pdf", linecol="black", gap.lim=0.30){
  
  # Import fasta file
  seqs.fasta <- readAAStringSet(seqs.file)
  
  # Generate MSA
  seqs.aln <- msa(seqs.fasta, "ClustalOmega") # can also use "Muscle" etc
  
  # Create empty dataframe template
  rows = c("-","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  df <- data.frame(matrix(nrow=length(rows), ncol=1))
  rownames(df) = rows
  
  # Make consensus matrix + sum num AAs at each position (except gap) - WORKAROUND FOR IF NO GAPS
  cons.mat <- consensusMatrix(seqs.aln)
  cons.mat_v2 <- merge(df, cons.mat, by="row.names", all=TRUE)
  rownames(cons.mat_v2) <- cons.mat_v2[,1]
  cons.mat_v2[,1:2] <- NULL
  cons.mat_v2[is.na(cons.mat_v2)] <- 0
  cons.mat2 <- as.data.frame(t(rbind(cons.mat_v2, colSums(cons.mat_v2 != 0)))) 
  
  cons.mat2 <- cons.mat2 %>% rename(gap = "-") %>%
    tibble::rownames_to_column("pos")
  
  # Filter output  
  cons.mat.fil <- cons.mat2 %>% filter(gap < length(seqs.fasta)*gap.lim) %>% # Filter to exclude positions where >30% sequences have a gap
    select(-pos) %>%
    tibble::rownames_to_column("pos")
  cons.fil.cons <- cons.mat.fil %>% filter(cons.mat.fil[,ncol(cons.mat.fil)] == 1 & gap == 0) %>% select(-gap) 
  
  # Save the plot as a function so can save it >> little hacky
  plot <- function(){conposer_geneplot(cons.mat.fil, cons.fil.cons, filename, linecol)}
  
  return(plot())
}

conposer_geneplot <- function(all, cons, filename, linecol){
  pdf(file=filename, height=2.5, width=8)
  par(mfrow=c(2,1), oma= c(2,2,0,0) + 0.1, mar = c(0, 0, 0.1, 0.1) + 0.1) 
  barplot(all[,ncol(all)], xaxs="i", yaxs="i", xlab="", ylab="count", col="black")
  plot(c(0, nrow(all)), c(0,0), type="n", xlab="", ylab="", yaxt="n", yaxs="i", xaxs="i") + 
    axis(1, lwd.ticks=2) +
    abline(v=cons$pos, col=linecol, lwd=2) +
    box(lwd=2)
  
  dev.off()
}
