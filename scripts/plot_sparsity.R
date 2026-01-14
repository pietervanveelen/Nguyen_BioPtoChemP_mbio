# plot sparsity
require(tidyverse)
require(phyloseq)

# sample_stats function
sample_stats<-function (physeq,...){
  m <- otu_table(physeq)
  marg<-ifelse(taxa_are_rows(m), yes = 2, no = 1)
  return(apply(m, marg, ...))
}

# sample_prev function
sample_prev<-function(physeq){
  m <- otu_table(physeq)
  if (taxa_are_rows(m)) {
    colSums(m > 0)
  }
  else {
    rowSums(m > 0)
  }
}

# taxa_stats function
taxa_stats<-function (physeq,...){
  m <- otu_table(physeq)
  marg<-ifelse(taxa_are_rows(m), yes = 1, no = 2)
  return(apply(m, marg, ...))
}

# taxa_prev function
taxa_prev<-function(physeq){
  m <- otu_table(physeq)
  if (taxa_are_rows(m)) {
    rowSums(m > 0)
  }
  else {
    colSums(m > 0)
  }
}

# plot_sparsity function
plot_sparsity<-function(physeq, title = NULL){
  # Melt the OTU table and merge associated metadata
  df<-psmelt(physeq) %>% as_tibble() %>% select(OTU, Sample, Abundance) # keep OTU, Sample, Abundance columns
  # Order OTU by their occurrence/prevalence
  df$OTU<-factor(x = df$OTU,
                 levels = names(sort(taxa_prev(physeq))) )
  # Order Samples by their richness (or completedness)
  df$Sample<-factor(x = df$Sample,
                    levels = names(sort(sample_prev(physeq))) )
  # Transform Abundance numeric vector into boolean presence/absence
  df$Abundance<-df$Abundance > 0
  # Plot OTU table
  p<-ggplot(df[df$Abundance,],
            aes(x=Sample, y=OTU))+
    geom_tile()+
    scale_x_discrete(drop = F) + 
    scale_y_discrete(drop=F) +
    theme_minimal()+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
          )+
    labs(x= paste(length(levels(df$Sample)),"samples - from low to high richness"),
         y= paste(length(levels(df$OTU)),"OTU -- from rare to dominant"))
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

