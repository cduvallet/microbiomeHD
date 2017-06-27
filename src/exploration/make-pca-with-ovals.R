library(RColorBrewer)
library(ggplot2)

## Make PCA of just healthy with ovals
wdir <- '/Users/claire/github/microbiomeHD/src/exploration'
setwd(wdir)

abun.fn <- 'mds.healthy_samples.rel_abun.csv'
combat.fn <- 'mds.healthy_samples.combat_rel_abun.csv'

abun <- read.table(abun.fn, header=TRUE, sep=',')
combat <- read.table(combat.fn, header=TRUE, sep=',')

# Try using manual color palette from this stackoverflow: http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Plot all points
p <- ggplot(abun, aes(x=PC1, y=PC2, color=dataset)) +
     geom_point(size=1) +
     stat_ellipse(aes(x=PC1, y=PC2, color=dataset), inherit.aes=FALSE, type='norm', linetype='dotted') +
     theme_bw() + scale_colour_manual(values=col_vector)

# Subsample 35 points and re-plot. 
# Note: studies with fewer than 35 controls will have duplicate points...
subabun <- abun %>% group_by(dataset) %>% sample_n(35, replace=TRUE)
p2 <- ggplot(subabun, aes(x=PC1, y=PC2, color=dataset)) +
      geom_point(size=2) +
      #stat_ellipse(aes(x=PC1, y=PC2, color=dataset), inherit.aes=FALSE, type='norm', linetype='dotted') +
      theme_bw() + scale_colour_manual(values=col_vector)

ggsave("mds_abundance_subsampled35.png", plot=p2)

## Combated data
subcombat <- combat %>% group_by(dataset) %>% sample_n(35, replace=TRUE)
p3 <- ggplot(subcombat, aes(x=PC1, y=PC2, color=dataset)) +
      geom_point(size=2) +
      #stat_ellipse(aes(x=PC1, y=PC2, color=dataset), inherit.aes=FALSE, type='norm', linetype='dotted') +
      theme_bw() + scale_colour_manual(values=col_vector)

ggsave("mds_combat_subsampled35.png", plot=p3)

