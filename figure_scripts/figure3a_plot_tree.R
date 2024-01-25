library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(tidyr)

##########
### Code to create draft figure of 24 subline mutation tree and growth rates.
### Modifications made in Adobe Illustrator.
##########

tree_file <- "tree_files/sc-bwes-cons-resolved-10.tree"
tree <- read.tree(tree_file)

growth_data <- read.csv("data/mouse_treatment/24_subclone_growth_dynamics.csv")
growth_data <- na.omit(gather(growth_data, key="subline", value="time"))

p <- ggtree(tree) + geom_tiplab(as_ylab = TRUE) + 
    geom_text(aes(x=branch, label=branch.length))
p <- p + geom_fruit(
    data = growth_data,
    geom = geom_boxplot,
    mapping = aes(
        x = time,
        y = subline
    ),
    axis.params=list(axis="xy", text.size=2)
)
p <- ggplotify::as.ggplot(p, angle=270)

svg("figures/figure1a_24subline_tree_code_output.svg")
p
dev.off()