library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(tidyr)
library(dplyr)

##########
### Code to create draft figure of 24 subline mutation tree and growth rates.
### Modifications made in Adobe Illustrator.
##########

tree_file <- "tree_files/sc-bwes-cons-unresolved-10.tree"
tree <- read.tree(tree_file)

growth_data <- read.csv("data/mouse_treatment/24_subclone_growth_dynamics.csv")
growth_data <- na.omit(gather(growth_data, key="subline", value="time"))

pcna_data <- data.frame()
for (file in list.files(path = "data/tpm/", full.names = TRUE)) {
    temp <- read.csv(file)
    temp <- filter(temp, temp$gene == "ENSMUSG00000027342.14_Pcna")
    pcna_data <- rbind(pcna_data, temp)
}
m <- mean(pcna_data$exprval)
s <- sd(pcna_data$exprval)
pcna_data <- pcna_data %>% mutate(zscore = (exprval - m) / s)

pcna_data <- pcna_data %>% mutate(logfc = log2(exprval/mean(exprval)))
tree <- ape::rotateConstr(tree, c("C7", "C20", "C8", "C16", "C11", "C15", "C18", "C13", "C6", "C24", "C19", "C9", "C21", "C4", "C22", "C1", "C12", "C14", "C3", "C10", "C5", "C17", "C23"))
colors <- c("white", "white", "white", "cadetblue2", "cadetblue2", "cadetblue2", "cadetblue2", "white", "white", "white", "white", "white", "white", "orange", "orange", "orange", "white", "bisque1", "bisque1", "bisque1", "white", "white", "white")

p <- ggtree(tree, ladderize=FALSE) + geom_tiplab(as_ylab = TRUE) + 
    geom_text(aes(x=branch, label=branch.length, angle=90))
p <- p + geom_fruit(
    data = growth_data,
    geom = geom_boxplot,
    mapping = aes(
        x = time,
        y = subline
    ),
    axis.params=list(axis="xy", text.size=2)
)
p <- p + geom_fruit(
    data = pcna_data,
    geom = geom_boxplot,
    offset = .35,
    mapping = aes(
        x = zscore,
        y = species
    ),
    axis.params=list(axis="xy", text.size=2)
)

p <- ggplotify::as.ggplot(p, angle=270)

svg("figures/figure2b_24subline_tree_code_output.svg")
p
dev.off()