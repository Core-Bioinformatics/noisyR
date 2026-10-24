## Figures from the CSVs written by run_benchmark.R and run_scaling.R.
## Run from the repository root: Rscript analysis/grn_noise/plot_results.R
library(ggplot2)
res.dir <- "analysis/grn_noise/results"
pal <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
theme_set(theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
                  legend.position = "top"))

bench <- read.csv(file.path(res.dir, "benchmark.csv"))
prs <- read.csv(file.path(res.dir, "pr_curves.csv"))
G.show <- max(bench$G[bench$G <= 2000])

## 1. PR curves of the full-data and consensus rankings, with each method's operating point
pr <- prs[prs$G == G.show & prs$inference != "direct" &
          (is.na(prs$scheme) | prs$scheme == "split"), ]
pr$panel <- paste(pr$inference, "network /", pr$truth, "truth")
pts <- bench[bench$G == G.show & bench$inference != "direct" &
             (is.na(bench$scheme) | bench$scheme == "split") &
             !grepl("^oracle", bench$method), ]
pts$panel <- paste(pts$inference, "network /", pts$truth, "truth")
p1 <- ggplot(pr, aes(recall, precision)) +
  geom_line(aes(linetype = ranking), colour = "grey45", linewidth = 0.6) +
  geom_point(data = pts, aes(colour = method, shape = method), size = 3, stroke = 1) +
  facet_wrap(~ panel) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = c(16, 17, 15, 18, 1, 2)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = paste0("Edge recovery, G = ", G.show, " genes, 5 disjoint cell splits"),
       linetype = "ranking", colour = NULL, shape = NULL) +
  guides(colour = guide_legend(nrow = 2))
ggsave(file.path(res.dir, "pr_operating_points.png"), p1, width = 10, height = 4.8, dpi = 150)

## 2. precision / recall / F1 against the similarity threshold (concordance, split)
sw <- read.csv(file.path(res.dir, "threshold_sweep.csv"))
sw <- sw[sw$method == "noisyR-concordance" & sw$scheme %in% c("split", "independent"), ]
sw$panel <- paste0(sw$inference, " / ", sw$truth, ", G = ", sw$G)
long <- rbind(data.frame(sw[, c("panel", "similarity.threshold")], metric = "precision", value = sw$precision),
              data.frame(sw[, c("panel", "similarity.threshold")], metric = "recall", value = sw$recall),
              data.frame(sw[, c("panel", "similarity.threshold")], metric = "F1", value = sw$F1))
p2 <- ggplot(long, aes(similarity.threshold, value, colour = metric)) +
  geom_vline(xintercept = 0.3, colour = "grey60", linetype = 2) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ panel, ncol = 4) +
  scale_colour_manual(values = pal[1:3]) +
  labs(title = "Concordance similarity: edge recovery vs similarity threshold (dashed: default 0.3)",
       x = "similarity threshold", y = NULL, colour = NULL)
ggsave(file.path(res.dir, "threshold_sweep.png"), p2, width = 11, height = 7, dpi = 150)

## 3. runtime by stage and number of genes
if(file.exists(file.path(res.dir, "scaling.csv"))){
  sc <- read.csv(file.path(res.dir, "scaling.csv"))
  sc <- sc[sc$stage != "prepare_edges, top.m = 10 G", ]
  p3 <- ggplot(sc, aes(G, pmax(seconds, 0.01), colour = stage, shape = stage)) +
    geom_line(linewidth = 0.6) + geom_point(size = 2.5) +
    facet_wrap(~ inference) +
    scale_x_log10(breaks = unique(sc$G)) + scale_y_log10() +
    scale_colour_manual(values = c(pal, "grey40", "black")) +
    scale_shape_manual(values = c(16, 17, 15, 18, 1, 2, 0, 5)) +
    labs(title = "Runtime by stage (5 replicates, 1000 cells)", x = "genes (G)",
         y = "seconds (log scale)", colour = NULL, shape = NULL) +
    theme(legend.position = "right")
  ggsave(file.path(res.dir, "scaling.png"), p3, width = 11, height = 4.5, dpi = 150)
}
