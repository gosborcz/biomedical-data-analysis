library("ggpubr")
library("dplyr")
library("magrittr")

df <- read.table(
  file       = "./two-way-examples/two-way-data-second-example.tsv",
  header     = TRUE,
  sep        = " ",              
  stringsAsFactors = FALSE
)

anova_results <- aov(behaviour ~ sex * age, data = df)

print(summary(anova_results))

df$group <- interaction(df$sex, df$age, sep = "_")

pw_all <- pairwise.t.test(df$behaviour, df$group, p.adjust.method = "fdr")

get_p <- function(mat, g1, g2) {
  rn <- rownames(mat);  cn <- colnames(mat)
  if (g1 %in% rn && g2 %in% cn)       return(mat[g1, g2])
  if (g2 %in% rn && g1 %in% cn)       return(mat[g2, g1])
  stop(sprintf("not found", g1, g2))
}

p_fem_vs_fem  <- get_p(pw_all$p.value, "female_old",  "female_young")
p_fem_vs_male <- get_p(pw_all$p.value, "female_young",    "male_old")
p_male_vs_fem <- get_p(pw_all$p.value, "male_young",  "female_young")


y_max  <- max(df$behaviour)
ann_df <- data.frame(
  group1     = rep("female_young", 3),
  group2     = c("female_old", "male_old", "male_young"),
  p          = c(p_fem_vs_fem, p_fem_vs_male, p_male_vs_fem),
  y.position = y_max + c(0.5, 1.0, 1.5)          
) %>% 
  dplyr::mutate(
    p.signif = dplyr::case_when(
      p <= 0.0001 ~ "****",
      p <= 0.001  ~ "***",
      p <= 0.01   ~ "**",
      p <= 0.05   ~ "*",
      TRUE        ~ "ns"
    )
  )


p <- ggboxplot(
  df,
  x        = "group",
  y        = "behaviour",
  color    = "group",
  palette  = "jco",
  add      = "jitter"
) +
  stat_pvalue_manual(
    data       = ann_df,
    label      = "p.signif",
    tip.length = 0,
    hide.ns    = FALSE
  ) +
  labs(
    title = "behaviour by sex and age group ",
    x     = "sex_age group",
    y     = "behaviour"
  ) +
  theme_pubr()

print(p)



