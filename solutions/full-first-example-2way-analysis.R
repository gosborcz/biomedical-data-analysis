library("ggpubr")
library("dplyr")
library("magrittr")

df <- read.table(
  file       = "./two-way-examples/two-way-data-first-example.tsv",
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

p_female <- get_p(pw_all$p.value, "female_old", "female_young")
p_male   <- get_p(pw_all$p.value, "male_old",   "male_young")


y_max  <- max(df$behaviour)
ann_df <- data.frame(
  sex        = c("female", "male"),
  group1     = "young",
  group2     = "old",
  p          = c(p_female, p_male),
  y.position = y_max + 0.5           
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
  x        = "age",
  y        = "behaviour",
  color    = "age",
  palette  = "jco",
  add      = "jitter",
  facet.by = "sex"
) +
  stat_pvalue_manual(
    data       = ann_df,
    label      = "p.signif",
    tip.length = 0,
    hide.ns    = FALSE
  ) +
  labs(
    title = "behaviour by age and sex",
    x     = "age group",
    y     = "behaviour"
  ) +
  theme_pubr()

print(p)


