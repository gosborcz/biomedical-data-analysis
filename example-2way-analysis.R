library("ggpubr")
library("dplyr")
library("magrittr")


df <- read.table(
  file       = "./two-way-examples/two-way-data-first-example.tsv",
  header     = TRUE,
  sep        = " ",              
  stringsAsFactors = FALSE
)

#run anova
anova_results <- aov(behaviour ~ sex * age, data = df)

print(summary(anova_results))



df$group <- interaction(df$sex, df$age, sep = "_") #this is a factor which we want for the t.test

# (A) run pairwise.t.test for posthoc comparisons; use the data in the df, experiment with the different p-value adjustment methods
# hint - google how to get help on a function
# save your result as a variable named pw_all

# (B) inspect how pw_all looks like, what is its class?
# the function below is a helper function to extract p-values

get_p <- function(mat, g1, g2) {
  rn <- rownames(mat)
  cn <- colnames(mat)
  if (g1 %in% rn && g2 %in% cn)
    return(mat[g1, g2])
  if (g2 %in% rn && g1 %in% cn)
    return(mat[g2, g1])
  stop(sprintf("not found", g1, g2)) # (C) what dose this line do?
}

p_female <- get_p(pw_all$p.value, "female_old", "female_young")
p_male   <- get_p(pw_all$p.value, "male_old",   "male_young")


y_max  <- max(df$behaviour) # (D) we will need this for plotting. why? (you can google the documentation for stat_pvalue_manual in ggpubr)

ann_df <- data.frame(
  sex        = c("female", "male"),
  group1     = "young",
  group2     = "old",
  p          = c(p_female, p_male),
  y.position = y_max + 0.5           
) 

# (E) inspect ann_df, ad a column p.signif that will have stars instead of p-values (hint mutate the p-value column)

p <- ggboxplot(
  df,
  x        = "age",
  y        = "behaviour",
  color    = "age",
  palette  = "jco", # (F) can you find other possible palletes? can you define your own palette?
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


