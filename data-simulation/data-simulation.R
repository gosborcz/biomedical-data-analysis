
set.seed(20250507)          

n_per_group <- 30           
levels_sex  <- c("male", "female")
levels_age  <- c("young", "old")

expand_design <- function() {
  expand.grid(sex = levels_sex,
              age = levels_age,
              rep = seq_len(n_per_group),
              KEEP.OUT.ATTRS = FALSE,
              stringsAsFactors = FALSE)[, 1:2]
}


dat_age_only <- expand_design()

dat_age_only$behaviour <- with(dat_age_only,
                               rnorm(n      = nrow(dat_age_only),
                                     mean   = ifelse(age == "young", 8, 5),
                                     sd     = 3)
)


anova_age_only <- aov(behaviour ~ sex * age, data = dat_age_only)
summary(anova_age_only)

write.table(dat_age_only, row.names = FALSE, quote = FALSE, './two-way-examples/two-way-data-first-example.tsv')


dat_interaction <- expand_design()

dat_interaction$behaviour <- with(dat_interaction, {
  mu <- ifelse(sex == "female" & age == "young", 10, 7)
  rnorm(n = length(mu), mean = mu, sd = 2)
})


anova_interaction <- aov(behaviour ~ sex * age, data = dat_interaction)
summary(anova_interaction)

write.table(dat_interaction, row.names = FALSE, quote = FALSE, './two-way-examples/two-way-data-second-example.tsv')

