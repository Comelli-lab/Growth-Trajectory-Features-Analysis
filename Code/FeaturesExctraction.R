# Author: Paraskevi Massara
# Created: April 4th, 2019
# First update: Jan 20th, 2020
# Second major update:  February 8th, 2021
#
# ========= Table of contents ==================================================
#TOC> ==========================================================================
#TOC>
#TOC> Section    Title                                            Line
#TOC> --------------------------------------------------------------------------
#TOC>   1        Dependencies
#TOC>   2        Catch22
#TOC>   3        Clinical Features
#TOC> ==========================================================================

# ====  1. Dependencies ========================================================


# Packages' installation

if (!require(dplyr, quietly=TRUE)) {
  install.packages("dplyr")
  library(dplyr)
}


# ========== 2. Clinical features ==============================================

# ---- Measure: zlen -------

m <- "zlen"
print(m)
wide <- pelotas_zlen_wide
subject = as.numeric(wide[, 1])
ageinmonths <- colnames(wide)[-1]
colnames(wide)[1] <- "subject"
wide = wide[, -c(1)]

# Intercepts
intercepts <- data.frame(subject, stringsAsFactors = F)
intercepts <-
  mutate(intercepts, intercepts_lm_zlen = numeric(length(subject)))
intercepts <-
  mutate(intercepts, intercepts_q_zlen = numeric(length(subject)))

# Slopes
slopes <- data.frame(subject, stringsAsFactors = FALSE)
slopes <- mutate(slopes, slopes_lm_zlen = numeric(length(subject)))
slopes <- mutate(slopes, slopes_q_zlen = numeric(length(subject)))

# AUC
auc <- data.frame(subject, stringsAsFactors = FALSE)
auc <- mutate(auc, auc_lm_zlen = numeric(length(subject)))
auc <- mutate(auc, auc_q_zlen = numeric(length(subject)))

# Tempo
tempo <- data.frame(subject, stringsAsFactors = FALSE)
tempo <- mutate(tempo, tempo_lm_zlen = numeric(length(subject)))
tempo <- mutate(tempo, tempo_q_zlen = numeric(length(subject)))

# Conditional growth is the difference in RSS between the child's linear model and the population's
# linear model

lm_total<-lm(zlen~ageinmonths, data=pelotas)
pop_rss_zlen<-anova(lm_total)["Residuals", "Sum Sq"]

q_total<-lm(zlen~poly(ageinmonths,2), data=pelotas)
pop_rss_zlen_q<-anova(q_total)["Residuals", "Sum Sq"]

conditional <- data.frame(subject, stringsAsFactors = FALSE)
conditional <- mutate(conditional, conditional_lm_zlen = numeric(length(subject)))
conditional <- mutate(conditional, conditional_q_zlen = numeric(length(subject)))

for (i in 1:nrow(wide)) {
  print(i)
  row.df <- cbind(as.numeric(ageinmonths), as.numeric(t(wide)[, i]))
  row.df <- as.data.frame(row.df)
  colnames(row.df) <- c("ageinmonths", m)
  row.df <- na.omit(row.df)
  
# Linear
  if (nrow(row.df) > 1) {
    row_lm = lm(row.df[[m]] ~ row.df$ageinmonths)
    conditional[which(conditional$subject == subject[i]),]$conditional_lm_zlen <-
        abs(pop_rss_zlen - anova(row_lm)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_zlen <-
      row_lm$coefficients[1]
    slopes[which(slopes$subject == subject[i]), ]$slopes_lm_zlen <-
      row_lm$coefficients[2]
    integrand <-
      function(x) {
        row_lm$coefficients[2] * x + row_lm$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_lm_zlen <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_lm_zlen <-
      auc[which(auc$subject == subject[i]), ]$auc_lm_zlen / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_lm_zlen <- NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_zlen <-
  #     NA
  #   slopes[which(slopes$subject == subject[i]), ]$slopes_lm_zlen <- NA
  #   auc[which(auc$subject == subject[i]), ]$auc_lm_zlen <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_lm_zlen <- NA
  # }

# Quadratic
  if (nrow(row.df) > 2) {
    row_q = lm(row.df[[m]] ~ poly(row.df$ageinmonths, 2))
    conditional[which(conditional$subject == subject[i]),]$conditional_q_zlen <-
        abs(pop_rss_zlen_q - anova(row_lm)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_zlen <-
      row_q$coefficients[1]
    slopes[which(slopes$subject == subject[i]), ]$slopes_q_zlen <-
      row_q$coefficients[3]
    integrand <-
      function(x) {
        row_q$coefficients[3] * I(x) ^ 2 + row_q$coefficients[2] * I(x) + row_q$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_q_zlen <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_q_zlen <-
      auc[which(auc$subject == subject[i]), ]$auc_q_zlen / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_q_zlen <- NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_zlen <-
  #     NA
  #   slopes[which(slopes$subject == subject[i]), ]$slopes_q_zlen <- NA
  #   auc[which(auc$subject == subject[i]), ]$auc_q_zlen <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_q_zlen <- NA
  # }
}

pelotas_features <-
  full_join(pelotas_features, intercepts, by = "subject")
pelotas_features <- full_join(pelotas_features, slopes, by = "subject")
pelotas_features <- full_join(pelotas_features, auc, by = "subject")
pelotas_features <- full_join(pelotas_features, tempo, by = "subject")
pelotas_features <- full_join(pelotas_features, conditional, by = "subject")

# ---- Measure: Height -------
m <- "height"
print(m)
wide <- pelotas_height_wide
subject = as.numeric(wide[, 1])
ageinmonths <- colnames(wide)[-1]
colnames(wide)[1] <- "subject"
wide = wide[, -c(1)]

# Intercepts
intercepts <- data.frame(subject, stringsAsFactors = F)
intercepts <-
  mutate(intercepts, intercepts_lm_height = numeric(length(subject)))
intercepts <-
  mutate(intercepts, intercepts_q_height = numeric(length(subject)))

# Slopes
slopes <- data.frame(subject, stringsAsFactors = FALSE)
slopes <- mutate(slopes, slopes_lm_height = numeric(length(subject)))
slopes <- mutate(slopes, slopes_q_height = numeric(length(subject)))

# AUC
auc <- data.frame(subject, stringsAsFactors = FALSE)
auc <- mutate(auc, auc_lm_height = numeric(length(subject)))
auc <- mutate(auc, auc_q_height = numeric(length(subject)))

# Tempo
tempo <- data.frame(subject, stringsAsFactors = FALSE)
tempo <- mutate(tempo, tempo_lm_height = numeric(length(subject)))
tempo <- mutate(tempo, tempo_q_height = numeric(length(subject)))

# Conditional growth is the difference in RSS between the child's linear model and the population's
# linear model

lm_total<-lm(raw_height~ageinmonths, data=pelotas)
pop_rss_height<-anova(lm_total)["Residuals", "Sum Sq"]

q_total<-lm(raw_height~poly(ageinmonths,2), data=pelotas)
pop_rss_height_q<-anova(q_total)["Residuals", "Sum Sq"]

conditional <- data.frame(subject, stringsAsFactors = FALSE)
conditional <- mutate(conditional, conditional_lm_height = numeric(length(subject)))
conditional <- mutate(conditional, conditional_q_height = numeric(length(subject)))

for (i in 1:nrow(wide)) {
  print(i)
  row.df <- cbind(as.numeric(ageinmonths), as.numeric(t(wide)[, i]))
  row.df <- as.data.frame(row.df)
  colnames(row.df) <- c("ageinmonths", m)
  row.df <- na.omit(row.df)
  
# Linear
  if (nrow(row.df) > 1) {
    row_lm = lm(row.df[[m]] ~ row.df$ageinmonths)
    conditional[which(conditional$subject == subject[i]),]$conditional_lm_height <-
        abs(pop_rss_height - anova(row_lm)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_height <-
      row_lm$coefficients[1]
    slopes[which(slopes$subject == subject[i]), ]$slopes_lm_height <-
      row_lm$coefficients[2]
    integrand <-
      function(x) {
        row_lm$coefficients[2] * x + row_lm$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_lm_height <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_lm_height <-
      auc[which(auc$subject == subject[i]), ]$auc_lm_height / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_lm_height <-
  #         NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_height <-
  #     NA
  #   slopes[which(slopes$subject == subject[i]), ]$slopes_lm_height <-
  #     NA
  #   auc[which(auc$subject == subject[i]), ]$auc_lm_height <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_lm_height <- NA
  # }

# Quadratic
  if (nrow(row.df) > 2) {
    row_q = lm(row.df[[m]] ~ poly(row.df$ageinmonths, 2))
    conditional[which(conditional$subject == subject[i]),]$conditional_q_height <-
        abs(pop_rss_height_q - anova(row_lm)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_height <-
      row_q$coefficients[1]
    slopes[which(slopes$subject == subject[i]), ]$slopes_q_height <-
      row_q$coefficients[3]
    integrand <-
      function(x) {
        row_q$coefficients[3] * I(x) ^ 2 + row_q$coefficients[2] * I(x) + row_q$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_q_height <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_q_height <-
      auc[which(auc$subject == subject[i]), ]$auc_q_height / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_q_height <- NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_height <-
  #     NA
  #   slopes[which(slopes$subject == subject[i]), ]$slopes_q_height <-
  #     NA
  #   auc[which(auc$subject == subject[i]), ]$auc_q_height <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_q_height <- NA
  # }
}

pelotas_features <-
  full_join(pelotas_features, intercepts, by = "subject")
pelotas_features <- full_join(pelotas_features, slopes, by = "subject")
pelotas_features <- full_join(pelotas_features, auc, by = "subject")
pelotas_features <- full_join(pelotas_features, tempo, by = "subject")
pelotas_features <- full_join(pelotas_features, conditional, by = "subject")

# ---- Measure: Weight -------

m <- "weight"
print(m)
wide <- pelotas_weight_wide
subject = as.numeric(wide[, 1])
ageinmonths <- colnames(wide)[-1]
colnames(wide)[1] <- "subject"
wide = wide[, -c(1)]

# Intercepts
intercepts <- data.frame(subject, stringsAsFactors = F)
intercepts <-
  mutate(intercepts, intercepts_lm_weight = numeric(length(subject)))
intercepts <-
  mutate(intercepts, intercepts_q_weight = numeric(length(subject)))

# Slopes
slopes <- data.frame(subject, stringsAsFactors = FALSE)
slopes <- mutate(slopes, slopes_lm_weight = numeric(length(subject)))
slopes <- mutate(slopes, slopes_q_weight = numeric(length(subject)))

# AUC
auc <- data.frame(subject, stringsAsFactors = FALSE)
auc <- mutate(auc, auc_lm_weight = numeric(length(subject)))
auc <- mutate(auc, auc_q_weight = numeric(length(subject)))

# Tempo
tempo <- data.frame(subject, stringsAsFactors = FALSE)
tempo <- mutate(tempo, tempo_lm_weight = numeric(length(subject)))
tempo <- mutate(tempo, tempo_q_weight = numeric(length(subject)))

# Conditional growth is the difference in RSS between the child's linear model and the population's
# linear model

lm_total<-lm(raw_weight~ageinmonths, data=pelotas)
pop_rss_weight<-anova(lm_total)["Residuals", "Sum Sq"]

q_total<-lm(raw_weight~poly(ageinmonths,2), data=pelotas)
pop_rss_weight_q<-anova(q_total)["Residuals", "Sum Sq"]

conditional <- data.frame(subject, stringsAsFactors = FALSE)
conditional <- mutate(conditional, conditional_lm_weight = numeric(length(subject)))
conditional <- mutate(conditional, conditional_q_weight = numeric(length(subject)))

for (i in 1:nrow(wide)) {
  print(i)
  row.df <- cbind(as.numeric(ageinmonths), as.numeric(t(wide)[, i]))
  row.df <- as.data.frame(row.df)
  colnames(row.df) <- c("ageinmonths", m)
  row.df <- na.omit(row.df)
  
# Linear
  if (nrow(row.df) > 1) {
    row_lm = lm(row.df[[m]] ~ row.df$ageinmonths)
    conditional[which(conditional$subject == subject[i]),]$conditional_lm_weight <-
        abs(pop_rss_weight - anova(row_lm)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_weight <-
      row_lm$coefficients[1]
    slopes[which(slopes$subject == subject[i]), ]$slopes_lm_weight <-
      row_lm$coefficients[2]
    integrand <-
      function(x) {
        row_lm$coefficients[2] * x + row_lm$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_lm_weight <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_lm_weight <-
      auc[which(auc$subject == subject[i]), ]$auc_lm_weight / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_lm_weight <- NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_weight <-
  #     NA
  #   slopes[which(slopes$subject == subject[i]), ]$slopes_lm_weight <-
  #     NA
  #   auc[which(auc$subject == subject[i]), ]$auc_lm_weight <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_lm_weight <- NA
  # }

# Quadratic
  if (nrow(row.df) > 2) {
    row_q = lm(row.df[[m]] ~ poly(row.df$ageinmonths, 2))
    conditional[which(conditional$subject == subject[i]),]$conditional_q_weight <-
        abs(pop_rss_weight_q - anova(row_q)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_weight <-
      row_q$coefficients[1]
    slopes[which(slopes$subject==subject[i]),]$slopes_q_weight<-row_q$coefficients[3]
    integrand <-
      function(x) {
        row_q$coefficients[3] * I(x) ^ 2 + row_q$coefficients[2] * I(x) + row_q$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_q_weight <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_q_weight <-
      auc[which(auc$subject == subject[i]), ]$auc_q_weight / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_q_weight <- NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_weight <-
  #     NA
  #   slopes[which(slopes$subject==subject[i]),]$slopes_q_weight<-NA
  #   auc[which(auc$subject == subject[i]), ]$auc_q_weight <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_q_weight <- NA
  # }
}

pelotas_features <-
  full_join(pelotas_features, intercepts, by = "subject")
pelotas_features <- full_join(pelotas_features, slopes, by = "subject")
pelotas_features <- full_join(pelotas_features, auc, by = "subject")
pelotas_features <- full_join(pelotas_features, tempo, by = "subject")
pelotas_features <- full_join(pelotas_features, conditional, by = "subject")

# ---- Measure: zWeight -------

m <- "zwei"
print(m)
wide <- pelotas_zwei_wide
subject = as.numeric(wide[, 1])
ageinmonths <- colnames(wide)[-1]
colnames(wide)[1] <- "subject"
wide = wide[, -c(1)]

# Intercepts
intercepts <- data.frame(subject, stringsAsFactors = F)
intercepts <-
  mutate(intercepts, intercepts_lm_zwei = numeric(length(subject)))
intercepts <-
  mutate(intercepts, intercepts_q_zwei = numeric(length(subject)))

# Slopes
slopes <- data.frame(subject, stringsAsFactors = FALSE)
slopes <- mutate(slopes, slopes_lm_zwei = numeric(length(subject)))
slopes <- mutate(slopes, slopes_q_zwei = numeric(length(subject)))

# AUC
auc <- data.frame(subject, stringsAsFactors = FALSE)
auc <- mutate(auc, auc_lm_zwei = numeric(length(subject)))
auc <- mutate(auc, auc_q_zwei = numeric(length(subject)))

# Tempo
tempo <- data.frame(subject, stringsAsFactors = FALSE)
tempo <- mutate(tempo, tempo_lm_zwei = numeric(length(subject)))
tempo <- mutate(tempo, tempo_q_zwei = numeric(length(subject)))

# Conditional growth is the difference in RSS between the child's linear model and the population's
# linear model

lm_total<-lm(zwei~ageinmonths, data=pelotas)
pop_rss_zwei<-anova(lm_total)["Residuals", "Sum Sq"]

q_total<-lm(zwei~poly(ageinmonths,2), data=pelotas)
pop_rss_zwei_q<-anova(q_total)["Residuals", "Sum Sq"]

conditional <- data.frame(subject, stringsAsFactors = FALSE)
conditional <- mutate(conditional, conditional_lm_zwei = numeric(length(subject)))
conditional <- mutate(conditional, conditional_q_zwei = numeric(length(subject)))

for (i in 1:nrow(wide)) {
  print(i)
  row.df <- cbind(as.numeric(ageinmonths), as.numeric(t(wide)[, i]))
  row.df <- as.data.frame(row.df)
  colnames(row.df) <- c("ageinmonths", m)
  row.df <- na.omit(row.df)
  
# Linear
  if (nrow(row.df) > 1) {
    row_lm = lm(row.df[[m]] ~ row.df$ageinmonths)
    conditional[which(conditional$subject == subject[i]),]$conditional_lm_zwei <-
        abs(pop_rss_zwei - anova(row_lm)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_zwei <-
      row_lm$coefficients[1]
    slopes[which(slopes$subject == subject[i]), ]$slopes_lm_zwei <-
      row_lm$coefficients[2]
    integrand <-
      function(x) {
        row_lm$coefficients[2] * x + row_lm$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_lm_zwei <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_lm_zwei <-
      auc[which(auc$subject == subject[i]), ]$auc_lm_zwei / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_lm_zwei <- NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_zwei <-
  #     NA
  #   slopes[which(slopes$subject == subject[i]), ]$slopes_lm_zwei <- NA
  #   auc[which(auc$subject == subject[i]), ]$auc_lm_zwei <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_lm_zwei <- NA
  # }

# Quadratic
  if (nrow(row.df) > 2) {
    row_q = lm(row.df[[m]] ~ poly(row.df$ageinmonths, 2))
    conditional[which(conditional$subject == subject[i]),]$conditional_q_zwei <-
        abs(pop_rss_zwei_q - anova(row_q)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_zwei <-
      row_q$coefficients[1]
    slopes[which(slopes$subject == subject[i]), ]$slopes_q_zwei <-
      row_q$coefficients[3]
    integrand <-
      function(x) {
        row_q$coefficients[3] * I(x) ^ 2 + row_q$coefficients[2] * I(x) + row_q$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_q_zwei <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_q_zwei <-
      auc[which(auc$subject == subject[i]), ]$auc_q_zwei / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_q_zwei <- NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_zwei <-
  #     NA
  #   slopes[which(slopes$subject == subject[i]), ]$slopes_q_zwei <- NA
  #   auc[which(auc$subject == subject[i]), ]$auc_q_zwei <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_q_zwei <- NA
  # }
}

pelotas_features <-
  full_join(pelotas_features, intercepts, by = "subject")
pelotas_features <- full_join(pelotas_features, slopes, by = "subject")
pelotas_features <- full_join(pelotas_features, auc, by = "subject")
pelotas_features <- full_join(pelotas_features, tempo, by = "subject")
pelotas_features <- full_join(pelotas_features, conditional, by = "subject")

# ---- Measure: zBMI -------

m <- "zbmi"
print(m)
wide <- pelotas_zbmi_wide
subject = as.numeric(wide[, 1])
ageinmonths <- colnames(wide)[-1]
colnames(wide)[1] <- "subject"
wide = wide[, -c(1)]

# Intercepts
intercepts <- data.frame(subject, stringsAsFactors = F)
intercepts <-
  mutate(intercepts, intercepts_lm_zbmi = numeric(length(subject)))
intercepts <-
  mutate(intercepts, intercepts_q_zbmi = numeric(length(subject)))

#Slopes
slopes <- data.frame(subject, stringsAsFactors = FALSE)
slopes <- mutate(slopes, slopes_lm_zbmi = numeric(length(subject)))
slopes <- mutate(slopes, slopes_q_zbmi = numeric(length(subject)))

#AUC
auc <- data.frame(subject, stringsAsFactors = FALSE)
auc <- mutate(auc, auc_lm_zbmi = numeric(length(subject)))
auc <- mutate(auc, auc_q_zbmi = numeric(length(subject)))

#Tempo
tempo <- data.frame(subject, stringsAsFactors = FALSE)
tempo <- mutate(tempo, tempo_lm_zbmi = numeric(length(subject)))
tempo <- mutate(tempo, tempo_q_zbmi = numeric(length(subject)))

# Conditional growth is the difference in RSS between the child's linear model and the population's
# linear model

lm_total<-lm(zbmi~ageinmonths, data=pelotas)
pop_rss_zbmi<-anova(lm_total)["Residuals", "Sum Sq"]

q_total<-lm(zbmi~poly(ageinmonths,2), data=pelotas)
pop_rss_zbmi_q<-anova(q_total)["Residuals", "Sum Sq"]

conditional <- data.frame(subject, stringsAsFactors = FALSE)
conditional <- mutate(conditional, conditional_lm_zbmi = numeric(length(subject)))
conditional <- mutate(conditional, conditional_q_zbmi = numeric(length(subject)))

for (i in 1:nrow(wide)) {
  print(i)
  row.df <- cbind(as.numeric(ageinmonths), as.numeric(t(wide)[, i]))
  row.df <- as.data.frame(row.df)
  colnames(row.df) <- c("ageinmonths", m)
  row.df <- na.omit(row.df)
  
# Linear
  if (nrow(row.df) > 1) {
    row_lm = lm(row.df[[m]] ~ row.df$ageinmonths)
    conditional[which(conditional$subject == subject[i]),]$conditional_lm_zbmi <-
        abs(pop_rss_zbmi - anova(row_lm)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_zbmi <-
      row_lm$coefficients[1]
    slopes[which(slopes$subject == subject[i]), ]$slopes_lm_zbmi <-
      row_lm$coefficients[2]
    integrand <-
      function(x) {
        row_lm$coefficients[2] * x + row_lm$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_lm_zbmi <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_lm_zbmi <-
      auc[which(auc$subject == subject[i]), ]$auc_lm_zbmi / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_lm_zbmi <- NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_lm_zbmi <-
  #     NA
  #   slopes[which(slopes$subject == subject[i]), ]$slopes_lm_zbmi <- NA
  #   auc[which(auc$subject == subject[i]), ]$auc_lm_zbmi <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_lm_zbmi <- NA
  # }

# Quadratic
  if (nrow(row.df) > 2) {
    row_q = lm(row.df[[m]] ~ poly(row.df$ageinmonths, 2))
    conditional[which(conditional$subject == subject[i]),]$conditional_q_zbmi <-
        abs(pop_rss_zbmi_q - anova(row_q)["Residuals", "Sum Sq"])
    intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_zbmi <-
      row_q$coefficients[1]
    slopes[which(slopes$subject == subject[i]), ]$slopes_q_zbmi <-
      row_q$coefficients[3]
    integrand <-
      function(x) {
        row_q$coefficients[3] * I(x) ^ 2 + row_q$coefficients[2] * I(x) + row_q$coefficients[1]
      }
    auc[which(auc$subject == subject[i]), ]$auc_q_zbmi <-
      integrate(integrand, row.df$ageinmonths[1], row.df$ageinmonths[length(row.df$ageinmonths)])$value
    tempo[which(tempo$subject == subject[i]), ]$tempo_q_zbmi <-
      auc[which(auc$subject == subject[i]), ]$auc_q_zbmi / (row.df[nrow(row.df), ]$ageinmonths - row.df[1, ]$ageinmonths)
  }
  # else {
  #     conditional[which(conditional$subject == subject[i]),]$conditional_q_zbmi <- NA
  #   intercepts[which(intercepts$subject == subject[i]), ]$intercepts_q_zbmi <-
  #     NA
  #   slopes[which(slopes$subject == subject[i]), ]$slopes_q_zbmi <- NA
  #   auc[which(auc$subject == subject[i]), ]$auc_q_zbmi <- NA
  #   tempo[which(tempo$subject == subject[i]), ]$tempo_q_zbmi <- NA
  # }
}

pelotas_features <-
  full_join(pelotas_features, intercepts, by = "subject")
pelotas_features <- full_join(pelotas_features, slopes, by = "subject")
pelotas_features <- full_join(pelotas_features, auc, by = "subject")
pelotas_features <- full_join(pelotas_features, tempo, by = "subject")
pelotas_features <- full_join(pelotas_features, conditional, by = "subject")

#----  Start-End 

pelotas_features_border <-
  data.frame(
    subject = subject,
    start_zbmi = apply(pelotas_zbmi_wide[,-1], 1,  function(x) min(x, na.rm = T)),
    end_zbmi = apply(pelotas_zbmi_wide[,-1], 1,  function(x) max(x, na.rm = T)),
    change_zbmi = apply(pelotas_zbmi_wide[,-1], 1,  function(x) max(x, na.rm = T)) - apply(pelotas_zbmi_wide[,-1], 1,  function(x) min(x, na.rm = T)),
    start_zlen = apply(pelotas_zlen_wide[,-1], 1,  function(x) min(x, na.rm = T)),
    end_zlen = apply(pelotas_zlen_wide[,-1], 1,  function(x) max(x, na.rm = T)),
    change_zlen = apply(pelotas_zlen_wide[,-1], 1,  function(x) max(x, na.rm = T)) - apply(pelotas_zlen_wide[,-1], 1,  function(x) min(x, na.rm = T)),
    start_zwei = apply(pelotas_zwei_wide[,-1], 1,  function(x) min(x, na.rm = T)),
    end_zwei = apply(pelotas_zwei_wide[,-1], 1,  function(x) max(x, na.rm = T)),
    change_zwei = apply(pelotas_zwei_wide[,-1], 1,  function(x) max(x, na.rm = T)) - apply(pelotas_zwei_wide[,-1], 1,  function(x) min(x, na.rm = T)),
    stringsAsFactors = F
  )
pelotas_features <-
  left_join(pelotas_features, pelotas_features_border, by = "subject")

#  Age at peak 

print("Age at peak")
peak_ages <-
  data.frame(
    subject = character(),
    peak_age_zbmi = numeric(),
    peak_zbmi = numeric(),
    peak_age_zlen = numeric(),
    peak_zlen = numeric(),
    peak_age_zwei = numeric(),
    peak_zwei = numeric(),
    stringsAsFactors = F
  )
peak_ages <- {
}

for (s in unique(pelotas$subject)) {
  subject_data <- pelotas[which(pelotas$subject == s), ]
  #avg_BMI<-mean(subject_data$PARENTBMI, na.rm = T)
  #peak_ages[nrow(peak_ages)+1,]<-c(subject=s, min_age_zbmi=pelotas[which.min(subject_data$zbmi),]$ageinmonths, peak_age_zbmi=pelotas[which.max(subject_data$zbmi),]$ageinmonths, peak_zbmi=pelotas[which.max(subject_data$zbmi),]$zbmi, peak_age_zlen=pelotas[which.max(subject_data$zlen),]$ageinmonths, peak_zlen=pelotas[which.max(subject_data$zlen),]$zlen, peak_age_zwei=pelotas[which.max(subject_data$zwei),]$ageinmonths, peak_zwei=pelotas[which.max(subject_data$zwei),]$zwei)
  peak_ages <-
    rbind(
      peak_ages,
      c(
        subject = s,
        peak_age_zbmi = subject_data[which.max(subject_data$zbmi), ]$ageinmonths,
        peak_zbmi = subject_data[which.max(subject_data$zbmi), ]$zbmi,
        peak_age_zlen = subject_data[which.max(subject_data$zlen), ]$ageinmonths,
        peak_zlen = subject_data[which.max(subject_data$zlen), ]$zlen,
        peak_age_zwei = subject_data[which.max(subject_data$zwei), ]$ageinmonths,
        peak_zwei = subject_data[which.max(subject_data$zwei), ]$zwei
      )
    )
}

peak_ages <- as.data.frame(peak_ages, stringsAsFactors = F)
peak_ages$peak_age_zbmi <- as.numeric(peak_ages$peak_age_zbmi)
peak_ages$peak_age_zlen <- as.numeric(peak_ages$peak_age_zlen)
peak_ages$peak_age_zwei <- as.numeric(peak_ages$peak_age_zwei)
peak_ages$peak_zbmi <- as.numeric(peak_ages$peak_zbmi)
peak_ages$peak_zlen <- as.numeric(peak_ages$peak_zlen)
peak_ages$peak_zwei <- as.numeric(peak_ages$peak_zwei)
peak_ages$subject <- as.numeric(peak_ages$subject)
pelotas_features <-
  full_join(pelotas_features, peak_ages, by = "subject")

# ---- Change in anthropometry >2 z scores -----------

print("Change in anthropometry")
change_anthropometry <- {
}

for (s in unique(pelotas$subject)) {
  subject_data <-
    arrange(pelotas[which(pelotas$subject == s), ], ageinmonths)
  rapid_growth_zbmi <- 0
  rapid_decrease_zbmi <- 0
  change_zbmi_YN <- 0
  rapid_growth_zlen <- 0
  rapid_decrease_zlen <- 0
  change_zlen_YN <- 0
  rapid_growth_zwei <- 0
  rapid_decrease_zwei <- 0
  change_zwei_YN <- 0
  for (i in 1:(nrow(subject_data) - 1)) {
    if (!(is.na(subject_data[i, ]$zbmi) |
          is.na(subject_data[i + 1, ]$zbmi))) {
      diff <- subject_data[i, ]$zbmi - subject_data[i + 1, ]$zbmi
      if (diff > 2) {
        rapid_growth_zbmi <- 1
        change_zbmi_YN <- 1
      }
      else if (diff < -2) {
        rapid_decrease_zbmi <- 1
        change_zbmi_YN <- 1
      }
    }
    if (!(is.na(subject_data[i, ]$zlen) |
          is.na(subject_data[i + 1, ]$zlen))) {
      diff <- subject_data[i, ]$zlen - subject_data[i + 1, ]$zlen
      if (diff > 2) {
        rapid_growth_zlen <- 1
        change_zlen_YN <- 1
      }
      else if (diff < -2) {
        rapid_decrease_zlen <- 1
        change_zlen_YN <- 1
      }
    }
    if (!(is.na(subject_data[i, ]$zwei) |
          is.na(subject_data[i + 1, ]$zwei))) {
      diff <- subject_data[i, ]$zwei - subject_data[i + 1, ]$zwei
      if (diff > 2) {
        rapid_growth_zwei <- 1
        change_zwei_YN <- 1
      }
      else if (diff < -2) {
        rapid_decrease_zwei <- 1
        change_zwei_YN <- 1
      }
    }
  }
  change_anthropometry <-
    rbind(
      change_anthropometry,
      c(
        s,
        rapid_growth_zbmi,
        rapid_decrease_zbmi,
        change_zbmi_YN,
        rapid_growth_zlen,
        rapid_decrease_zlen,
        change_zlen_YN,
        rapid_growth_zwei,
        rapid_decrease_zwei,
        change_zwei_YN
      )
    )
  
}
change_anthropometry <-
  as.data.frame(change_anthropometry, stringsAsFactors = F)
colnames(change_anthropometry) <-
  c(
    "subject",
    "rapid_growth_zbmi",
    "rapid_decrease_zbmi",
    "change_zbmi_YN",
    "rapid_growth_zlen",
    "rapid_decrease_zlen",
    "change_zlen_YN",
    "rapid_growth_zwei",
    "rapid_decrease_zwei",
    "change_zwei_YN"
  )
change_anthropometry$subject <-
  as.numeric(change_anthropometry$subject)
pelotas_features <-
  full_join(pelotas_features, change_anthropometry, by = "subject")



for (j in 1:ncol(pelotas_features)) {
    print(paste0(colnames(pelotas_features)[j], ": ", sum(is.na(pelotas_features[, j]))))
}

summary_table(pelotas_features, summaries = qsummary(pelotas_features), by = NULL)



