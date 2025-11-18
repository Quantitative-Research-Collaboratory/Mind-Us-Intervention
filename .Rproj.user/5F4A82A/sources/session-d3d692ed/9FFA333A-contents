library(mice)
library(psych)
library(dplyr)
library(ggplot2)
library(jtools)
library(svglite)

# Load in Data
imputeddat <- mar.9shared@imputations
names.dat <- names(imputeddat[[1]])


# Setup Matrices to Save Means and SDs
means <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))
means0 <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))
means1 <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))
sds <- matrix (nrow = 20, ncol = ncol(imputeddat[[1]]))
sds0 <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))
sds1 <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))

# Calculate Means and SDs
for (i in 1:20){
  temp <- dplyr::filter(imputeddat[[i]], Group == 0)
  means0[i,] <- t(describe(temp)$mean)
  sds0[i,] <- t(describe(temp)$sd)
  temp <- dplyr::filter(imputeddat[[i]], Group == 1)
  means1[i,] <- t(describe(temp)$mean)
  sds1[i,] <- t(describe(temp)$sd)
  means[i,] <- t(describe(imputeddat[[i]])$mean)
  sds[i,] <- t(describe(imputeddat[[i]])$sd)
}


# Change means and SDs to dataframes
means <- data.frame(means)
names(means) <- names.dat
means0 <- data.frame(means0)
names(means0) <- names.dat
means1 <- data.frame(means1)
names(means1) <- names.dat
sds <- data.frame(sds)
names(sds) <- names.dat
sds0 <- data.frame(sds0)
names(sds0) <- names.dat
sds1 <- data.frame(sds1)
names(sds1) <- names.dat

# mean of means for whole sample
describe(means)
#mean of means for control
describe(means0)
# mean of means for treatment
describe(means1)
#mean of sds for whole sample
describe(sds)
#mean of sds for control
describe(sds0)
# mean of sds for treatment
describe(sds1)

# Calculate Cohen's ds
cohensds <- (describe(means0)$mean - describe(means1)$mean)/sqrt((describe(sds0)$mean^2+describe(sds1)$mean^2)/2)
names(cohensds) <- row.names(describe(means0))

## Graphs

# Depression
dep.graph <- data.frame(matrix(nrow = 6))
dep.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
dep.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# PHQ9 Baseline Treatment
dep_pre_trt <- lapply(imputeddat, function(x) lm(PHQ9_0~factor(Group, levels = c("1", "0")),x))
dep.pre.trt.sum <- summary(dep_pre_trt_pool <- pool(dep_pre_trt))
dep.graph$Est[1] <- dep.pre.trt.sum$estimate[1]
dep.graph$SE[1] <- dep.pre.trt.sum$std.error[1]
dep.graph$df[1] <- dep.pre.trt.sum$df[1]

# PHQ9 Baseline Control
dep_pre_cont <- lapply(imputeddat, function(x) lm(PHQ9_0~Group,x))
dep.pre.cont.sum <- summary(dep_pre_cont_pool <- pool(dep_pre_cont))
dep.graph$Est[2] <- dep.pre.cont.sum$estimate[1]
dep.graph$SE[2] <- dep.pre.cont.sum$std.error[1]
dep.graph$df[2] <- dep.pre.cont.sum$df[1]

# PHQ9 Midtreatment Treatment
dep_mid_trt <- lapply(imputeddat, function(x) lm(PHQ9_1~factor(Group, levels = c("1", "0")),x))
dep.mid.trt.sum <- summary(dep_mid_trt_pool <- pool(dep_mid_trt))
dep.graph$Est[3] <- dep.mid.trt.sum$estimate[1]
dep.graph$SE[3] <- dep.mid.trt.sum$std.error[1]
dep.graph$df[3] <- dep.mid.trt.sum$df[1]

# PHQ9 Midtreatment Control
dep_mid_cont <- lapply(imputeddat, function(x) lm(PHQ9_1~Group,x))
dep.mid.trt.cont <- summary(dep_mid_cont_pool <- pool(dep_mid_cont))
dep.graph$Est[4] <- dep.mid.trt.cont$estimate[1]
dep.graph$SE[4] <- dep.mid.trt.cont$std.error[1]
dep.graph$df[4] <- dep.mid.trt.cont$df[1]

# PHQ9 Posttreatment Treatment
dep_post_trt <- lapply(imputeddat, function(x) lm(PHQ9_2~factor(Group, levels = c("1", "0")),x))
dep.post.trt.sum <- summary(dep_post_trt_pool <- pool(dep_post_trt))
dep.graph$Est[5] <- dep.post.trt.sum$estimate[1]
dep.graph$SE[5] <- dep.post.trt.sum$std.error[1]
dep.graph$df[5] <- dep.post.trt.sum$df[1]

# PHQ9 Posttreatment Control
dep_post_cont <- lapply(imputeddat, function(x) lm(PHQ9_2~Group,x))
dep.post.cont.sum <- summary(dep_post_cont_pool <- pool(dep_post_cont))
dep.graph$Est[6] <- dep.post.cont.sum$estimate[1]
dep.graph$SE[6] <- dep.post.cont.sum$std.error[1]
dep.graph$df[6] <- dep.post.cont.sum$df[1]

dep.graph$CI.low <- dep.graph$Est - dep.graph$SE*qt(0.975,dep.graph$df)
dep.graph$CI.high <- dep.graph$Est + dep.graph$SE*qt(0.975,dep.graph$df)

## ANXIETY
anx.graph <- data.frame(matrix(nrow = 6))
anx.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
anx.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# GAD7 Baseline Treatment
anx_pre_trt <- lapply(imputeddat, function(x) lm(GAD7_0~factor(Group, levels = c("1", "0")),x))
anx.pre.trt.sum <- summary(anx_pre_trt_pool <- pool(anx_pre_trt))
anx.graph$Est[1] <- anx.pre.trt.sum$estimate[1]
anx.graph$SE[1] <- anx.pre.trt.sum$std.error[1]
anx.graph$df[1] <- anx.pre.trt.sum$df[1]

# GAD7 Baseline Control
anx_pre_cont <- lapply(imputeddat, function(x) lm(GAD7_0~Group,x))
anx.pre.cont.sum <- summary(anx_pre_cont_pool <- pool(anx_pre_cont))
anx.graph$Est[2] <- anx.pre.cont.sum$estimate[1]
anx.graph$SE[2] <- anx.pre.cont.sum$std.error[1]
anx.graph$df[2] <- anx.pre.cont.sum$df[1]

# GAD7 Midtreatment Treatment
anx_mid_trt <- lapply(imputeddat, function(x) lm(GAD7_1~factor(Group, levels = c("1", "0")),x))
anx.mid.trt.sum <- summary(anx_mid_trt_pool <- pool(anx_mid_trt))
anx.graph$Est[3] <- anx.mid.trt.sum$estimate[1]
anx.graph$SE[3] <- anx.mid.trt.sum$std.error[1]
anx.graph$df[3] <- anx.mid.trt.sum$df[1]

# GAD7 Midtreatment Control
anx_mid_cont <- lapply(imputeddat, function(x) lm(GAD7_1~Group,x))
anx.mid.cont.sum <- summary(anx_mid_cont_pool <- pool(anx_mid_cont))
anx.graph$Est[4] <- anx.mid.cont.sum$estimate[1]
anx.graph$SE[4] <- anx.mid.cont.sum$std.error[1]
anx.graph$df[4] <- anx.mid.cont.sum$df[1]

# GAD7 Posttreatment Treatment
anx_post_trt <- lapply(imputeddat, function(x) lm(GAD7_2~factor(Group, levels = c("1", "0")),x))
anx.post.trt.sum <- summary(anx_post_trt_pool <- pool(anx_post_trt))
anx.graph$Est[5] <- anx.post.trt.sum$estimate[1]
anx.graph$SE[5] <- anx.post.trt.sum$std.error[1]
anx.graph$df[5] <- anx.post.trt.sum$df[1]

# GAD7 Posttreatment Control
anx_post_cont <- lapply(imputeddat, function(x) lm(GAD7_2~Group,x))
anx.post.cont.sum <- summary(anx_post_cont_pool <- pool(anx_post_cont))
anx.graph$Est[6] <- anx.post.cont.sum$estimate[1]
anx.graph$SE[6] <- anx.post.cont.sum$std.error[1]
anx.graph$df[6] <- anx.post.cont.sum$df[1]

anx.graph$CI.low <- anx.graph$Est - anx.graph$SE*qt(0.975,anx.graph$df)
anx.graph$CI.high <- anx.graph$Est + anx.graph$SE*qt(0.975,anx.graph$df)

## STRESS
str.graph <- data.frame(matrix(nrow = 6))
str.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
str.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# PSS Baseline Treatment
str_pre_trt <- lapply(imputeddat, function(x) lm(PSS_0~factor(Group, levels = c("1", "0")),x))
str.pre.trt.sum <- summary(str_pre_trt_pool <- pool(str_pre_trt))
str.graph$Est[1] <- str.pre.trt.sum$estimate[1]
str.graph$SE[1] <- str.pre.trt.sum$std.error[1]
str.graph$df[1] <- str.pre.trt.sum$df[1]

# PSS Baseline Control
str_pre_cont <- lapply(imputeddat, function(x) lm(PSS_0~Group,x))
str.pre.cont.sum <- summary(str_pre_cont_pool <- pool(str_pre_cont))
str.graph$Est[2] <- str.pre.cont.sum$estimate[1]
str.graph$SE[2] <- str.pre.cont.sum$std.error[1]
str.graph$df[2] <- str.pre.cont.sum$df[1]

# PSS Midtreatment Treatment
str_mid_trt <- lapply(imputeddat, function(x) lm(PSS_1~factor(Group, levels = c("1", "0")),x))
str.mid.trt.sum <- summary(str_mid_trt_pool <- pool(str_mid_trt))
str.graph$Est[3] <- str.mid.trt.sum$estimate[1]
str.graph$SE[3] <- str.mid.trt.sum$std.error[1]
str.graph$df[3] <- str.mid.trt.sum$df[1]

# PSS Midtreatment Control
str_mid_cont <- lapply(imputeddat, function(x) lm(PSS_1~Group,x))
str.mid.cont.sum <- summary(str_mid_cont_pool <- pool(str_mid_cont))
str.graph$Est[4] <- str.mid.cont.sum$estimate[1]
str.graph$SE[4] <- str.mid.cont.sum$std.error[1]
str.graph$df[4] <- str.mid.cont.sum$df[1]

# PSS Posttreatment Treatment
str_post_trt <- lapply(imputeddat, function(x) lm(PSS_2~factor(Group, levels = c("1", "0")),x))
str.post.trt.sum <- summary(str_post_trt_pool <- pool(str_post_trt))
str.graph$Est[5] <- str.post.trt.sum$estimate[1]
str.graph$SE[5] <- str.post.trt.sum$std.error[1]
str.graph$df[5] <- str.post.trt.sum$df[1]

# PSS Posttreatment Control
str_post_cont <- lapply(imputeddat, function(x) lm(PSS_2~Group,x))
str.post.cont.sum <- summary(str_post_cont_pool <- pool(str_post_cont))
str.graph$Est[6] <- str.post.cont.sum $estimate[1]
str.graph$SE[6] <- str.post.cont.sum $std.error[1]
str.graph$df[6] <- str.post.cont.sum $df[1]

str.graph$CI.low <- str.graph$Est - str.graph$SE*qt(0.975,str.graph$df)
str.graph$CI.high <- str.graph$Est + str.graph$SE*qt(0.975,str.graph$df)

# BEAQ
beaq.graph <- data.frame(matrix(nrow = 6))
beaq.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
beaq.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))


# BEAQ Baseline Treatment
beaq_pre_trt <- lapply(imputeddat, function(x) lm(BEAQ_0~factor(Group, levels = c("1", "0")),x))
beaq.pre.trt.sum <- summary(beaq_pre_trt_pool <- pool(beaq_pre_trt))
beaq.graph$Est[1] <- beaq.pre.trt.sum$estimate[1]
beaq.graph$SE[1] <- beaq.pre.trt.sum$std.error[1]
beaq.graph$df[1] <- beaq.pre.trt.sum$df[1]

# BEAQ Baseline Control
beaq_pre_cont <- lapply(imputeddat, function(x) lm(BEAQ_0~Group,x))
beaq.pre.cont.sum <- summary(beaq_pre_cont_pool <- pool(beaq_pre_cont))
beaq.graph$Est[2] <- beaq.pre.cont.sum$estimate[1]
beaq.graph$SE[2] <- beaq.pre.cont.sum$std.error[1]
beaq.graph$df[2] <- beaq.pre.cont.sum$df[1]

# BEAQ Midtreatment Treatment
beaq_mid_trt <- lapply(imputeddat, function(x) lm(BEAQ_1~factor(Group, levels = c("1", "0")),x))
beaq.mid.trt.sum <- summary(beaq_mid_trt_pool <- pool(beaq_mid_trt))
beaq.graph$Est[3] <- beaq.mid.trt.sum $estimate[1]
beaq.graph$SE[3] <- beaq.mid.trt.sum $std.error[1]
beaq.graph$df[3] <- beaq.mid.trt.sum $df[1]

# BEAQ Midtreatment Control
beaq_mid_cont <- lapply(imputeddat, function(x) lm(BEAQ_1~Group,x))
beaq.mid.cont.sum <- summary(beaq_mid_cont_pool <- pool(beaq_mid_cont))
beaq.graph$Est[4] <- beaq.mid.cont.sum$estimate[1]
beaq.graph$SE[4] <- beaq.mid.cont.sum$std.error[1]
beaq.graph$df[4] <- beaq.mid.cont.sum$df[1]

# BEAQ Posttreatment Treatment
beaq_post_trt <- lapply(imputeddat, function(x) lm(BEAQ_2~factor(Group, levels = c("1", "0")),x))
beaq.post.trt.sum <- summary(beaq_post_trt_pool <- pool(beaq_post_trt))
beaq.graph$Est[5] <- beaq.post.trt.sum$estimate[1]
beaq.graph$SE[5] <- beaq.post.trt.sum$std.error[1]
beaq.graph$df[5] <- beaq.post.trt.sum$df[1]

# BEAQ Posttreatment Control
beaq_post_cont <- lapply(imputeddat, function(x) lm(BEAQ_2~Group,x))
beaq.post.cont.sum <- summary(beaq_post_cont_pool <- pool(beaq_post_cont))
beaq.graph$Est[6] <- beaq.post.cont.sum$estimate[1]
beaq.graph$SE[6] <- beaq.post.cont.sum$std.error[1]
beaq.graph$df[6] <- beaq.post.cont.sum$df[1]

beaq.graph$CI.low <- beaq.graph$Est - beaq.graph$SE*qt(0.975,beaq.graph$df)
beaq.graph$CI.high <- beaq.graph$Est + beaq.graph$SE*qt(0.975,beaq.graph$df)


# MAAS (Mindfulness) Create data for plot
maas.graph <- data.frame(matrix(nrow = 6))
maas.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
maas.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# maas Baseline Treatment
maas_pre_trt <- lapply(imputeddat, function(x) lm(MAAS_0~factor(Group, levels = c("1", "0")),x))
maas.pre.trt.sum <- summary(maas_pre_trt_pool <- pool(maas_pre_trt))
maas.graph$Est[1] <- maas.pre.trt.sum $estimate[1]
maas.graph$SE[1] <- maas.pre.trt.sum $std.error[1]
maas.graph$df[1] <- maas.pre.trt.sum $df[1]

# maas Baseline Control
maas_pre_cont <- lapply(imputeddat, function(x) lm(MAAS_0~Group,x))
maas.pre.cont.sum <- summary(maas_pre_cont_pool <- pool(maas_pre_cont))
maas.graph$Est[2] <- maas.pre.cont.sum$estimate[1]
maas.graph$SE[2] <- maas.pre.cont.sum$std.error[1]
maas.graph$df[2] <- maas.pre.cont.sum$df[1]

# maas Midtreatment Treatment
maas_mid_trt <- lapply(imputeddat, function(x) lm(MAAS_1~factor(Group, levels = c("1", "0")),x))
maas.mid.trt.sum <- summary(maas_mid_trt_pool <- pool(maas_mid_trt))
maas.graph$Est[3] <- maas.mid.trt.sum$estimate[1]
maas.graph$SE[3] <- maas.mid.trt.sum$std.error[1]
maas.graph$df[3] <- maas.mid.trt.sum$df[1]

# maas Midtreatment Control
maas_mid_cont <- lapply(imputeddat, function(x) lm(MAAS_1~Group,x))
maas.mid.cont.sum <- summary(maas_mid_cont_pool <- pool(maas_mid_cont))
maas.graph$Est[4] <- maas.mid.cont.sum$estimate[1]
maas.graph$SE[4] <- maas.mid.cont.sum$std.error[1]
maas.graph$df[4] <- maas.mid.cont.sum$df[1]

# maas Posttreatment Treatment
maas_post_trt <- lapply(imputeddat, function(x) lm(MAAS_2~factor(Group, levels = c("1", "0")),x))
maas.post.trt.sum <- summary(maas_post_trt_pool <- pool(maas_post_trt))
maas.graph$Est[5] <- maas.post.trt.sum$estimate[1]
maas.graph$SE[5] <- maas.post.trt.sum$std.error[1]
maas.graph$df[5] <- maas.post.trt.sum$df[1]


# maas Posttreatment Control
maas_post_cont <- lapply(imputeddat, function(x) lm(MAAS_2~Group,x))
maas.post.cont.sum <- summary(maas_post_cont_pool <- pool(maas_post_cont))
maas.graph$Est[6] <- maas.post.cont.sum$estimate[1]
maas.graph$SE[6] <- maas.post.cont.sum$std.error[1]
maas.graph$df[6] <- maas.post.cont.sum$df[1]

maas.graph$CI.low <- maas.graph$Est - maas.graph$SE*qt(0.975,maas.graph$df)
maas.graph$CI.high <- maas.graph$Est + maas.graph$SE*qt(0.975,maas.graph$df)

# SCSSF (Self-Compassion) Create data for plot
scssf.graph <- data.frame(matrix(nrow = 6))
scssf.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
scssf.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# scssf Baseline Treatment
scssf_pre_trt <- lapply(imputeddat, function(x) lm(SCSSF_0~factor(Group, levels = c("1", "0")),x))
scssf.pre.trt.sum <- summary(scssf_pre_trt_pool <- pool(scssf_pre_trt))
scssf.graph$Est[1] <- scssf.pre.trt.sum$estimate[1]
scssf.graph$SE[1] <- scssf.pre.trt.sum$std.error[1]
scssf.graph$df[1] <- scssf.pre.trt.sum$df[1]

# scssf Baseline Control
scssf_pre_cont <- lapply(imputeddat, function(x) lm(SCSSF_0~Group,x))
scssf.pre.cont.sum <- summary(scssf_pre_cont_pool <- pool(scssf_pre_cont))
scssf.graph$Est[2] <- scssf.pre.cont.sum$estimate[1]
scssf.graph$SE[2] <- scssf.pre.cont.sum$std.error[1]
scssf.graph$df[2] <- scssf.pre.cont.sum$df[1]

# scssf Midtreatment Treatment
scssf_mid_trt <- lapply(imputeddat, function(x) lm(SCSSF_1~factor(Group, levels = c("1", "0")),x))
scssf.mid.trt.sum <- summary(scssf_mid_trt_pool <- pool(scssf_mid_trt))
scssf.graph$Est[3] <- scssf.mid.trt.sum$estimate[1]
scssf.graph$SE[3] <- scssf.mid.trt.sum$std.error[1]
scssf.graph$df[3] <- scssf.mid.trt.sum$df[1]

# scssf Midtreatment Control
scssf_mid_cont <- lapply(imputeddat, function(x) lm(SCSSF_1~Group,x))
scssf.mid.cont.sum <- summary(scssf_mid_cont_pool <- pool(scssf_mid_cont))
scssf.graph$Est[4] <- scssf.mid.cont.sum$estimate[1]
scssf.graph$SE[4] <- scssf.mid.cont.sum$std.error[1]
scssf.graph$df[4] <- scssf.mid.cont.sum$df[1]

# scssf Posttreatment Treatment
scssf_post_trt <- lapply(imputeddat, function(x) lm(SCSSF_2~factor(Group, levels = c("1", "0")),x))
scssf.post.trt.sum <- summary(scssf_post_trt_pool <- pool(scssf_post_trt))
scssf.graph$Est[5] <- scssf.post.trt.sum$estimate[1]
scssf.graph$SE[5] <- scssf.post.trt.sum$std.error[1]
scssf.graph$df[5] <- scssf.post.trt.sum$df[1]


# scssf Posttreatment Control
scssf_post_cont <- lapply(imputeddat, function(x) lm(SCSSF_2~Group,x))
scssf.post.cont.sum <- summary(scssf_post_cont_pool <- pool(scssf_post_cont))
scssf.graph$Est[6] <- scssf.post.cont.sum$estimate[1]
scssf.graph$SE[6] <- scssf.post.cont.sum$std.error[1]
scssf.graph$df[6] <- scssf.post.cont.sum$df[1]

scssf.graph$CI.low <- scssf.graph$Est - scssf.graph$SE*qt(0.975,scssf.graph$df)
scssf.graph$CI.high <- scssf.graph$Est + scssf.graph$SE*qt(0.975,scssf.graph$df)

# RRS (Rumination) Create data for plot
rrs.graph <- data.frame(matrix(nrow = 6))
rrs.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
rrs.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# rrs Baseline Treatment
rrs_pre_trt <- lapply(imputeddat, function(x) lm(RRS_0~factor(Group, levels = c("1", "0")),x))
rrs.pre.trt.sum <- summary(rrs_pre_trt_pool <- pool(rrs_pre_trt))
rrs.graph$Est[1] <- rrs.pre.trt.sum$estimate[1]
rrs.graph$SE[1] <- rrs.pre.trt.sum$std.error[1]
rrs.graph$df[1] <- rrs.pre.trt.sum$df[1]

# rrs Baseline Control
rrs_pre_cont <- lapply(imputeddat, function(x) lm(RRS_0~Group,x))
rrs.pre.cont.sum <- summary(rrs_pre_cont_pool <- pool(rrs_pre_cont))
rrs.graph$Est[2] <- rrs.pre.cont.sum$estimate[1]
rrs.graph$SE[2] <- rrs.pre.cont.sum$std.error[1]
rrs.graph$df[2] <- rrs.pre.cont.sum$df[1]

# rrs Midtreatment Treatment
rrs_mid_trt <- lapply(imputeddat, function(x) lm(RRS_1~factor(Group, levels = c("1", "0")),x))
rrs.mid.trt.sum <- summary(rrs_mid_trt_pool <- pool(rrs_mid_trt))
rrs.graph$Est[3] <- rrs.mid.trt.sum$estimate[1]
rrs.graph$SE[3] <- rrs.mid.trt.sum$std.error[1]
rrs.graph$df[3] <- rrs.mid.trt.sum$df[1]

# rrs Midtreatment Control
rrs_mid_cont <- lapply(imputeddat, function(x) lm(RRS_1~Group,x))
rrs.mid.cont.sum <- summary(rrs_mid_cont_pool <- pool(rrs_mid_cont))
rrs.graph$Est[4] <- rrs.mid.cont.sum$estimate[1]
rrs.graph$SE[4] <- rrs.mid.cont.sum$std.error[1]
rrs.graph$df[4] <- rrs.mid.cont.sum$df[1]

# rrs Posttreatment Treatment
rrs_post_trt <- lapply(imputeddat, function(x) lm(RRS_2~factor(Group, levels = c("1", "0")),x))
rrs.post.trt.sum <- summary(rrs_post_trt_pool <- pool(rrs_post_trt))
rrs.graph$Est[5] <- rrs.post.trt.sum$estimate[1]
rrs.graph$SE[5] <- rrs.post.trt.sum$std.error[1]
rrs.graph$df[5] <- rrs.post.trt.sum$df[1]


# rrs Posttreatment Control
rrs_post_cont <- lapply(imputeddat, function(x) lm(RRS_2~Group,x))
rrs.post.cont.sum <- summary(rrs_post_cont_pool <- pool(rrs_post_cont))
rrs.graph$Est[6] <- rrs.post.cont.sum$estimate[1]
rrs.graph$SE[6] <- rrs.post.cont.sum$std.error[1]
rrs.graph$df[6] <- rrs.post.cont.sum$df[1]

rrs.graph$CI.low <- rrs.graph$Est - rrs.graph$SE*qt(0.975,rrs.graph$df)
rrs.graph$CI.high <- rrs.graph$Est + rrs.graph$SE*qt(0.975,rrs.graph$df)

# ERQ_ES (Emotional Suppression) Create data for plot
erq.es.graph <- data.frame(matrix(nrow = 6))
erq.es.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
erq.es.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# erq.es Baseline Treatment
erq.es_pre_trt <- lapply(imputeddat, function(x) lm(ERQ_ES_0~factor(Group, levels = c("1", "0")),x))
erq.es.pre.trt.sum <- summary(erq.es_pre_trt_pool <- pool(erq.es_pre_trt))
erq.es.graph$Est[1] <- erq.es.pre.trt.sum$estimate[1]
erq.es.graph$SE[1] <- erq.es.pre.trt.sum$std.error[1]
erq.es.graph$df[1] <- erq.es.pre.trt.sum$df[1]

# erq.es Baseline Control
erq.es_pre_cont <- lapply(imputeddat, function(x) lm(ERQ_ES_0~Group,x))
erq.es.pre.cont <- summary(erq.es_pre_cont_pool <- pool(erq.es_pre_cont))
erq.es.graph$Est[2] <- erq.es.pre.cont$estimate[1]
erq.es.graph$SE[2] <- erq.es.pre.cont$std.error[1]
erq.es.graph$df[2] <- erq.es.pre.cont$df[1]

# erq.es Midtreatment Treatment
erq.es_mid_trt <- lapply(imputeddat, function(x) lm(ERQ_ES_1~factor(Group, levels = c("1", "0")),x))
erq.es.mid.trt.sum <- summary(erq.es_mid_trt_pool <- pool(erq.es_mid_trt))
erq.es.graph$Est[3] <- erq.es.mid.trt.sum$estimate[1]
erq.es.graph$SE[3] <- erq.es.mid.trt.sum$std.error[1]
erq.es.graph$df[3] <- erq.es.mid.trt.sum$df[1]

# erq.es Midtreatment Control
erq.es_mid_cont <- lapply(imputeddat, function(x) lm(ERQ_ES_1~Group,x))
erq.es.mid.cont.sum <- summary(erq.es_mid_cont_pool <- pool(erq.es_mid_cont))
erq.es.graph$Est[4] <- erq.es.mid.cont.sum$estimate[1]
erq.es.graph$SE[4] <- erq.es.mid.cont.sum$std.error[1]
erq.es.graph$df[4] <- erq.es.mid.cont.sum$df[1]

# erq.es Posttreatment Treatment
erq.es_post_trt <- lapply(imputeddat, function(x) lm(ERQ_ES_2~factor(Group, levels = c("1", "0")),x))
erq.es.post.trt <- summary(erq.es_post_trt_pool <- pool(erq.es_post_trt))
erq.es.graph$Est[5] <- erq.es.post.trt$estimate[1]
erq.es.graph$SE[5] <- erq.es.post.trt$std.error[1]
erq.es.graph$df[5] <- erq.es.post.trt$df[1]

# erq.es Posttreatment Control
erq.es_post_cont <- lapply(imputeddat, function(x) lm(ERQ_ES_2~Group,x))
erq.es.post.cont.sum <- summary(erq.es_post_cont_pool <- pool(erq.es_post_cont))
erq.es.graph$Est[6] <- erq.es.post.cont.sum$estimate[1]
erq.es.graph$SE[6] <- erq.es.post.cont.sum$std.error[1]
erq.es.graph$df[6] <- erq.es.post.cont.sum$df[1]

erq.es.graph$CI.low <- erq.es.graph$Est - erq.es.graph$SE*qt(0.975,erq.es.graph$df)
erq.es.graph$CI.high <- erq.es.graph$Est + erq.es.graph$SE*qt(0.975,erq.es.graph$df)

# ERQ_CR (Cognitive Reappraisal) Create data for plot
erq.cr.graph <- data.frame(matrix(nrow = 6))
erq.cr.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
erq.cr.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# erq.cr Baseline Treatment
erq.cr_pre_trt <- lapply(imputeddat, function(x) lm(ERQ_CR_0~factor(Group, levels = c("1", "0")),x))
erq.cr.pre.trt.sum <- summary(erq.cr_pre_trt_pool <- pool(erq.cr_pre_trt))
erq.cr.graph$Est[1] <- erq.cr.pre.trt.sum$estimate[1]
erq.cr.graph$SE[1] <- erq.cr.pre.trt.sum$std.error[1]
erq.cr.graph$df[1] <- erq.cr.pre.trt.sum$df[1]

# erq.cr Baseline Control
erq.cr_pre_cont <- lapply(imputeddat, function(x) lm(ERQ_CR_0~Group,x))
erq.cr.pre.cont <- summary(erq.cr_pre_cont_pool <- pool(erq.cr_pre_cont))
erq.cr.graph$Est[2] <- erq.cr.pre.cont$estimate[1]
erq.cr.graph$SE[2] <- erq.cr.pre.cont$std.error[1]
erq.cr.graph$df[2] <- erq.cr.pre.cont$df[1]

# erq.cr Midtreatment Treatment
erq.cr_mid_trt <- lapply(imputeddat, function(x) lm(ERQ_CR_1~factor(Group, levels = c("1", "0")),x))
erq.cr.mid.trt.sum <- summary(erq.cr_mid_trt_pool <- pool(erq.cr_mid_trt))
erq.cr.graph$Est[3] <- erq.cr.mid.trt.sum$estimate[1]
erq.cr.graph$SE[3] <- erq.cr.mid.trt.sum$std.error[1]
erq.cr.graph$df[3] <- erq.cr.mid.trt.sum$df[1]

# erq.cr Midtreatment Control
erq.cr_mid_cont <- lapply(imputeddat, function(x) lm(ERQ_CR_1~Group,x))
erq.cr.mid.cont.sum <- summary(erq.cr_mid_cont_pool <- pool(erq.cr_mid_cont))
erq.cr.graph$Est[4] <- erq.cr.mid.cont.sum$estimate[1]
erq.cr.graph$SE[4] <- erq.cr.mid.cont.sum$std.error[1]
erq.cr.graph$df[4] <- erq.cr.mid.cont.sum$df[1]

# erq.cr Posttreatment Treatment
erq.cr_post_trt <- lapply(imputeddat, function(x) lm(ERQ_CR_2~factor(Group, levels = c("1", "0")),x))
erq.cr.post.trt <- summary(erq.cr_post_trt_pool <- pool(erq.cr_post_trt))
erq.cr.graph$Est[5] <- erq.cr.post.trt$estimate[1]
erq.cr.graph$SE[5] <- erq.cr.post.trt$std.error[1]
erq.cr.graph$df[5] <- erq.cr.post.trt$df[1]

# erq.cr Posttreatment Control
erq.cr_post_cont <- lapply(imputeddat, function(x) lm(ERQ_CR_2~Group,x))
erq.cr.post.cont.sum <- summary(erq.cr_post_cont_pool <- pool(erq.cr_post_cont))
erq.cr.graph$Est[6] <- erq.cr.post.cont.sum$estimate[1]
erq.cr.graph$SE[6] <- erq.cr.post.cont.sum$std.error[1]
erq.cr.graph$df[6] <- erq.cr.post.cont.sum$df[1]

erq.cr.graph$CI.low <- erq.cr.graph$Est - erq.cr.graph$SE*qt(0.975,erq.cr.graph$df)
erq.cr.graph$CI.high <- erq.cr.graph$Est + erq.cr.graph$SE*qt(0.975,erq.cr.graph$df)


# Depression: Create the plot
dep.plot <- ggplot(dep.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(size = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(4, 12)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "Depression (PHQ-9)") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

ggsave("dep_changeplot_mi.png", plot = dep.plot)
ggsave("dep_changeplot_mi.svg", plot = dep.plot)


# Anxiety: Create the plot
anx.plot <- ggplot(anx.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(size = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(4, 12)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "Anxiety (GAD-7)") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

ggsave("anx_changeplot_mi.png", plot = anx.plot)
ggsave("anx_changeplot_mi.svg", plot = anx.plot)

# STRESS: Create the plot
str.plot <- ggplot(str.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(size = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(13, 25)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "Stress (PSS)") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

ggsave("str_changeplot_mi.png", plot = str.plot)
ggsave("str_changeplot_mi.svg", plot = str.plot)

# Create BEAQ (Experiential Avoidance) Plot
beaq.plot <- ggplot(beaq.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(linewidth = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(40, 55)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "Experiential Avoidance (BEAQ)") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

beaq.plot

ggsave("beaq_changeplot_mi.png", plot = beaq.plot)
ggsave("beaq_changeplot_mi.svg", plot = beaq.plot)

# Create MAAS (Mindfulness) Plot
maas.plot <- ggplot(maas.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(linewidth = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(45, 65)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "Mindfulness (MAAS)") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

maas.plot

ggsave("maas_changeplot_mi.png", plot = maas.plot)
ggsave("maas_changeplot_mi.svg", plot = maas.plot)

# Create SCSSF (Self-Compassion) Plot
scssf.plot <- ggplot(scssf.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(linewidth = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(28, 40)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "Self-Compassion (SCSSF)") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

scssf.plot

ggsave("scssf_changeplot_mi.png", plot = scssf.plot)
ggsave("scssf_changeplot_mi.svg", plot = scssf.plot)

# Create RRS (Rumination) Plot
rrs.plot <- ggplot(rrs.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(linewidth = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(8, 15)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "Rumination (RRS)") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

rrs.plot

ggsave("rrs_changeplot_mi.png", plot = rrs.plot)
ggsave("rrs_changeplot_mi.svg", plot = rrs.plot)


# Create ERQ_ES (Emotion Suppression) Plot
erq.es.plot <- ggplot(erq.es.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(linewidth = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(12, 18)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "Emotional Suppression (ERQ-ES)") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

erq.es.plot

ggsave("erq.es_changeplot_mi.png", plot = erq.es.plot)
ggsave("erq.es_changeplot_mi.svg", plot = erq.es.plot)

# Create ERQ_CR (Cognitive Reappraisal) Plot
erq.cr.plot <- ggplot(erq.cr.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(linewidth = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(24,32)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "Cognitive Reappraisal (ERQ-CR)") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

erq.cr.plot

ggsave("erq.cr_changeplot_mi.png", plot = erq.cr.plot)
ggsave("erq.cr_changeplot_mi.svg", plot = erq.cr.plot)



