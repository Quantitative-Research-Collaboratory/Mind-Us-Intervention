library(mice)
library(psych)
library(dplyr)
library(ggplot2)
library(jtools)
library(svglite)


imputeddat <- mar.9shared@imputations
names.dat <- names(imputeddat[[1]])

means <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))
means0 <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))
means1 <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))
sds <- matrix (nrow = 20, ncol = ncol(imputeddat[[1]]))
sds0 <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))
sds1 <- matrix(nrow = 20, ncol = ncol(imputeddat[[1]]))
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

cohensds <- (describe(means0)$mean - describe(means1)$mean)/sqrt((describe(sds0)$mean^2+describe(sds1)$mean^2)/2)
names(cohensds) <- row.names(describe(means0))

## Graphs

# Depression Baseline Control
dep_pre_cont <- lapply(imputeddat, function(x) lm(PHQ9_0~Group,x))
summary(dep_pre_cont_pool <- pool(dep_pre_cont))

# Depression Baseline Treatment
dep_pre_trt <- lapply(imputeddat, function(x) lm(PHQ9_0~factor(Group, levels = c("1", "0")),x))
summary(dep_pre_trt_pool <- pool(dep_pre_trt))

# Depression Midtreatment Control
dep_mid_cont <- lapply(imputeddat, function(x) lm(PHQ9_1~Group,x))
summary(dep_mid_cont_pool <- pool(dep_mid_cont))

# Depression Midtreatment Treatment
dep_mid_trt <- lapply(imputeddat, function(x) lm(PHQ9_1~factor(Group, levels = c("1", "0")),x))
summary(dep_mid_trt_pool <- pool(dep_mid_trt))

# Depression Posttreatment Control
dep_post_cont <- lapply(imputeddat, function(x) lm(PHQ9_2~Group,x))
summary(dep_post_cont_pool <- pool(dep_post_cont))

# Depression Posttreatment Treatment
dep_post_trt <- lapply(imputeddat, function(x) lm(PHQ9_2~factor(Group, levels = c("1", "0")),x))
summary(dep_post_trt_pool <- pool(dep_post_trt))

### ANXIETY

# Anxiety Baseline Control
anx_pre_cont <- lapply(imputeddat, function(x) lm(GAD7_0~Group,x))
summary(anx_pre_cont_pool <- pool(anx_pre_cont))

# Anxiety Baseline Treatment
anx_pre_trt <- lapply(imputeddat, function(x) lm(GAD7_0~factor(Group, levels = c("1", "0")),x))
summary(anx_pre_trt_pool <- pool(anx_pre_trt))

# Anxiety Midtreatment Control
anx_mid_cont <- lapply(imputeddat, function(x) lm(GAD7_1~Group,x))
summary(anx_mid_cont_pool <- pool(anx_mid_cont))

# Anxiety Midtreatment Treatment
anx_mid_trt <- lapply(imputeddat, function(x) lm(GAD7_1~factor(Group, levels = c("1", "0")),x))
summary(anx_mid_trt_pool <- pool(anx_mid_trt))

# Anxiety Posttreatment Control
anx_post_cont <- lapply(imputeddat, function(x) lm(GAD7_2~Group,x))
summary(anx_post_cont_pool <- pool(anx_post_cont))

# Anxiety Posttreatment Treatment
anx_post_trt <- lapply(imputeddat, function(x) lm(GAD7_2~factor(Group, levels = c("1", "0")),x))
summary(anx_post_trt_pool <- pool(anx_post_trt))


### STRESS

# Stress Baseline Control
str_pre_cont <- lapply(imputeddat, function(x) lm(PSS_0~Group,x))
summary(str_pre_cont_pool <- pool(str_pre_cont))

# Stress Baseline Treatment
str_pre_trt <- lapply(imputeddat, function(x) lm(PSS_0~factor(Group, levels = c("1", "0")),x))
summary(str_pre_trt_pool <- pool(str_pre_trt))

# Stress Midtreatment Control
str_mid_cont <- lapply(imputeddat, function(x) lm(PSS_1~Group,x))
summary(str_mid_cont_pool <- pool(str_mid_cont))

# Stress Midtreatment Treatment
str_mid_trt <- lapply(imputeddat, function(x) lm(PSS_1~factor(Group, levels = c("1", "0")),x))
summary(str_mid_trt_pool <- pool(str_mid_trt))

# Stress Posttreatment Control
str_post_cont <- lapply(imputeddat, function(x) lm(PSS_2~Group,x))
summary(str_post_cont_pool <- pool(str_post_cont))

# Stress Posttreatment Treatment
str_post_trt <- lapply(imputeddat, function(x) lm(PSS_2~factor(Group, levels = c("1", "0")),x))
summary(str_post_trt_pool <- pool(str_post_trt))


# BEAQ
beaq.graph <- data.frame(matrix(nrow = 6))
beaq.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
beaq.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))


# BEAQ Baseline Treatment
beaq_pre_trt <- lapply(imputeddat, function(x) lm(BEAQ_0~factor(Group, levels = c("1", "0")),x))
beaq.graph$Est[1] <- summary(beaq_pre_trt_pool <- pool(beaq_pre_trt))$estimate[1]
beaq.graph$SE[1] <- summary(beaq_pre_trt_pool <- pool(beaq_pre_trt))$std.error[1]
beaq.graph$df[1] <- summary(beaq_pre_trt_pool <- pool(beaq_pre_trt))$df[1]

# BEAQ Baseline Control
beaq_pre_cont <- lapply(imputeddat, function(x) lm(BEAQ_0~Group,x))
beaq.graph$Est[2] <- summary(beaq_pre_cont_pool <- pool(beaq_pre_cont))$estimate[1]
beaq.graph$SE[2] <- summary(beaq_pre_cont_pool <- pool(beaq_pre_cont))$std.error[1]
beaq.graph$df[2] <- summary(beaq_pre_cont_pool <- pool(beaq_pre_cont))$df[1]

# BEAQ Midtreatment Treatment
beaq_mid_trt <- lapply(imputeddat, function(x) lm(BEAQ_1~factor(Group, levels = c("1", "0")),x))
beaq.graph$Est[3] <- summary(beaq_mid_trt_pool <- pool(beaq_mid_trt))$estimate[1]
beaq.graph$SE[3] <- summary(beaq_mid_trt_pool <- pool(beaq_mid_trt))$std.error[1]
beaq.graph$df[3] <- summary(beaq_mid_trt_pool <- pool(beaq_mid_trt))$df[1]

# BEAQ Midtreatment Control
beaq_mid_cont <- lapply(imputeddat, function(x) lm(BEAQ_1~Group,x))
beaq.graph$Est[4] <- summary(beaq_mid_cont_pool <- pool(beaq_mid_cont))$estimate[1]
beaq.graph$SE[4] <- summary(beaq_mid_cont_pool <- pool(beaq_mid_cont))$std.error[1]
beaq.graph$df[4] <- summary(beaq_mid_cont_pool <- pool(beaq_mid_cont))$df[1]

# BEAQ Posttreatment Treatment
beaq_post_trt <- lapply(imputeddat, function(x) lm(BEAQ_2~factor(Group, levels = c("1", "0")),x))
beaq.graph$Est[5] <- summary(beaq_post_trt_pool <- pool(beaq_post_trt))$estimate[1]
beaq.graph$SE[5] <- summary(beaq_post_trt_pool <- pool(beaq_post_trt))$std.error[1]
beaq.graph$df[5] <- summary(beaq_post_trt_pool <- pool(beaq_post_trt))$df[1]


# BEAQ Posttreatment Control
beaq_post_cont <- lapply(imputeddat, function(x) lm(BEAQ_2~Group,x))
beaq.graph$Est[6] <- summary(beaq_post_cont_pool <- pool(beaq_post_cont))$estimate[1]
beaq.graph$SE[6] <- summary(beaq_post_cont_pool <- pool(beaq_post_cont))$std.error[1]
beaq.graph$df[6] <- summary(beaq_post_cont_pool <- pool(beaq_post_cont))$df[1]

beaq.graph$CI.low <- beaq.graph$Est - beaq.graph$SE*qt(0.975,beaq.graph$df)
beaq.graph$CI.high <- beaq.graph$Est + beaq.graph$SE*qt(0.975,beaq.graph$df)


# MAAS (Mindfulness) Create data for plot
maas.graph <- data.frame(matrix(nrow = 6))
maas.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
maas.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))


# maas Baseline Treatment
maas_pre_trt <- lapply(imputeddat, function(x) lm(MAAS_0~factor(Group, levels = c("1", "0")),x))
maas.graph$Est[1] <- summary(maas_pre_trt_pool <- pool(maas_pre_trt))$estimate[1]
maas.graph$SE[1] <- summary(maas_pre_trt_pool <- pool(maas_pre_trt))$std.error[1]
maas.graph$df[1] <- summary(maas_pre_trt_pool <- pool(maas_pre_trt))$df[1]

# maas Baseline Control
maas_pre_cont <- lapply(imputeddat, function(x) lm(MAAS_0~Group,x))
maas.graph$Est[2] <- summary(maas_pre_cont_pool <- pool(maas_pre_cont))$estimate[1]
maas.graph$SE[2] <- summary(maas_pre_cont_pool <- pool(maas_pre_cont))$std.error[1]
maas.graph$df[2] <- summary(maas_pre_cont_pool <- pool(maas_pre_cont))$df[1]

# maas Midtreatment Treatment
maas_mid_trt <- lapply(imputeddat, function(x) lm(MAAS_1~factor(Group, levels = c("1", "0")),x))
maas.graph$Est[3] <- summary(maas_mid_trt_pool <- pool(maas_mid_trt))$estimate[1]
maas.graph$SE[3] <- summary(maas_mid_trt_pool <- pool(maas_mid_trt))$std.error[1]
maas.graph$df[3] <- summary(maas_mid_trt_pool <- pool(maas_mid_trt))$df[1]

# maas Midtreatment Control
maas_mid_cont <- lapply(imputeddat, function(x) lm(MAAS_1~Group,x))
maas.graph$Est[4] <- summary(maas_mid_cont_pool <- pool(maas_mid_cont))$estimate[1]
maas.graph$SE[4] <- summary(maas_mid_cont_pool <- pool(maas_mid_cont))$std.error[1]
maas.graph$df[4] <- summary(maas_mid_cont_pool <- pool(maas_mid_cont))$df[1]

# maas Posttreatment Treatment
maas_post_trt <- lapply(imputeddat, function(x) lm(MAAS_2~factor(Group, levels = c("1", "0")),x))
maas.graph$Est[5] <- summary(maas_post_trt_pool <- pool(maas_post_trt))$estimate[1]
maas.graph$SE[5] <- summary(maas_post_trt_pool <- pool(maas_post_trt))$std.error[1]
maas.graph$df[5] <- summary(maas_post_trt_pool <- pool(maas_post_trt))$df[1]


# maas Posttreatment Control
maas_post_cont <- lapply(imputeddat, function(x) lm(MAAS_2~Group,x))
maas.graph$Est[6] <- summary(maas_post_cont_pool <- pool(maas_post_cont))$estimate[1]
maas.graph$SE[6] <- summary(maas_post_cont_pool <- pool(maas_post_cont))$std.error[1]
maas.graph$df[6] <- summary(maas_post_cont_pool <- pool(maas_post_cont))$df[1]

maas.graph$CI.low <- maas.graph$Est - maas.graph$SE*qt(0.975,maas.graph$df)
maas.graph$CI.high <- maas.graph$Est + maas.graph$SE*qt(0.975,maas.graph$df)

# SCSSF (Self-Compassion) Create data for plot
scssf.graph <- data.frame(matrix(nrow = 6))
scssf.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
scssf.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# scssf Baseline Treatment
scssf_pre_trt <- lapply(imputeddat, function(x) lm(SCSSF_0~factor(Group, levels = c("1", "0")),x))
scssf.graph$Est[1] <- summary(scssf_pre_trt_pool <- pool(scssf_pre_trt))$estimate[1]
scssf.graph$SE[1] <- summary(scssf_pre_trt_pool <- pool(scssf_pre_trt))$std.error[1]
scssf.graph$df[1] <- summary(scssf_pre_trt_pool <- pool(scssf_pre_trt))$df[1]

# scssf Baseline Control
scssf_pre_cont <- lapply(imputeddat, function(x) lm(SCSSF_0~Group,x))
scssf.graph$Est[2] <- summary(scssf_pre_cont_pool <- pool(scssf_pre_cont))$estimate[1]
scssf.graph$SE[2] <- summary(scssf_pre_cont_pool <- pool(scssf_pre_cont))$std.error[1]
scssf.graph$df[2] <- summary(scssf_pre_cont_pool <- pool(scssf_pre_cont))$df[1]

# scssf Midtreatment Treatment
scssf_mid_trt <- lapply(imputeddat, function(x) lm(SCSSF_1~factor(Group, levels = c("1", "0")),x))
scssf.graph$Est[3] <- summary(scssf_mid_trt_pool <- pool(scssf_mid_trt))$estimate[1]
scssf.graph$SE[3] <- summary(scssf_mid_trt_pool <- pool(scssf_mid_trt))$std.error[1]
scssf.graph$df[3] <- summary(scssf_mid_trt_pool <- pool(scssf_mid_trt))$df[1]

# scssf Midtreatment Control
scssf_mid_cont <- lapply(imputeddat, function(x) lm(SCSSF_1~Group,x))
scssf.graph$Est[4] <- summary(scssf_mid_cont_pool <- pool(scssf_mid_cont))$estimate[1]
scssf.graph$SE[4] <- summary(scssf_mid_cont_pool <- pool(scssf_mid_cont))$std.error[1]
scssf.graph$df[4] <- summary(scssf_mid_cont_pool <- pool(scssf_mid_cont))$df[1]

# scssf Posttreatment Treatment
scssf_post_trt <- lapply(imputeddat, function(x) lm(SCSSF_2~factor(Group, levels = c("1", "0")),x))
scssf.graph$Est[5] <- summary(scssf_post_trt_pool <- pool(scssf_post_trt))$estimate[1]
scssf.graph$SE[5] <- summary(scssf_post_trt_pool <- pool(scssf_post_trt))$std.error[1]
scssf.graph$df[5] <- summary(scssf_post_trt_pool <- pool(scssf_post_trt))$df[1]


# scssf Posttreatment Control
scssf_post_cont <- lapply(imputeddat, function(x) lm(SCSSF_2~Group,x))
scssf.graph$Est[6] <- summary(scssf_post_cont_pool <- pool(scssf_post_cont))$estimate[1]
scssf.graph$SE[6] <- summary(scssf_post_cont_pool <- pool(scssf_post_cont))$std.error[1]
scssf.graph$df[6] <- summary(scssf_post_cont_pool <- pool(scssf_post_cont))$df[1]

scssf.graph$CI.low <- scssf.graph$Est - scssf.graph$SE*qt(0.975,scssf.graph$df)
scssf.graph$CI.high <- scssf.graph$Est + scssf.graph$SE*qt(0.975,scssf.graph$df)

# RRS (Rumination) Create data for plot
rrs.graph <- data.frame(matrix(nrow = 6))
rrs.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
rrs.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# rrs Baseline Treatment
rrs_pre_trt <- lapply(imputeddat, function(x) lm(RRS_0~factor(Group, levels = c("1", "0")),x))
rrs.graph$Est[1] <- summary(rrs_pre_trt_pool <- pool(rrs_pre_trt))$estimate[1]
rrs.graph$SE[1] <- summary(rrs_pre_trt_pool <- pool(rrs_pre_trt))$std.error[1]
rrs.graph$df[1] <- summary(rrs_pre_trt_pool <- pool(rrs_pre_trt))$df[1]

# rrs Baseline Control
rrs_pre_cont <- lapply(imputeddat, function(x) lm(RRS_0~Group,x))
rrs.graph$Est[2] <- summary(rrs_pre_cont_pool <- pool(rrs_pre_cont))$estimate[1]
rrs.graph$SE[2] <- summary(rrs_pre_cont_pool <- pool(rrs_pre_cont))$std.error[1]
rrs.graph$df[2] <- summary(rrs_pre_cont_pool <- pool(rrs_pre_cont))$df[1]

# rrs Midtreatment Treatment
rrs_mid_trt <- lapply(imputeddat, function(x) lm(RRS_1~factor(Group, levels = c("1", "0")),x))
rrs.graph$Est[3] <- summary(rrs_mid_trt_pool <- pool(rrs_mid_trt))$estimate[1]
rrs.graph$SE[3] <- summary(rrs_mid_trt_pool <- pool(rrs_mid_trt))$std.error[1]
rrs.graph$df[3] <- summary(rrs_mid_trt_pool <- pool(rrs_mid_trt))$df[1]

# rrs Midtreatment Control
rrs_mid_cont <- lapply(imputeddat, function(x) lm(RRS_1~Group,x))
rrs.graph$Est[4] <- summary(rrs_mid_cont_pool <- pool(rrs_mid_cont))$estimate[1]
rrs.graph$SE[4] <- summary(rrs_mid_cont_pool <- pool(rrs_mid_cont))$std.error[1]
rrs.graph$df[4] <- summary(rrs_mid_cont_pool <- pool(rrs_mid_cont))$df[1]

# rrs Posttreatment Treatment
rrs_post_trt <- lapply(imputeddat, function(x) lm(RRS_2~factor(Group, levels = c("1", "0")),x))
rrs.graph$Est[5] <- summary(rrs_post_trt_pool <- pool(rrs_post_trt))$estimate[1]
rrs.graph$SE[5] <- summary(rrs_post_trt_pool <- pool(rrs_post_trt))$std.error[1]
rrs.graph$df[5] <- summary(rrs_post_trt_pool <- pool(rrs_post_trt))$df[1]


# rrs Posttreatment Control
rrs_post_cont <- lapply(imputeddat, function(x) lm(RRS_2~Group,x))
rrs.graph$Est[6] <- summary(rrs_post_cont_pool <- pool(rrs_post_cont))$estimate[1]
rrs.graph$SE[6] <- summary(rrs_post_cont_pool <- pool(rrs_post_cont))$std.error[1]
rrs.graph$df[6] <- summary(rrs_post_cont_pool <- pool(rrs_post_cont))$df[1]

rrs.graph$CI.low <- rrs.graph$Est - rrs.graph$SE*qt(0.975,rrs.graph$df)
rrs.graph$CI.high <- rrs.graph$Est + rrs.graph$SE*qt(0.975,rrs.graph$df)

# ERQ_ES (Emotional Suppression) Create data for plot
erq.es.graph <- data.frame(matrix(nrow = 6))
erq.es.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
erq.es.graph$Time <- factor(c("Pre", "Pre", "Mid", "Mid", "Post", "Post"), levels = c("Pre", "Mid", "Post"))

# erq.es Baseline Treatment
erq.es_pre_trt <- lapply(imputeddat, function(x) lm(ERQ_ES_0~factor(Group, levels = c("1", "0")),x))
erq.es.graph$Est[1] <- summary(erq.es_pre_trt_pool <- pool(erq.es_pre_trt))$estimate[1]
erq.es.graph$SE[1] <- summary(erq.es_pre_trt_pool <- pool(erq.es_pre_trt))$std.error[1]
erq.es.graph$df[1] <- summary(erq.es_pre_trt_pool <- pool(erq.es_pre_trt))$df[1]

# erq.es Baseline Control
erq.es_pre_cont <- lapply(imputeddat, function(x) lm(ERQ_ES_0~Group,x))
erq.es.graph$Est[2] <- summary(erq.es_pre_cont_pool <- pool(erq.es_pre_cont))$estimate[1]
erq.es.graph$SE[2] <- summary(erq.es_pre_cont_pool <- pool(erq.es_pre_cont))$std.error[1]
erq.es.graph$df[2] <- summary(erq.es_pre_cont_pool <- pool(erq.es_pre_cont))$df[1]

# erq.es Midtreatment Treatment
erq.es_mid_trt <- lapply(imputeddat, function(x) lm(ERQ_ES_1~factor(Group, levels = c("1", "0")),x))
erq.es.graph$Est[3] <- summary(erq.es_mid_trt_pool <- pool(erq.es_mid_trt))$estimate[1]
erq.es.graph$SE[3] <- summary(erq.es_mid_trt_pool <- pool(erq.es_mid_trt))$std.error[1]
erq.es.graph$df[3] <- summary(erq.es_mid_trt_pool <- pool(erq.es_mid_trt))$df[1]

# erq.es Midtreatment Control
erq.es_mid_cont <- lapply(imputeddat, function(x) lm(ERQ_ES_1~Group,x))
erq.es.graph$Est[4] <- summary(erq.es_mid_cont_pool <- pool(erq.es_mid_cont))$estimate[1]
erq.es.graph$SE[4] <- summary(erq.es_mid_cont_pool <- pool(erq.es_mid_cont))$std.error[1]
erq.es.graph$df[4] <- summary(erq.es_mid_cont_pool <- pool(erq.es_mid_cont))$df[1]

# erq.es Posttreatment Treatment
erq.es_post_trt <- lapply(imputeddat, function(x) lm(ERQ_ES_2~factor(Group, levels = c("1", "0")),x))
erq.es.graph$Est[5] <- summary(erq.es_post_trt_pool <- pool(erq.es_post_trt))$estimate[1]
erq.es.graph$SE[5] <- summary(erq.es_post_trt_pool <- pool(erq.es_post_trt))$std.error[1]
erq.es.graph$df[5] <- summary(erq.es_post_trt_pool <- pool(erq.es_post_trt))$df[1]


# erq.es Posttreatment Control
erq.es_post_cont <- lapply(imputeddat, function(x) lm(ERQ_ES_2~Group,x))
erq.es.graph$Est[6] <- summary(erq.es_post_cont_pool <- pool(erq.es_post_cont))$estimate[1]
erq.es.graph$SE[6] <- summary(erq.es_post_cont_pool <- pool(erq.es_post_cont))$std.error[1]
erq.es.graph$df[6] <- summary(erq.es_post_cont_pool <- pool(erq.es_post_cont))$df[1]

erq.es.graph$CI.low <- erq.es.graph$Est - erq.es.graph$SE*qt(0.975,erq.es.graph$df)
erq.es.graph$CI.high <- erq.es.graph$Est + erq.es.graph$SE*qt(0.975,erq.es.graph$df)


# Depression Graph Create Dataframe
dep.graph <- data.frame(matrix(nrow = 6))
dep.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
dep.graph$Time <- c("Pre", "Pre", "Mid", "Mid", "Post", "Post")
dep.graph$Est <- c(9.96, 8.96, 7.86, 8.71, 6.31, 8.15)
dep.graph$SE <- c(0.62, 0.64, 0.66, 0.66, 0.65, 0.65)
dep.graph$df <- c(151, 151, 140, 150, 134, 150)
dep.graph$CI.low <- dep.graph$Est - dep.graph$SE*qt(0.975,dep.graph$df)
dep.graph$CI.high <- dep.graph$Est + dep.graph$SE*qt(0.975,dep.graph$df)

# Create the plot
dep.plot <- ggplot(dep.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(size = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(4, 12)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "PHQ-9") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

ggsave("dep_changeplot_mi.svg", plot = dep.plot)


anx.graph <- data.frame(matrix(nrow = 6))
anx.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
anx.graph$Time <- c("Pre", "Pre", "Mid", "Mid", "Post", "Post")
anx.graph$Est <- c(9.12, 7.91, 6.88, 7.97, 5.62, 7.72)
anx.graph$SE <- c(0.52, 0.54, 0.59, 0.59, 0.53, 0.53)
anx.graph$df <- c(151, 151, 134, 150, 141, 150)
anx.graph$CI.low <- anx.graph$Est - anx.graph$SE*qt(0.975,anx.graph$df)
anx.graph$CI.high <- anx.graph$Est + anx.graph$SE*qt(0.975,anx.graph$df)


# Create the plot
anx.plot <- ggplot(anx.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(size = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(4, 12)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "GAD-7") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

ggsave("anx_changeplot_mi.svg", plot = anx.plot)

str.graph <- data.frame(matrix(nrow = 6))
str.graph$Cond <- factor(c("Treatment","Control","Treatment","Control","Treatment","Control"))
str.graph$Time <- c("Pre", "Pre", "Mid", "Mid", "Post", "Post")
str.graph$Est <- c(21.55, 20.15, 18.38, 19.38, 16.22, 19.34)
str.graph$SE <- c(0.62, 0.64, 0.76, 0.77, 0.69, 0.68)
str.graph$df <- c(151, 151, 140, 150, 132, 150)
str.graph$CI.low <- str.graph$Est - str.graph$SE*qt(0.975,str.graph$df)
str.graph$CI.high <- str.graph$Est + str.graph$SE*qt(0.975,str.graph$df)

# Create the plot
str.plot <- ggplot(str.graph, aes(x = Time, y = Est, color = Cond, group = Cond)) +
  geom_line(size = 1) +  # Line plot
  geom_point(size = 3) +  # Points
  geom_errorbar(aes(ymin = CI.low, ymax = CI.high), width = 0.5) + # Error bars
  #scale_x_continuous(breaks = c(0, 4, 8, 12, 16), labels = c("Baseline", "4", "8", "12", "16")) +
  scale_y_continuous(limits = c(13, 25)) +  # Adjust limits for y-axis
  labs(x = "Time", y = "PSS") +  # Axis labels
  #scale_color_manual(values = c("0 (Control)" = "#1f77b4", "1 (Treatment)" = "#ff7f0e")) +  # Colors
  theme_apa() +  # Clean theme
  theme(legend.position = "right")

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


