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
  temp <- filter(imputeddat[[i]], Group == 0)
  means0[i,] <- t(describe(temp)$mean)
  sds0[i,] <- t(describe(temp)$sd)
  temp <- filter(imputeddat[[i]], Group == 1)
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
