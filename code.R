
library(nlme)
library(lme4)
library(optimx)

rm(list = ls())

# Include data
# source("progs/03-prep-registry-data.R")


library(glmmTMB)
library(metafor)

dat <- dat[dat$bearingsurface != "Ceramicised Metal/XLPE", ]

# Get variance-covariance matrix
dat$within <- paste0(dat$study, dat$bearingsurface)
# Use the observed as predicted at the first run
dat$os.pred <- dat$os
# Calculate the SE of log(OS)
dat$lnse <- exp(-exp(dat$se/(dat$os*log(dat$os))))

cor <- diag(nrow(dat))

for(i in 1:(nrow(dat) - 1)){
  for(j in (i+1):nrow(dat)){
    if(dat$within[i] == dat$within[j]){
      os.k <- dat$os.pred[i]
      os.kt <- dat$os.pred[j]
      cor.val <- sqrt(os.k * (1 - os.kt)/(1 - os.k) * os.kt)
      cor[i, j] <- cor[j, i] <- exp(-exp(cor.val))
    }
  }
}

vcov <- dat$lnse %*% t(dat$lnse) * cor

# One random effect
fit2.0 <- rma.mv(y~ bearingsurface * poly(ln.year, 2),
                 V = vcov,
                 random = ~ 1 | study,
                 data = dat)

fit2.1 <- rma.mv(y~ bearingsurface * ln.year,
                 V = vcov,
                 random = ~ 1 | study,
                 data = dat)
fit2.2 <- rma.mv(y ~ bearingsurface + ln.year,
                 V = vcov,
                 random = ~ 1 | study,
                 data = dat)

fit2.year <- rma.mv(y~ bearingsurface + ln.year,
                    V = vcov,
                    random = ~ ln.year | study,
                    data = dat)

fit2.bs <- rma.mv(y~ bearingsurface + ln.year,
                  V = vcov,
                  random = ~ bearingsurface | study,
                  data = dat)

fit2.bo <- rma.mv(y~ bearingsurface + ln.year,
                  V = vcov,
                  random = list( ~ ln.year | study,
                                 ~ bearingsurface | study),
                  data = dat)


fit2.year.int <- rma.mv(y~ bearingsurface * ln.year,
                        V = vcov,
                        random = ~ ln.year | study,
                        data = dat)

fit2.year.polyint <- rma.mv(y~ bearingsurface * poly(ln.year, 2),
                            V = vcov,
                            random = ~ ln.year | study,
                            data = dat)

perf_out <- performance::compare_performance(
  fit2.year, fit2.year.int, fit2.year.polyint,
  fit2.year, fit2.bs, fit2.bo,
  fit2.0, fit2.1,
  fit2.2
)

writexl::write_xlsx(perf_out, path = "Output/model-performance.xlsx")

fit <- fit2.2


# Plot ==============
library(dplyr)
library(ggplot2)
library(patchwork)

newDat <- expand.grid(year = 1:30,
                      bearingsurface = unique(dat$bearingsurface))
newDat$ln.year <- log(newDat$year)

trans.fn <- function(x){
  exp(-exp(x))
}

pred <- predict(fit,
                newmods = model.matrix(formula(fit), newDat)[,-1],
                transf = trans.fn)
pred <- as.data.frame(pred)

pd <- position_dodge(1)

p1 <- dat |>
  mutate(lci = 1 - as.numeric(lci)/100,
         ucl = 1 - as.numeric(uci)/100) |>
  ggplot(aes(year, os, col = bearingsurface, shape = study)) +
  geom_point(position = pd) +
  geom_pointrange(aes(ymin = lci, ymax = ucl), position = pd) +
  scale_shape_manual(values = c(3, 4, 8, 13, 16, 17, 18, 15))+
  scale_color_brewer(palette = "Set1") +
  labs(x = "Year", y = "Survival probability",
       shape = 'Data source',
       col = "Bearing materials") +
  theme_bw(base_size = 16)

predicted_data <- bind_cols(newDat, pred)

p2 <- rbind(predicted_data[predicted_data$year <= 20,] |>
              mutate(obs_pred = "obs"),
            predicted_data[predicted_data$year >= 20,] |>
              mutate(obs_pred = "pred")) |>
  ggplot(aes(year, pred, col = bearingsurface)) +
  geom_line(aes(linetype = obs_pred), show.legend = FALSE) +
  # geom_ribbon(aes(ymin = pi.lb, ymax = pi.ub, fill = bearingsurface),
  #             alpha = 0.1, colour = NA) +
  geom_line(aes(y = pi.lb, col = bearingsurface), linetype = "dotted") +
  geom_line(aes(y = pi.ub, col = bearingsurface), linetype = "dotted") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_linetype_manual(values = c("obs" = "solid", pred = "dashed")) +
  geom_vline(xintercept = 20, linetype = "dotted", color = "black") +
  labs(x = "Year", y = "Survival probability",
       fill = "Bearing materials",
       col = "Bearing materials") +
  theme_bw(base_size = 16)

p <- p1/p2 +
  plot_annotation(tag_levels = 'a')

ggsave("Output/Figure 4.pdf", p,
       width= 10, height = 10,
       device = cairo_pdf)


openxlsx2::write_xlsx(predicted_data,
                      file = "output/predicted.xlsx")




