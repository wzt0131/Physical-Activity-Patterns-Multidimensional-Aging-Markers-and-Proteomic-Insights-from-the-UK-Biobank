# Linear regression
# (only MRI brain-age outcome is shown as sample; the rest of the code remains the same)

# Weekend warrior pattern - Brain age
library(forestplot)
sports_MRI$smoke<-as.factor(sports_MRI$smoke)
sports_MRI$alcohol<-as.factor(sports_MRI$alcohol)
sports_MRI$Sex<-as.factor(sports_MRI$Sex)
sports_MRI$diabetes<-as.factor(sports_MRI$diabetes)
sports_MRI$high_bloodpressure<-as.factor(sports_MRI$high_bloodpressure)
sports_MRI$ethnic<-as.factor(sports_MRI$ethnic)
sports_MRI$education<-as.factor(sports_MRI$education)
sports_MRI$employment<-as.factor(sports_MRI$employment)
sports_MRI$Townsend.deprivation.index.at.recruitment<-as.numeric(sports_MRI$Townsend.deprivation.index.at.recruitment)
sports_MRI$Body.mass.index..BMI....Instance.0<-as.numeric(sports_MRI$Body.mass.index..BMI....Instance.0)
sports_MRI$weekly<-as.numeric(sports_MRI$weekly)
sports_MRI$BAG<-as.numeric(sports_MRI$BAG)
sports_MRI$weekend <- relevel(sports_MRI$weekend, ref = "inactive")
model <- lm(BAG ~ weekend+smoke+alcohol+Sex+high_bloodpressure
            +ethnic+education+employment+Townsend.deprivation.index.at.recruitment
            +Body.mass.index..BMI....Instance.0, data = sports_MRI)
model_summary <- tidy(model, conf.int = TRUE)
model_summary$term <- recode(model_summary$term,
                             "weekendactivenormal" = "activenormal",
                             "weekendactiveww" = "activeww")
ref_row <- data.frame(
  term = "inactive",  
  estimate = 0,       
  std.error = 0,
  statistic = 0,
  p.value = NA,
  conf.low = 0,      
  conf.high = 0     
)
model_summary <- rbind(ref_row, model_summary)
model_summary <- model_summary[model_summary$term %in% c("inactive", "activenormal", "activeww"), ]

format_p_value <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
}
model_summary$p.value<- format_p_value(model_summary$p.value)
labeltext <- sapply(1:nrow(model_summary), function(i) {
  if (model_summary$term[i] == "inactive") {
    sprintf("%s\n%.2f", model_summary$term[i], model_summary$estimate[i])
  } else {
    sprintf("%s\n%.2f (%.2f, %.2f)\n%s", 
            model_summary$term[i], model_summary$estimate[i], model_summary$conf.low[i], model_summary$conf.high[i], model_summary$p.value[i])
  }
})
labeltext <- matrix(labeltext, ncol = 1)
forestplot(
  labeltext = labeltext,
  mean = model_summary$estimate,
  lower = model_summary$conf.low,
  upper = model_summary$conf.high,
  zero = 0,  
  clip = c(-0.5, 1, 0.5) , 
  lineheight = "auto",
  boxsize = 0.15,
  col = fpColors(box = '#1f78b4', line = "darkblue"),
  txt_gp = fpTxtGp(label = gpar(cex = 1.0),
                   ticks = gpar(cex = 0.5)),
  xlab = "Coefficient",
  title = "Weekend Warrior")

# Daily exercise timing - Brain age
sportsday_MRI$smoke<-as.factor(sportsday_MRI$smoke)
sportsday_MRI$alcohol<-as.factor(sportsday_MRI$alcohol)
sportsday_MRI$Sex<-as.factor(sportsday_MRI$Sex)
sportsday_MRI$diabetes<-as.factor(sportsday_MRI$diabetes)
sportsday_MRI$high_bloodpressure<-as.factor(sportsday_MRI$high_bloodpressure)
sportsday_MRI$ethnic<-as.factor(sportsday_MRI$ethnic)
sportsday_MRI$education<-as.factor(sportsday_MRI$education)
sportsday_MRI$employment<-as.factor(sportsday_MRI$employment)
sportsday_MRI$Townsend.deprivation.index.at.recruitment<-as.numeric(sportsday_MRI$Townsend.deprivation.index.at.recruitment)
sportsday_MRI$Body.mass.index..BMI....Instance.0<-as.numeric(sportsday_MRI$Body.mass.index..BMI....Instance.0)
sportsday_MRI$BAG<-as.numeric(sportsday_MRI$BAG)
sportsday_MRI$sports_day <- relevel(sportsday_MRI$sports_day, ref = "inactive")
model <- lm(BAG ~ sports_day+smoke+alcohol+Sex+high_bloodpressure
            +ethnic+education+employment+Townsend.deprivation.index.at.recruitment
            +Body.mass.index..BMI....Instance.0, data = sportsday_MRI)
model_summary <- tidy(model, conf.int = TRUE)
model_summary$term <- recode(model_summary$term,
                             "sports_daymorning" = "morning",
                             "sports_dayafternoon" = "afternoon",
                             "sports_dayevening" = "evening",
                             "sports_daymixed" = "mixed")
ref_row <- data.frame(
  term = "inactive",  
  estimate = 0,      
  std.error = 0,
  statistic = 0,
  p.value = NA,
  conf.low = 0,      
  conf.high = 0     
)
model_summary <- rbind(ref_row, model_summary)
model_summary <- model_summary[model_summary$term %in% c("inactive", "morning","afternoon", "evening","mixed"), ]
new_order <- c(1,3,4,5,2)
model_summary <- model_summary[order(new_order), ]

format_p_value <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
}
model_summary$p.value<- format_p_value(model_summary$p.value)
labeltext <- sapply(1:nrow(model_summary), function(i) {
  if (model_summary$term[i] == "inactive") {
    sprintf("%s\n%.2f", model_summary$term[i], model_summary$estimate[i])
  } else {
    sprintf("%s\n%.2f (%.2f, %.2f)\n%s", 
            model_summary$term[i], model_summary$estimate[i], model_summary$conf.low[i], model_summary$conf.high[i], model_summary$p.value[i])
  }
})
labeltext <- matrix(labeltext, ncol = 1)
forestplot(
  labeltext = labeltext,
  mean = model_summary$estimate,
  lower = model_summary$conf.low,
  upper = model_summary$conf.high,
  zero = 0,  
  clip = c(-1, 0, 0.5) ,  
  lineheight = "auto",
  boxsize = 0.2,
  col = fpColors(box = '#1f78b4', line = "darkblue"),
  txt_gp = fpTxtGp(label = gpar(cex = 0.75),
                   ticks = gpar(cex = 0.5)),
  xlab = "Coefficient",
  title = "Physical Activity Timing")

# Weekly MVPA - Brain age
library(rms)
library(ggplot2)

ddist <- datadist(sports_MRI)  
options(datadist = "ddist")
model <- ols(BAG ~ rcs(weekly, 4)+smoke+alcohol+Sex+Townsend.deprivation.index.at.recruitment+Body.mass.index..BMI....Instance.0+high_bloodpressure+ethnic+education+employment, data = sports_MRI)
pred <- Predict(model, weekly = seq(0, 750, length = 10000)
                ,Townsend.deprivation.index.at.recruitment = mean(sports_MRI$Townsend.deprivation.index.at.recruitment), 
                Body.mass.index..BMI....Instance.0 = mean(sports_MRI$Body.mass.index..BMI....Instance.0), 
                smoke = "0",alcohol = "1",Sex = "1",high_bloodpressure = "0",ethnic = "1",education = "1",employment = "1")
pred_df <- data.frame(weekly= pred$weekly, yhat = pred$yhat,lower=pred$lower,upper=pred$upper)
y_start <- pred_df$yhat[1]
pred_df$yhat <- pred_df$yhat - y_start
pred_df$lower <- pred_df$lower - y_start
pred_df$upper <- pred_df$upper - y_start
min_point <- pred_df$weekly[which.min(pred_df$yhat)]
min_y <- pred_df$yhat[which.min(pred_df$yhat)]
ggplot(pred_df, aes(x = weekly, y = yhat)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "blue") +
  geom_vline(xintercept = min_point, color = "red", linetype = "dashed") +
  geom_text(aes(x = min_point,  y = min_y, label = paste("x =", round(min_point, 2))),color = "red", hjust = -0.1, vjust = 1.5) +
  labs(title = "MVPA (minutes/week)",x = "", y = "BAG") +
  theme_minimal()
