# Mediation analysis
# (only phenotypic age outcome is shown as sample; the rest of the code remains the same)

# Metabolite mediation
library(mediation)
library(glmnet)
library(broom)
library(dplyr)
sports_phenoage_try<-merge(sports_phenoage,mid,by='Participant.ID')
sports_phenoage_try<-subset(sports_phenoage_try,select= -Age.at.recruitment.y)
metabolites <- names(sports_phenoage_try)[34:284]
results <- data.frame(
  Metabolite = character(),
  Estimate = numeric(),
  `Std. Error` = numeric(),
  `t value` = numeric(),
  `Pr(>|t|)` = numeric(),
  stringsAsFactors = FALSE
)
for (metabolite in metabolites) {
  formula <- as.formula(paste0("weekly ~ ", metabolite, " + ", paste(names(sports_phenoage_try)[21:29], collapse = "+")))
  model <- lm(formula, data = sports_phenoage_try)
  model_tidy <- tidy(model)
  print(names(model_tidy))
  relevant_row <- model_tidy %>%
    filter(term == metabolite) %>%
    select(term, estimate, std.error, statistic, p.value)
  colnames(relevant_row) <- c("Metabolite", "Estimate", "Std. Error", "t value", "Pr(>|t|)")
  results <- bind_rows(results, relevant_row)
}
P_B<-results %>%
  filter(`Pr(>|t|)` < 0.05/251)
cols_to_select <- P_B$Metabolite
filtered_data <- mid %>%
  select(1:2, all_of(cols_to_select))
if(!is.data.frame(filtered_data)) {
  filtered_data <- as.data.frame(filtered_data)
}
y_try<-sports_phenoage_try%>% select(Participant.ID,weekly)
filtered_data<-merge(filtered_data,y_try,by='Participant.ID')
names(filtered_data)
set.seed(123)
x <- filtered_data[,3:34] 
y <- filtered_data[, 35]        
set.seed(123)
train_size <- floor(0.7 * length(y))
trainIndex <- sample(seq_len(length(y)), size = train_size)
x_train <- x[trainIndex, ]
x_val <- x[-trainIndex, ]
y_train <- y[trainIndex]
y_val <- y[-trainIndex]
library(caret)
model <- train(
  x = x_train, y = y_train,
  method = 'glmnet', 
  trControl = trainControl('cv', number = 10), 
  tuneLength = 10
)
sum(is.na(x_train))
sum(is.na(y_train))
model$bestTune
coef(model$finalModel, model$bestTune$lambda) 
coefficients <- coef(model$finalModel, model$bestTune$lambda) 
coefficients_df <- as.data.frame(as.matrix(coefficients))
coefficients <- coefficients[-1] 
new_data<-sports_phenoage_try%>%dplyr::select(Apolipoprotein.A1...Instance.0,
                                              Average.Diameter.for.HDL.Particles...Instance.0, 
                                              Cholesterol.in.Large.HDL...Instance.0,                                          
                                              Cholesterol.in.Medium.HDL...Instance.0,                                    
                                              Cholesterol.in.Very.Large.HDL...Instance.0,                                
                                              Cholesteryl.Esters.in.HDL...Instance.0,                                     
                                              Cholesteryl.Esters.in.Large.HDL...Instance.0,                               
                                              Cholesteryl.Esters.in.Medium.HDL...Instance.0,                              
                                              Cholesteryl.Esters.in.Very.Large.HDL...Instance.0,                            
                                              Concentration.of.HDL.Particles...Instance.0,                                  
                                              Concentration.of.Large.HDL.Particles...Instance.0,                            
                                              Concentration.of.Medium.HDL.Particles...Instance.0,                           
                                              Concentration.of.Very.Large.HDL.Particles...Instance.0,                        
                                              Free.Cholesterol.in.HDL...Instance.0,                                         
                                              Free.Cholesterol.in.Large.HDL...Instance.0,                                    
                                              Free.Cholesterol.in.Medium.HDL...Instance.0,                                  
                                              Free.Cholesterol.to.Total.Lipids.in.Medium.HDL.percentage...Instance.0,     
                                              Free.Cholesterol.to.Total.Lipids.in.Very.Large.HDL.percentage...Instance.0,   
                                              HDL.Cholesterol...Instance.0,                                               
                                              Phosphatidylcholines...Instance.0,                                           
                                              Phosphoglycerides...Instance.0,                                             
                                              Phospholipids.in.HDL...Instance.0,                                         
                                              Phospholipids.in.Large.HDL...Instance.0,                                     
                                              Phospholipids.in.Medium.HDL...Instance.0,                                    
                                              Total.Cholines...Instance.0,                                                 
                                              Total.Concentration.of.Lipoprotein.Particles...Instance.0,                  
                                              Total.Lipids.in.HDL...Instance.0,                                             
                                              Total.Lipids.in.Large.HDL...Instance.0,                                     
                                              Total.Lipids.in.Medium.HDL...Instance.0,                                   
                                              Total.Lipids.in.Very.Large.HDL...Instance.0,                               
                                              Total.Phospholipids.in.Lipoprotein.Particles...Instance.0,                  
                                              Triglycerides.to.Phosphoglycerides.ratio...Instance.0)
metabolite_scores <- as.matrix(new_data) %*% coefficients
metabolite_scores <- as.data.frame(as.matrix(metabolite_scores))
metabolite_scores$metabolite_scores<-metabolite_scores$V1+13536.36870
sports_phenoage_try<-cbind(sports_phenoage_try,metabolite_scores)
sports_phenoage_try <- as.data.frame(sports_phenoage_try)
sports_phenoage_try$smoke<-as.factor(sports_phenoage_try$smoke)
sports_phenoage_try$alcohol<-as.factor(sports_phenoage_try$alcohol)
sports_phenoage_try$Sex<-as.factor(sports_phenoage_try$Sex)
sports_phenoage_try$diabetes<-as.factor(sports_phenoage_try$diabetes)
sports_phenoage_try$high_bloodpressure<-as.factor(sports_phenoage_try$high_bloodpressure)
sports_phenoage_try$ethnic<-as.factor(sports_phenoage_try$ethnic)
sports_phenoage_try$education<-as.factor(sports_phenoage_try$education)
sports_phenoage_try$employment<-as.factor(sports_phenoage_try$employment)
sports_phenoage_try$Townsend.deprivation.index.at.recruitment<-as.numeric(sports_phenoage_try$Townsend.deprivation.index.at.recruitment)
sports_phenoage_try$Body.mass.index..BMI....Instance.0<-as.numeric(sports_phenoage_try$Body.mass.index..BMI....Instance.0)
model_m <- lm(metabolite_scores~ weekly+smoke+alcohol+Sex+high_bloodpressure+education+employment+ethnic+Townsend.deprivation.index.at.recruitment+Body.mass.index..BMI....Instance.0, data = sports_phenoage_try)
model_y <- lm(PhenoAgeAccel ~ weekly+metabolite_scores+smoke+alcohol+Sex+high_bloodpressure+education+employment+ethnic+Townsend.deprivation.index.at.recruitment+Body.mass.index..BMI....Instance.0, data = sports_phenoage_try)
mediation_result <- mediate(
  model.m = model_m, 
  model.y = model_y,  
  treat = "weekly",       
  mediator = "metabolite_scores",  
  boot = TRUE,        
  sims = 1000        
)
summary(mediation_result)

# Proteome mediation
library(tidyverse)
library(broom)
library(sandwich)
library(mediation) 
library(ggplot2)
sports_phenoage_pro<-merge(sports_phenoage,protein_filtered,BY="Participant.ID")
protein_cols <- 36:2955
sports_phenoage_pro[, protein_cols] <- scale(sports_phenoage_pro[, protein_cols])
cat("\Standardization completed", length(protein_cols), "column\n")
library(dplyr)
library(broom)
protein_results <- data.frame()
protein_names <- names(sports_phenoage_pro)[36:2955]  
cat("find", length(protein_names))
sports_phenoage_pro <- sports_phenoage_pro %>%
  mutate(weekly_std = scale(weekly))
for(protein in protein_names) {
  formula <- as.formula(paste(
    protein, 
    "~weekly_std"
  ))
  fit <- lm(formula, data = sports_phenoage_pro)
  result <- tidy(fit) %>%
    filter(term == "weekly_std") %>%  
    mutate(protein = protein, .before = term)
  protein_results <- rbind(protein_results, result)
  if(which(protein == protein_names) %% 100 == 0) {
    cat("Completed", which(protein == protein_names), "/", length(protein_names), "protein analyses")
  }
}
protein_results$p_adjusted <- p.adjust(protein_results$p.value, method = "bonferroni")
significant_proteins <- protein_results %>%
  filter(p_adjusted < 0.05 & abs(estimate) > 0.25)
cat("find", nrow(significant_proteins), "significant proteins\n")
mediation_results <- data.frame()
sports_phenoage_pro$weekly_std <- as.numeric(sports_phenoage_pro$weekly_std)
sports_phenoage_pro$PhenoAgeAccel_std <- as.numeric(scale(sports_phenoage_pro$PhenoAgeAccel))
library(mediation) 
results_list <- list()
for(i in seq_along(significant_proteins$protein)) {
  protein <- significant_proteins$protein[i]
  cat("analyze (", i, "/", length(significant_proteins$protein), "):", protein, "\n")
  tryCatch({
    model_mediator <- lm(as.formula(paste(protein, "~ weekly_std+ smoke + alcohol + Sex + Townsend.deprivation.index.at.recruitment + Body.mass.index..BMI....Instance.0 + high_bloodpressure + ethnic + education + employment")), 
                         data = sports_phenoage_pro)
    model_outcome <- lm(as.formula(paste("PhenoAgeAccel_std ~ weekly_std +", protein,"+smoke+alcohol+Sex+Townsend.deprivation.index.at.recruitment+Body.mass.index..BMI....Instance.0+high_bloodpressure+ethnic+education+employment")), 
                        data = sports_phenoage_pro)
    med_fit <- mediate(model_mediator, model_outcome, 
                       treat = "weekly_std", mediator = protein,
                       boot = TRUE, sims = 1000)
    result_summary <- summary(med_fit)
    results_list[[protein]] <- data.frame(
      protein = protein,
      ACME = result_summary$d0,
      ACME_p = result_summary$d0.p,
      ADE = result_summary$z0,
      ADE_p = result_summary$z0.p,
      Total_Effect = result_summary$tau.coef,
      Prop_Mediated = result_summary$n0,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("  false:", e$message, "\n")
    results_list[[protein]] <- data.frame(
      protein = protein,
      ACME = NA,
      ACME_p = NA,
      ADE = NA,
      ADE_p = NA,
      Total_Effect = NA,
      Prop_Mediated = NA,
      stringsAsFactors = FALSE
    )
  })
}
mediation_results <- do.call(rbind, results_list)
mediation_results$ACME_p_adjusted <- p.adjust(mediation_results$ACME_p, method = "bonferroni")
significant_mediation <- mediation_results %>%
  filter(ACME_p_adjusted < 0.05)
write.csv(significant_mediation, "significant_proteins_mediation.csv", row.names = FALSE)
