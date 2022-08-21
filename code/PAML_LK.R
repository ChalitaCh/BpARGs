#Clean the working environment
rm(list = ls())
graphics.off()

#load the dataset(s)
genes_LK <- read.csv("../data/genes_likelihood.csv", header = TRUE)

#Calculate the log-likelihood ratio and their p values

#using equation log-likelihood = 2*(L2 - L1)
#where L2 = log-likelihood of model 2 (complicated one)
#and L1 = log-likelihood of model 1 (less complicated)

#parameters for each model for significant testing df = P2 - P1
#where P2 is no. of free parameters in model 2 and P1 is no. of free parameters in model 1
#model0 para = 1
#model1 para = 2
#model2 para = 4

genes_LK$M2vsM1 <- 2 * (genes_LK$M2 - genes_LK$M1)

genes_LK$p.valueM2vsM1 <- pchisq(genes_LK$M2vsM1, df = 4-2, lower.tail = FALSE)

genes_LK$M2vsM0 <- 2 * (genes_LK$M2 - genes_LK$M0)

genes_LK$p.valueM2vsM0 <- pchisq(genes_LK$M2vsM0, df = 4-1, lower.tail = FALSE)

genes_LK$M1vsM0 <- 2 * (genes_LK$M1 - genes_LK$M0)

genes_LK$p.valueM1vsM0 <- pchisq(genes_LK$M1vsM0, df = 2-1, lower.tail = FALSE)

genes_LK$M0vsM2 <- 2 * (genes_LK$M0 - genes_LK$M2)

genes_LK$p.valueM0vsM2 <- pchisq(genes_LK$M0vsM2, df = 4-2, lower.tail = FALSE)

#Calculate the AIC of each model
#using  (-2 * log(L)) + 2 * p
#L = the likelihood of the model
#p = no. of free parameters in the model
genes_LK$AIC0 <- (-2*genes_LK$M0) + 2*1
genes_LK$AIC1 <- (-2*genes_LK$M1) + 2*2
genes_LK$AIC2 <- (-2*genes_LK$M2) + 2*4

#find which model has the lowest AIC score
genes_LK$model <- apply(genes_LK[,11:13], 1,which.min)
genes_LK$model <- ifelse(genes_LK$model == "c(AIC0 = 1)", "Model0", genes_LK$model)
genes_LK$model <- ifelse(genes_LK$model == "c(AIC1 = 2)", "Model1", genes_LK$model)
genes_LK$model <- ifelse(genes_LK$model == "c(AIC2 = 3)", "Model2", genes_LK$model)

write.csv(genes_LK, file = "../results/ARGs_selection_pressure.csv", quote = FALSE, row.names = FALSE)