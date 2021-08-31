library(tidyverse)

# plots in supplementary material that compare p-values of 
# FCoxPH models with false shape information with CoxPH model

## lung cancer
load("./data/lung/lung_loocv_result.Rdata")
resdf = data.frame(permres)
ggplot(resdf, aes(x=permres))+
  geom_histogram(binwidth=0.000005, fill="white",color="black") +
  geom_vline(xintercept=0.000012, color="red",linetype="dashed") + # p-value of the CoxPH model
  labs(x="P-values") +
  theme_classic() +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"))

## brain tumor
load("./data/brain/brain_loocv_result.Rdata")
resdf = data.frame(permres)
ggplot(resdf, aes(x=permres))+
  geom_histogram(binwidth=0.001, fill="white",color="black") +
  geom_vline(xintercept=0.0017, color="red",linetype="dashed") + # p-value of the CoxPH model
  labs(x="P-values") +
  theme_classic() +
  theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "cm"))
