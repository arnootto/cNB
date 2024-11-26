# The cNB package
An R package to implement the methodology described in *The Contamianted Negegative Binomial Regression Model* (2024).

To install the package, use the following code in R
```{r}
#install.packages("devtools")
library(devtools)
install_github("arnootto/cNB-RM")
```
## Example
Code to reproduce Example 5.1: Number of visits to a doctor in a year
```{r}
library(cNB)
library(COUNT)
data("badhealth")
est=ml.cnb(formula = numvisit ~ badh+age,data=badhealth)
```
