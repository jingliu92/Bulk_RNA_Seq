# Detection of Deferentially Expressed Gene

Transcript expression was exported as read count matrices from GTF fles.

<img width="614" alt="Screen Shot 2022-10-12 at 12 39 38 PM" src="https://user-images.githubusercontent.com/100873921/195411249-f98bf2bf-0899-4d78-a310-adb43de7af2e.png">
## Pre-processing Data

### Identify Isoform Switches Using R Package `IsoformSwitchAnalyzeR`

Installation of `IsoformSwitchAnalyzeR`
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IsoformSwitchAnalyzeR")
```
loading the R package
```{r}
library(IsoformSwitchAnalyzeR)
```
