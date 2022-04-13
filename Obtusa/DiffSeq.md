#### Differential statistics

```
https://popgen.nescent.org/PopDiffSequenceData.html

library("apex")
library("adegenet")
library("pegas")
library("mmod")
library("poppr")
library("radiator")

#Import multifasta 
data <- read.multiFASTA(c("COI_var_st.fas", "18S_var_st.fas"))
````


````
getLocusNames(data)
(setLocusNames(data) <- gsub("_HL_ST.fas", "", getLocusNames(data)))
```

