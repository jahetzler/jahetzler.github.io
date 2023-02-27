#### Differential statistics in R

Differential statistic test included here are F-statistics based on different measures of Fst and AMOVA.

[article](https://popgen.nescent.org/PopDiffSequenceData.html)

R packages used in this analysis 
```
library("apex")
library("adegenet")
library("pegas")
library("mmod")
library("poppr")
library("radiator")
```

Import Multifasta files

```
#Import multifasta 
data <- read.multiFASTA(c("COI_var_st.fas", "18S_var_st.fas"))
```

Rename locus names to COI and 18S. 
```
getLocusNames(data)
(setLocusNames(data) <- gsub("_var_st.fas", "", getLocusNames(data)))
```

#Create genind object
```
data.gid <- multidna2genind(data, mlst = TRUE)
ploidy(data.gid)<-1 # haploid mtDNA
data.gid

/// GENIND OBJECT /////////

 // 200 individuals; 2 loci; 77 alleles; size: 95.6 Kb

 // Basic content
   @tab:  200 x 77 matrix of allele counts
   @loc.n.all: number of alleles per locus (range: 7-70)
   @loc.fac: locus factor for the 77 columns of @tab
   @all.names: list of allele names for each locus
   @ploidy: ploidy of each individual  (range: 1-1)
   @type:  codom
   @call: df2genind(X = xdfnum, ind.names = x@labels, ploidy = 1)

 // Optional content
   - empty -
```

Genind object include the 200 sampled individuals, 2 loci (COI and 18S), and 77 variable sites (COI: 70, 18S: 7). <br/>
We will add sampling sites and fjords to the optional content in the genind object.

```
my_strata <- data.frame(regions = rep(c("SAL", "SKJ"), each = 100), 
                        populations = rep(c("A", "B", "C", "D"), each = 50))
strata(data.gid) <- my_strata
setPop(data.gid) <- ~populations
data.gid

// Optional content
   @pop: population of each individual (group size range: 50-50)
   @strata: a data frame with 2 columns ( regions, populations )

```
With this we now have the additional information from where the samples are from and we can run some overall differential statistics test on the data. 

```
diff_stats(data.gid) 

$per.locus
            Hs        Ht         Gst  Gprime_st          D
COI  0.8802020 0.8832005 0.003395022 0.03774336 0.03337268
X18S 0.5951515 0.6125879 0.028463449 0.09286085 0.05742515

$global
        Hs         Ht    Gst_est  Gprime_st      D_het     D_mean 
0.73767677 0.74789419 0.01366159 0.06912418 0.05193300 0.04221315 
```

we are seeing for each locus and global calculations for within-group heterozygosity (HS); total heterozygosity (Ht); Nei's Gst (Gst);
Hedrick's G”st (Gprime_st), Jost's D (D). 

A more approprate test for mtDNA would be the AMOVA based meirmans Φst

```
Phi_st_Meirmans(data.gid) 

$per.locus
        COI        X18S 
-0.05086222  0.06964849 

$global
[1] 0.04581864

```

Pairwise global Fst tests for Nei's Gst, Hedrick's Gst, and Joost's D

```
# Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- Gst)  
pairwise_Gst_Nei(data.gid, linearized = FALSE) 

A           B           C
B 0.007977253                        
C 0.010071518 0.010380902            
D 0.012177607 0.007049356 0.007292493

# Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- Gst')  
pairwise_Gst_Hedrick(data.gid, linearized = FALSE)
A          B          C
B 0.05237285                      
C 0.07427677 0.07977650           
D 0.09009671 0.05465308 0.06509810

# Calculates pairwise Gst. If linearized = TRUE, it calculates 1/(1- D)  
pairwise_D(data.gid, linearized = FALSE, hsht_mean = "harmonic") 

A          B          C
B 0.04885409                      
C 0.06599322 0.07305464           
D 0.08164839 0.04458731 0.05565753
```
*Global pairwise test for rDNA and mtDNA might not be the best idea, so do individual runs with COI and 18S*


