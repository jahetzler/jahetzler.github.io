## Welcome to my GitHub Reposatory

I'm a biologist, with a master of biosciences from Nord Universitet with an keen interest in population genetics, bioinformatics, map making and datavisualisation.

## xkcd samling map - Norwegian coast

This is my attempt at making a sampling map with data from my thesis in the style of Randall Munroe's xkcd webcomic maps. 
This map is made in R and utilises data from offical Norwegian map database: [Kartkatalogen](https://www.geonorge.no/).

- [xkcd sampling map](xkcd sampling map/xkcd_map.md) 

## Master thesis: population genetics on _Orchomenella obtusa_

This is a step-by-step guide for the bioinformatics for my master thesis. A population genetics study on the scavenging amphipod species _Orchomenella obtusa_ in fjords surrounding Saltstraumen (Skjerstadfjorden and Saltenfjorden). Here I'm trying to determine if the sill separating the fjords are acting as a genetic barrier for the selected species (spoiler alert: it's not), by using amplicon sequences of the mitochondrial gene COI and rRNA gene 18S.

Thesis link (expected to be out by late 2024): [Population genetic structure of the scavenging amphipod O. obtusa in a deep fjord system.](https://nordopen.nord.no/nord-xmlui/handle/11250/2731119) 

Pipeline for COI and 18S amplicon sequence variant calling: 
Quality control, adapter trimming, mapping, quality filtering, variant calling and creating multifasta files.

- [O. obtusa Variant Calling](Obtusa/obtusapopgen.md)

Genetic variability and differential statistics
- [O. obtusa Statistics](Obtusa/DiffSeq.md)

Haplotype network with PopART
- [O. obtusa PopART](Obtusa/PopART.md)

