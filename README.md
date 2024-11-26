## Welcome to my GitHub Reposatory

I'm a biologist, with a master in biosciences from Nord University, and an keen interest in population genetics, bioinformatics, data analysis, map making and datavisualisation.

## Python - Pydeck: Fangstdata

Testing out pydeck to plot large datasets (700K points) on a map. Map showing where fishing boats have undocked their load and from where it was caught. 

- [FDIR](FDir/pydeck/fangstdata.md)

## R: Tolkien map from open source data

Trying to make a map over different parts of norway in the style of J.R.R Tolkien using map data from geonorge.
Example image, code and data for vestvågøy included here.

- [Tolkien map](TolkienMap/tolkienMap.md)

## Shiny R: SHUe - Statistikk for Høyere Utdanning'e

- [SHUe](https://jhetzler.shinyapps.io/SHUe)

Project of mine for visualizing openly available data for higher education in Norway. For simplicity sake the data is subset to include only one univeristy (NORD), but if there is any interest or I have time I might add the other universities in Norway.

General information about universities are gathered from [Database for statistikk om høyere utdanning](https://dbh.hkdir.no/), for applications related to Norwegian Resaerch Council I've used the data available at [data.norge.no](https://data.norge.no/datasets/d23bbbfa-fdad-4dab-b31b-4daadbfa3299), and for EU Horizon projects I've gathered data from [data.europa.eu](https://data.europa.eu/data/datasets/cordis-eu-research-projects-under-horizon-europe-2021-2027?locale=en) 

## R: xkcd samling map - Norwegian coast

This is my attempt at making a sampling map with data from my thesis in the style of Randall Munroe's xkcd webcomic maps. 
This map is made in R and utilises data from offical Norwegian map database: [Kartkatalogen](https://www.geonorge.no/).

- [xkcd sampling map](xkcdSamplingMap/xkcd_map.md) 

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

