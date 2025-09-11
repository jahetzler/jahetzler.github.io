## Welcome to my GitHub Reposatory

I'm a biologist, with a master in biosciences from Nord University, and an keen interest in population genetics, bioinformatics, data analysis, map making and datavisualisation.

## Shiny R - Visualisering av partifinansiering (updated autumn 2025)
- [Partifinansiering](https://jhetzler.shinyapps.io/Partifinansiering)

This interactive dashboard visualises the financial records of Norwegian political parties from 2017 to 2024 using dynamic Sankey diagrams that make money flows easy to read. The underlying data is drawn from the official Party Funding portal (partifinansiering.no), which publishes reported donations, accounting figures and metadata required under the Political Parties Act.


## Shiny R: SHUe - Statistikk for Høyere Utdanning (updated for spring 2025)

- [SHUe](https://jhetzler.shinyapps.io/SHUe)

Dashboard project for visualizing openly available data for higher education in Norway. For simplicity sake the data is subset to include only one univeristy (NORD), but if there is any interest or I have time I might add other universities in Norway.

General information about universities are gathered from [Database for statistikk om høyere utdanning](https://dbh.hkdir.no/), for applications related to Norwegian Resaerch Council I've used the data available at [data.norge.no](https://data.norge.no/datasets/d23bbbfa-fdad-4dab-b31b-4daadbfa3299), and for EU Horizon projects I've gathered data from [data.europa.eu](https://data.europa.eu/data/datasets/cordis-eu-research-projects-under-horizon-europe-2021-2027?locale=en) 

## Python - Pydeck: Fangstdata

Pydeck test to plot large datasets (700K points) on a map. Map showing where fishing boats have unloaded and from where it was caught. Example map with fewer points uploaded as example.

- [FDIR](FDir/pydeck/fangstdata.md)

## R/Python: Fiskeridirektoratet

Ongoing project from fiskeridirektoratets on ERS raports with records of caught fish. Mainly used for exploring different graphical modules in R and python, but examples could be of use to others. 

- [FDIR](FDir/FDir.md)

## R: Tolkien map from open source data

Trying to make a map over different parts of norway in the style of J.R.R Tolkien utilizing data from the offical Norwegian map database: [Kartkatalogen](https://www.geonorge.no/).

Example image, code and data for vestvågøy included here.

- [Tolkien map](TolkienMap/tolkienMap.md)

## R: xkcd samling map - Norwegian coast

This is my attempt at making a sampling map with data from my thesis in the style of Randall Munroe's xkcd webcomic maps. 
This map is made in R and utilises data from the offical Norwegian map database: [Kartkatalogen](https://www.geonorge.no/).

- [xkcd sampling map](xkcdSamplingMap/xkcd_map.md) 

## Master thesis: population genetics on _Orchomenella obtusa_

This is a step-by-step guide for the bioinformatics for my master thesis. A population genetics study on the scavenging amphipod species _Orchomenella obtusa_ in fjords surrounding Saltstraumen (Skjerstadfjorden and Saltenfjorden). Here I'm trying to determine if the sill separating the fjords are acting as a genetic barrier for the selected species (spoiler alert: it's not), by using amplicon sequences of the mitochondrial gene COI and rRNA gene 18S.

Pipeline for COI and 18S amplicon sequence variant calling: 
Quality control, adapter trimming, mapping, quality filtering, variant calling and creating multifasta files.

- [O. obtusa Variant Calling](Obtusa/obtusapopgen.md)

Genetic variability and differential statistics
- [O. obtusa Statistics](Obtusa/DiffSeq.md)

Haplotype network with PopART
- [O. obtusa PopART](Obtusa/PopART.md)

