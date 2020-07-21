# Citation
Constanze Ciavarella, Neil Ferguson, (2020), Deriving fine-scale models of human mobility from aggregated origin-destination flow data, Imperial College London, under review

---
# Description

## Problem this paper seeks to address
Flow data is often only available at administrative unit level; however, IBMs and other spatial simulations are implemented at a much finer scale (down to 100m).

## Aim
Fit mathematical models of human mobility to flow data, allowing the spatial resolution of the mobility model to be higher than that of the data itself

## Data
- CDR-based admin-level flow data from Kenya and Namibia
- population density per square km

## Models
Mathematical models of human mobility (adapted to trips between overnight locations)
- gravity model (one version)
- radiation models (several versions)

## Methods
Implement human mobility models at different spatial scales (administartive unit level or square patches created by grids of different sizes): if the resolution of the model is higher than that of the data, aggregated model estimates up in order to fit them to the available data

## Terminology
- admin: administrative unit of the level present in the flow data
- pixel: ~1km^2 area for which we have LandScan (ambient) population data
- cell:  area that result from joining pixels into square shapes or admin shapes

---
# Requirements

## R
see `code/R/main.R` for packages to install. Install `geos` on Mac using `homebrew` as follows
```bash
brew install geos
```

## C
Install `gcc` and `OpenMP` on Mac using `homebrew` as follows
```bash
brew install gcc
brew install libomp
```

---
# Repo structure

## code/
- `main.R`:  run complete analysis
- `C/`:      C code for simulation
- `R/`:      R functions for analyis
- `Rutils/`: R utility functions

## data/
- `raw/`:        original unmodified data
- `input/`:      data cleaned for input into C simulation
- `processing/`: output from C simulation
- `results/`:    cleaned results from C simulation
- `figures/`:    figures for reference and publication
- `tables/`:     tables for reference and publication

---
# Data sources

**N.B. No data is provided as part of this repository!**

Shapefiles are available from GADM at:
- https://gadm.org/download_country_v3.html (Kenya)
- https://gadm.org/download_country_v3.html (Namibia)

Flow data is available from:
- Wesolowski, 2012, Quantifying the Impact of Human Mobility on Malaria, doi:10.1126/science.1223467 (Kenya)
- Ruktanonchai, 2016, Identifying Malaria Transmission Foci for Elimination Using Human Mobility Data, doi:10.1371/journal.pcbi.1004846 (Namibia)

Ambient population was retrieved from Landscan (https://landscan.ornl.gov).

---
# Data format used in code
Below a description of the format of the input data to the code in this repository. This should make it easier to understand the workings of the raw data cleaning functions.

## CDR flow data
### Kenya
* File:
`data/raw/Kenya_flow_data.csv`
* Source:
Wesolowski, 2012, Quantifying the Impact of Human Mobility on Malaria
* Retrieval date:
18.10.2016 (attachment to above paper)
* Description:
Administrative unit 1 levels, refers to period from June 2008 to June 2009 but is temporally aggregated,
14'816'512 subscribers (32% of population),
data is given as counts between administrative units,
no within-unit travel, lacks directionality (we don't know the subscribers' unit of residence)
* Format:
[1] "From"    : string,  admin 1 of origin (for a single trip, not residency)
[2] "To"      : string,  admin 1 of destination (for a single trip)
[3] "Amt"     : integer, observed count of trips from said origin to said destination
[4] "PopFrom" : integer, population resident at the origin of the trip
[5] "PopTo"   : integer, population resident at the destination of the trip
[6] "Dist"    : double,  distance between origin and destination
[7] "xFrom"   : double,  x-coordinate of the origin
[8] "yFrom"   : double,  y-coordinate of the origin
[9] "xTo"     : double,  x-coordinate of the destination
[10] "yTo"    : double,  y-coordinate of the destination
* Notes:
N.B. coordinates refer to lower left corner of grid square
### Namibia
* File:
`data/raw/Namibia_flow_data.csv`
* Source:
Ruktanonchai, 2016, Identifying Malaria Transmission Foci for Elimination Using Human Mobility Data
* Retrieval date:
18.10.2016 (attachment to above paper)
* Description:
Administrative unit 2 levels, refers to period from October 2010 to September 2011 but is temporally aggregated,
1.19 million subscribers (54% of population),
data is given as proportion of nights residents of a specific admin 2 spent in each unit (including their residential unit),
records only trips with overnight stays (authors were interested in malaria transmission)
* Format:
  * matrix,
  * column [1] (rownames) : string, admin 2 of residency
  * columnnames not given, but same as rownames: admin 2 of destination
  * entry [i][j]: double, proportion of nights residents of admin i spent in admin j

## Administrative unit codes
### Kenya
* File:
`data/raw/Kenya_names.txt`
* Source:
GADM
* Retrieval date:
02.02.2017
* Description:
administrative unit 1 level codes together with country and admin 1 names
* Format:
  * [1] : character, admin unit 1 level code
  * [2] : character, country name
  * [3] : character, admin unit 1 level name
### Namibia
* File:
`data/raw/Namibia_names.txt`
* Source:
GADM
* Retrieval date:
29.11.2017
* Description:
administrative unit 2 level codes together with country, admin 1 names and admin 2 names
* Format:
  * [1] : character, admin unit 1 level code
  * [2] : character, country name
  * [3] : character, admin unit 1 level name
  * [4] : character, admin unit 2 level name

## Ambient population by square km (approx)
### Kenya
* File:
`data/raw/Kenya_pop.txt`
* Source:
LandScan and GADM
* Retrieval date:
02.02.2017
* Description:
LandScan dataset has a resolution of 30 arc seconds (approximately 1km at the equator) and counts
ambient rather than resident population
* Format:
  * [1] : double, longitude
  * [2] : double, latitude
  * [3] : integer, population of 2005
  * [4] : -- (ignore) --
  * [5] : character, administrative unit 1 level code
* Notes:

### Namibia
File:
`data/raw/Namibia_pop.txt`
* Source:
LandScan and GADM
* Retrieval date:
29.11.2016
* Description:
LandScan dataset has a resolution of 30 arc seconds (approximately 1km at the equator) and counts
ambient rather than resident population
* Format:
  * [1] : double, longitude
  * [2] : double, latitude
  * [3] : integer, population of 2010*
  * [4] : character, administrative unit 2 level code
