# Chestnut-breasted Whiteface — Occupancy & Climate Projection Analysis

Supplementary code repository to reproduce analysis in the manuscript:

Ongoing climate change drives a distributional contraction in a range-restricted Australian dryland-specialist bird

DOI: tbd

This analysis looks at chestnut-breasted whiteface (*Aphelocephala pectoralis*) occupancy dynamics and projected distribution shifts under climate change. The analysis combines Bayesian single-season occupancy models (fitted in JAGS) with Boosted Regression Trees to model historical occupancy and project future distributions under three SSP scenarios using NARCliM2.0 regional climate projections.

------------------------------------------------------------------------

## Repository structure

```         
cbwf-distribution/
├── scripts/
│   ├── 1.1_CBWF_cov_proccess_MODIS_FC.R
│   ├── 1.2_CBWF_cov_proccess_temperature.R
│   ├── 1.3_CBWF_cov_proccess_rainfall.R
│   ├── 1.4_CBWF_cov_proccess_final_outputs.R
│   ├── 2.1_occu_null.jags
│   ├── 2.2_occu_global.jags
│   ├── 2.3_CBWF_occupancy.R
│   ├── 3.1_CBWF_BRT_data_AGCD.R
│   ├── 3.2_CBWF_BRT_data_narclim.R
│   ├── 3.3_CBWF_BRT_narclim_ensemble_averaging.R
│   ├── 3.4_CBWF_BRT_analysis.R
    └── functions/  #helper functions used throughout the scripts
│
├── data/
│   ├── table/
│   │   ├── CBWF_surveys.csv
│   │   ├── CBWF_site_locations.csv
│   │   └── CBWF_records_sdm.csv
│   ├── raster/
│   │   ├── agcd/
│   │   │   ├── max_temp/
│   │   │   └── rainfall/
│   │   ├── Guerschman_FC_Monthly/
│   │   ├── narclim/
│   │   │   ├── ACCESS-ESM1-5/
│   │   │   ├── EC-Earth3-Veg/
│   │   │   ├── MPI-ESM1-2-HR/
│   │   │   ├── NorESM2-MM/
│   │   │   └── UKESM1-0-LL/
│   │   └── NVIS/
│   │       └── NVIS6_0_AUST_EXT_MVG_ALB.tif
│   └── vector/
│       └── ibra_subregions.gpkg
│
├── outputs/
│   ├── raster/
│   │   ├── agcd/
│   │   ├── guerschman_fc_monthly/
│   │   └── narclim/
│   │       ├── ensemble/
│   │       ├── predictions/
│   │       └── rainfall/
│   ├── table/
│   └── vector/
│
└── figures/
```

## Data requirements

Most of the data needed to run this analysis is already included in the repository under `data/`. The exceptions are the large climate raster datasets (AGCD, MODIS Fractional Cover, and NARCliM2.0), which are downloaded automatically by the scripts themselves the first time they are run — there is nothing to download manually. Simply run the scripts in order (see [Analysis workflow](#analysis-workflow) below) and the repository will populate itself as it goes.

### Included in the repository

**Survey data** (`data/table/`):

| File                      | Description                                  |
|---------------------------|----------------------------------------------|
| `CBWF_surveys.csv`        | Field survey detection/non-detection records |
| `CBWF_site_locations.csv` | Survey site coordinates                      |
| `CBWF_records_sdm.csv`    | Species presence records for SDM             |

**Spatial vector data** (`data/vector/`):

| File | Description | Source |
|------------------------|------------------------|------------------------|
| `ibra_subregions.gpkg` | IBRA bioregional boundaries | [DCCEEW](https://www.dcceew.gov.au/environment/land/nrs/science/ibra) |

**NVIS** (`data/raster/NVIS/`):

| File | Description | Source |
|------------------------|------------------------|------------------------|
| `NVIS6_0_AUST_EXT_MVG_ALB.tif` | National Vegetation Information System — Major Vegetation Groups | [DCCEEW NVIS data portal](https://www.dcceew.gov.au/environment/land/vegetation/national-vegetation-information-system) |

### Downloaded automatically by the scripts

These climate raster datasets are large and are **not** stored in the repository. Scripts `1.1`–`3.2` download them on demand into the corresponding `data/raster/` subfolders (skipping files that already exist, so re-running is safe):

**AGCD** (`data/raster/agcd/`) — Australian Gridded Climate Data, downloaded from the [NCI THREDDS server](https://thredds.nci.org.au/thredds/catalog/zv2/agcd/v1/catalog.html) by scripts `1.2, 1.3` and `3.1`:

- `max_temp/` — Daily maximum temperature (`agcd_v1_tmax_mean_r005_daily_YYYY.nc`)
- `rainfall/` — Monthly total precipitation (`agcd_v2_precip_total_r001_monthly_YYYY.nc`)

**MODIS Fractional Cover** (`data/raster/Guerschman_FC_Monthly/`) — Monthly PV/NPV/bare cover, downloaded from the [NCI THREDDS server](https://thredds.nci.org.au/thredds/catalog/ub8/au/FractCover/FC.v3.2.0/catalog.html) and processed in script `1.1`.

**NARCliM2.0** (`data/raster/narclim/`) — Regional climate model projections for five global climate models (ACCESS-ESM1-5, EC-Earth3-Veg, MPI-ESM1-2-HR, NorESM2-MM, UKESM1-0-LL) across three SSP scenarios (ssp126, ssp245, ssp370) and two time periods (2040–2059, 2080–2099), downloaded from the NCI THREDDS server (NSW Government NARCliM2.0 products) by script `3.2`. Helper functions for working with these files are in `scripts/functions/narclim_functions.R`.

------------------------------------------------------------------------

## Analysis workflow {#analysis-workflow}

Scripts are numbered to indicate execution order. Running them in sequence downloads the required climate rasters and populates `data/` and `outputs/` as it goes:

| Script | Purpose |
|------------------------------------|------------------------------------|
| `1.1` | Process MODIS Fractional Cover: extract PV, NPV, and bare cover change (2001–2024) |
| `1.2` | Process AGCD temperature: compute extreme heat day trends (1990–2024) |
| `1.3` | Process AGCD rainfall for occupancy modelling |
| `1.4` | Consolidate and finalise all covariates |
| `2.1` / `2.2` | JAGS model definitions (null and global occupancy models) |
| `2.3` | Fit Bayesian single-season occupancy models; posterior predictive checks |
| `3.1` | Prepare historical climate covariates from AGCD (three time periods) |
| `3.2` | Prepare NARCliM2.0 future climate covariates |
| `3.3` | Ensemble-average NARCliM2.0 predictions across the five climate models |
| `3.4` | Fit BRT; predict to historical periods and future SSP scenarios; produce figures |
