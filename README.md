# Fire is Essential for Enhanced Productivity in Fire-Adapted Ecosystems

Fire plays a key role in the health of fire-adapted ecosystems. But in areas where natural fire patterns are disrupted, it’s hard to know how ecosystems will respond to fire. Our research shows that regular fire is essential for Everglades pineland resilience. Contemporary productivity, measured by NDVI, was closely linked to fire history. Interestingly, how the land was doing before a fire—its pre-fire NDVI—was the strongest determinant of how quickly ecosystems recovered from fire, whether it occurred as a prescribed or wild fire. Regular fires enhance productivity and promote rapid post-fire recovery.

## Scripts are ordered by major workflow steps:

### 1 = Data Preparation and Integration (fire perimeter files, Ladsat data, vegetation layers)

-   **1.1_AOI.R**: Produces the raster grid (Raster30x30.tif) and Everglades boundary files (EVG_bound.shp) used throughout this project.

-   **1.2_VegLayers.R**: Generates the upland vegetation layer for Everglades National Park and Big Cypress National Preserve (Uplands_raster.tif). Generates upland sample points (Sample_pts_upland_052925.shp).

-   **1.4_FireHistoryLayers.R**: Extracts raster data to upland sample pts to generate fire history dataframes (FireHistory_df.csv, FireHistory_df.RDATA, FireYears_df.csv, FireYears_df.RDATA)

-   **1.5_SamplingR.R**: Removes unburned locations (Burned_smpl_pts_updated.shp), selects recovery sample points based on fire history (Recov_smpl_pts_updated.shp), and selects baseline sample points based on fire history (BL_smpl_pts_updated.shp).

-   **1.6a_ImageProcessing_MGM.R**: $${\color{red}RUN \space ONCE}$$. Opens tar files and saves Landsat .tif image. Images were downloaded manually from Earth Explorer: https://earthexplorer.usgs.gov/. 

-   **1.6b_ImageProcessing_MGM.R**: Extracts spectral data from tile-specific stacks to upland sample points and filters data for clouds and QAQC. Merges all tile-specific dataframes to make a Spectral Master dataframe and calculates spectral indices (Spec_Master_updated.RDATA).

-   **1.7a_DAYMET_MGM.R**: $${\color{red}RUN \space ONCE}$$. Downloads monthly climate data from DAYMET.

-   **1.7b_DAYMET_MGM.R**: Stacks monthly layers and masks stacks to the Everglades boundary (precip_EVG_updated.tif, tmax_EVG_updated.tif, tmin_EVG_updated.tif)

-   **1.8_Severity_MGM.R**: Calculates fire severity as the change in NBR pre- and post-fire for the Recovery Sample Points (Severity_updated.csv, Sev_df_updated.RDATA)

-   **1.9_LandscapeSummary_MGM.R**: Summarizes the fire history (Landscape_Summary_updated.RDATA) and climate patterns (Upland_DAYMET_updated.RDATA, Uplands_PDSI1_updated.RDATA, Uplands_PDSI2_updated.RDATA, Uplands_PDSI3_updated.RDATA, mean_pdsi_combo_updated.RDATA) for the landscape across pineland communities during the study period.

### 2 = The NDVI Model Workflow (Develop the dataframe used for the baseline NDVI model and make the model)

### 3 = Recovery Dataset Development (Develop the dataframe used for the driver analysis)

### 4 = Analysis Workflow (Driver Analysis for the recovery dataframes)

5 = Data visualization
