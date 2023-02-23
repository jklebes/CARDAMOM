## CARDAMOM Coupled Models

This files provides information on the currently coupled models within CARDAMOM. Table 1 provides a summary of the model name, alternate names used in publications and / or source code, publication containing the first or most complete description and status (e.g. published, under development, abandoned). The model name is the combination of the nominal model family (e.g. DALEC) and classifiers indicating the specific sub-models used. If an available classifier is not used in a given models name, please assume that option "0" is implicit. Sub-model specific information in summaries in subsequent tables.

The current classifiers used are: 
A = Assimilation (Table 2)
C = C-cycle (Table 3)
D = Manual disturbance (Table 4)
F = Fire (Table 5)
H = Hydrology (Table 6)
N = N-cycle (Table 7)
P = Phenology (Table 8)
R = Respiration (Table 9)
M = Managed Grassland / arable crops (Table 10)

#### TABLE 1. CARDAMOM models.

|Ref| Model (Code) Name        | Alternate Names     | Short Description                                   | Details (POC)                 | Status    |
|---|--------------------------|---------------------|-----------------------------------------------------|-------------------------------|-----------|
| 1 | DALEC_A0C0D1F2H0M0N0P0R0 | DALEC_CDEA_LU_FIRES, C1, M1           | Fire model updated from F1 to F2| Bloom & Williams et al., 2015| Published |
| 2 | DALEC_A1D2F2H1        | DALEC_CDEA_ACM2, C6, M2                  | Fire model updated from F1 to F2| Smallman et al., 2021        | Published |
| 3 | DALEC_A1D2F2H2        | DALEC_CDEA_ACM2_BUCKET, C7, M3           | Fire model updated from F1 to F2| Smallman et al., 2021        | Published |
| 4 | DALEC_A1D2F2H2R1      | DALEC_CDEA_ACM2_BUCKET_RmRg, C10, M4     | Fire model updated from F1 to F2| Smallman et al., 2021        | Published |
| 5 | DALEC_A1C1D2F2H2R1    | DALEC_CDEA_ACM2_BUCKET_RmRg_CWD, C11, M5 | Fire model updated from F1 to F2| Smallman et al., 2021        | Published |
| 6 | DALEC_A1C1D2F2H2P1R1  | DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT, C12| Fire model updated from F1 to F2| Smallman et al., 2021        | Published |
| 7 | DALEC_A1C1D2F2H1P2    | DALEC_GSI_DFOL_CWD_FR, G1                | Fire model updated from F1 to F2| Smallman et al., 2017        | Published |
| 8 | DALEC_A1C1D2F2H2P2    | DALEC_GSI_BUCKET, G2                     | Fire model updated from F1 to F2| Smallman et al., 2017        | Published |


#### TABLE 2. Assimilation due to photosynthic activity sub-models

| Model Name | Short Description                                                          | Details                   |
|------------|----------------------------------------------------------------------------|---------------------------|
| A0         | Aggregated canopy model for canopy scale photosynthetic activity           | Published (Williams et al., 1997)  |
| A1         | Aggregated canopy model for canopy scale photosynthesis and plant ~ soil water cycle | Published (Smallman & Williams 2019)  |
| A2         | As A2 but with Farquhar equations | Unpublished |

#### TABLE 3. Carbon pools represented and connected structure sub-models

| Model Name | Short Description                                                          | Details                   |
|------------|----------------------------------------------------------------------------|---------------------------|
| C0         | 6-pool representation (labile, foliage, fine root, wood, litter, som) where labile supplies foliage only | Published (Bloom & Williams 2015)  |
| C1         | 7-pool representation (labile, foliage, fine root, wood, fol+root litter, wood litter, som) where labile supplies foliage only | Published (Smallman et al., 2021)  |

#### TABLE 4. Disturbance due to direct human mechanical intervention sub-models

| Model Name | Short Description                                                          | Details                   |
|------------|----------------------------------------------------------------------------|---------------------------|
| D0         | No disturbance model | Published (Bloom & Williams 2015)  |
| D1         | Foliage and wood removed as determined by input fractional loss, no litter residues | Published (Bloom & Williams 2016)  |
| D2         | All C pools undergo removal and litter residues based on a fixed number of scenarios but driven by fractional cover loss | Published (Smallman et al., 2021)  |

#### TABLE 5. Fire sub-models

| Model Name | Short Description                                                          | Details                   |
|------------|----------------------------------------------------------------------------|---------------------------|
| F1                               | Fire resiliance and tissue specific combustion factors are hardcoded       | Published (Bloom et al., 2016)  |
| F2                               | Fire resiliance and tissue specific combustion factors are parameterised   | Published (Smallman et al., 2021)  |

#### TABLE 6. Hydrology sub-models

| Model Name | Short Description                                                          | Details                   |
|------------|----------------------------------------------------------------------------|---------------------------|
| H0                               | No water cycle representation                                        |                                  |
| H1                               | 3-pool water cycle, 30cm top layer, variable rooting depth layer and remainder of soil. Soil coupled to C-cycle via fine root biomass and supply ~ demand model of stomatal conductance (A2). NOTE: that this model version assumes soil water content remains at field capacity but allows supply demand mechanisms. | Published (Smallman & Williams, 2019)  |
| H2                               | 3-pool water cycle, 30cm top layer, variable rooting depth layer and remainder of soil. Soil coupled to C-cycle via fine root biomass and supply ~ demand model of stomatal conductance (A2) | Published (Smallman & Williams, 2019)  |

#### TABLE 7. Nitrogen pools represented and connected structure sub-models

| Model Name | Short Description                                                          | Details                   |
|------------|----------------------------------------------------------------------------|---------------------------|
| N0         | No nitrogen cycle representation                                     |                                  |


#### TABLE 8. Plant Phenology (inc. foliage, wood and fine roots) sub-models

| Model Name | Short Description                                                          | Details                   |
|------------|----------------------------------------------------------------------------|---------------------------|
| P0         | Combined Deciduous Evergreen Analytical model. Uses a day of year based approach to canopy growth and loss. Fixed fraction for NPP allocation. 1st order kinetic for non-foliage turnovers. | Published (Bloom & Williams et al., 2015) |
| P1         | As P0 except 1st order kinetics modified by a Michaelis-Menten function on wood turnover. | Published (Famiglietti et al., 2021) |
| P2         | As P0 except CDEA replaced by Growing Season Index (GSI) Canopy growth and mortality a linear function of daylength, temperature and vapour pressure deficit. | Published (Smallman et al., 2017) |

#### TABLE 9. Respiration (either autotrophic or heterotrophic) sub-models

| Model Name | Short Description                                                          | Details                   |
|------------|----------------------------------------------------------------------------|---------------------------|
| R0         | Autotrophic respiration estimated as a parameterisable fixed fraction of photosynthesis, while heterotrophic respiration follows 1st order kinetics with exponential temperature modification | Published (Bloom & Williams et al., 2015)  |
| R1         | Maintenance respiration estimated as a parameterisable fixed fraction of photosynthesis, growth respiration is a hardcoded fraction of NPP. Heterotrophic respiration follows 1st order kinetics with exponential temperature modification | Published (Bloom & Williams et al., 2015)  |


