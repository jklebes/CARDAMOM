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


| Model Name                       | Alternate Names     | Short Description                                   | Details (POC)                   | Status      |
|----------------------------------|---------------------|-----------------------------------------------------|-------------------------------|----------|
| DALEC_A1C1D1F1P1R1               | DALEC_CDEA_LU_FIRES, C1, M1 |                                             | Bloom & Williams et al., 2015 | Published   |
| ...                              | ...                 |                                                     | Bloom et al., 2016            | Published   |

#### TABLE 2. Assimilation due to photosynthic activity sub-models

| Model Name                       | Short Description                                                          | Details                   |
|----------------------------------|----------------------------------------------------------------------------|---------------------------|
| A1                               | Aggregated canopy model for canopy scale photosynthetic activity           | Published (Williams et al., 1997)  |
| A2                               | Aggregated canopy model for canopy scale photosynthesis and plant ~ soil water cycle | Published (Smallman & Williams 2019)  |

#### TABLE 3. Carbon pools represented and connected structure sub-models

| Model Name                       | Short Description                                                          | Details                   |
|----------------------------------|----------------------------------------------------------------------------|---------------------------|
| C1                               | 6-pool representation (labile, foliage, fine root, wood, litter, som)| Published (Bloom & Williams 2015)  |
| C2                               | 7-pool representation (labile, foliage, fine root, wood, fol+root litter, wood litter, som) | Published (Smallman et al., 2021)  |

#### TABLE 4. Disturbance due to direct human mechanical intervention sub-models

| Model Name                       | Short Description                                                          | Details                   |
|----------------------------------|----------------------------------------------------------------------------|---------------------------|
| D1                               | 6-pool representation (labile, foliage, fine root, wood, litter, som)| Published (Bloom & Williams 2015)  |
| D2                               | 7-pool representation (labile, foliage, fine root, wood, fol+root litter, wood litter, som) | Published (Smallman et al., 2021)  |

#### TABLE 5. Fire sub-models

| Model Name                       | Short Description                                                          | Details                   |
|----------------------------------|----------------------------------------------------------------------------|---------------------------|
| F1                               |  | Published (Bloom et al., 2016)  |
| F2                               |   | Published (Smallman et al., 2021)  |

#### TABLE 6. Hydrology sub-models

| Model Name                       | Short Description                                                          | Details                   |
|----------------------------------|----------------------------------------------------------------------------|---------------------------|
| H0                               | No water cycle representation                                        |                                  |
| H1                               | 3-pool water cycle, 30cm top layer, variable rooting depth layer and remainder of soil. Soil coupled to C-cycle via fine root biomass and supply ~ demand model of stomatal conductance (A2). NOTE: that this model version assumes soil water content remains at field capacity but allows supply demand mechanisms. | Published (Smallman & Williams, 2019)  |
| H1                               | 3-pool water cycle, 30cm top layer, variable rooting depth layer and remainder of soil. Soil coupled to C-cycle via fine root biomass and supply ~ demand model of stomatal conductance (A2) | Published (Smallman & Williams, 2019)  |


#### TABLE 7. Nitrogen pools represented and connected structure sub-models

| Model Name                       | Short Description                                                          | Details                   |
|----------------------------------|----------------------------------------------------------------------------|---------------------------|
| N0                               | No nitrogen cycle representation                                     |                                  |


#### TABLE 8. Phenology sub-models

| Model Name                       | Short Description                                                          | Details                   |
|----------------------------------|----------------------------------------------------------------------------|---------------------------|
| P1                               | Combined Deciduous Evergreen Analytical model. Uses a day of year based approach to canopy growth and loss| Published (Bloom & Williams et al., 2015)  |
| P2                               | Growing Season Index (GSI) Canopy growth and mortality a linear function of daylength, temperature and vapour pressure deficit. | Published (Smallman et al., 2017)  |

#### TABLE 9. Respiration (either autotrophic or heterotrophic) sub-models

| Model Name                       | Short Description                                                          | Details                   |
|----------------------------------|----------------------------------------------------------------------------|---------------------------|
| R1                               | Autotrophic respiration estimated as a parameterisable fixed fraction of photosynthesis, while heterotrophic respiration follows 1st order kinetics with exponential temperature modification | Published (Bloom & Williams et al., 2015)  |

