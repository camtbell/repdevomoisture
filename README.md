This is the information file (README) for the data and code used in the paper "The effect of moisture during development on phenotypes of egg-laying reptiles: A systematic review and meta-analysis" that is provided in this OSF repository. This paper was published in TBD (DOI). The following explains each file organized into a series of folders:

1_Overall Database
-------------------
This folder contains the file named ***database.xlsx***, which is a XLS spreadsheet containing all the original, unpaired data that was extracted from 34 papers included in the meta-analysis. This XLS spreadsheet contains two worksheets. The first contains the raw, extracted data and is called 'database'. The second sheet is named 'metadata' and holds an explanation for all the columns and variables in the database.


2_Code 
--------
This folder contains all our code that we used to run the MLMA and meta-regression for each phenotypic trait included in our meta-analysis. The files are:

- incubation_model.R contains the code for our models of the incubation duration trait and refers to the incubation_data.csv (in the '3_trait data' folder) and incubation_species.csv (in the '4_species data' folder) datasets

- length_model.R contains the code for our models of the length trait and refers to the length_data.csv (in the '3_trait data' folder) and length_species.csv (in the '4_species data' folder) datasets 

- mass_model.R contains the code for our models of the mass trait and refers to the mass_data.csv (in the '3_trait data' folder) and mass_species.csv (in the '4_species data' folder) datasets 

- sex_model.R contains the code for our models of the sex ratio trait and refers to the sex_data.csv (in the '3_trait data' folder) and sex_species.csv (in the '4_species data' folder) datasets

- survival_model.R contains the code for our models of the hatching success (survival) trait and refers to the survival_data.csv (in the '3_trait data' folder) and survival_species.csv (in the '4_species data' folder) datasets'

There are also a series of script files for our figures in R. These are: 
- figures 2 and 4_modelOR.r that contains the code for Figures 2 and 4
- figure 5_modelRR.r that contains the code for Figure 5
- figure 4_contour plot that contains the code for Figure 4

There is also a file called "conversion_example.r" which we made use of for our interpretation of the moderators on effect sizes in our study.


3_Trait Data 
-------------
This folder contains all the trait datasets we used in our meta-analyses. Each trait was analysed seperately. All this data has already been paired as described in the manuscript. This folder contains the following files:

- incubation_data.csv that contains all paired trait data for incubation duration

- length_data.csv that contains all paired trait data for length

- mass_data.csv that contains all paired trait data for mass

- sex_data.csv that contains all paired trait data for sex

- survival_data.csv that contains all paired trait data for hatching success (e.g., survival)

- trait_metadata.csv a metadata file explaining all columns and variables for the trait datasets. These are consistent across all 5 datasets.


4_Species Data
---------------
This folder contains datasets with genus and species names of all reptiles present within the appropriate trait datasets. It contains the following files:

- incubation_species.csv contains genus and species names of all reptiles present in the incubation duration dataset (incubation_data.csv in the 3_trait data folder)

- length_species.csv contains genus and species names of all reptiles present in the length dataset (length_data.csv in the 3_trait data folder)

- mass_species.csv contains genus and species names of all reptiles present in the mass dataset (mass_data.csv in the 3_trait data folder)

- sex_species.csv contains genus and species names of all reptiles present in the sex dataset (sex_data.csv in the 3_trait data folder)

- survival_species.csv contains genus and species names of all reptiles present in the survival dataset (survival_data.csv in the 3_trait data folder)


5_Rep Moist Meta.Rproj
-----------------------
This is an R Project file that can be opened with R Studio and used to link the code and data files together for review. 
