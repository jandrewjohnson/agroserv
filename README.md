# Agroserv implementation of SEALS model

## Quick Start:
Creates 12 geotiffs of SA that use the following categorization:

    Forest = 1
    Savanna = 2
    EXISTING Pasture = 3
    EXISTING Cropland = 4
    Urban = 5
    Other = 6
    Water = 7
    NEW Pasture 8
    NEW Cropland 9
    ndv = 255

The 12 scenarios are named as follows. foronly means expansion happens only in forests. 

    'amazon_foronly_cerrado_nosavfor_fraglocal',
    'amazon_foronly_cerrado_savonly_fraglocal',
    'amazon_foronly_cerrado_bau_fraglocal',
    'amazon_nosavfor_cerrado_nosavfor_fraglocal',
    'amazon_nosavfor_cerrado_savonly_fraglocal',
    'amazon_nosavfor_cerrado_bau_fraglocal',
    'amazon_foronly_cerrado_nosavfor_fragregional',
    'amazon_foronly_cerrado_savonly_fragregional',
    'amazon_foronly_cerrado_bau_fragregional',
    'amazon_nosavfor_cerrado_nosavfor_fragregional',
    'amazon_nosavfor_cerrado_savonly_fragregional',
    'amazon_nosavfor_cerrado_bau_fragregional',
    
In these, the keywords mean that in the preceeding biomoe, there is expansion according to:

    foronly: expansion only in forests
    nosavfor: no expansion happens in savanna or forests
    savonly: expansion only in savanna
    bau: expansion happens without constraint on LULC type and at the historical rates.
    
The final keyword, fraglocal or fagregional, denotes how broadly distributed the projected changes are placed on the landscape. Fraglocal means that all expansion happens in the 0.25 deg grid-cell (native resolution to IPBES) while fragregional means that all expansion happens in a 6.0 deg grid-cell (created by upsampling the 0.25) while preserving spatial sums.

## Introduction
Because SEALS is still a relatively new software tool for NatCap, a specialized implementation, provided here, was required for this project. To install, just clone the repo and then run seals_model_agroserv.py, fix any import statements (there will be several, check out the requirements.txt for faster fixing) iteratively until the code runs. The only line of code that needs to be modified to make it work on the virtual machine is the one that sets the project dir:

p = hb.ProjectFlow(DIR_FOR_SCENARIOS)

Set the argument to this object to be the location of where the runtime data is (see the results I sent for an example organization).

## Data Input and Ingestion

All data was provided in agroserv_base_data.zip sent by justin on 8/15/2018. For the script to run, this folder must be extracted into the same parent directory as the code repo as follows

- parent_dir
    - agroserv_0-0-1 (code dir)
    - agroserv_base_data

To see how the data were generated, open and/or run process_agroserv_inputs.py. Note that this will require updating file paths to where the input data is stored on the local computer, but given the very detailed nature of the code, it should be very clear which file(s) from the input that Avery sent me should go where.

## Generating outputs and defining scenarios

seals_model_agroserv.py generates 10 geotiffs for Cerrado and Amazon. Once SEALS has run, you may want to combine the different zones into a single file that represents the scenarios as defined by Avery. This can be done using postprocess_agroserve_outputs.py, which generates 12 secenarios by taking permutations of the 10 maps according to the factors below. Note that you will have to update the pathnames to be specific to where you ran put the inputs and outputs.

### Factors

Cerrado deforestation scenario a. 1.5 mha per year w/ savannah-forest split based on SSP 3 b. 1.5 mha per year savannah only conversion c. zero conversion of savannah or forest

Amazon deforestation scenario a. 1.5 mha per year all forest b. zero conversion of savannah or forest

Fragmentation a. with constraint that total change correctly sums up within an LUH SSP 3 grid-cell after harmonization of quantities withing grid cell of each land use and land cover class as per our call. You wrote that this will occur at 15 minutes, but that was a typo right? The resolution is actually 0.25 degrees for these. b. with constraint that total change correctly sums up within biome

Because there are 3 cerrado deforestation scenarios but only two for the Amazon, we get 6 scenanarios by mixing and matching the two regions. Then, permuting the two fragmentation scenarios gives us 12.

This is why the files were provided separately for Amazon vs Cerrado so you can correctly mix-n-match them. This should be just a quick gdal_translate to merge the two extents, but let me know if this is a challenge and i could also provide the joined results in the next run.

## Apendices
### Conversion tables for different LULC maps

p.esa_names_to_mbs_names = OrderedDict()
p.esa_names_to_mbs_names['ndv'] = 'ndv'
p.esa_names_to_mbs_names['crop_rainfed'] = 'Cropland'
p.esa_names_to_mbs_names['crop_rainfed_herb'] = 'Cropland'
p.esa_names_to_mbs_names['crop_rainfed_tree'] = 'Cropland'
p.esa_names_to_mbs_names['crop_irrigated'] = 'Cropland'
p.esa_names_to_mbs_names['crop_natural_mosaic'] = 'Cropland'
p.esa_names_to_mbs_names['natural_crop_mosaic'] = 'Cropland'
p.esa_names_to_mbs_names['tree_broadleaved_evergreen'] = 'Forest'
p.esa_names_to_mbs_names['tree_broadleaved_deciduous_closed_to_open_15'] = 'Forest'
p.esa_names_to_mbs_names['tree_broadleaved_deciduous_closed_40'] = 'Forest'
p.esa_names_to_mbs_names['tree_broadleaved_deciduous_open_15_40'] = 'Forest'
p.esa_names_to_mbs_names['tree_needleleaved_evergreen_closed_to_open_15'] = 'Forest'
p.esa_names_to_mbs_names['tree_needleleaved_evergreen_closed_to_open_15'] = 'Forest'
p.esa_names_to_mbs_names['tree_needleleaved_evergreen_open_15_40'] = 'Forest'
p.esa_names_to_mbs_names['tree_needleleaved_deciduous_closed_to_open_15'] = 'Forest'
p.esa_names_to_mbs_names['tree_needleleaved_deciduous_closed_40'] = 'Forest'
p.esa_names_to_mbs_names['tree_needleleaved_deciduous_open_15_40'] = 'Forest'
p.esa_names_to_mbs_names['tree_mixed_type'] = 'Forest'
p.esa_names_to_mbs_names['Mosaic_tree_and_shrub_50_herbaceous_cover_50'] = 'Savanna'
p.esa_names_to_mbs_names['Mosaic_herbaceous_cover_50_tree_and_shrub_50'] = 'Savanna'
p.esa_names_to_mbs_names['Shrubland'] = 'Savanna'
p.esa_names_to_mbs_names['Evergreen_shrubland'] = 'Savanna'
p.esa_names_to_mbs_names['Deciduous_shrubland_'] = 'Savanna'
p.esa_names_to_mbs_names['Grassland'] = 'Pasture'
p.esa_names_to_mbs_names['Lichens_and_mosses'] = 'Other'
p.esa_names_to_mbs_names['Sparse_vegetation_tree_shrub_herbaceous_cover_15'] = 'Savanna'
p.esa_names_to_mbs_names['Sparse_tree_15'] = 'Other'
p.esa_names_to_mbs_names['Sparse_shrub_15'] = 'Savanna'
p.esa_names_to_mbs_names['Sparse_herbaceous_cover_15'] = 'Savanna'
p.esa_names_to_mbs_names['Tree_cover_flooded_fresh_or_brakish_water'] = 'Other'
p.esa_names_to_mbs_names['Tree_cover_flooded_saline_water'] = 'Other'
p.esa_names_to_mbs_names['Shrub_or_herbaceous_cover_flooded_fresh_saline_brakish_water'] = 'Other'
p.esa_names_to_mbs_names['Urban_areas'] = 'Urban'
p.esa_names_to_mbs_names['Bare_areas'] = 'Other'
p.esa_names_to_mbs_names['Consolidated_bare_areas'] = 'Other'
p.esa_names_to_mbs_names['Unconsolidated_bare_areas'] = 'Other'
p.esa_names_to_mbs_names['Water_bodies'] = 'Water'
p.esa_names_to_mbs_names['Permanent_snow_and_ice'] = 'Other'

p.mbs_names_to_codes = OrderedDict()
p.mbs_names_to_codes['Forest'] = 1
p.mbs_names_to_codes['Savanna'] = 2
p.mbs_names_to_codes['Pasture'] = 3
p.mbs_names_to_codes['Cropland'] = 4
p.mbs_names_to_codes['Urban'] = 5
p.mbs_names_to_codes['Other'] = 6
p.mbs_names_to_codes['Water'] = 7
p.mbs_names_to_codes['split_between_crops_and_pastures'] = 8
p.mbs_names_to_codes['ndv'] = 255


p.luh_names_to_luh_codes = OrderedDict()
p.luh_names_to_luh_codes['range'] = 1
p.luh_names_to_luh_codes['secdf'] = 2 # NOTE: Half of this will be allocated to forest, half to Savanna
p.luh_names_to_luh_codes['secdn'] = 3
p.luh_names_to_luh_codes['urban'] = 4
p.luh_names_to_luh_codes['c3ann'] = 5
p.luh_names_to_luh_codes['c3nfx'] = 6
p.luh_names_to_luh_codes['c3per'] = 7
p.luh_names_to_luh_codes['c4ann'] = 8
p.luh_names_to_luh_codes['c4per'] = 9
p.luh_names_to_luh_codes['pastr'] = 10
p.luh_names_to_luh_codes['primf'] = 11
p.luh_names_to_luh_codes['primn'] = 12


p.luh_names_to_mbs_names = OrderedDict()
p.luh_names_to_mbs_names['range'] = 'Pasture'
p.luh_names_to_mbs_names['secdf'] = 'Forest'
p.luh_names_to_mbs_names['secdn'] = 'Savanna'
p.luh_names_to_mbs_names['urban'] = 'Urban'
p.luh_names_to_mbs_names['c3ann'] = 'Cropland'
p.luh_names_to_mbs_names['c3nfx'] = 'Cropland'
p.luh_names_to_mbs_names['c3per'] = 'Cropland'
p.luh_names_to_mbs_names['c4ann'] = 'Cropland'
p.luh_names_to_mbs_names['c4per'] = 'Cropland'
p.luh_names_to_mbs_names['pastr'] = 'Pasture'
p.luh_names_to_mbs_names['primf'] = 'Forest'
p.luh_names_to_mbs_names['primn'] = 'Savanna'


p.luh_names_to_filenames = OrderedDict()
p.luh_names_to_filenames['range'] = 'range ^ area_fraction ^ rangeland.tif'
p.luh_names_to_filenames['secdf'] = 'secdf ^ area_fraction ^ potentially forested secondary land.tif'
p.luh_names_to_filenames['secdn'] = 'secdn ^ area_fraction ^ potentially non-forested secondary land.tif'
p.luh_names_to_filenames['urban'] = 'urban ^ area_fraction ^ urban land.tif'
p.luh_names_to_filenames['c3ann'] = 'c3ann ^ area_fraction ^ C3 annual crops.tif'
p.luh_names_to_filenames['c3nfx'] = 'c3nfx ^ area_fraction ^ C3 nitrogen-fixing crops.tif'
p.luh_names_to_filenames['c3per'] = 'c3per ^ area_fraction ^ C3 perennial crops.tif'
p.luh_names_to_filenames['c4ann'] = 'c4ann ^ area_fraction ^ C4 annual crops.tif'
p.luh_names_to_filenames['c4per'] = 'c4per ^ area_fraction ^ C4 perennial crops.tif'
p.luh_names_to_filenames['pastr'] = 'pastr ^ area_fraction ^ managed pasture.tif'
p.luh_names_to_filenames['primf'] = 'primf ^ area_fraction ^ forested primary land.tif'
p.luh_names_to_filenames['primn'] = 'primn ^ area_fraction ^ non-forested primary land.tif'




p.mb_names_to_mb_codes = OrderedDict()
p.mb_names_to_mb_codes['Forest'] = 1
p.mb_names_to_mb_codes['Natural_Forest_Formations'] = 2
p.mb_names_to_mb_codes['Dense_Forest'] = 3
p.mb_names_to_mb_codes['Open_Forest'] = 4
p.mb_names_to_mb_codes['Mangrove'] = 5
p.mb_names_to_mb_codes['Flooded_forest'] = 6
p.mb_names_to_mb_codes['Degraded_Forest'] = 7
p.mb_names_to_mb_codes['Secondary_Forest'] = 8
p.mb_names_to_mb_codes['Forestry'] = 9
p.mb_names_to_mb_codes['Non_Forest_Natural_Formations'] = 10
p.mb_names_to_mb_codes['Non_forest_Natural_Wetlands'] = 11
p.mb_names_to_mb_codes['Vegetation_Country_'] = 12
p.mb_names_to_mb_codes['Other_non_forest_formations'] = 13
p.mb_names_to_mb_codes['Agricultural_Use'] = 14
p.mb_names_to_mb_codes['Pasture'] = 15
p.mb_names_to_mb_codes['Pasture_in_Natural_Fields'] = 16
p.mb_names_to_mb_codes['Other_Pastures'] = 17
p.mb_names_to_mb_codes['Agriculture'] = 18
p.mb_names_to_mb_codes['Annual_Cultures'] = 19
p.mb_names_to_mb_codes['Semi_Perennial_Cultures_'] = 20
p.mb_names_to_mb_codes['Crop_Mosaic'] = 28
p.mb_names_to_mb_codes['Agriculture_or_Grassland'] = 21
p.mb_names_to_mb_codes['not_vegetated'] = 22
p.mb_names_to_mb_codes['Beaches_and_dunes'] = 23
p.mb_names_to_mb_codes['Urban_infrastructure'] = 24
p.mb_names_to_mb_codes['Other_non_vegetated_areas'] = 25
p.mb_names_to_mb_codes['Water_bodies'] = 26
p.mb_names_to_mb_codes['ndv'] = 27

p.mbs_names_to_codes = OrderedDict()
p.mbs_names_to_codes['Forest'] = 1
p.mbs_names_to_codes['Savanna'] = 2
p.mbs_names_to_codes['Pasture'] = 3
p.mbs_names_to_codes['Cropland'] = 4
p.mbs_names_to_codes['Urban'] = 5
p.mbs_names_to_codes['Other'] = 6
p.mbs_names_to_codes['Water'] = 7
p.mbs_names_to_codes['split_between_crops_and_pastures'] = 8
p.mbs_names_to_codes['ndv'] = 255

p.mb_names_to_mbs_names = OrderedDict()
p.mb_names_to_mbs_names['Forest'] = 'Other' # NOTE that THIS categorization of Forest is because it shouldnt be in the data and all should be as Dense Forest. Included for completeness wrt mapbiomas classification
p.mb_names_to_mbs_names['Natural_Forest_Formations'] = 'Other'
p.mb_names_to_mbs_names['Dense_Forest'] = 'Forest'
p.mb_names_to_mbs_names['Open_Forest'] = 'Savanna'
p.mb_names_to_mbs_names['Mangrove'] = 'Other'
p.mb_names_to_mbs_names['Flooded_forest'] = 'Other'
p.mb_names_to_mbs_names['Degraded_Forest'] = 'Other'
p.mb_names_to_mbs_names['Secondary_Forest'] = 'Other'
p.mb_names_to_mbs_names['Forestry'] = 'Other'
p.mb_names_to_mbs_names['Non_Forest_Natural_Formations'] = 'Other'
p.mb_names_to_mbs_names['Non_forest_Natural_Wetlands'] = 'Other'
p.mb_names_to_mbs_names['Vegetation_Country_'] = 'Other'
p.mb_names_to_mbs_names['Other_non_forest_formations'] = 'Other'
p.mb_names_to_mbs_names['Agricultural_Use'] = 'Other'
p.mb_names_to_mbs_names['Pasture'] = 'Pasture'
p.mb_names_to_mbs_names['Pasture_in_Natural_Fields'] = 'Other'
p.mb_names_to_mbs_names['Other_Pastures'] = 'Other'
p.mb_names_to_mbs_names['Agriculture'] = 'Cropland'
p.mb_names_to_mbs_names['Annual_Cultures'] = 'Other'
p.mb_names_to_mbs_names['Semi_Perennial_Cultures_'] = 'Other'
p.mb_names_to_mbs_names['Crop_Mosaic'] = 'Other'
p.mb_names_to_mbs_names['Agriculture_or_Grassland'] = 'split_between_crops_and_pastures' # NOTE: split according to existing split, will be recategorized into 3, 4
p.mb_names_to_mbs_names['not_vegetated'] = 'Other'
p.mb_names_to_mbs_names['Beaches_and_dunes'] = 'Other'
p.mb_names_to_mbs_names['Urban_infrastructure'] = 'Urban'
p.mb_names_to_mbs_names['Other_non_vegetated_areas'] = 'Other'
p.mb_names_to_mbs_names['Water_bodies'] = 'Water'
p.mb_names_to_mbs_names['ndv'] = 'ndv' # NDV






















