# Agroserv implementation of SEALS model

Because SEALS is still a relatively new software tool for NatCap, a specialized implementation, provided here, was required for this project. To install, just clone the repo and then run seals_model_agroserv.py, fix any import statements (there will be several, check out the requirements.txt for faster fixing) iteratively until the code runs. The input data is not in version control, so when this code gets running to the point of not finding the data, contact Justin.

## Results
The code above generates 10 geotiffs separately for Cerrado and Amazon. This generates 12 secenarios by taking permutations of the 10 maps according to the following:

Factors

1. Cerrado deforestation scenario
a. 1.5 mha per year w/ savannah-forest split based on SSP 3
b. 1.5 mha per year savannah only conversion
c. zero conversion of savannah or forest

2. Amazon deforestation scenario
a.  1.5 mha per year all forest
b. zero conversion of savannah or forest

3. Fragmentation
a. with constraint that total change correctly sums up within an LUH SSP 3 grid-cell after harmonization of quantities withing grid cell of each land use and land cover class  as per our call. You wrote that this will occur at 15 minutes, but that was a typo right? The resolution is actually 0.25 degrees for these.
b. with constraint that total change correctly sums up within biome

Because there are 3 cerrado deforestation scenarios but only two for the Amazon, we get 6 scenanarios by mixing and matching the two regions. Then, permuting the two fragmentation scenarios gives us 12.

This is why the files were provided separately for Amazon vs Cerrado so you can correctly mix-n-match them. This should be just a quick gdal_translate to merge the two extents, but let me know if this is a challenge and i could also provide the joined results in the next run.
