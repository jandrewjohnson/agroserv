import os
import shutil
import sys

sys.path.extend(['C:\OneDrive\Projects'])
import netCDF4
import numpy as np
from collections import OrderedDict
import logging

import hazelbean as hb

L = hb.get_logger('process_ssp_scenarios')

p = hb.ProjectFlow('../projects/input_processing')

L.setLevel(logging.DEBUG)


p.scenarios_data_dir = os.path.join(hb.EXTERNAL_BULK_DATA_DIR, 'lulc\\luh2')
p.scenario_names = [
    'IMAGE-ssp126',
    'AIM-ssp370',
    # 'MAGPIE-ssp585',
]
p.esa_lulc_base_data_path = os.path.join(hb.BASE_DATA_DIR, 'lulc', 'esacci',
                                         'ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif')
p.years = [2015, 2030]

def extract_luh_netcdf_to_geotiffs(input_nc_uri, output_dir, dim_selection_indices):
    """
    WARNING: currently only works given the nc structure of LUH2

    dim_selection_indices allows you to select all the data specific to an asigned set of indices.
    For example, with a standard time, lat, lon dimension structure of the .nc, assigning
    dim_selection_indices = [7, None, None] would grab time=7, all lats, all lons. If dim_selection_indices
    is a single value, it assumes that will be the selection for the first dimension..
    """

    if type(dim_selection_indices) in [float, int]:
        dim_selection_indices = [dim_selection_indices, None, None]

    print('input_nc_uri', input_nc_uri)
    ncfile = netCDF4.Dataset(input_nc_uri, 'r')
    # print('ncfile', ncfile)
    dim_names = ncfile.dimensions.keys()
    dims_data = OrderedDict()
    dim_lengths = []
    for dim_name in dim_names:
        if dim_name != 'bounds':
            dims_data[dim_name] = ncfile.variables[dim_name][:]
            dim_lengths.append(len(dims_data[dim_name]))

    # ASSUME THIS IS GLOBAL to derive the geotransform
    res = (360.0) / len(dims_data['lon'])
    geotransform = [-180.0, res, 0.0, 90.0, 0.0, -1 * res]
    projection = 'wgs84'

    for var_name, var in ncfile.variables.items():
        if 'standard_name' in var.__dict__:
            var_string = str(var_name) + ' ^ ' + var.__dict__['standard_name'] + ' ^ ' + var.__dict__['long_name']

            L.info('Processing ' + var_string)

            if var_name not in dim_names:
                # BROKEN, this wouldn't work with other dim_selection_indices schemes.
                array = var[:][dim_selection_indices[0]]

                no_data_value = -9999
                array = np.where(array > 1E19, no_data_value, array)

                output_geotiff_uri = os.path.join(output_dir, var_string + '.tif')
                hb.save_array_as_geotiff(array, output_geotiff_uri, ndv=no_data_value,
                                         geotransform_override=geotransform, projection_override=projection)

                # output_geotiff_rp_uri = hb.suri(output_geotiff_uri, 'eck')
                # wkt = hb.get_wkt_from_epsg_code(54012)
                # hb.reproject_dataset_uri(output_geotiff_uri, 0.25, wkt, 'bilinear', output_geotiff_rp_uri)



def extract_lulc(p):
    if p.run_this:
        for scenario_name in p.scenario_names:
            scenario_dir_name = 'RCP' + str(scenario_name.split('-')[1][4]) + '.' + str(scenario_name.split('-')[1][5]) + ' ' + str(scenario_name.split('-')[1][0:4])
            scenario_dir = os.path.join(p.cur_dir, scenario_dir_name)
            p.scenario_input_data_dir = os.path.join(p.scenarios_data_dir, scenario_dir_name)
            hb.create_dirs(scenario_dir)
            filename = 'multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-' + scenario_name + '-2-1-f_gn_2015-2100.nc'
            states_path = os.path.join(p.scenario_input_data_dir, filename)

            for year in p.years:
                year_dir = os.path.join(scenario_dir, str(year))
                hb.create_dirs(year_dir)
                L.info('Extracting from ' + states_path)
                extract_luh_netcdf_to_geotiffs(states_path, year_dir, year - 2015) #0 = 2015, last year is 85=2100
def convert_states_to_cropland_extent(p):
    if p.run_this:

        def add_crop_layers_from_dir(input_dir):

            crop_layer_names = [
                "c4per ^ area_fraction ^ C4 perennial crops.tif",
                "c4ann ^ area_fraction ^ C4 annual crops.tif",
                "c3per ^ area_fraction ^ C3 perennial crops.tif",
                "c3nfx ^ area_fraction ^ C3 nitrogen-fixing crops.tif",
                "c3ann ^ area_fraction ^ C3 annual crops.tif",
            ]
            uris_to_combine = [os.path.join(input_dir, i) for i in crop_layer_names]
            match_af = hb.ArrayFrame(uris_to_combine[0])
            proportion_cultivated = np.zeros(match_af.shape)
            mask = np.where((match_af.data >= 0.0) & (match_af.data <= 1.0))
            for uri in uris_to_combine:
                proportion_cultivated[mask] += hb.ArrayFrame(uri).data[mask]

            return proportion_cultivated

        scenario_dir_name = 'RCP' + str(p.scenario_names[0].split('-')[1][4]) + '.' + str(
            p.scenario_names[0].split('-')[1][5]) + ' ' + str(p.scenario_names[0].split('-')[1][0:4])
        match_path = os.path.join(p.extract_lulc_dir, scenario_dir_name, str(p.years[0]), "c4per ^ area_fraction ^ C4 perennial crops.tif")


        for scenario_name in p.scenario_names:
            scenario_dir = os.path.join(p.extract_lulc_dir, scenario_name)
            for year in p.years:
                input_dir = os.path.join(scenario_dir, str(year))

                scenario_dir_name = 'RCP' + str(scenario_name.split('-')[1][4]) + '.' + str(
                    scenario_name.split('-')[1][5]) + ' ' + str(scenario_name.split('-')[1][0:4])
                scenario_dir = os.path.join(p.extract_lulc_dir, scenario_dir_name, str(year))
                array = add_crop_layers_from_dir(scenario_dir)

                output_dir = os.path.join(p.cur_dir, scenario_dir_name, str(year))
                hb.create_dirs(output_dir)
                output_path = os.path.join(output_dir, 'proportion_cropland.tif')

                hb.save_array_as_geotiff(array, output_path, match_path)
        else:
            pass
def raster_stats_on_proportion_cropland(p):
    if p.run_this:

        scenario_dir_name = 'RCP' + str(p.scenario_names[0].split('-')[1][4]) + '.' + str(
            p.scenario_names[0].split('-')[1][5]) + ' ' + str(p.scenario_names[0].split('-')[1][0:4])
        match_path = os.path.join(p.extract_lulc_dir, scenario_dir_name, str(p.years[0]), "c4per ^ area_fraction ^ C4 perennial crops.tif")

        p.biomes_shapefile_path = os.path.join(p.input_dir, r"extent_shapefiles\Brazil_biomes.shp")
        for scenario_name in p.scenario_names:

            scenario_dir_name = 'RCP' + str(scenario_name.split('-')[1][4]) + '.' + str(
                scenario_name.split('-')[1][5]) + ' ' + str(scenario_name.split('-')[1][0:4])

            scenario_dir = os.path.join(p.convert_states_to_cropland_extent_dir, scenario_dir_name)
            for year in p.years:
                input_dir = os.path.join(scenario_dir, str(year))
                input_raster_path = os.path.join(input_dir, 'proportion_cropland.tif')

                base_raster_path_band = (input_raster_path, 1)
                aggregate_vector_path = 5
                aggregate_field_name = 5
                stats = hb.zonal_statistics(base_raster_path_band, p.biomes_shapefile_path,  'objectid')

                for k, v in stats.items():
                    print(input_raster_path,'#',k,'#',v['sum'] / v['count'] * 100)
def raster_stats_on_luh(p):
    if p.run_this:

        scenario_dir_name = 'RCP' + str(p.scenario_names[0].split('-')[1][4]) + '.' + str(
            p.scenario_names[0].split('-')[1][5]) + ' ' + str(p.scenario_names[0].split('-')[1][0:4])
        match_path = os.path.join(p.extract_lulc_dir, scenario_dir_name, str(p.years[0]), "c4per ^ area_fraction ^ C4 perennial crops.tif")

        p.biomes_shapefile_path = os.path.join(p.input_dir, r"extent_shapefiles\Brazil_biomes.shp")

        luh_filenames = [
            "pastr ^ area_fraction ^ managed pasture.tif",
            "primf ^ area_fraction ^ forested primary land.tif",
            "primn ^ area_fraction ^ non-forested primary land.tif",
            "range ^ area_fraction ^ rangeland.tif",
            "secdf ^ area_fraction ^ potentially forested secondary land.tif",
            "secdn ^ area_fraction ^ potentially non-forested secondary land.tif",
            "secma ^  ^ secondary mean age.tif",
            "secmb ^ vegetation_carbon_content ^ secondary mean biomass carbon density.tif",
            "urban ^ area_fraction ^ urban land.tif",
            "c3ann ^ area_fraction ^ C3 annual crops.tif",
            "c3nfx ^ area_fraction ^ C3 nitrogen-fixing crops.tif",
            "c3per ^ area_fraction ^ C3 perennial crops.tif",
            "c4ann ^ area_fraction ^ C4 annual crops.tif",
            "c4per ^ area_fraction ^ C4 perennial crops.tif",
        ]


        for scenario_name in p.scenario_names:
            scenario_dir_name = 'RCP' + str(scenario_name.split('-')[1][4]) + '.' + str(scenario_name.split('-')[1][5]) + ' ' + str(scenario_name.split('-')[1][0:4])
            scenario_input_data_dir = os.path.join(p.extract_lulc_dir, scenario_dir_name)


            for year in p.years:
                year_dir = os.path.join(scenario_input_data_dir, str(year))

                for luh_filename in luh_filenames:
                    luh_path = os.path.join(year_dir, luh_filename)
                    luh_path_tuple = (luh_path, 1)

                    stats = hb.zonal_statistics(luh_path_tuple, p.biomes_shapefile_path,  'objectid')

                    for k, v in stats.items():
                        print(luh_path, '#', k, '#', v['sum'] / v['count'] * 100)
def resample_mapbiomas(p):
    # WARNING, i corrupted amazonia by accident and haven't rerun (5 hr) yet
    p.mapbiomas_filenames = [
        "AMAZONIA.tif",
        # "CAATINGA.tif",
        # "CERRADO.tif",
        # "MATAATLANTICA.tif",
        # "PAMPA.tif",
        # "PANTANAL.tif",
    ]

    p.fine_pixel_size = hb.get_raster_info(p.esa_lulc_base_data_path)['pixel_size']
    p.mapbiomas_pixel_size = hb.get_raster_info(os.path.join(p.input_dir, 'fine_lulc', p.mapbiomas_filenames[0]))['pixel_size']

    if p.run_this:
        for filename in p.mapbiomas_filenames:
            path = os.path.join(p.input_dir, 'fine_lulc', filename)
            resampled_path = os.path.join(p.cur_dir, filename)

            # NOTE that all years are in separate bands (1=2000, 17=2016)
            options = list(hb.DEFAULT_GTIFF_CREATION_OPTIONS) + ['NUM_THREADS=ALL_CPUS']
            hb.warp_raster(path, p.fine_pixel_size, resampled_path, 'mode', gtiff_creation_options=options)
p.base_year_band_number = 16  # 2015
def extract_base_year_and_stitch_biomes(p):
    p.all_years_lulc_path = os.path.join(p.cur_dir, 'lulc_mapbiomas_300m.tif')
    p.base_year_lulc_path = os.path.join(p.cur_dir, 'lulc_2015_mapbiomas_300m.tif')
    if p.run_this:
        paths = []
        for filename in p.mapbiomas_filenames:
            path = os.path.join(p.resample_mapbiomas_dir, filename)
            paths.append(path)
        hb.create_gdal_virtual_raster_using_file(paths, p.all_years_lulc_path, remove_generator_files=False, srcnodata=0, dstnodata=255)
        p.base_year_lulc_array = hb.as_array(p.all_years_lulc_path, band_number=p.base_year_band_number)
        hb.save_array_as_geotiff(p.base_year_lulc_array, p.base_year_lulc_path, p.all_years_lulc_path, ndv=0)

def clip_luh_to_region(p):
    p.luh_input_dir = os.path.join(p.extract_lulc_dir, 'RCP7.0 ssp3')

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

    p.years_to_process = ['2015', '2030']

    p.biomes_shapefile = os.path.join(p.input_dir, r"extent_shapefiles\Brazil_biomes.shp")

    if p.run_this:
        for year in p.years_to_process:
            for name, filename in p.luh_names_to_filenames.items():
                path = os.path.join(p.luh_input_dir, year, filename)
                L.info(path)

                output_path = os.path.join(p.cur_dir, 'luh_' + str(year) + '_' + filename)
                hb.clip_raster_by_vector(path, output_path, p.biomes_shapefile, all_touched=True)
            for filename in ['proportion_cropland.tif']:
                path = os.path.join(p.convert_states_to_cropland_extent_dir, 'RCP7.0 ssp3', year, filename)
                output_path = os.path.join(p.cur_dir, 'luh_' + str(year) + '_' + filename)
                hb.clip_raster_by_vector(path, output_path, p.biomes_shapefile, all_touched=True)



def recategorize_luh(p):
    if p.run_this:
        for year in p.years_to_process:
            # NOTE cropland is aggregated above, just gonna use that but copy it to preserve folder structure.
            proportion_cropland_old_path = os.path.join(p.clip_luh_to_region_dir, 'luh_' + year + '_proportion_cropland.tif')
            proportion_cropland_current_path = os.path.join(p.cur_dir, 'luh_' + year + '_proportion_cropland.tif')
            shutil.copy(proportion_cropland_old_path, proportion_cropland_current_path)

            # Pasture is rangeland + pasture to align with mapbiomas
            luh_rangeland_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_' + year + '_range ^ area_fraction ^ rangeland.tif'))
            luh_pasture_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_' + year + '_pastr ^ area_fraction ^ managed pasture.tif'))
            proportion_pasture_path = os.path.join(p.cur_dir, 'luh_' + year + '_proportion_pasture.tif')

            # TODOO started good idea about having valid_masks saved as af properties but ran out of time.
            # proportion_pasture_current_af = hb.add_smart(luh_rangeland_af, luh_pasture_af, proportion_pasture_path, p.brazil_valid_mask_path, -9999.0)

            proportion_pasture_array = np.where((luh_pasture_af.valid_mask == 1) & (luh_rangeland_af.valid_mask == 1), luh_pasture_af.data + luh_rangeland_af.data, luh_pasture_af.ndv)
            hb.save_array_as_geotiff(proportion_pasture_array, proportion_pasture_path, luh_pasture_af.path)

            # Forest is primary forest plus half of secondary forst
            luh_primary_forest_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_' + year + '_primf ^ area_fraction ^ forested primary land.tif'))
            luh_secondary_forest_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_' + year + '_secdf ^ area_fraction ^ potentially forested secondary land.tif'))
            luh_primary_savanna_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_' + year + '_primn ^ area_fraction ^ non-forested primary land.tif'))
            luh_secondary_savanna_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_' + year + '_secdn ^ area_fraction ^ potentially non-forested secondary land.tif'))

            proportion_forest_path = os.path.join(p.cur_dir, 'luh_' + year + '_proportion_forest.tif')
            proportion_forest_array = np.where((luh_primary_forest_af.valid_mask == 1) & (luh_secondary_forest_af.valid_mask == 1), luh_primary_forest_af.data + 0.5 * luh_secondary_forest_af.data, luh_secondary_forest_af.ndv)
            hb.save_array_as_geotiff(proportion_forest_array, proportion_forest_path, luh_pasture_af.path)

            proportion_savanna_path = os.path.join(p.cur_dir, 'luh_' + year + '_proportion_savanna.tif')
            proportion_savanna_array = np.where((luh_primary_savanna_af.valid_mask == 1) & (luh_secondary_savanna_af.valid_mask == 1), luh_primary_savanna_af.data + 0.5 * luh_secondary_savanna_af.data, luh_secondary_savanna_af.ndv)
            hb.save_array_as_geotiff(proportion_savanna_array, proportion_savanna_path, luh_pasture_af.path)

            # NOTE urban is aggregated above, just gonna use that but copy it to preserve folder structure.
            proportion_urban_old_path = os.path.join(p.clip_luh_to_region_dir, 'luh_' + year + '_urban ^ area_fraction ^ urban land.tif')
            proportion_urban_current_path = os.path.join(p.cur_dir, 'luh_' + year + '_proportion_urban.tif')
            shutil.copy(proportion_urban_old_path, proportion_urban_current_path)




def reclassify_mb_to_mbs(p):
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

    p.mb_codes_to_mbs_codes = OrderedDict()
    for k, v in p.mb_names_to_mbs_names.items():
        p.mb_codes_to_mbs_codes[int(p.mb_names_to_mb_codes[k])] = p.mbs_names_to_codes[v]

    if p.run_this:
        input_tuple = (os.path.join(p.extract_base_year_and_stitch_biomes_dir, 'lulc_2015_mapbiomas_300m.tif'), 1)
        temp_path = hb.temp(remove_at_exit=True)
        hb.reclassify_raster(input_tuple, p.mb_codes_to_mbs_codes, temp_path, 1, 255)

        af = hb.ArrayFrame(temp_path)
        rand_array = np.random.randint(3, 5, af.shape) # 5 because randint exclusive on left end.
        a = np.where(af.data == 8, rand_array, af.data)
        output_path = os.path.join(p.cur_dir, 'lulc_2015_mbs.tif')
        hb.save_array_as_geotiff(a, output_path, temp_path, ndv=255)

def calculate_coarse_change_maps(p):
    """1) in net, the remaining 0.4 mha /yr should go to pasture. 2) as we discussed, there should also be a far amount of crop expansion into pasture and a fair amount of pasture expansion into savannah and forest. typically something like 65% of all deforestation goes to pasture."""

    # NOTE, this is the main output of the preprocessing step, creating 2 change maps prepended with the int of the class it will replace.
    p.cropland_change_mha_per_year = 1.1
    p.pasture_change_mha_per_year = 0.4


    p.global_ha_per_cell_15m_path = os.path.join(hb.BASE_DATA_DIR, "misc\ha_per_cell_15m_BUT_STILL_AT_5M_DATA.tif")
    p.ha_per_cell_15m_path = os.path.join(p.cur_dir, "ha_per_cell_15m.tif")
    hb.clip_raster_by_vector(p.global_ha_per_cell_15m_path, p.ha_per_cell_15m_path, p.biomes_shapefile, all_touched=True)

    if p.run_this:
        ha_per_cell_15m = hb.as_array(p.ha_per_cell_15m_path) * 9.0 #adjusting for 5 to 15 min size dif

        # PASTURE
        proportion_pasture_2015_af = hb.ArrayFrame(os.path.join(p.recategorize_luh_dir, 'luh_2015_proportion_pasture.tif'))
        proportion_pasture_2030_af = hb.ArrayFrame(os.path.join(p.recategorize_luh_dir, 'luh_2030_proportion_pasture.tif'))

        pasture_dif = (proportion_pasture_2030_af.data - proportion_pasture_2015_af.data) * ha_per_cell_15m
        pasture_dif_path = os.path.join(p.cur_dir, 'pasture_change.tif')
        hb.save_array_as_geotiff(pasture_dif, pasture_dif_path, proportion_pasture_2015_af.path)
        pasture_dif_af = hb.ArrayFrame(pasture_dif_path)

        pasture_increase = np.sum(pasture_dif[(pasture_dif > 0) & (pasture_dif != pasture_dif_af.ndv)])
        pasture_decrease = np.sum(pasture_dif[(pasture_dif < 0) & (pasture_dif != pasture_dif_af.ndv)])

        print('pasture_increase', pasture_increase, 'per year', pasture_increase/15)
        print('pasture_decrease', pasture_decrease, 'per year', pasture_decrease/15)

        # CROPLAND
        proportion_cropland_2015_af = hb.ArrayFrame(os.path.join(p.recategorize_luh_dir, 'luh_2015_proportion_cropland.tif'))
        proportion_cropland_2030_af = hb.ArrayFrame(os.path.join(p.recategorize_luh_dir, 'luh_2030_proportion_cropland.tif'))

        cropland_dif = (proportion_cropland_2030_af.data - proportion_cropland_2015_af.data) * ha_per_cell_15m
        cropland_dif_path = os.path.join(p.cur_dir, 'cropland_change.tif')
        hb.save_array_as_geotiff(cropland_dif, cropland_dif_path, proportion_cropland_2015_af.path)
        cropland_dif_af = hb.ArrayFrame(cropland_dif_path)

        cropland_increase = np.sum(cropland_dif[(cropland_dif > 0) & (cropland_dif != cropland_dif_af.ndv)])
        cropland_decrease = np.sum(cropland_dif[(cropland_dif < 0) & (cropland_dif != cropland_dif_af.ndv)])

        print('cropland_increase', cropland_increase, 'per year', cropland_increase/15)
        print('cropland_decrease', cropland_decrease, 'per year', cropland_decrease/15)


        # PHASE 2 NOTE: Our assumption that pasture stays about the same conflicts with luh. For phase 1 i scaled this away but bring this up to the group.
        ## Treatment that dealt with manageed pasture and rangeland separately. ignored for now, but may be useful later.
        # # MANAGED PASTURE
        # proportion_managed_pasture_2015_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_2015_pastr ^ area_fraction ^ managed pasture.tif'))
        # proportion_managed_pasture_2030_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_2030_pastr ^ area_fraction ^ managed pasture.tif'))
        #
        # managed_pasture_dif = (proportion_managed_pasture_2030_af.data - proportion_managed_pasture_2015_af.data) * ha_per_cell_15m
        # managed_pasture_dif_path = os.path.join(p.cur_dir, 'managed_pasture_change.tif')
        # hb.save_array_as_geotiff(managed_pasture_dif, managed_pasture_dif_path, proportion_managed_pasture_2015_af.path)
        # managed_pasture_dif_af = hb.ArrayFrame(managed_pasture_dif_path)
        #
        # managed_pasture_increase = np.sum(managed_pasture_dif[(managed_pasture_dif > 0) & (managed_pasture_dif != managed_pasture_dif_af.ndv)])
        # managed_pasture_decrease = np.sum(managed_pasture_dif[(managed_pasture_dif < 0) & (managed_pasture_dif != managed_pasture_dif_af.ndv)])
        #
        # print('managed_pasture_increase', managed_pasture_increase, 'per year', managed_pasture_increase/15)
        # print('managed_pasture_decrease', managed_pasture_decrease, 'per year', managed_pasture_decrease/15)
        #
        # # RANGELAND
        # proportion_rangeland_2015_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_2015_range ^ area_fraction ^ rangeland.tif'))
        # proportion_rangeland_2030_af = hb.ArrayFrame(os.path.join(p.clip_luh_to_region_dir, 'luh_2030_range ^ area_fraction ^ rangeland.tif'))
        #
        # rangeland_dif = (proportion_rangeland_2030_af.data - proportion_rangeland_2015_af.data) * ha_per_cell_15m
        # rangeland_dif_path = os.path.join(p.cur_dir, 'rangeland_change.tif')
        # hb.save_array_as_geotiff(rangeland_dif, rangeland_dif_path, proportion_rangeland_2015_af.path)
        # rangeland_dif_af = hb.ArrayFrame(rangeland_dif_path)
        #
        # rangeland_increase = np.sum(rangeland_dif[(rangeland_dif > 0) & (rangeland_dif != rangeland_dif_af.ndv)])
        # rangeland_decrease = np.sum(rangeland_dif[(rangeland_dif < 0) & (rangeland_dif != rangeland_dif_af.ndv)])
        #
        # print('rangeland_increase', rangeland_increase, 'per year', rangeland_increase/15)
        # print('rangeland_decrease', rangeland_decrease, 'per year', rangeland_decrease/15)

        cropland_scalar = (p.cropland_change_mha_per_year * 15 * 1000000) / cropland_increase # years and conversion to ha
        pasture_scalar = (p.pasture_change_mha_per_year * 15 * 1000000) / pasture_increase
        print('cropland_scalar', cropland_scalar)
        print('pasture_scalar', pasture_scalar)

        # NOW DO THE SCALING
        pasture_scaled_change = pasture_scalar * pasture_dif
        pasture_scaled_change[pasture_scaled_change<0] = 0
        pasture_scaled_change_path = os.path.join(p.cur_dir, '3_pasture_scaled_change.tif') # NOTE, prepended int corresponds with mbs class
        hb.save_array_as_geotiff(pasture_scaled_change, pasture_scaled_change_path, proportion_pasture_2015_af.path)
        pasture_scaled_change_af = hb.ArrayFrame(pasture_scaled_change_path)

        cropland_scaled_change = cropland_scalar * cropland_dif
        cropland_scaled_change[cropland_scaled_change<0] = 0
        cropland_scaled_change_path = os.path.join(p.cur_dir, '4_cropland_scaled_change.tif')
        hb.save_array_as_geotiff(cropland_scaled_change, cropland_scaled_change_path, proportion_cropland_2015_af.path)
        cropland_scaled_change_af = hb.ArrayFrame(cropland_scaled_change_path)


        # NOW do an mean preserving average resample to generate regional fragmentaiton inputs.
        pixel_scalar = 24
        coarser_pixel_size = cropland_scaled_change_af.geotransform[1] * pixel_scalar
        options = hb.DEFAULT_GTIFF_CREATION_OPTIONS + ['COMPRESS=LZW']

        coarser_pixel_pre_path = os.path.join(p.cur_dir, '3_pasture_scaled_change_scaled_' + str(pixel_scalar) + '_pre.tif')
        hb.warp_raster(pasture_scaled_change_af.path, (coarser_pixel_size, -coarser_pixel_size), coarser_pixel_pre_path, 'average', gtiff_creation_options=options)
        coarser_pixel_pre_af = hb.ArrayFrame(coarser_pixel_pre_path)
        new_array = coarser_pixel_pre_af.data
        coarser_pixel_path = os.path.join(p.cur_dir, '3_pasture_scaled_change_scaled_' + str(pixel_scalar) + '.tif')
        scaled_array = np.where(new_array != coarser_pixel_pre_af.ndv, new_array * pixel_scalar ** 2, coarser_pixel_pre_af.ndv)
        hb.save_array_as_geotiff(scaled_array, coarser_pixel_path, coarser_pixel_pre_path)
        coarser_pixel_af = hb.ArrayFrame(coarser_pixel_path)

        coarser_pixel_pre_path = os.path.join(p.cur_dir, '4_cropland_scaled_change_scaled_' + str(pixel_scalar) + '_pre.tif')
        hb.warp_raster(cropland_scaled_change_af.path, (coarser_pixel_size, -coarser_pixel_size), coarser_pixel_pre_path, 'average', gtiff_creation_options=options)
        coarser_pixel_pre_af = hb.ArrayFrame(coarser_pixel_pre_path)
        new_array = coarser_pixel_pre_af.data
        coarser_pixel_path = os.path.join(p.cur_dir, '4_cropland_scaled_change_scaled_' + str(pixel_scalar) + '.tif')
        scaled_array = np.where(new_array != coarser_pixel_pre_af.ndv, new_array * pixel_scalar ** 2, coarser_pixel_pre_af.ndv)
        hb.save_array_as_geotiff(scaled_array, coarser_pixel_path, coarser_pixel_pre_path)
        coarser_pixel_af = hb.ArrayFrame(coarser_pixel_path)





# START HERE, I think i'm ready to downscale in seals? Main thing to do is figure out how to deal with 2 types of change at once.
def extract_individual_region_shapefiles(p):

    if p.run_this:
        hb.convert_shapefile_to_multiple_shapefiles_by_id(p.biomes_shapefile, 'objectid', p.cur_dir)

def downscale_simple_luh(p):

    if p.run_this:

        output_path = os.path.join(p.cur_dir, 'luh_300m.tif')


if __name__ == '__main__':

    extract_lulc_task = p.add_task(extract_lulc)
    convert_states_to_cropland_extent_task = p.add_task(convert_states_to_cropland_extent)
    raster_stats_on_proportion_cropland_task = p.add_task(raster_stats_on_proportion_cropland)
    raster_stats_on_luh_task = p.add_task(raster_stats_on_luh)
    resample_mapbiomas_task = p.add_task(resample_mapbiomas)
    extract_base_year_and_stitch_biomes_task = p.add_task(extract_base_year_and_stitch_biomes)
    clip_luh_to_region_task = p.add_task(clip_luh_to_region)
    recategorize_luh_task = p.add_task(recategorize_luh)
    reclassify_mb_to_mbs_task = p.add_task(reclassify_mb_to_mbs)
    calculate_coarse_change_maps_task = p.add_task(calculate_coarse_change_maps)
    extract_individual_region_shapefiles_task = p.add_task(extract_individual_region_shapefiles)

    extract_lulc_task.skip_existing = 1
    convert_states_to_cropland_extent_task.skip_existing = 1
    raster_stats_on_proportion_cropland_task.run = 0
    raster_stats_on_luh_task.run = 0
    resample_mapbiomas_task.skip_existing = 1
    extract_base_year_and_stitch_biomes_task.skip_existing = 1
    clip_luh_to_region_task.skip_existing = 1
    recategorize_luh_task.skip_existing = 1
    reclassify_mb_to_mbs_task.skip_existing = 1
    calculate_coarse_change_maps_task.skip_existing = 0
    extract_individual_region_shapefiles_task.skip_existing = 1


    notes_from_avery = """
    
    # START HERE: Working for 1 full run. Need to create a region ID map so that I can modify the coarse change maps per the other scenarios (e.g., no deforestation in amazonia).
    
# Scenario methods
2015 is base year
deforestation in [baudef, nosfdef, onlysavadef*] *only run in cerrado
fragmentation in [fraglocal, fragregional]
forestcode in [fc, nofc]

Filenames (all are for 2030): 
cerrado_baudef_fragregional_nofc.tif
cerrado_nosfdef_fragregional_nofc.tif
cerrado_onlysavadef_fragregional_nofc.tif
    cerrado_baudef_fraglocal_nofc.tif
cerrado_nosfdef_fraglocal_nofc.tif
cerrado_onlysavadef_fraglocal_nofc.tif
cerrado_baudef_fragregional_fc.tif
cerrado_nosfdef_fragregional_fc.tif
cerrado_onlysavadef_fragregional_fc.tif
cerrado_baudef_fraglocal_fc.tif
cerrado_nosfdef_fraglocal_fc.tif
cerrado_onlysavadef_fraglocal_fc.tif
amazon_baudef_fragregional_nofc.tif
amazon_nosfdef_fragregional_nofc.tif
amazon_onlysavadef_fragregional_nofc.tif
amazon_baudef_fraglocal_nofc.tif
    amazon_nosfdef_fraglocal_nofc.tif
amazon_onlysavadef_fraglocal_nofc.tif
amazon_baudef_fragregional_fc.tif
amazon_nosfdef_fragregional_fc.tif
amazon_onlysavadef_fragregional_fc.tif
amazon_baudef_fraglocal_fc.tif
amazon_nosfdef_fraglocal_fc.tif
amazon_onlysavadef_fraglocal_fc.tif

# Prepping LUH projections
First remap LUH2015 to mapbiomas classes
    only primary forest is forest. + half of the other
    in cerrado, scale bau to 10000-15000 km2 1.5 mil ha per year of combined forest and savannah loss, do this by scaling SSP3
    rangeland and pasture in luh combine to equal mb data
    
Crops in bau have 1.1mha/yr expansion.
Pasture has 0.0 NET mha/yr expansion, but shifts as crop replaces it but it replaces natural land.
    
alter fractions of soy, annual agriculture, forest, and savannah from the LUH SSP 3 to harmonize with local knowledge of BAU rates of change
    First adjust the 2015 LUH quantities by biome according to the areas statistics we went over in the call.
    Then, you will use rates of growth for soy, annual agriculture determine agricultural areas.
    
    
    
# Fragmentation note:
    constraint that total change correctly sums up within an LUH SSP 3 grid-cell after harmonization of quantities withing grid cell of each land use and land cover class  as per our call


# Final product, stamp the downscaled data onto a trivially resampled 300m map of ssp3    

    """

    p.execute()
