import os
import shutil
import sys

sys.path.extend(['C:\OneDrive\Projects'])
from collections import OrderedDict
import logging
import hazelbean as hb

L = hb.get_logger('process_ssp_scenarios')

p = hb.ProjectFlow('../projects/output_processing')

p.scenario_and_region_paths = [
    r"C:\OneDrive\Projects\cge\agroserv\projects\amazon_foronly_fragregional",
    r"C:\OneDrive\Projects\cge\agroserv\projects\amazon_nosavfor_fragregional",
    r"C:\OneDrive\Projects\cge\agroserv\projects\cerrado_nosavfor_fragregional",
    r"C:\OneDrive\Projects\cge\agroserv\projects\cerrado_savonly_fragregional",
    r"C:\OneDrive\Projects\cge\agroserv\projects\amazon_nosavfor_fraglocal",
    r"C:\OneDrive\Projects\cge\agroserv\projects\cerrado_nosavfor_fraglocal",
    r"C:\OneDrive\Projects\cge\agroserv\projects\cerrado_savonly_fraglocal",
    r"C:\OneDrive\Projects\cge\agroserv\projects\amazon_foronly_fraglocal",
    r"C:\OneDrive\Projects\cge\agroserv\projects\cerrado_bau_fraglocal",
    r"C:\OneDrive\Projects\cge\agroserv\projects\cerrado_bau_fragregional",
]

def extract_and_rename_outputs(p):

    if p.run_this:
        output_tifs = []
        for path in p.scenario_and_region_paths:

            output_path = os.path.join(path, 'output')
            scenario_tifs = hb.get_most_recent_timestamped_file_in_dir(output_path, include_extensions='.tif')
            print('scenario_tifs', scenario_tifs)
            output_tifs.append(scenario_tifs)

        hb.pp(output_tifs)
        for path in output_tifs:
            dst = os.path.join(p.cur_dir, '_'.join(hb.explode_path(path)['file_root'].split('_')[0:-3]) + '.tif')
            print('path, dst', path, dst)
            shutil.copy(path, dst)

def reclass_esa_to_mbs(p):



    p.esa_names_to_esa_codes = OrderedDict()
    for k, v in hb.esacci_extended_short_class_descriptions.items():
        p.esa_names_to_esa_codes[v] = k

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

    p.esa_reclassed_path = os.path.join(p.cur_dir, 'esa_2015_sa_mbs.tif')

    p.esa_codes_to_mbs_codes = OrderedDict()
    for k, v in p.esa_names_to_mbs_names.items():
        p.esa_codes_to_mbs_codes[int(p.esa_names_to_esa_codes[k])] = p.mbs_names_to_codes[v]

    if p.run_this:
        input_tuple = (r"C:\OneDrive\Projects\cge\agroserv\projects\input_processing\input\ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif", 1)
        esa_clipped_path = os.path.join(p.cur_dir, 'esa_2015_sa.tif')
        aoi_path = r"C:\OneDrive\Projects\cge\agroserv\projects\sa_bb.shp"
        hb.clip_raster_by_vector(input_tuple[0], esa_clipped_path, aoi_path, output_data_type=1, nodata_target=255)


        input_tuple = (esa_clipped_path, 1)
        temp_path = hb.temp(remove_at_exit=True)

        hb.reclassify_raster(input_tuple, p.esa_codes_to_mbs_codes, p.esa_reclassed_path, 1, 255)

        # af = hb.ArrayFrame(temp_path)
        # rand_array = np.random.randint(3, 5, af.shape)  # 5 because randint exclusive on left end.
        # a = np.where(af.data == 8, rand_array, af.data)
        # output_path = os.path.join(p.cur_dir, 'lulc_2015_mbs.tif')
        # hb.save_array_as_geotiff(a, output_path, temp_path, ndv=255)

def stitch_files_by_scenario(p):
    p.scenario_names = [
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

    ]


    for scenario_name in p.scenario_names:
        amazon_string = '_'.join(scenario_name.split('_')[0:2] + [scenario_name.split('_')[-1]])
        cerrado_string = '_'.join(scenario_name.split('_')[2:4] + [scenario_name.split('_')[-1]])

        to_stitch = [
            p.esa_reclassed_path,
            hb.get_most_recent_timestamped_file_in_dir('C:/OneDrive/Projects/cge/agroserv/projects/' + amazon_string + '/output', amazon_string, '.tif'),
            hb.get_most_recent_timestamped_file_in_dir('C:/OneDrive/Projects/cge/agroserv/projects/' + cerrado_string + '/output', cerrado_string, '.tif'),
        ]
        output_tif_path = os.path.join(p.cur_dir, scenario_name + '.tif')

        print('scenario_name', scenario_name, to_stitch)
        hb.pp(to_stitch)

        hb.create_gdal_virtual_raster_using_file(to_stitch, output_tif_path, remove_generator_files=True, srcnodata=255, dstnodata=255)

L.setLevel(logging.DEBUG)

if __name__ == '__main__':
    extract_and_rename_outputs_task = p.add_task(extract_and_rename_outputs)
    extract_and_rename_outputs_task.run = 1
    extract_and_rename_outputs_task.skip_existing = 1

    reclass_esa_to_mbs_task = p.add_task(reclass_esa_to_mbs)
    reclass_esa_to_mbs_task.run = 1
    reclass_esa_to_mbs_task.skip_existing = 1

    stitch_files_by_scenario_task = p.add_task(stitch_files_by_scenario)
    stitch_files_by_scenario_task.run = 1
    stitch_files_by_scenario_task.skip_existing = 1

    p.execute()