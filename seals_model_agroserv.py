import collections
import logging
import os
import numpy as np
from osgeo import gdal

import hazelbean as hb
from hazelbean.ui import model, inputs

import seals_utils
import seals_ui

logging.basicConfig(level=logging.WARNING)
hb.ui.model.LOGGER.setLevel(logging.WARNING)
hb.ui.inputs.LOGGER.setLevel(logging.WARNING)

L = hb.get_logger('seals', logging_level='warning')
L.setLevel(logging.INFO)

logging.getLogger('Fiona').setLevel(logging.WARNING)
logging.getLogger('fiona.collection').setLevel(logging.WARNING)

np.seterr(divide='ignore', invalid='ignore')

dev_mode = True


def generate_batch_zones(p):
    if p.enable_batch_mode:
        if not p.use_existing_batch:
            hb.convert_shapefile_to_multiple_shapefiles_by_id(p.area_of_interest_path, p.batch_id, p.cur_dir)

        list_of_aoi_paths = hb.list_filtered_paths_nonrecursively(p.cur_dir, include_extensions='.shp')

        # This is the one part of the model code that needs to be aware of the structure of project_flow
        # When a batch is designed, it first does something (like creating the shapefiles above), and then it defines
        # p.iterator_replacements, which is a dict of name, value pairs that will modify the p object when the iterated task is run
        # according to it's position among the iterating tasks.
        p.iterator_replacements = collections.OrderedDict()

        # Simple replacement of the aoi to use
        p.iterator_replacements['area_of_interest_path'] = list_of_aoi_paths

        # Trickier replacement that will redefine the parent dir for each task so that it also WRITES in the correct output location
        p.iterator_replacements['cur_dir_parent_dir'] = [os.path.splitext(i)[0] for i in list_of_aoi_paths]
    else:
        print('Not running in batch mode.')

        # Because the seals model is enherently batchable, i chose to still define a single iteration into replacements:
        p.iterator_replacements = collections.OrderedDict()
        p.iterator_replacements['area_of_interest_path'] = [p.area_of_interest_path]
        p.iterator_replacements['cur_dir_parent_dir'] = [hb.file_root(p.area_of_interest_path)]


def process_coarse_change_maps(p):
    L.info('process_coarse_change_maps.')

    # Change maps are in this directory and must be of the format [CLASS_ID_INT]_[someting, but anything else].tif
    if not os.path.isdir(p.coarse_change_maps_dir):
        p.coarse_change_maps_dir = os.path.split(p.coarse_change_maps_dir)[0]
        if not os.path.isdir(p.coarse_change_maps_dir):
            raise NameError('Unable to parse coarse_change_maps_dir.')
    tifs_in_dir = hb.list_filtered_paths_nonrecursively(p.coarse_change_maps_dir, include_extensions='.tif')
    p.change_map_paths = []
    for path in tifs_in_dir:
        try:
            rendered_int = int(hb.file_root(path).split('_')[0])
        except:
            rendered_int = None
        if isinstance(rendered_int, int):
            p.change_map_paths.append(path)
    p.change_map_raster_infos = [hb.get_raster_info(i) for i in p.change_map_paths]

    # Test that all the change maps are the same properties.
    if len(set([i['geotransform'] for i in p.change_map_raster_infos])) != 1:
        raise NameError('The maps in coarse change maps dir are not all the same shape, projection, etc, or they have been improperly named/formatted.')

    # p.current_change_in_crop_extent_path = os.path.join(p.cur_dir, 'change_in_crop_extent.tif')
    p.current_change_map_paths = []
    for path in p.change_map_paths:
        new_path = os.path.join(p.cur_dir, os.path.split(path)[1])
        p.current_change_map_paths.append(new_path)
        if p.run_this:  # NOTE NONSTANDARD placement of run_this
            hb.clip_raster_by_vector(path, new_path, p.area_of_interest_path,
                                     resample_method='nearest',
                                     all_touched=True, verbose=True,
                                     ensure_fits=True, gtiff_creation_options=hb.DEFAULT_GTIFF_CREATION_OPTIONS)

    L.info('current_change_map_paths' + str(p.current_change_map_paths))


def create_lulc(p):
    L.info('Creating class-types lulc.')

    p.name_from_iterator_replacements = hb.file_root(p.area_of_interest_path)
    p.base_year_current_zone_lulc_path = os.path.join(p.cur_dir, 'base_year_' + p.name_from_iterator_replacements + '.tif')
    # Create match paths of both data types
    p.match_int_path = p.base_year_current_zone_lulc_path

    # p.lulc_simplified_path = os.path.join(p.cur_dir, 'lulc_simplified.tif')
    p.lulc_simplified_path = p.base_year_current_zone_lulc_path
    p.valid_mask_path = os.path.join(p.cur_dir, 'valid_mask_high_res.tif')
    if p.run_this:

        hb.clip_while_aligning_to_coarser(p.base_year_lulc_path, p.base_year_current_zone_lulc_path, p.area_of_interest_path,
                                          p.current_change_map_paths[0], resample_method='nearest',
                                          output_data_type=1, nodata_target=255,
                                          all_touched=True, verbose=True,
                                          ensure_fits=True, gtiff_creation_options=hb.DEFAULT_GTIFF_CREATION_OPTIONS)

        # Set NDV masking based on AOI of current zone.

        hb.create_valid_mask_from_vector_path(p.area_of_interest_path, p.base_year_current_zone_lulc_path, p.valid_mask_path)
        p.valid_mask = hb.as_array(p.valid_mask_path)
        hb.set_ndv_by_mask_path(p.base_year_current_zone_lulc_path, p.valid_mask_path)

        lulc_ds = gdal.Open(p.base_year_current_zone_lulc_path)
        lulc_band = lulc_ds.GetRasterBand(1)
        lulc_array = lulc_band.ReadAsArray().astype(np.int)

        # AGROSERVE CHANGE: In the event that there isn't a csv, assume that it is because the task has been done manually prior. Could make this programatic with conditionals.
        if os.path.exists(p.lulc_class_types_path):
            # load the simplified class correspondnce as a nested dictionary.
            lulc_class_types_odict = hb.file_to_python_object(p.lulc_class_types_path, declare_type='DD')

            # For cythonization reasons, I need to ensure this comes in as ints
            lulc_class_types_ints_dict = dict()

            for row_name in lulc_class_types_odict.keys():
                lulc_class_types_ints_dict[int(row_name)] = int(lulc_class_types_odict[row_name]['lulc_class_type'])

            # # 1 is agriculture, 2 is mixed ag/natural, 3 is natural, 4 is urban, 0 is no data
            lulc_simplified_array = hb.reclassify_int_array_by_dict_to_ints(lulc_array, lulc_class_types_ints_dict)
            no_data_value_override = hb.get_nodata_from_uri(p.base_year_current_zone_lulc_path)
            hb.save_array_as_geotiff(lulc_simplified_array, p.lulc_simplified_path, p.base_year_current_zone_lulc_path, data_type_override=1, set_inf_to_no_data_value=False,
                                     no_data_value_override=no_data_value_override, compress=True)
    else:
        p.lulc_simplified_path = p.base_year_current_zone_lulc_path


def create_physical_suitability(p):
    L.info('Creating physical suitability layer from base data.')
    # THIS OLD CODE documents how physical suitability was calculated, though now it's included as a base datum.
    # dem_unaligned_path = hb.temp('.tif', folder=p.workspace_dir, remove_at_exit=True) #hb.temp('.tif', remove_at_exit=True)
    # stats_to_calculate = ['TRI']
    # hb.clip_hydrosheds_dem_from_aoi(dem_unaligned_path, p.area_of_interest_path, p.match_float_path)
    # hb.calculate_topographic_stats_from_dem(dem_unaligned_path, p.physical_suitability_dir, stats_to_calculate=stats_to_calculate, output_suffix='unaligned')
    # dem_path = os.path.join(p.physical_suitability_dir, 'dem.tif')
    # hb.align_dataset_to_match(dem_unaligned_path, p.match_float_path, dem_path, aoi_uri=p.area_of_interest_path)
    # for stat in stats_to_calculate:
    #     stat_unaligned_path = os.path.join(p.physical_suitability_dir, stat + '_unaligned.tif')
    #     hb.delete_path_at_exit(stat_unaligned_path)
    #     stat_path = os.path.join(p.physical_suitability_dir, stat + '.tif')
    #     hb.align_dataset_to_match(stat_unaligned_path, p.match_float_path, stat_path, resample_method='bilinear',
    #                               align_to_match=True, aoi_uri=p.area_of_interest_path)
    # soc_path = os.path.join(p.physical_suitability_dir, 'soc.tif')
    # hb.align_dataset_to_match(p.base_data_soc_path, p.match_int_path, soc_path, aoi_uri=p.area_of_interest_path, output_data_type=7)
    # tri_path = os.path.join(p.physical_suitability_dir, 'tri.tif')
    # hb.align_dataset_to_match(p.base_data_tri_path, p.match_int_path, tri_path, aoi_uri=p.area_of_interest_path, output_data_type=7)
    # # TODOO Create cythonized array_sum_product()
    # p.physical_suitability_path = os.path.join(p.physical_suitability_dir, 'physical_suitability.tif')
    # soc_array = hb.as_array(soc_path)
    # tri_array = hb.as_array(tri_path)
    # physical_suitability_array = np.log(soc_array) - np.log(tri_array)

    p.global_physical_suitability_path = os.path.join(p.model_base_data_dir, 'physical_suitability_compressed.tif')
    p.physical_suitability_path = os.path.join(p.cur_dir, 'physical_suitability.tif')
    if p.run_this:
        # hb.clip_raster_by_vector(p.global_physical_suitability_path, p.physical_suitability_path, p.coarse_res_aoi_path, all_touched=True)
        hb.clip_while_aligning_to_coarser(p.global_physical_suitability_path, p.physical_suitability_path, p.area_of_interest_path,
                                          p.current_change_map_paths[0], resample_method='nearest',
                                          all_touched=True, verbose=True,
                                          ensure_fits=True, gtiff_creation_options=hb.DEFAULT_GTIFF_CREATION_OPTIONS)

        # hb.clip_dataset_uri(p.global_physical_suitability_path, p.coarse_res_aoi_path, p.physical_suitability_path, False, all_touched=False)
        physical_suitability_array = hb.as_array(p.physical_suitability_path)

        # Mask out nans etc
        # with np.warnings.catch_warnings():
        #     np.warnings.filterwarnings('ignore', r'RuntimeWarning: invalid value encountered')
        # np.warnings.filterwarnings('ignore', r'RuntimeWarning: invalid value encountered')

        np.seterr(divide='ignore', invalid='ignore')
        physical_suitability_array = np.where(physical_suitability_array > -1000, physical_suitability_array, 0)
        physical_suitability_array = np.where(physical_suitability_array < 100000000, physical_suitability_array, 0)

        hb.save_array_as_geotiff(physical_suitability_array, p.physical_suitability_path, p.match_int_path, data_type_override=6, compress=True)


def create_convolution_inputs(p):
    lulc_array = hb.as_array(p.lulc_simplified_path)

    ndv = hb.get_nodata_from_uri(p.lulc_simplified_path)

    # Get which values exist in simplified_lulc
    unique_values = list(hb.enumerate_array_as_odict(lulc_array).keys())
    unique_values = [int(i) for i in unique_values]
    ignore_values = [ndv]
    unique_values = [i for i in unique_values if i not in ignore_values]
    p.simplified_lulc_classes = unique_values

    p.classes_with_change = [int(os.path.split(i)[1].split('_')[0]) for i in p.current_change_map_paths]

    p.convolution_inputs_dir = p.cur_dir
    if p.run_this:
        L.info('Creating binaries for classes ' + str(p.simplified_lulc_classes))
        binary_paths = []
        for unique_value in p.simplified_lulc_classes:
            # binary_array = np.zeros(lulc_array.shape)
            binary_array = np.where(lulc_array == unique_value, 1, 0).astype(np.uint8)
            binary_path = os.path.join(p.convolution_inputs_dir, 'class_' + str(unique_value) + '_binary.tif')
            binary_paths.append(binary_path)
            hb.save_array_as_geotiff(binary_array, binary_path, p.lulc_simplified_path, compress=True)

        convolution_params = hb.file_to_python_object(p.class_proximity_parameters_path, declare_type='DD', output_key_data_type=str, output_value_data_type=float)
        convolution_paths = []
        for i, v in enumerate(p.simplified_lulc_classes):
            L.info('Calculating convolution for class ' + str(v))
            binary_array = hb.as_array(binary_paths[i])
            convolution_metric = seals_utils.distance_from_blurred_threshold(binary_array, convolution_params[str(v)]['clustering'], 0.5, convolution_params[str(v)]['decay'])
            convolution_path = os.path.join(p.convolution_inputs_dir, 'class_' + str(p.simplified_lulc_classes[i]) + '_convolution.tif')
            convolution_paths.append(convolution_path)
            hb.save_array_as_geotiff(convolution_metric, convolution_path, p.lulc_simplified_path, compress=True)

        pairwise_params = hb.file_to_python_object(p.pairwise_class_relationships_path, declare_type='DD', output_key_data_type=str, output_value_data_type=float)
        for i in p.simplified_lulc_classes:
            i_convolution_path = os.path.join(p.convolution_inputs_dir, 'class_' + str(i) + '_convolution.tif')
            i_convolution_array = hb.as_array(i_convolution_path)
            for j in p.classes_with_change:
                L.info('Processing effect of ' + str(i) + ' on ' + str(j))
                adjacency_effect_path = os.path.join(p.convolution_inputs_dir, 'adjacency_effect_of_' + str(i) + '_on_' + str(j) + '.tif')
                adjacency_effect_array = i_convolution_array * pairwise_params[str(i)][str(j)]
                hb.save_array_as_geotiff(adjacency_effect_array, adjacency_effect_path, p.lulc_simplified_path, compress=True)

        for i in p.classes_with_change:
            L.info('Combining adjacency effects for class ' + str(i))
            combined_adjacency_effect_array = np.ones(lulc_array.shape)
            combined_adjacency_effect_path = os.path.join(p.convolution_inputs_dir, 'combined_adjacency_effect_' + str(i) + '.tif')
            for j in p.simplified_lulc_classes:
                current_uri = os.path.join(p.convolution_inputs_dir, 'adjacency_effect_of_' + str(j) + '_on_' + str(i) + '.tif')  # NOTICE SWITCHED I and J
                current_effect = hb.as_array(current_uri)
                combined_adjacency_effect_array *= current_effect + 1.0  # Center on 1 so that 0.0 has no effect
            hb.save_array_as_geotiff(combined_adjacency_effect_array, combined_adjacency_effect_path, p.lulc_simplified_path, compress=True)


def create_conversion_eligibility(p):
    p.conversion_eligibility_dir = p.cur_dir

    if p.run_this:
        # Prevent illogical conversion eg new ag onto existing ag, or new ag onto urban
        conversion_eligibility_params = hb.file_to_python_object(p.conversion_eligibility_path, declare_type='DD', output_key_data_type=str, output_value_data_type=int)
        simplified_lulc_array = hb.as_array(p.lulc_simplified_path)
        for i in p.classes_with_change:

            conversion_eligibility_raster_path = os.path.join(p.conversion_eligibility_dir, str(i) + '_conversion_eligibility.tif')
            conversion_eligibility_array = np.zeros(simplified_lulc_array.shape).astype(np.float64)
            for j in p.simplified_lulc_classes:
                conversion_eligibility_array = np.where(simplified_lulc_array == j, conversion_eligibility_params[str(j)][str(i)], conversion_eligibility_array)
            hb.save_array_as_geotiff(conversion_eligibility_array, conversion_eligibility_raster_path, p.match_int_path, data_type_override=6, compress=True)


def create_overall_suitability(p):
    ## FOR REFERENCE ipbes_global_change_in_crop_extent_calculated_here but I should redo it for all the years/rcps/ssps
    # # Generate global cahnge map
    # proportion_ag_2000_path = r"C:\OneDrive\Projects\cge\seals\projects\gluc_mn\input\proportion_ag_2000.tif"
    # proportion_ag_2100_path = r"C:\OneDrive\Projects\cge\seals\projects\gluc_mn\input\proportion_ag_2100.tif"
    # proportion_ag_2100_aligned_path = r"C:\OneDrive\Projects\cge\seals\projects\gluc_mn\input\proportion_ag_2100_aligned.tif"
    # global_change_in_crop_extent_path = r"C:\OneDrive\Projects\cge\seals\projects\gluc_mn\input\global_change_in_crop_extent.tif"
    # hb.align_dataset_to_match(proportion_ag_2100_path, proportion_ag_2000_path, proportion_ag_2100_aligned_path, force_to_match=True, all_touched=True)
    # a = hb.as_array(proportion_ag_2100_aligned_path).astype(np.float64)
    # b = hb.as_array(proportion_ag_2000_path).astype(np.float64)
    # change_in_crop_extent_array = a - b
    # change_in_crop_extent_ha_array = change_in_crop_extent_array * 10000
    # hb.save_array_as_geotiff(change_in_crop_extent_ha_array, global_change_in_crop_extent_path, proportion_ag_2100_aligned_path, data_type_override=7)

    # TODOO NOTE went crop and pasture specific here.
    p.overall_suitability_dir = p.cur_dir
    hb.create_directories(p.overall_suitability_dir)
    p.overall_suitability_paths = []
    # p.overall_suitability_path = os.path.join(p.overall_suitability_dir, 'overall_suitability.tif')
    # p.weighted_overall_suitability_path = os.path.join(p.overall_suitability_dir, 'weighted_overall_suitability.tif')

    if p.run_this:
        # NOTE, here the methods assume ONLY crop will be changing insofar as the physical suitability is defined wrt crops; 0 is 1 becasue already  got rid of 0 in unique values
        physical_suitability_array = hb.as_array(p.physical_suitability_path)

        for i in p.classes_with_change:
            suitability_path = os.path.join(p.overall_suitability_dir, 'overall_suitability_' + str(i) + '.tif')
            p.overall_suitability_paths.append(suitability_path)
            combined_adjacency_effect_path = os.path.join(p.convolution_inputs_dir, 'combined_adjacency_effect_' + str(i) + '.tif')

            adjacency_effect_array = hb.as_array(combined_adjacency_effect_path)
            print('adjacency_effect_array mean', np.mean(adjacency_effect_array))
            adjacency_effect_array = seals_utils.normalize_array(adjacency_effect_array)  # Didn't put this in HB because didn't want to redo the 0.4.0 release.
            conversion_eligibility_raster_path = os.path.join(p.create_conversion_eligibility_dir, str(i) + '_conversion_eligibility.tif')
            conversion_eligibility_array = hb.as_array(conversion_eligibility_raster_path)
            print('conversion_eligibility_array mean', np.mean(conversion_eligibility_array))
            try:
                physical_suitability_importance = float(p.physical_suitability_importance)
            except:
                physical_suitability_importance = 0.5
                L.warning('Could not interpret physical suitability importance. Using default of 0.5')

            physical_suitability_array = seals_utils.normalize_array(physical_suitability_array)
            print('physical_suitability_array mean', np.mean(physical_suitability_array))

            overall_suitability_array = (adjacency_effect_array + (physical_suitability_importance * physical_suitability_array)) * conversion_eligibility_array
            print('overall_suitability_array mean', np.mean(overall_suitability_array))
            overall_suitability_array = np.where(np.isnan(overall_suitability_array), 0, overall_suitability_array)
            overall_suitability_array = np.where(overall_suitability_array < 0, 0, overall_suitability_array)
            print('overall_suitability_array mean', np.mean(overall_suitability_array))
            hb.save_array_as_geotiff(overall_suitability_array, suitability_path, p.match_int_path, data_type_override=6, compress=True)

        # # NOTE, i force a resample here, thus making the ipbes also 5min. consider making exclipcit.
        # temp1 = hb.temp('.tif', folder=p.workspace_dir, remove_at_exit=True)
        # hb.align_dataset_to_match(p.current_change_in_crop_extent_path, p.overall_suitability_path, temp1, resample_method='bilinear')
        # change_to_allocate_array_5min = hb.as_array(temp1)

        # weighted_overall_suitability_array = change_to_allocate_array_5min * overall_suitability_array
        # hb.save_array_as_geotiff(weighted_overall_suitability_array, p.weighted_overall_suitability_path, p.match_int_path, data_type_override=6, compress=True)


def create_allocation_from_change_map(p):
    # AGROSERVE shortcut note: assumed that it happens in SEQUENCE first cropland then pasture.
    if p.run_this:

        lulc_array = hb.as_array(p.lulc_simplified_path)
        new_lulc_array = np.copy(lulc_array)
        for change_map_path in p.current_change_map_paths:

            change_to_allocate_array = hb.as_array(change_map_path)

            # Often it is the case that the number of cells that will be allocated is greater than the amount of high-res cells actually available for conversion. This happens only if the
            # conversion_elligibility.csv rules out cells (it will not happen if only adjacency and physical suitability is done, as there will be SOME places allbethem terrible.
            num_cells_skipped = np.zeros(change_to_allocate_array.shape)

            class_to_allocate = int(os.path.split(change_map_path)[1].split('_')[0])
            # combined_adjacency_effect_path = os.path.join(p.convolution_inputs_dir, 'combined_adjacency_effect_' + str(class_to_allocate) + '.tif')

            current_overall_suitability_path = os.path.join(p.create_overall_suitability_dir, 'overall_suitability_' + str(class_to_allocate) + '.tif')

            overall_suitability_array = hb.as_array(current_overall_suitability_path)

            # Test that map resolutions are workable multiples of each other
            assert int(round(overall_suitability_array.shape[0] / change_to_allocate_array.shape[0])) == int(
                round(overall_suitability_array.shape[1] / change_to_allocate_array.shape[1]))
            aspect_ratio = int(round(overall_suitability_array.shape[0] / change_to_allocate_array.shape[0]))

            L.info('Beginning allocation using allocation ratio of ' + str(aspect_ratio))
            L.info('Sizes involved: overall_suitability_array, ' + str(overall_suitability_array.shape) + ' change_to_allocate_array, ' + str(change_to_allocate_array.shape))

            # TODOO SHORTCUT, used unprojected units
            ha_per_source_cell = 300 ** 2 / 100 ** 2
            change_array = np.zeros(lulc_array.shape)
            combined_rank_array = np.zeros(lulc_array.shape).astype(np.int)

            # TODOO Note that i ignored smaller-than-chunk shards.
            for change_map_region_row in range(change_to_allocate_array.shape[0]):
                L.info('Starting horizontal row ' + str(change_map_region_row))
                for change_map_region_col in range(change_to_allocate_array.shape[1]):
                    if not change_to_allocate_array[change_map_region_row, change_map_region_col] > 0:
                        num_cells_to_allocate = 0
                    else:
                        num_cells_to_allocate = int(round(change_to_allocate_array[change_map_region_row, change_map_region_col] / ha_per_source_cell))

                    if num_cells_to_allocate > 0:
                        source_map_starting_row = change_map_region_row * aspect_ratio
                        source_map_starting_col = change_map_region_col * aspect_ratio
                        combined_adjacency_effect_chunk = overall_suitability_array[source_map_starting_row: source_map_starting_row + aspect_ratio,
                                                          source_map_starting_col: source_map_starting_col + aspect_ratio]

                        ranked_chunk, sorted_keys = hb.get_rank_array_and_keys(combined_adjacency_effect_chunk, ndv=0)

                        if num_cells_to_allocate > len(sorted_keys[0]):
                            previous_num_cells_to_allocate = num_cells_to_allocate
                            num_skipped = num_cells_to_allocate - len(sorted_keys[0])
                            num_cells_to_allocate = len(sorted_keys[0])
                            L.warning(
                                'Allocation algorithm requested to allocate more cells than were available for transition given the suitability constraints. Num requested: ' + str(
                                    previous_num_cells_to_allocate) + ', Num allocated: ' + str(len(sorted_keys[0])) + ', Num skipped ' + str(num_skipped))
                            num_cells_skipped[change_map_region_row, change_map_region_col] = num_skipped

                        sorted_keys_array = np.array(sorted_keys)

                        # Create a tuple (ready for use as a numpy key) of the top allocation_amoutn keys
                        keys_to_change = (sorted_keys_array[0][0:num_cells_to_allocate], sorted_keys_array[1][0:num_cells_to_allocate])

                        change_chunk = np.zeros(ranked_chunk.shape)
                        change_chunk[keys_to_change] = 1

                        ## TODOO this was useful but there was a 29x29 vs 30x30 error. Renable after fix.
                        # Just for visualization purposes, who what all the ranked zones look like together when mosaiced.
                        combined_rank_array[source_map_starting_row: source_map_starting_row + aspect_ratio,
                        source_map_starting_col: source_map_starting_col + aspect_ratio] = ranked_chunk

                        # TODOO BUG, there's a slight shift to the right that comes in here.
                        change_array[source_map_starting_row: source_map_starting_row + aspect_ratio,
                        source_map_starting_col: source_map_starting_col + aspect_ratio] = change_chunk

            L.info('Processing outputted results.')
            # NOTE: Shortcut here that i just didn't allow it to overwrite if there wasn't any left. Better way would be precalibration or spreading elsewhere?
            # TODOO shortcut, should have the new classes be smartly chosen and or based on max present value in lulc.
            new_lulc_array = np.where((change_array == 1), class_to_allocate + 5, new_lulc_array)  # NOTE, pasture will be 8 thus, crops 9

            p.change_array_path = hb.ruri(os.path.join(p.cur_dir, str(class_to_allocate) + '_change_array.tif'))
            hb.save_array_as_geotiff(change_array, p.change_array_path, p.match_int_path, compress=True)

            p.num_cells_skipped_path = hb.ruri(os.path.join(p.cur_dir, str(class_to_allocate) + '_num_cells_skipped.tif'))
            hb.save_array_as_geotiff(num_cells_skipped, p.num_cells_skipped_path, change_map_path, compress=True)

            p.combined_rank_array_path = hb.ruri(os.path.join(p.cur_dir, str(class_to_allocate) + '_combined_rank_array.tif'))
            hb.save_array_as_geotiff(combined_rank_array, p.combined_rank_array_path, p.match_int_path, compress=True)

        p.projected_lulc_path = hb.ruri(os.path.join(p.cur_dir, 'projected_lulc.tif'))
        hb.save_array_as_geotiff(new_lulc_array, p.projected_lulc_path, p.match_int_path, no_data_value_override=255, data_type_override=1, compress=True)

        p.projected_lulc_masked_path = hb.ruri(os.path.join(p.cur_dir, 'projected_lulc_masked.tif'))
        hb.set_ndv_by_mask_path(p.projected_lulc_path, p.valid_mask_path, p.projected_lulc_masked_path)


def stitch_projections(p):
    if p.run_this:
        files_to_stitch = hb.list_filtered_paths_recursively(p.generate_batch_zones_dir, include_extensions='.tif', include_strings='projected_lulc_masked_20', depth=None)
        scenario_name = os.path.split(p.workspace_dir)[1]
        p.stitched_path = hb.ruri(os.path.join(p.output_dir, scenario_name + '.tif'))

        L.info('Stitching together all of the generated LULCs.')
        hb.create_gdal_virtual_raster_using_file(files_to_stitch, p.stitched_path, dstnodata=255)

        if p.output_base_map_path:
            L.info('Stamping generated lulcs with extent_shift_match_path of output_base_map_path ' + str(p.output_base_map_path))
            ndv = hb.get_datatype_from_uri(p.stitched_path)
            # extent_shift_match_path can be used to make the size of the output file exactly match something else.
            hb.create_gdal_virtual_raster_using_file(files_to_stitch, p.stitched_path, p.model_base_data_lulc_path, srcnodata=ndv, dstnodata=255)

            base_raster_path_band_list = [(p.stitched_path, 1),
                                          (p.model_base_data_lulc_path, 1)]

            def fill_where_missing(a, b):
                # NOTE == not work here because a.any() or a.all() error. Crappy workaround is inequalities.
                return np.where((a >= 255) & (b <= 255), b, a)

            # Because SEALS doesn't run for small islands, we fill in any missing values based on the base data input lulc.
            target_raster_path = r"C:\OneDrive\Projects\cge\seals\projects\unilever" + '/lulc_scenario1.tif'
            datatype_target = 1
            nodata_target = 255
            opts = ['TILED=YES', 'BIGTIFF=IF_SAFER', 'COMPRESS=lzw']
            hb.raster_calculator(base_raster_path_band_list, fill_where_missing, target_raster_path,
                                 datatype_target, nodata_target, gtiff_creation_options=opts)

            try:
                import geoecon as ge
                ge.add_geotiff_overview_file(target_raster_path)
            except:
                pass


main = ''
if __name__ == '__main__':
    from hazelbean.ui import model, inputs

    # TODOO Full solution would remove this hard-code and instead have it be from the user input. Note that this means the tree and it's task dirs are set only upon execute and nothing can be created in the tasknode creation algorithm besides the raw relationships..
    p = hb.ProjectFlow('../projects/russia_example')

    # pre_batch_task = p.add_task(validate_project_input_data)
    # pre_batch_task.creates_dir = False

    define_zones_iterator = p.add_iterator(generate_batch_zones)

    process_coarse_change_maps_task = p.add_task(process_coarse_change_maps, parent=define_zones_iterator)
    create_lulc_task = p.add_task(create_lulc, parent=define_zones_iterator)
    create_physical_suitability_task = p.add_task(create_physical_suitability, parent=define_zones_iterator)
    create_convolution_inputs_task = p.add_task(create_convolution_inputs, parent=define_zones_iterator)
    create_conversion_eligibility_task = p.add_task(create_conversion_eligibility, parent=define_zones_iterator)
    create_overall_suitability_task = p.add_task(create_overall_suitability, parent=define_zones_iterator)
    create_allocation_from_change_map_task = p.add_task(create_allocation_from_change_map, parent=define_zones_iterator)

    process_coarse_change_maps_task.run = 1
    create_lulc_task.run = 1
    create_physical_suitability_task.run = 1
    create_convolution_inputs_task.run = 1
    create_conversion_eligibility_task.run = 1
    create_overall_suitability_task.run = 1
    create_allocation_from_change_map_task.run = 1

    process_coarse_change_maps_task.skip_existing = 0
    create_lulc_task.skip_existing = 0
    create_physical_suitability_task.skip_existing = 0
    create_convolution_inputs_task.skip_existing = 0
    create_conversion_eligibility_task.skip_existing = 0
    create_overall_suitability_task.skip_existing = 0
    create_allocation_from_change_map_task.skip_existing = 0

    p.add_task(stitch_projections)

    ui = seals_ui.SealsUI(p)
    ui.run()
    EXITCODE = inputs.QT_APP.exec_()  # Enter the Qt application event loop. Without this line the UI will launch and then close.
