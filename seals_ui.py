import logging
from collections import OrderedDict

import numpy as np

import hazelbean as hb
from hazelbean.ui import model, inputs
from hazelbean.ui import validation

logging.basicConfig(level=logging.WARNING)
hb.ui.model.LOGGER.setLevel(logging.WARNING)
hb.ui.inputs.LOGGER.setLevel(logging.WARNING)

L = hb.get_logger('seals')
L.setLevel(logging.INFO)

logging.getLogger('Fiona').setLevel(logging.WARNING)
logging.getLogger('fiona.collection').setLevel(logging.WARNING)

np.seterr(divide='ignore', invalid='ignore')

dev_mode = True
# TODOO NOTE This funtion must be here as it is automatically called by the invest model. Consider making this more flexible in the next release.
@validation.invest_validator
def validate(args, limit_to=None):
    validation_error_list = []
    return validation_error_list

class SealsUI(model.HazelbeanModel):
    def __init__(self, p):
        self.p = p
        model.HazelbeanModel.__init__(self,
                                   # label=u'seals',
                                   label=u'SEALS: Spatial Economic Allocation Land-change Simulator',
                                   target=p.execute,
                                   validator=validate,
                                   localdoc='../documentation')

        self.area_of_interest_path = inputs.File(args_key='area_of_interest_path',
            # default_path=p.default_paths['area_of_interest_path'],
            helptext="A shapefile with a single polygon that will be used to clip any other datasets.",
            label='Area of interest',
            validator=None)
        self.add_input(self.area_of_interest_path)

        self.base_year_lulc_path = inputs.File(args_key='base_year_lulc_path',
            helptext=("File path for the land-use, land-cover map onto which expansion will be added."),
            label='Land-use, land-cover (raster)',
            validator=None)
        self.add_input(self.base_year_lulc_path)

        # MODE SWITCH WOULD GO HERE

        self.coarse_change_maps_dir = inputs.Folder(args_key='coarse_change_maps_dir',
            helptext=("Dir to several rasters or a single raster that defines how many new hectares of each class will be in each grid-cell. By definition, this must be coarser resolution than the LULC map (otherwise you alread know all you need to know). Running the model will allocate the hectarage changes in this coarser map to the best locations on the higher-resolution LULC map."),
            label='Coarse change maps directory',
            validator=None)
        self.add_input(self.coarse_change_maps_dir)

        self.physical_suitability_path = inputs.File(args_key='physical_suitability_path',
            helptext=("File path to a raster that defines physical suitability (e.g. from soils or slope, NOT from adjacency) on a 0-1 scale (1 is most suitable). Where this raster is zero, no expansion can happen."),
            label='Physical suitability (raster)',
            validator=None)
        self.add_input(self.physical_suitability_path)

        self.physical_suitability_importance = inputs.Text(
            args_key='physical_suitability_importance',
            helptext=("File path to a raster that defines physical suitability (e.g. from soils or slope, NOT from adjacency) on a 0-1 scale (1 is most suitable). Where this raster is zero, no expansion can happen."),
            label='Physical suitability importance',
            validator=None)
        self.add_input(self.physical_suitability_importance)

        self.lulc_class_types_path = inputs.File(
            args_key='lulc_class_types_path',
            helptext=('CSV that defines which LULC classes below to which LULC class types. Has two columns, "lulc_id" for the original class IDs in the input LULC map, and "lulc_class_type", which is an index for each of the class types. This step typically simplifies an LULC map with many similar classes to one with fewer classes (that might allow for easier understanding of expansion relationships). For example forest and grass both might be simplified to \"natural\".'),
            label='LULC class types (CSV)',
            validator=None)
        self.add_input(self.lulc_class_types_path)

        self.class_proximity_parameters_path = inputs.File(
            args_key='class_proximity_parameters_path',
            helptext=("A csv with 3 columns, lulc_class_type, clustering, decay and 1 row for each class that has any effect on  agriculture. Clustering defines how much fragmentation affects the definition of class proximity ( lower value, like 1, means that pixels must be quite close to be considered part of a cluster while higher values, like 8, mean that pixels count as a cluster even if spread out considerably."),
            label='Class proximity parameters (CSV)',
            validator=None)
        self.add_input(self.class_proximity_parameters_path)

        self.pairwise_class_relationships_path = inputs.File(
            args_key='pairwise_class_relationships_path',
            helptext=("A csv that defines the relationships between all N classes in the simplified LULC. The first row and first column contain the class-ids from the simplified LULC while the interior values define how the class in the row_id attracts or repulses (1, -1 respectively) the col_id class;"),
            label='Pairwise class relationships (CSV)',
            validator=None)
        self.add_input(self.pairwise_class_relationships_path)

        self.conversion_eligibility_path = inputs.File(
            args_key='conversion_eligibility_path',
            helptext=("A csv that defines the relationships between all N classes in the simplified LULC. The first row and first column contain the class-ids from the simplified LULC while the interior values define how the class in the row_id attracts or repulses (1, -1 respectively) the col_id class;"),
            label='Conversion eligibility (CSV)',
            validator=None)
        self.add_input(self.conversion_eligibility_path)

        self.output_base_map_path = inputs.File(
            args_key='output_base_map_path',
            helptext=("If given, results will be placed on top of this map, inheriting the base map's final size and LULC where no data."),
            label='Output base map (optional)',
            validator=None)
        self.add_input(self.output_base_map_path)

        self.intermediate_dir = inputs.Folder('intermediate_dir', helptext='help', args_key='intermediate_dir')
        self.add_input(self.intermediate_dir)

        # # NOTE, containers dont need a seperate interactivity slot. has it  by default it seems
        # self.advanced_options_container = inputs.Container(
        #     args_key='advanced_options_container',
        #     expandable=True,
        #     expanded=False,
        #     interactive=True,
        #     label='Show advanced options')
        # self.add_input(self.advanced_options_container)

        self.enable_batch_mode = inputs.Checkbox('Enable batch mode', helptext='help', args_key='enable_batch_mode')
        self.enable_batch_mode.checkbox.setChecked(False)

        self.use_existing_batch = inputs.Checkbox('Use existing batch directory', helptext='help', args_key='use_existing_batch')
        self.use_existing_batch.checkbox.setChecked(True)

        self.skip_existing_batch_components = inputs.Checkbox('Skip existing batch components', helptext='help', args_key='skip_existing_batch_components')
        self.skip_existing_batch_components.checkbox.setChecked(True)

        self.batch_id = inputs.Text(
            args_key='batch_id',
            helptext=("Name of the column within the batch shapefile that identifies the name of each batch region."),
            label='Batch ID',
            validator=None)

        if dev_mode:
            self.add_input(self.enable_batch_mode)
            self.add_input(self.use_existing_batch)
            self.add_input(self.skip_existing_batch_components)
            self.add_input(self.batch_id)

    def generate_args_from_inputs(self):
        """Used to geenrate args automatically rather than manuually adding them.
         e.g., args[self.create_simplied_lulc.args_key] = self.create_simplied_lulc.value(),
         Note that this then means that the args_key must be exactly correct."""
        args = OrderedDict()
        input_types_to_read = [
            inputs.Text,
            inputs.Checkbox,
            inputs.Container,
            inputs.Dropdown,
            inputs.File,
            inputs.FileButton,
            inputs.FolderButton,
            inputs.Folder,
            inputs.InVESTModelInput,
            inputs.Multi,
            # inputs.FileSystemRunDialog,
            # inputs.FileDialog,
        ]

        for k,v in self.__dict__.items():
            if type(v) in input_types_to_read:
                args[v.args_key] = v.value()
        return args

    def assemble_args(self):
        return self.generate_args_from_inputs()
