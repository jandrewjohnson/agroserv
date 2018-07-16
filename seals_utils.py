import logging

import numpy as np
import scipy

import hazelbean as hb
from hazelbean.ui import model, inputs

logging.basicConfig(level=logging.WARNING)
hb.ui.model.LOGGER.setLevel(logging.WARNING)
hb.ui.inputs.LOGGER.setLevel(logging.WARNING)

L = hb.get_logger('seals', logging_level='warning')
L.setLevel(logging.WARNING)

logging.getLogger('Fiona').setLevel(logging.WARNING)
logging.getLogger('fiona.collection').setLevel(logging.WARNING)

np.seterr(divide='ignore', invalid='ignore')

dev_mode = True


def normalize_array(array, low=0, high=1, log_transform=True):
    if log_transform:
        min = np.min(array)
        max = np.max(array)
        to_add = float(min * -1.0 + 1.0)
        array = array + to_add

        array = np.log(array)

    min = np.min(array)
    max = np.max(array)

    normalizer = (high - low) / (max - min)

    output_array = (array - min) *  normalizer

    return output_array

def distance_from_blurred_threshold(input_array, sigma, threshold, decay):
    """
    Blur the input with a guassian using sigma (higher sigma means more blur). In areas of the blur above the threshold,
    return 1 - blurred so thtat values near zero indicate strong presence. In areas below the threshold, return


    The positive attribute of this func is it gives an s-curve relationship between 0 and 1 with a smoothed discontinuity
    around the threshold while never converging too close to 1 even at extreme distances without requiring slow calculation
    of a large footrint convolution.
    """

    blurred = scipy.ndimage.filters.gaussian_filter(input_array, sigma).astype(np.float32)

    blurred_below = np.where(blurred < threshold, 1, 0).astype(np.float32)  # NOTE opposite logic because EDT works only where true.

    if np.any(blurred_below == 0):
        distance = scipy.ndimage.morphology.distance_transform_edt(blurred_below).astype(np.float32)
    else:
        # Interesting choice here that I wasn't sure about how to address:
        # In the event that there are NO significant shapes, and thus blurred_below is all ones, what should be the distance?
        L.warning('WARNING NOTHING came up as above the blurring threshold!')
        metric = np.ones(blurred.shape)
        return metric

    outside = 1.0 - (1.0 / ((1 + float(decay)) ** (distance) + (1.0 / float(threshold) - 1)))  # 1 -  eponential distance decay from blur above threshold minus scalar that makes it match the level of
    inside = np.ones(blurred.shape).astype(np.float32) - blurred  # lol. just use 1 - the blurred value when above the threshold.

    metric = np.where(blurred_below == 1, outside, inside).astype(np.float32)

    metric = np.where(metric > 0.9999999, 1, metric)
    metric = np.where(metric < 0.0000001, 0, metric)
    metric = 1 - metric

    return metric