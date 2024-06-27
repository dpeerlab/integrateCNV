import scanpy as sc
import os

from . import prepare_regions
from . import score


def load_data():
    # Get current directory
    current_dir = os.path.dirname(os.path.realpath(__file__))
    # Load data
    ad = sc.read(os.path.join(current_dir, 'resources', 'ad.h5ad'))
    return ad
