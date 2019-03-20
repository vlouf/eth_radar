
"""
Determination of the echo top height from radar PPI data using the
Lakshmanan et al. (2013). method. This script format the input radar data file
into what is expected for the echotop.py script to run. It also saves the data.

@title: eth_radar
@author: Valentin Louf <valentin.louf@monash.edu>
@institution: Monash University
@date: 20/03/2019
@version: 1
.. autosummary::
    :toctree: generated/
    save_data
    main
"""
# Python Standard Library
import os
import argparse
import warnings
import traceback

# Other Libraries
import xarray as xr
import numpy as np
import crayons
import echotop


def save_data(radar, cth_grid):
    """
    TODO: save output data function.
    """

    return None


def main():
    """
    Read input radar file, extract the reflectivity, compute the echo top height
    (0 dB by default) and save the result.

    Arguments:
    ==========
    INFILE: global str
        Input file name
    ETH_THLD: global float
        Echo-top-height threshold (default 0 dB).
    REFL_NAME: global str
        Name of the reflectivity field used to compute the echo top height
        ('reflectivity' by default) in the radar data.
    """
    # Load data
    radar = xr.open_dataset(INFILE)
    r = radar.range.values
    azimuth = radar.azimuth.values
    elevation = radar.elevation.values
    try:
        refl = radar[REFL_NAME].values
    except KeyError:
        traceback.print_exc()
        print(crayons.red(f'Error: wrong reflectivity field name. You gave {REFL_NAME}. ' +
                          f'The possible field names are: {radar.keys()}'))

    st_sweep = radar.sweep_start_ray_index.values
    ed_sweep = radar.sweep_end_ray_index.values
    print(crayons.green(f'{os.path.basename(INFILE)} data loaded.'))

    # Compute ETH
    cth = echotop.compute_cloud_top(r, azimuth, elevation, st_sweep, ed_sweep, refl, eth_thld=ETH_THLD)
    print(crayons.green(f'{ETH_THLD}-dB echo top height computed on polar coordinates.'))

    # Grid data
    th = 450 - azimuth[slice(st_sweep[0], ed_sweep[0] + 1)]
    th[th < 0] += 360

    R, A = np.meshgrid(r, th)
    x = R * np.cos(np.pi * A / 180)
    y = R * np.sin(np.pi * A / 180)

    xgrid = np.arange(-145e3, 146e3, 2500)
    cth_grid = echotop.grid_radar(cth, x, y, xgrid, xgrid)
    print(crayons.green(f'Data gridded.'))

    save_data(radar, cth_grid)

    return None


if __name__ == "__main__":
    # Numba problem workaround
    os.environ["NUMBA_DISABLE_INTEL_SVML"] = "1"

    # Parse arguments
    parser_description = """Raw radar PPIs processing. It provides Quality
control, filtering, attenuation correction, dealiasing, unfolding, hydrometeors
calculation, and rainfall rate estimation."""
    parser = argparse.ArgumentParser(description=parser_description)
    parser.add_argument(
        '-i',
        '--input',
        dest='infile',
        type=str,
        help='Input file',
        required=True)
    parser.add_argument(
        '-o',
        '--output',
        dest='outdir',
        type=str,
        help='Output directory.',
        required=True)
    parser.add_argument(
        '-e',
        '--eth-thld',
        dest='eth_thld',
        type=float,
        help='Echo top height threshold (e.g. 0, 5, 17, ...). 0 dB by default.',
        default=0.0)
    parser.add_argument(
        '-f',
        '--dbz-name',
        dest='db_name',
        type=str,
        help='Radar reflectivity field name.',
        default='reflectivity')


    args = parser.parse_args()
    INFILE = args.infile
    OUTPATH = args.outdir
    ETH_THLD = args.eth_thld
    REFL_NAME = args.db_name

    main()