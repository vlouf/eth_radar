
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
import uuid
import datetime
import argparse
import warnings
import traceback

# Other Libraries
import xarray as xr
import numpy as np
import crayons
import echotop


def save_data(radar, x, y, cth_grid):
    """
    Save data to netCDF4 format using xarray.

    Parameters:
    ===========
    radar:
        The input radar dataset.
    x:
        x-dimension of output grid
    y:
        y-dimension of output grid
    cth_grid:
        Echo top height output grid.

    """
    dtime = radar.time.to_pandas()
    date = dtime[0].strftime("%Y%m%d.%H%M")
    outfilename = f"twp10cpolgrid.c1.eth{ETH_THLD}.{date}.nc"
    outfilename = os.path.join(OUTPATH, outfilename)

    metadata = radar.attrs.copy()

    metadata['creator_name'] = 'Valentin Louf'
    metadata['creator_email'] = 'valentin.louf@bom.gov.au'
    metadata['creator_url'] = 'github.com/vlouf'
    metadata['institution'] = 'Monash University'
    metadata['acknowledgement'] = 'This work has been supported by the U.S. Department of Energy Atmospheric Systems' +\
                                  ' Research Program through the grant DE-SC0014063. Data may be freely distributed.'
    metadata['created'] = datetime.datetime.utcnow().isoformat()
    metadata['uuid'] = str(uuid.uuid4())
    metadata['processing_level'] = 'c1'
    metadata['source'] = os.path.basename(INFILE)
    metadata['references'] = 'cf. 10.1175/WAF-D-12-00084.1'

    obsolete_keys = ['Conventions', 'version', 'history', 'field_names']
    for key in obsolete_keys:
        try:
            metadata.pop(key)
        except KeyError:
            pass

    dataset = xr.Dataset({
        'x': (('x'), x.astype(np.int32)),
        'y': (('y'), y.astype(np.int32)),
        'echo_top_height': (('y', 'x'), cth_grid),
    })

    dataset.x.attrs['standard_name'] = 'projection_x_coordinate'
    dataset.x.attrs['long_name'] = 'X distance on the projection plane from the origin'
    dataset.x.attrs['units'] = 'm'

    dataset.y.attrs['standard_name'] = 'projection_y_coordinate'
    dataset.y.attrs['long_name'] = 'Y distance on the projection plane from the origin'
    dataset.y.attrs['units'] = 'm'

    dataset.echo_top_height.attrs['long_name'] = f"{ETH_THLD}dB_echo_top_height"
    dataset.echo_top_height.attrs['description'] = '{} dBZ radar echo top height'.format(ETH_THLD)
    dataset.echo_top_height.attrs['units'] = 'm'
    dataset.echo_top_height.attrs['_FillValue'] = FILLVALUE
    dataset.echo_top_height.attrs['valid_min'] = np.int32(0)
    dataset.echo_top_height.attrs['valid_max'] = np.int32(25000)

    dataset.to_netcdf(outfilename, encoding=args)
    print(crayons.green(os.path.basename(outfilename) + " written."))

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
    cth_grid = echotop.KD_nn_interp(cth, x, y, xgrid, xgrid, nnearest = 24, maxdist = 2500) #nearest=24 should be enough to sample out to 2500m on a 1000m grid
    cth_grid = np.ma.masked_invalid(cth_grid).astype(np.int32).filled(FILLVALUE)
    print(crayons.green(f'Data gridded.'))

    save_data(radar, x=xgrid, y=xgrid, cth_grid=cth_grid)

    return None


if __name__ == "__main__":
    FILLVALUE = -9999

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

    if not os.path.isfile(INFILE):
        parser.error(f"Input file {INFILE} does not exists.")

    if not os.path.isdir(OUTPATH):
        parser.error(f"Output directory {OUTPATH} does not exists.")

    main()