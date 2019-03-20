
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

# Other Libraries
import crayons
import echotop


def save_data():
    """
    TODO: save output data function.
    """

    return None


def main():
    """
    TODO: Read data, format them into the appropriate args needed to call echotop.py.
    """

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

    args = parser.parse_args()
    INFILE = args.infile
    OUTPATH = args.outdir
    ETH_THLD = args.eth_thld

    main()