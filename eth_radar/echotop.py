"""
Determination of the echo top height from radar PPI data using the
Lakshmanan et al. (2013). method
@title: echotop
@author: Valentin Louf <valentin.louf@monash.edu>
@copyright: Valentin Louf (2018-2019)
@institution: Monash University
@reference: Lakshmanan et al. (2013), "An Improved Method for Estimating Radar
            Echo-Top Height". Weather Forecast. 28, 481–488,
            doi:10.1175/WAF-D-12-00084.1.

.. autosummary::
    :toctree: generated/
    cloud_top_height
    grid_cloud_top
"""
import numpy as np
from numba import jit


@jit
def cloud_top_height(
    r,
    azimuth,
    elevation,
    st_sweep,
    ed_sweep,
    refl,
    eth_thld=0,
    noise_thld=-2,
    min_range=15e3,
    verbose=False,
):
    """
    Estimating Radar Echo-Top Height using the improved method from Lakshmanan
    et al. (2013).

    Parameters:
    ===========
    r: <nr>
        Radar range.
    azimuth: <time>
        Radar azimuth.
    elevation: <time>
        Radar elevation.
    st_sweep: <nsweep>
        Radar sweep start ray index.
    ed_sweep: <nsweep>
        Radar sweep end ray index.
    refl: <time, nr>
        Radar reflectivity.
    eth_thld: float
        Threshold value (e.g., 0 dBZ, 18 dBZ, ...) used to compute the echo top
    noise_thld: float
        Signal to noise cutoff threshold value.
    min_range: float
        Minimum range in meter at which the echo top height are computed to
        avoid the cone of silence, generally 15 km.
    verbose: bool
        Print debug messages

    Returns:
    ========
    cloudtop: <na, nr>
        Cloud top height in meters, dimensions are na: length of the azimuth
        array of the first sweep, and nr: length of the input 'r' array.
    """
    earth_radius = 6371000

    na0 = st_sweep[1]
    nsweep = len(st_sweep)
    cloudtop = np.zeros((na0, len(r))) + np.NaN
    ground_range = np.zeros((nsweep, len(r)))
    elev_ref = elevation[0]

    for i, st in enumerate(st_sweep):
        ground_range[i, :] = r * np.cos(np.pi * elevation[st + 1] / 180)

    for i in range(1, len(st_sweep)):
        st = st_sweep[i]
        ed = ed_sweep[i]
        elev_ref = elevation[st_sweep[i - 1]]
        elev_iter = elevation[st]

        st_ref = st_sweep[i - 1]
        ed_ref = ed_sweep[i - 1]

        if verbose:
            print(i, st, ed, elev_iter, elev_ref)

        for j in range(na0):
            nazi_ref = np.argmin(np.abs(azimuth[st_ref:ed_ref] - azimuth[j])) + st_ref
            nazi_iter = np.argmin(np.abs(azimuth[st:ed] - azimuth[j])) + st

            if np.abs(azimuth[nazi_ref] - azimuth[nazi_iter]) > 5:
                continue

            for k in range(len(r)):
                if r[k] < min_range:
                    continue

                gr_ref = ground_range[i - 1, k]
                npos = np.argmin(np.abs(ground_range[i, :] - gr_ref))

                if np.abs(ground_range[i, npos] - ground_range[0, k] > 1000):
                    continue

                refb = refl[nazi_ref, k]
                refa = refl[nazi_iter, npos]

                if refb < noise_thld:
                    continue

                if refa < noise_thld or np.isnan(refa):
                    refa = noise_thld

                height = np.NaN
                if refb > eth_thld and refa < eth_thld:
                    if elev_iter == 90:
                        theta_total = elev_iter + 0.5
                    else:
                        theta_total = (eth_thld - refa) * (elev_ref - elev_iter) / (refb - refa) + elev_ref

                    # Include correction for Earth sphericity.
                    height = (r[k] * np.sin(np.pi * theta_total / 180) +
                              np.sqrt(r[k] ** 2 + earth_radius ** 2) - earth_radius)

                    if np.isnan(height):
                        continue

                    if height > cloudtop[j, k] or np.isnan(cloudtop[j, k]):
                        cloudtop[j, k] = height

    return cloudtop


@jit
def grid_cloud_top(data, xradar, yradar, xgrid, ygrid, theta_3db=1.5, rmax=150e3, gatespacing=250):
    """
    This function grid the cloud top height data (which are in polar
    coordinates) onto a Cartesian grid. This gridding technique is made to
    properly handle the absence of data (i.e. absence of clouds) while other
    gridding techniques tend to propagate NaN values.

    Parameters:
    ===========
    data: <ny, nx>
        Data to grid, ideally cloud top heights.
    xradar: <ny, nx>
        x-axis Cartesian coordinates array of the input data
    yradar: <ny, nx>
        y-axis Cartesian coordinates array of the input data
    xgrid: <ny_out, nx_out>
        x-axis Cartesian coordinates array for the output data
    ygrid: <ny_out, nx_out>
        y-axis Cartesian coordinates array for the output data
    theta_3db: float
        Maximum resolution angle in degrees for polar coordinates.
    rmax: float
        Maximum range of the data (same unit as x/y).
    gatespacing: float
        Gate-to-gate resolution (same unit as x/y).

    Returns:
    ========
    eth_out: <ny_out, nx_out>
        Gridded data.
    """
    if xradar.shape != data.shape:
        raise IndexError("Bad dimensions")

    if len(xgrid.shape) < len(xradar.shape):
        xgrid, ygrid = np.meshgrid(xgrid, ygrid)

    eth_out = np.zeros(xgrid.shape) + np.NaN

    for i in range(len(xgrid)):
        for j in range(len(ygrid)):
            cnt = 0
            zmax = 0
            xi = xgrid[j, i]
            yi = ygrid[j, i]

            if xi ** 2 + yi ** 2 > rmax ** 2:
                continue

            width = 0.5 * (np.sqrt(xi ** 2 + yi ** 2) * theta_3db * np.pi / 180)
            if width < gatespacing:
                width = gatespacing

            for k in range(data.shape[1]):
                for l in range(data.shape[0]):
                    xr = xradar[l, k]
                    yr = yradar[l, k]

                    if (
                        xr >= xi - width
                        and xr < xi + width
                        and yr >= yi - width
                        and yr < yi + width
                    ):
                        if data[l, k] > 0:
                            zmax = zmax + data[l, k]
                            cnt = cnt + 1

            if cnt != 0:
                eth_out[j, i] = zmax / cnt

    return eth_out
