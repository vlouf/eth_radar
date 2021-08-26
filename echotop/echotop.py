"""
Determination of the echo top height from radar PPI data using the
Lakshmanan et al. (2013). method
@title: echotop
@author: Valentin Louf <valentin.louf@monash.edu>
@date: 26/08/2021
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
from scipy.spatial import cKDTree


@jit
def cloud_top_height(
    r, azimuth, elevation, st_sweep, ed_sweep, refl, eth_thld=0, noise_thld=-2, min_range=15e3, verbose=False,
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
            nbeam_ref = np.argmin(np.abs(azimuth[st_ref:ed_ref] - azimuth[j])) + st_ref
            nbeam_iter = np.argmin(np.abs(azimuth[st:ed] - azimuth[j])) + st

            if np.abs(azimuth[nbeam_ref] - azimuth[nbeam_iter]) > 5:
                continue

            for k in range(len(r)):
                if r[k] < min_range:
                    continue

                gr_ref = ground_range[i - 1, k]
                npos = np.argmin(np.abs(ground_range[i, :] - gr_ref))

                if np.abs(ground_range[i, npos] - ground_range[0, k] > 1000):
                    continue

                refb = refl[nbeam_ref, k]
                refa = refl[nbeam_iter, npos]

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
                    height = (
                        r[k] * np.sin(np.pi * theta_total / 180) + np.sqrt(r[k] ** 2 + earth_radius ** 2) - earth_radius
                    )

                    if np.isnan(height):
                        continue

                    if height > cloudtop[j, k] or np.isnan(cloudtop[j, k]):
                        cloudtop[j, k] = height

    return cloudtop


@jit
def column_max_reflectivity(r, azimuth, elevation, st_sweep, ed_sweep, refl):
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

    Returns:
    ========
    cloudtop: <na, nr>
        Cloud top height in meters, dimensions are na: length of the azimuth
        array of the first sweep, and nr: length of the input 'r' array.
    """    
    na0 = st_sweep[1]
    nsweep = len(st_sweep)
    cloudtop = np.zeros((na0, len(r))) + np.NaN
    ground_range = np.zeros((nsweep, len(r)))    

    for i, st in enumerate(st_sweep):
        ground_range[i, :] = r * np.cos(np.pi * elevation[st + 1] / 180)

    for i in range(1, len(st_sweep)):
        st = st_sweep[i]
        ed = ed_sweep[i]        

        st_ref = st_sweep[i - 1]
        ed_ref = ed_sweep[i - 1]

        for j in range(na0):
            nbeam_ref = np.argmin(np.abs(azimuth[st_ref:ed_ref] - azimuth[j])) + st_ref
            nbeam_iter = np.argmin(np.abs(azimuth[st:ed] - azimuth[j])) + st

            if np.abs(azimuth[nbeam_ref] - azimuth[nbeam_iter]) > 5:
                continue

            for k in range(len(r)):

                gr_ref = ground_range[i - 1, k]
                npos = np.argmin(np.abs(ground_range[i, :] - gr_ref))

                if np.abs(ground_range[i, npos] - ground_range[0, k] > 1000):
                    continue

                refb = refl[nbeam_ref, k]
                refa = refl[nbeam_iter, npos]                
                
                cloudtop[j, k] = np.nanmax([refa, refb, cloudtop[j, k]])

    return cloudtop


def grid_cloud_top(data_in, x_in, y_in, x_out, y_out, nnearest=1, maxdist=None):
    """
    Nearest neighbour interpolation using scipy KDTree.

    Parameters:
    ===========
    data_in: ndarray of float with shape (n1, n2)
        Data values to interpolate in input coordinate space
    x_in: ndarray of float with shape (n1, n2)
        x values of input coordinate space (e.g., require conversion from polar to Catesian first)
    y_in: ndarray of float with shape (n1, n2)
        y values of input coordinate space
    x_out: ndarray of float with shape (n1a, n2a)
        x values of output coordinate space
    y_out: ndarray of float with shape (n1a, n2a)
        x values of output coordinate space
    nnearest: int
        maximum number of nearest neighbours to consider when filling NaN values
    maxdist: float (in units of Cartesian space)
        maximum distance of nearest neighbours to consider when filling NaN values

    Returns:
    ========
    vals_out: ndarray of float with shape (n1a, n2a)
    """

    def _make_coord_arrays(x):
        """
        Make sure that the coordinates are provided as ndarray
        of shape (numpoints, ndim)
        Parameters
        ----------
        x : ndarray of float with shape (numpoints, ndim)
            OR a sequence of ndarrays of float with len(sequence)==ndim and
            the length of the ndarray corresponding to the number of points
        """
        if type(x) in [list, tuple]:
            x = [item.ravel() for item in x]
            x = np.array(x).transpose()
        elif type(x) == np.ndarray:
            if x.ndim == 1:
                x = x.reshape(-1, 1)
            elif x.ndim == 2:
                pass
            else:
                raise Exception("Cannot deal wih 3-d arrays, yet.")
        return x

    # transform output coordinates into pairs of coordinates
    coord_out = _make_coord_arrays([x_out.ravel(), y_out.ravel()])
    vals_in = data_in.ravel()

    # build KDTree
    tree = cKDTree(np.c_[x_in.ravel(), y_in.ravel()])

    # query tree using output coordinates
    dists, idx = tree.query(coord_out, k=nnearest + 1)
    # avoid bug, if there is only one neighbor at all
    if dists.ndim == 1:
        dists = dists[:, np.newaxis]
        idx = idx[:, np.newaxis]
    # get first neighbour

    vals_out = vals_in[idx[:, 0]]
    dists_cp = dists[..., 0].copy()

    # iteratively fill NaN with next neighbours
    isnan = np.isnan(vals_out)
    nanidx = np.argwhere(isnan)[..., 0]
    if nnearest > 1 & np.count_nonzero(isnan):
        for i in range(nnearest - 1):
            vals_out[isnan] = vals_in[idx[:, i + 1]][isnan]
            dists_cp[nanidx] = dists[..., i + 1][nanidx]
            isnan = np.isnan(vals_out)
            nanidx = np.argwhere(isnan)[..., 0]
            if not np.count_nonzero(isnan):
                break

    # apply max distance
    if maxdist is not None:
        vals_out = np.where(dists_cp > maxdist, np.nan, vals_out)

    return np.reshape(vals_out, x_out.shape)
