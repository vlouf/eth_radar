# eth_radar

Method to estimate radar echo top height. This technique uses a volume scan to determine the maximum elevation angle at which a certain reflectivity threshold is exceeded. This is based upon the work of [Lakshmanan et al. (2013)][1].

## Description

Instead of simply finding the highest elevation angle within a volume, or a virtual volume where reflectivity exceeds some threshold, the modified algorithm is as follows:

- Find the maximum elevation angle (θb) where reflectivity (Zb) exceeds the echo-top reflectivity threshold. 
- If θb is not the highest elevation scan in the virtual volume, obtain the reflectivity value (Za) at the next higher elevation angle (θa). Then, the echo-top height is given by the height of the radar beam at an elevation angle:

where ZT is the threshold value (e.g., 0 dBZ, 18 dBZ) used to compute the echo top.

- If θb is the highest elevation scan available, set . This condition is met far away from the radar if higher-elevation scans have shorter ranges than a base “surveillance” scan and very close to the radar if the highest-elevation scan does not sample the top of the cloud. Under these circumstances, θT is set to be the top of the beam containing dBZ ≥ ZT; that is, the traditional echo-top algorithm is followed when there are no data available from a higher-elevation scan.

## Dependecies

This module requires [numpy][2] and [numba][3].

## Reference

[1.][1] Lakshmanan, V., Hondl, K., Potvin, C. K. & Preignitz, D. An Improved Method for Estimating Radar Echo-Top Height. Weather Forecast. 28, 481–488 (2013).

[1]: https://journals.ametsoc.org/doi/10.1175/WAF-D-12-00084.1
[2]: http://www.numpy.org/
[3]: http://numba.pydata.org/