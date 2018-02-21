"""functions to re-generate Geisler 1984 figures
"""
import numpy as np
from scipy.signal import convolve2d


def gauss_1d(x, amplitude, sigma):
    """generate 1d Gaussian centered at 0
    """
    return amplitude*np.exp(-0.5*np.power(x/sigma, 2))/(2*sigma)


def linespread_function(x=np.linspace(-4.25, 4.25, 100), a1=.684, s1=.443, a2=.587, s2=2.035):
    """generate (normalized) linespread function for specified parameters (un-numbered eqt, pg 2)

    defaults are the values used in the paper (all in arc-minutes). normalized as described in the
    legend for figure 1
    """
    lsf = gauss_1d(x, a1, s1) + gauss_1d(x, a2, s2)
    # np.trapz integrates using the "composite trapezoidal rule"
    return lsf / np.trapz(lsf, x)


def gauss_2d(coord, amplitude, sigma):
    """generate 2d Gaussian centered at 0
    """
    dist = np.sqrt(np.power(coord[0], 2) + np.power(coord[1], 2))
    return amplitude*np.exp(-0.5*np.power(dist/sigma, 2))/(2*sigma)


def pointspread_function(coord=np.meshgrid(np.linspace(-4.25, 4.25, 100),
                                           np.linspace(-4.25, 4.25, 100)),
                         a1=.684, s1=.443, a2=.587, s2=2.035):
    """generate (normalized) pointspread function for specified parameters

    defaults are the values used in the paper (all in arc-minutes)

    coord: 2-tuple/list of x and y, as returned by np.meshgrid

    """

    g_sum = gauss_2d(coord, a1, s1) + gauss_2d(coord, a2, s2)
    return g_sum


def absorb(lum, psf, a=0.28, d=0.2, s=3.1416, t=0.68, e555=0.5):
    """calculate the mean number of quanta absorbed (eqt 2)

    lum: luminance distr (parameterized by x, y) in cd/m2

    psf: point-spread function (param. by x, y) [a.u.]

    a: cross-sectional area of the receptor in sq minutes

    d: duration of stimulus in seconds

    s: pupil area in sq mm

    t: transmittance of ocular media

    e555: quantum efficiency of photoreceptor at lambda=555 nm
    """
    # 2D convolution. mode/boundary/fillvalue are defaults, but we include them for clarity why go
    # with this? photorecepter, modeled by psf, is discrete in space so, lum go beyond edge
    # (mode=full), but nothing will be absorbed (fillvalue=0)
    l_p_conv = convolve2d(lum, psf, mode='full', boundary='fill', fillvalue=0)
    return a*d*s*t*e555 * 347.8 * np.sum(np.sum(l_p_conv)), l_p_conv
