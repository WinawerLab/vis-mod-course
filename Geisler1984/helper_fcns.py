"""functions to re-generate Geisler 1984 figures
"""
import warnings
import itertools
import numpy as np
import pandas as pd
from scipy.signal import convolve2d
from scipy import optimize
import matplotlib.pyplot as plt


def get_middle(x, dim=0):
    mid = x.shape[dim] / 2
    if mid == np.floor(mid):
        warnings.warn("x has no middle index, returning the floor instead")
    return int(np.floor(mid))


def gauss_1d(x, amplitude, sigma):
    """generate 1d Gaussian centered at 0
    """
    return amplitude*np.exp(-0.5*np.power(x/sigma, 2))/(2*sigma)


def linespread_function(x=np.linspace(-4.25, 4.25, 101), a1=.684, s1=.443, a2=.587, s2=2.035):
    """generate (normalized) linespread function for specified parameters (un-numbered eqt, pg 2)

    defaults are the values used in the paper (all in arc-minutes). normalized as described in the
    legend for figure 1
    """
    lsf = gauss_1d(x, a1, s1) + gauss_1d(x, a2, s2)
    # np.trapz integrates using the "composite trapezoidal rule"
    return lsf / np.trapz(lsf, x)


def gauss_2d(x, y, amplitude, sigma):
    """generate 2d Gaussian centered at 0

    will convert x and y into coord, by calling np.meshgrid(x, y)
    """
    coord = np.meshgrid(x, y)
    dist = np.sqrt(np.power(coord[0], 2) + np.power(coord[1], 2))
    return amplitude*np.exp(-0.5*np.power(dist/sigma, 2))/(2*sigma)


def pointspread_function(x=np.linspace(-4.25, 4.25, 101), y=np.linspace(-4.25, 4.25, 101), a1=.684,
                         s1=.443, a2=.587, s2=2.035):
    """generate (normalized) pointspread function for specified parameters

    defaults are the values used in the paper (all in arc-minutes)
    """
    psf = gauss_2d(x, y, a1, s1) + gauss_2d(x, y, a2, s2)
    # np.trapz integrates using the "composite trapezoidal rule"
    return psf / np.trapz(np.trapz(psf, x), y)


def construct_photoreceptor_lattice(n_receptors_per_side):
    """construct hexagonal photoreceptor lattice with specified number of receptors

    in the paper, the receptors are modeled as points in a hexagonal lattice, with each point
    representing the center of the receptor and the distance between them (.6 arc-minutes)
    corresponding to the diameter of the receptors. based on the number of receptors, we can then
    determine the size of receptor lattice in visual degrees and we sample it finely enough so that
    each pixel represents .02 arc-minutes, which allows us to approximate this hexagonal lattice
    reasonably well (since we're trying to create a hexagonal lattice on a square grid, there's no
    way to do it exactly).

    we construct a basically square lattice. for example, if n_receptors_per_side=10, our lattice
    will have 10 rows of 10 photoreceptors each. based on the fact that sizes of the receptor
    lattice given in the paper were all square numbers (400, 10000), it seems reasonable to assume
    Geisler did something similar.

    Since each receptor center is .6 minutes apart, if we have row with 5 receptors, we'll need a
    total of (5-1) * .6 = 2.4 arc-minutes for all receptors (the first and last receptors will then
    be in the first and last pixel). Since the next row will be offset by .3 minutes (see fig 2 for
    an example image), we'll need 2.7 arc-minutes in the x direction. finally, since each row is
    separated by exactly .3*sqrt(3)~.52 minutes (there's no offset; this follows from the
    trigonometry), we'll need (5-1) * .52 = 2.08 arc-minutes in the y direction. In order to make
    sure we end up with an origin (0, 0), our lattice will run from a negative value to a positive
    value (with absolute value equal to half the max), and this half-max value must be divisible by
    .02. therefore our lattice is will be based on:

    np.meshgrid(np.arange(-1.36, 1.36+.02, .02), np.arange(-1.04, 1.04+.02, .02))

    (the extra `+.02` is there because of how python's range / numpy.arange works: they go up to
    but don't include the max value, so we set it one step beyond the max value we actually want)

    returns lattice, x_minutes, y_minutes, x, y. x/y_minutes are the number of minutes on either
    side of 0 in the created lattice (so, from example above, 1.36 and 1.04, respectively). x and y
    are the arrays that go from -x/y_minutes to +x/y_minutes and are necessary when calling the
    other functions in this file
    """
    receptor_diameter = .6
    # this is approximately (.6/2)*tan(pi/3)
    receptor_height_offset = .52
    min_per_pix = .02
    x_minutes = ((n_receptors_per_side - 1) * receptor_diameter + (receptor_diameter / 2.)) / 2
    y_minutes = ((n_receptors_per_side - 1) * receptor_height_offset) / 2
    if int(x_minutes / min_per_pix) != (x_minutes / min_per_pix):
        x_minutes += min_per_pix / 2
    if int(y_minutes / min_per_pix) != (y_minutes / min_per_pix):
        y_minutes += min_per_pix / 2
    x = np.arange(-x_minutes, x_minutes + min_per_pix, min_per_pix)
    y = np.arange(-y_minutes, y_minutes + min_per_pix, min_per_pix)
    lattice = np.zeros((len(y), len(x)))
    # by construction, we know that all these divisions will return integers
    for row_num, y_idx in enumerate(range(len(y))[::int(receptor_height_offset/min_per_pix)]):
        if row_num % 2:
            x_indices = range(len(x))[::int(receptor_diameter/min_per_pix)]
        else:
            # we subtract 1 here because python is 0 indexed
            x_indices = range(len(x))[int((receptor_diameter/2)/min_per_pix)-1::int(receptor_diameter/min_per_pix)]
        for x_idx in x_indices:
            lattice[y_idx, x_idx] = 1
    return lattice, x_minutes, y_minutes, x, y


def minutes_to_n_receptors(arcmin_diam):
    """given the size of an image in minutes, return the number of receptors necessary to have a
    receptor lattice that is larger
    """
    # this is approximately (.6/2)*tan(pi/3)
    receptor_height_offset = .52
    return np.ceil(arcmin_diam / receptor_height_offset) + 1


def retinal_image_convolution(lum, psf):
    """calculate the retinal image by convolving the luminance distribution with the pointspread

    note that you shouldn't interpret the numbers too much in the output from this function,
    because we should probably do some scaling to account for the size of the pupil and the
    transmittance of the ocular media. we handle that in the mean_photons_absorbed function
    instead.
    """
    # 2D convolution. what we're creating is the image on the (specified patch of) retina and we're
    # not modeling anything beyond that. so we'll throw everything away. Ideally, we just want to
    # make sure that the image is big enough that we're not throwing anything away.
    retinal_image = convolve2d(lum, psf, mode='same')
    return retinal_image


def retinal_image(lum, psf, lum_distance=None, output_shape=None):
    """calculate retinal image from luminance points

    lum_distance: in pixel, only necessar if more than one lum value
    """
    if not hasattr(lum, '__len__') or len(lum) == 1:
        return lum * psf
    else:
        if output_shape is None:
            retinal_image = np.zeros((psf.shape[0]+lum_distance, psf.shape[1]+lum_distance))
        else:
            retinal_image = np.zeros(output_shape)
        center_0, center_1 = get_middle(retinal_image, 0), get_middle(retinal_image, 1)
        psf_center_0, psf_center_1 = get_middle(psf, 0), get_middle(psf, 1)
        lum_distance_half = int(np.floor(lum_distance / 2))
        if np.mod(psf.shape[0], 2) == 1:
            floor_factor_0 = 1
        if np.mod(psf.shape[0], 2) == 0:
            floor_factor_0 = 0
        if np.mod(psf.shape[1], 2) == 1:
            floor_factor_1 = 1
        if np.mod(psf.shape[1], 2) == 0:
            floor_factor_1 = 0
        retinal_image[center_0-psf_center_0:center_0+psf_center_0+floor_factor_0,
                      center_1-psf_center_1-lum_distance_half:center_1+psf_center_1-lum_distance_half+floor_factor_1] += lum[0] * psf
        retinal_image[center_0-psf_center_0:center_0+psf_center_0+floor_factor_0,
                      center_1-psf_center_1+lum_distance_half:center_1+psf_center_1+lum_distance_half+floor_factor_1] += lum[1] * psf
        return retinal_image


def retina_photons_absorbed(retinal_image, receptor_lattice):
    """return the receptor lattice, with values scaled to show the amount of photons absorbed
    """
    if retinal_image.shape != receptor_lattice.shape:
        dim0 = (receptor_lattice.shape[0] - retinal_image.shape[0]) / 2.
        dim1 = (receptor_lattice.shape[1] - retinal_image.shape[1]) / 2.
        retinal_image = np.pad(retinal_image, ((int(np.ceil(dim0)), int(np.floor(dim0))),
                                               (int(np.ceil(dim1)), int(np.floor(dim1)))),
                               'constant')
    return retinal_image * receptor_lattice


def mean_photons_absorbed(retinal_image, receptor_lattice, a=0.28, d=0.2, s=3.1416, t=0.68,
                          e555=0.5):
    """mean number of photons absorbed (eqt 2)

    lum: luminance distr (parameterized by x, y) in cd/m2

    psf: point-spread function (param. by x, y) [a.u.]

    a: cross-sectional area of the receptor in sq minutes

    d: duration of stimulus in seconds

    s: pupil area in sq mm

    t: transmittance of ocular media

    e555: quantum efficiency of photoreceptor at lambda=555 nm

    we first scale the retinal image by all the constants. we then sample the retinal image where
    there are photoreceptors (as shown by the receptor lattice) and return this list of
    values. this list contains the "mean photons absorbed" or alpha_i values for all photoreceptors
    i.
    """
    if retinal_image.shape != receptor_lattice.shape:
        dim0 = (receptor_lattice.shape[0] - retinal_image.shape[0]) / 2.
        dim1 = (receptor_lattice.shape[1] - retinal_image.shape[1]) / 2.
        retinal_image = np.pad(retinal_image, ((int(np.ceil(dim0)), int(np.floor(dim0))),
                                               (int(np.ceil(dim1)), int(np.floor(dim1)))),
                               'constant')
    photons_absorbed = a*d*s*t*e555 * 347.8 * retinal_image
    return photons_absorbed[receptor_lattice.astype(bool)]


def calc_d_prime(alpha, beta):
    """calculates d prime for two lists of photon absorptions, alpha and beta; equation 3
    """
    alpha = np.array(alpha)
    beta = np.array(beta)
    log_ratio = np.log(beta / alpha)
    # sometimes there will be nans, which we replace with 0
    log_ratio[np.isnan(log_ratio)] = 0
    log_ratio[np.isinf(log_ratio)] = 0
    numerator = np.sum((beta - alpha) * log_ratio)
    denominator = np.sum((beta + alpha) * log_ratio**2)
    return numerator / np.sqrt(.5 * denominator)


def calc_deltaN(alpha, beta):
    """calculate deltaN as used in equation 4

    deltaN: average difference in number of effectively absorbed quanta from stimuli alpha and
    beta, which I'm interpreting to be: $\sum_i(\beta_i - \alpha_i)$ (based on d-prime)
    """
    return np.sum(beta - alpha)


def calc_N(alpha, beta):
    """calculate N, as used throughout the text

    N: mean number of quanta per stimulus, which I'm interpreting to be:
    $\sum_i{\frac{\beta_i+alpha_i}{2}}$ (based on d-prime)
    """
    return np.sum((alpha + beta) / 2)


def check_intensity_discrimination(alpha, beta):
    """checks whether d prime reduces to equation 4

    deltaN and N are really confusingly defined in the text, but I think this is what it
    means.
    """
    return abs(calc_deltaN(alpha, beta)) / np.sqrt(calc_N(alpha, beta))


def intensity_discrimination_task(lum_a, lum_b):
    """run the intensity discrimination task
    """
    psf = pointspread_function()
    rec_lattice, _, _, _, _ = construct_photoreceptor_lattice(minutes_to_n_receptors(8.5))
    ret_im_a = retinal_image(lum_a, psf)
    ret_im_b = retinal_image(lum_b, psf)
    absorbed_a = mean_photons_absorbed(ret_im_a, rec_lattice)
    absorbed_b = mean_photons_absorbed(ret_im_b, rec_lattice)
    return (calc_d_prime(absorbed_a, absorbed_b), calc_N(absorbed_a, absorbed_b),
            calc_deltaN(absorbed_a, absorbed_b))


def optimize_intensity_discrimination_task(lum_a, d_prime=1.36):
    bounds = [(lum_a, np.inf)]
    def obj_func(lum_b):
        return np.square(intensity_discrimination_task(lum_a, lum_b)[0] - d_prime)
    return optimize.minimize(obj_func, lum_a+.01, bounds=bounds)


def figure4(lum_a=[.2, .5, .75, 1, 2, 3, 4, 5, 6, 10, 15, 20], d_prime=1.36):
    """recreate figure 4

    d_prime=1.15 seems to get the best match to the Banks data.
    """
    solts = []
    for a in lum_a:
        solt = optimize_intensity_discrimination_task(a, d_prime)
        solts.append(intensity_discrimination_task(a, solt.x))
    plot_solts = np.log(np.array(solts)[:, 1:])
    plt.plot(plot_solts[:, 0], plot_solts[:, 1], label='Our solution', zorder=3)
    plt.plot([1, 5], [.697, 2.68], label='Geisler, 1984')
    x = np.arange(0, np.exp(6))
    y = 1.36 * np.sqrt(x)
    plt.plot(np.log(x), np.log(y), 'k--', label='$\Delta N = 1.36\sqrt{N}$')
    plt.legend()
    plt.xlabel('LOG N')
    plt.ylabel('LOG $\Delta$N')
    plt.title('INTENSITY DISCRIMINATION')
    return solts


def resolution_task(deltaTheta, lum=4):
    """run the resolution task

    deltaTheta: in units of arc-minutes
    """
    rec_lattice, x_minutes, y_minutes, x, y = construct_photoreceptor_lattice(minutes_to_n_receptors(8.5 + deltaTheta))
    psf = pointspread_function()
    ret_im_a = retinal_image(lum, psf)
    min_per_pix = .02
    deltaThetaPix = int(deltaTheta / min_per_pix)
    ret_im_b = retinal_image([lum / 2., lum / 2.], psf, deltaThetaPix)
    absorbed_a = mean_photons_absorbed(ret_im_a, rec_lattice)
    absorbed_b = mean_photons_absorbed(ret_im_b, rec_lattice)
    return ret_im_a, ret_im_b, rec_lattice, x_minutes, y_minutes
    return (calc_d_prime(absorbed_a, absorbed_b), calc_N(absorbed_a, absorbed_b),
            deltaTheta)
