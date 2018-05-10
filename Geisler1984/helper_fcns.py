"""functions to re-generate Geisler 1984 figures
"""
import warnings
import itertools
import numpy as np
import pandas as pd
from scipy.signal import convolve2d
from scipy import optimize
from scipy import interpolate
import matplotlib.pyplot as plt


MIN_PER_PIX = .02

# To better estimate and sample photoreceptor absroptions, we need to arbitrarily sample the gaussian-based pointspread functions.
# To this end, we'll create a class (Our_Gauss) which can be used to estimate the PDF of a gaussian at an arbitrary location.
# We approximate the normalization with a trapezoidal sum over the 2D Gaussian function within +/- 4.25 arc-minutes, approximated with "norm_sampling" samples
class Our_Gauss():
    
    def _exp_func(self, x, y, mu_x, mu_y, sigma, amplitude):
        x, y = np.meshgrid(x, y)
        dist = np.sqrt((x - mu_x)**2 + (y - mu_y)**2)
        return amplitude * np.exp(-.5 * (dist/sigma)**2) / (2*sigma)
    
    def __init__(self, mu, amplitude, sigma, norm_sampling=1001):
        self.mu_x = mu[0]
        self.mu_y = mu[1]
        self.amplitude = amplitude
        self.sigma = sigma
        x = np.linspace(-4.25 + mu[0], 4.25 + mu[0], norm_sampling)
        y = np.linspace(-4.25 + mu[1], 4.25 + mu[1], norm_sampling)
        example_psf = self._exp_func(x, y, mu[0], mu[1], sigma, amplitude)
        self.norm_constant = np.trapz(np.trapz(example_psf, x), y)
    
    def pdf(self, x, y, norm=True, diag_only=True):
        value = self._exp_func(x, y, self.mu_x, self.mu_y, self.sigma, self.amplitude)
        if diag_only:
            value = np.diagonal(value)
        if norm:
            return value / self.norm_constant
        else:
            return value
        
# Given a 2D Gaussian, we can then create a pointspread function, which is comprised from the sum of two gaussians (see Geisler, 1984)
# Here, we perform the normalization only after the summation, so that the volume under the entire curve is 1
class Pointspread_Function(Our_Gauss):
    
    def __init__(self, mu, amplitude, sigma, norm_sampling=1001):
        mu = np.array(mu)
        self.gauss1 = Our_Gauss(mu, amplitude[0], sigma[0], norm_sampling)
        self.gauss2 = Our_Gauss(mu, amplitude[1], sigma[1], norm_sampling)
        x = np.linspace(-4.25 + mu[0], 4.25 + mu[0], norm_sampling)
        y = np.linspace(-4.25 + mu[1], 4.25 + mu[1], norm_sampling)
        example_psf = self.gauss1.pdf(x, y, False, False) + self.gauss2.pdf(x, y, False, False)
        self.norm_constant = np.trapz(np.trapz(example_psf, x), y)
        
    def pdf(self, x, y, norm=True):
        value = self.gauss1.pdf(x, y, False) + self.gauss2.pdf(x, y, False)
        if norm:
            return value / self.norm_constant
        else:
            return value

# This function creates (x, y) coordinate pairs for where each photoreceptor in our lattice is centered.
# The construction here is based off of the construct_photoreceptor_lattice below.
# We accomplish this by tiling two x-y row pairs: 
#   the first is centered on 0 in X and Y; 
#   the second is offset in X by receptor_diameter/2, and receptor_height_offset in Y
# In both cases, each x location is offset by receptor_diameter and each y location by 2*receptor_height_offset
def get_photoreceptor_locations(x_minutes, y_minutes):
    """go out x_minutes in each direction
    """
    receptor_diameter = .6
    # this is approximately (.6/2)*tan(pi/3)
    receptor_height_offset = .52
    x1 = np.union1d(np.arange(0, x_minutes, receptor_diameter), np.arange(0, -x_minutes, -receptor_diameter))
    y1 = np.union1d(np.arange(0, y_minutes, 2*receptor_height_offset), np.arange(0, -y_minutes, -2*receptor_height_offset))
    x2 = np.union1d(np.arange(receptor_diameter/2, x_minutes, receptor_diameter), np.arange(receptor_diameter/2, -x_minutes, -receptor_diameter))
    y2 = np.union1d(np.arange(receptor_height_offset, y_minutes, 2*receptor_height_offset), np.arange(receptor_height_offset, -y_minutes, -2*receptor_height_offset))
    return np.vstack([list(itertools.product(x1,y1)), list(itertools.product(x2,y2))]);


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


def construct_photoreceptor_lattice(n_receptors_per_side, scale_factor=0):
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

    scale_factor: int. power of 2 by which to scale up the receptor lattice. If the minutes per
    pixel is .02, we don't always have a fine enough sampling to find the correct dprime value. so
    we need to scale it up (changing min_per_pix to .02 / 2 = .01 or .02 / 4 = .005) in order to
    get a better sampling.

    returns lattice, x_minutes, y_minutes, x, y. x/y_minutes are the number of minutes on either
    side of 0 in the created lattice (so, from example above, 1.36 and 1.04, respectively). x and y
    are the arrays that go from -x/y_minutes to +x/y_minutes and are necessary when calling the
    other functions in this file
    """
    receptor_diameter = .6
    # this is approximately (.6/2)*tan(pi/3)
    receptor_height_offset = .52
    min_per_pix = MIN_PER_PIX / 2**scale_factor
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


def visualize_receptor_lattice(lattice, x, y, mode='continuous', **kwargs):
    """visualize the receptor lattice

    using imshow doesn't work very well, because our receptors are only single pixels, they easily
    get lost.

    x and y are the arrays returned by construct_photoreceptor_lattice

    mode: {'continuous', 'binary', 'categorical'}. only considered if the lattice we get has values
    other than 0 and 1. in that case, if binary we plot points anywhere there's a non-zero value,
    all the same color (set by the kwarg 'c', default black). if continuous, we use the color to
    show what value. if categorical, we need a cmap (dictionary) that maps between the non-zero
    values found in lattice and colors to plot with (by default we use {-1: 'red', 1: 'blue', 2:
    'green'}, but you probably want to specify a better one).
    """
    rec_pos = np.where(lattice)
    if not ((lattice == 0) | (lattice == 1)).all() and mode == 'continuous':
        # then this isn't just 0s and 1s, and we need to grab colors
        c = lattice[rec_pos]
    elif not ((lattice == 0) | (lattice == 1)).all() and mode == 'categorical':
        cmap = kwargs.pop('cmap', {-1: 'red', 1: 'blue', 2: 'green'})
        c = [cmap[i] for i in lattice[rec_pos]]
    else:
        c = kwargs.pop('c', 'k')
    ax = kwargs.pop('ax', plt.gca())
    plotted = ax.scatter(x[rec_pos[1]], y[rec_pos[0]], c=c, **kwargs)
    ax.set_xlim((x.min(), x.max()))
    ax.set_ylim((y.min(), y.max()))
    return plotted


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
        ymin = center_0 - psf_center_0
        ymax = center_0 + psf_center_0 + floor_factor_0
        xmin0 = center_1 - psf_center_1 - lum_distance_half
        xmax0 = center_1 + psf_center_1 - lum_distance_half + floor_factor_1
        xmin1 = center_1 - psf_center_1 + lum_distance_half
        xmax1 = center_1 + psf_center_1 + lum_distance_half + floor_factor_1
        lum0 = lum[0] * psf
        lum1 = lum[1] * psf
        if xmin0 < 0:
            lum0 = lum0[:, -xmin0:]
            lum0 = np.pad(lum0, ((0, 0), (0, -xmin0)), 'constant')
            xmax0 += -xmin0
            xmin0 = 0
        elif xmax0 > retinal_image.shape[1]:
            end_offset = retinal_image.shape[1]-(xmax0+1)
            lum0 = lum0[:, :end_offset]
            lum0 = np.pad(lum0, ((0, 0), (-end_offset, 0)), 'constant')
            xmin0 += end_offset + 1
            xmax0 = retinal_image.shape[1]
        if xmin1 < 0:
            lum1 = lum1[:, -xmin0:]
            lum1 = np.pad(lum1, ((0, 0), (0, -xmin1)), 'constant')
            xmax1 += -xmin1
            xmin1 = 0
        elif xmax1 > retinal_image.shape[1]:
            end_offset = retinal_image.shape[1]-(xmax1+1)
            lum1 = lum1[:, :end_offset]
            lum1 = np.pad(lum1, ((0, 0), (-end_offset, 0)), 'constant')
            xmin1 += end_offset + 1
            xmax1 = retinal_image.shape[1]
        retinal_image[ymin:ymax, xmin0:xmax0] += lum0
        retinal_image[ymin:ymax, xmin1:xmax1] += lum1
        return retinal_image


def retina_photons_absorbed(retinal_image, receptor_lattice):
    """return the receptor lattice, with values scaled to show the amount of photons absorbed
    """
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
    photons_absorbed = a*d*s*t*e555 * 347.8 * retinal_image
    return photons_absorbed[receptor_lattice.astype(bool)]


def calc_d_prime(alpha, beta):
    """calculates d prime for two lists of photon absorptions, alpha and beta; equation 3
    """
    alpha = np.array(alpha)
    beta = np.array(beta)
    log_ratio = np.log(beta / alpha)
    if np.all(log_ratio == 0):
        # in this case, all absorptions are identical and so there's no information to use
        return 0
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
    rec_lattice, _, _, x, y = construct_photoreceptor_lattice(minutes_to_n_receptors(8.5))
    psf = pointspread_function(x, y)
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
    plot_solts = np.log10(np.array(solts)[:, 1:])
    plt.plot(plot_solts[:, 0], plot_solts[:, 1], label='Our solution', zorder=3)
    plt.plot([1, 5], [.697, 2.68], label='Geisler, 1984')
    x = np.arange(0, np.exp(6))
    y = 1.36 * np.sqrt(x)
    plt.plot(np.log10(x), np.log10(y), 'k--', label='$\Delta N = 1.36\sqrt{N}$')
    plt.legend()
    plt.xlabel('LOG N')
    plt.ylabel('LOG $\Delta$N')
    plt.title('INTENSITY DISCRIMINATION')
    return solts


def resolution_task(deltaTheta, lum=4, scale_factor=0):
    """run the resolution task

    deltaTheta: in units of arc-minutes
    """
    rec_lattice, x_minutes, y_minutes, x, y = construct_photoreceptor_lattice(
        minutes_to_n_receptors(8.5 + deltaTheta), scale_factor)
    psf = pointspread_function(x, y)
    ret_im_a = retinal_image(lum, psf)
    min_per_pix = MIN_PER_PIX / 2**scale_factor
    deltaThetaPix = int(deltaTheta / min_per_pix)
    ret_im_b = retinal_image([lum / 2., lum / 2.], psf, deltaThetaPix, psf.shape)
    absorbed_a = mean_photons_absorbed(ret_im_a, rec_lattice)
    absorbed_b = mean_photons_absorbed(ret_im_b, rec_lattice)
    # return ret_im_a, ret_im_b, rec_lattice, x, y
    return (calc_d_prime(absorbed_a, absorbed_b), calc_N(absorbed_a, absorbed_b),
            deltaThetaPix * min_per_pix)


def optimize_resolution_task(lum, d_prime_target=1.36, init_deltaThetaPix=100, init_step_size=8,
                             scale_factor=0):
    """optimize resolution task.

    this uses a custom loop because we have a discrete, convex function

    init_deltaTheta is in units of pixels

    returns deltaTheta in units of arc-minutes
    """
    min_per_pix = MIN_PER_PIX / 2**scale_factor
    def res_task(deltaThetaPix):
        return resolution_task(deltaThetaPix * min_per_pix, lum, scale_factor)
    current_dprime = res_task(init_deltaThetaPix)[0]
    loss = current_dprime - d_prime_target
    current_loss_sign = {True: 1, False: -1}.get(loss > 0)
    current_deltaThetaPix = init_deltaThetaPix
    step_size = int(init_step_size)
    # print(current_deltaThetaPix, loss, res_task(current_deltaThetaPix))
    while step_size >= 1:
        while {True: 1, False: -1}.get(loss > 0) == current_loss_sign:
            current_deltaThetaPix -= step_size * current_loss_sign
            if current_deltaThetaPix <= 0:
                current_deltaThetaPix = 1
            solt = res_task(current_deltaThetaPix)[0]
            loss = solt - d_prime_target
            # print(current_deltaThetaPix, loss, res_task(current_deltaThetaPix))
        current_loss_sign = {True: 1, False: -1}.get(loss > 0)
        step_size = int(step_size / 2.)
    deltaThetasPix = [current_deltaThetaPix-1, current_deltaThetaPix, current_deltaThetaPix+1]
    dprimes = [res_task(i)[0] for i in deltaThetasPix]
    f = interpolate.interp1d(x=dprimes, y=deltaThetasPix)
    deltaThetaMin = f(d_prime_target) * min_per_pix
    print("Interpolated delta theta %s" % deltaThetaMin)
    return resolution_task(deltaThetaMin, lum, scale_factor)


def figure5(lum=[.01, .05, .1, .5, 1, 5, 10, 15, 20], d_prime=1.36, scale_factor=0):
    solts = []
    for l in lum:
        solts.append(optimize_resolution_task(l, d_prime, init_step_size=8*(2**scale_factor),
                                              scale_factor=scale_factor))
    solts = np.array(solts)
    # this returns deltaTheta in arc-minutes, we want it in arc-seconds
    plt.semilogy(np.log10(solts[:, 1]), solts[:, 2] * 60, label='Our solution', zorder=3)
    # plt.plot([1, 5], [.697, 2.68], label='Geisler, 1984')
    # x = np.arange(0, np.exp(6))
    # y = 1.36 * np.sqrt(x)
    # plt.plot(np.log10(x), np.log10(y), 'k--', label='$\Delta N = 1.36\sqrt{N}$')
    plt.legend()
    plt.xlabel('LOG N')
    plt.ylabel('$\Delta\Theta$ (arc-seconds)')
    plt.title('RESOLUTION TASK')
    return solts
