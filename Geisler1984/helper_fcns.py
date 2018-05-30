"""functions to re-generate Geisler 1984 figures
"""
import warnings
import itertools
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


class Our_Gauss():
    """create a 2d Gaussian that we can arbitrarily sample

    To better estimate and sample photoreceptor absroptions, we need to arbitrarily sample the
    gaussian-based pointspread functions. To this end, we'll create a class (Our_Gauss) which can
    be used to estimate the PDF of a gaussian at an arbitrary location. We approximate the
    normalization with a trapezoidal sum over the 2D Gaussian function within +/- 4.25 arc-minutes,
    approximated with "norm_sampling" samples
    """
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


class Pointspread_Function(Our_Gauss):
    """create a pointspread function from the sum of two gaussians (see Geisler, 1984)

    Here, we perform the normalization only after the summation, so that the volume under the
    entire curve is 1
    """
    def __init__(self, mu, amplitude=[.684, .587], sigma=[.443, 2.035], norm_sampling=1001):
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


def get_photoreceptor_locations(x_minutes, y_minutes):
    """Returns (x, y) coordinate pairs for the center of each photoreceptor in our lattice.

    We accomplish this by tiling two x-y row pairs:

    - the first is centered on 0 in X and Y

    - the second is offset in X by receptor_diameter/2, and receptor_height_offset in Y (both of
    which are determined by the physiology and taken from the Geisler paper)

    In both cases, each x location is offset by receptor_diameter and each y location by
    2*receptor_height_offset

    x_minutes and y_minutes specify the distance (in arc-minutes) that we go out in either
    direction (so the resulting lattice will go from -x_minutes to x_minutes in the horizontal
    direction and similarly in for the vertical)
    """
    receptor_diameter = .6
    # this is approximately (.6/2)*tan(pi/3)
    receptor_height_offset = .52
    x1 = np.union1d(np.arange(0, x_minutes, receptor_diameter),
                    np.arange(0, -x_minutes, -receptor_diameter))
    y1 = np.union1d(np.arange(0, y_minutes, 2*receptor_height_offset),
                    np.arange(0, -y_minutes, -2*receptor_height_offset))
    x2 = np.union1d(np.arange(receptor_diameter/2, x_minutes, receptor_diameter),
                    np.arange(receptor_diameter/2, -x_minutes, -receptor_diameter))
    y2 = np.union1d(np.arange(receptor_height_offset, y_minutes, 2*receptor_height_offset),
                    np.arange(receptor_height_offset, -y_minutes, -2*receptor_height_offset))
    return np.vstack([list(itertools.product(x1, y1)), list(itertools.product(x2, y2))])


def get_middle(x, dim=0):
    """get the middle index of x along dimension dim
    """
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


def mean_photons_absorbed(psf, photoreceptors, lum, a=.28, d=.2, s=3.1416, t=.68, e555=.5):
    """calculate the mean number of photons absorbed

    psf: the pointspread function, an instance of our Pointspread_Function class

    photoreceptors: list of (x, y) pairs specifying the centers of the photoreceptors in our
    lattice, as returned by our get_photoreceptor_locations function

    lum: float, the luminance of the light we're showing the subject.

    all other parameters are constants taken from the Geisler paper.
    """
    return a*d*s*t*e555*347.8*lum * psf.pdf(photoreceptors[:, 0], photoreceptors[:, 1])


def intensity_discrimination_task(lum_a, lum_b, lattice_x=4.25, lattice_y=4.25, psf=None,
                                  debug=False):
    """run the intensity discrimination task

    Given the luminance of two stimuli at the same location, return the d_prime and N for the ideal
    observer

    lum_a, lum_b: floats, luminances of the two different stimuli

    lattice_x, lattice_y: floats, distance in each direction to go from 0 (in arcminutes). the
    receptor lattice will therefore run from -x_minutes to x_minutes in the horizontal direction
    and -y_minutes to y_minutes in the vertical direction

    psf: instance of our Pointspread_Function class, optional. If not specified, will create the
    default one. The same pointspread function is used for both a and b. You should specify it if
    you are optimizing this function, because initializing the psf can take some time.

    debug: boolean. if True, also returns the photoreceptors coordinates and the lists showing the
    amount absorbed in tasks a and b at each of those locations
    """
    if psf is None:
        psf = Pointspread_Function((0, 0))
    photoreceptors = get_photoreceptor_locations(lattice_x, lattice_y)
    absorbed_a = mean_photons_absorbed(psf, photoreceptors, lum_a)
    absorbed_b = mean_photons_absorbed(psf, photoreceptors, lum_b)
    to_return = [calc_d_prime(absorbed_a, absorbed_b), calc_N(absorbed_a, absorbed_b)]
    if debug:
        to_return.extend([photoreceptors, absorbed_a, absorbed_b])
    return to_return


def optimize_intensity_discrimination_task(lum_a, d_prime_target=1.36):
    """Optimize the intensity discrimination task.

    Given the luminance of one stimulus (lum_a) and the specified d prime value (d_prime_target),
    we optimize to find the luminance of the second stimulus.

    lum_a: float. The luminance of the first stimulus
    """
    # We know that lum_b is always greater than lum_a
    bounds_lum = [(lum_a, np.inf)]
    # the psf center is  at (0, 0)
    psf_all = Pointspread_Function((0, 0))

    # the way optimize.minimize works, we optimize a single function, so this wrapper accomplishes
    # that
    def obj_func(lum_b):
        return np.square(intensity_discrimination_task(lum_a, lum_b, psf=psf_all)[0]
                         - d_prime_target)

    return optimize.minimize(obj_func, lum_a+0.01, bounds=bounds_lum)


def figure4(lum_a=[.2, .5, .75, 1, 2, 3, 4, 5, 6, 10, 15, 20], d_prime=1.36):
    """recreate figure 4

    This does not match the figure in the Geisler paper exactly for some reason; we have a small
    shift.
    """
    solts = []
    for a in lum_a:
        solt = optimize_intensity_discrimination_task(a, d_prime)
        # solt = optimize_intensity_discrimination_task(a, d_prime)
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


def resolution_task(deltaTheta, lum=4, debug=False):
    """run the resolution task

    In this task, the ideal observer is trying to differentiate between one stimulus of luminance
    `lum` at the center of its receptor lattice and two stimuli with half that luminance thta are
    separated by `deltaTheta` (and therefore lie `deltaTheta/2` from the center of the receptor
    lattice).

    returns the d-prime and N (average photons absorbed).

    deltaTheta: float. Space between in units of arc-minutes

    lum: float. Luminance of the first stimulus.

    debug: boolean. if True, also returns the photoreceptors coordinates and the lists showing
    the amount absorbed in tasks a and b at each of those locations.
    """
    single_psf = Pointspread_Function((0, 0), norm_sampling=301)
    double_psf_1 = Pointspread_Function((-deltaTheta/2., 0), norm_sampling=101)
    double_psf_2 = Pointspread_Function((deltaTheta/2., 0), norm_sampling=101)
    photoreceptors = get_photoreceptor_locations(np.maximum(4.25, deltaTheta),
                                                 np.maximum(4.25, deltaTheta))
    absorbed_a = mean_photons_absorbed(single_psf, photoreceptors, lum)
    absorbed_b = (mean_photons_absorbed(double_psf_1, photoreceptors, lum/2.) +
                  mean_photons_absorbed(double_psf_2, photoreceptors, lum/2.))
    to_return = [calc_d_prime(absorbed_a, absorbed_b), calc_N(absorbed_a, absorbed_b)]
    if debug:
        to_return.extend([photoreceptors, absorbed_a, absorbed_b])
    return to_return


def optimize_resolution_task(lum, d_prime_target=1.36):
    """optimize the resolution task.

    Given the luminance of the first stimulus and the specified d-prime value (d_prime_target), we
    optimize the find the separation between the two stimuli (deltaTheta).
    """
    # this is just a guess, but it seems to work fine.
    bounds = [(.01, np.inf)]

    def obj_func(deltaTheta):
        return np.square(resolution_task(deltaTheta, lum)[0] - d_prime_target)

    # this initial guess seems to be approximately correct
    return optimize.minimize(obj_func, 10/(lum/.01), bounds=bounds)


def figure5(lum=[.01, .05, .1, .5, 1, 5, 10, 15, 20], d_prime=1.36, scale_factor=0):
    """recreate figure 5

    there are several differences between this and the figure in the Geisler paper:

    1. We only attempt to recreate one line, the one with no blurring of the point source (the
    bottom line in figure 5, the one labeled with "0.0")

    2. At very low logN (corresponding to low values of lum), our line is no longer
    straight. Unclear why this happens.

    3. The part of our line that is straight is parallel to the corresponding line in Geisler's
    figure, but does not have the same values.
    """
    solts = []
    for l in lum:
        solt = optimize_resolution_task(l, d_prime)
        print("Found solution for luminance %s: %s" % (l, solt.x))
        solts.append([solt.x[0], resolution_task(solt.x, l)[1]])
    solts = np.array(solts)
    # this returns deltaTheta in arc-minutes, we want it in arc-seconds
    plt.semilogy(np.log10(solts[:, 1]), solts[:, 0] * 60, label='Our solution', zorder=3)
    plt.plot([.84025, 5.18603], [58.385, 4.47235], label='Geisler, 1984')
    plt.legend()
    plt.xlabel('LOG N')
    plt.ylabel('$\Delta\Theta$ (arc-seconds)')
    plt.title('RESOLUTION TASK')
    return solts
