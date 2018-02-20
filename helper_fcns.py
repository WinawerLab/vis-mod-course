import numpy as np
from scipy.signal import convolve2d

def gauss_1d(x, a, s):
	# assume center at 0
  return a*np.exp(-0.5*np.power(x/s, 2))/(2*s);

def gauss_2d(coord, a, s):
	# assume center at 0, 0
  
  dist = np.sqrt(np.power(coord[0], 2) + np.power(coord[1], 2));
  return a*np.exp(-0.5*np.power(dist/s, 2))/(2*s);

def lsf(x, a1, s1, a2, s2):

  return gauss_1d(x, a1, s1) + gauss_1d(x, a2, s2);

def psf(coord, a1, s1, a2, s2):
	# returns non-normalized 2d gaussian
	# to normalize, simply divide g_sum by total area (s.t. sum of area is 1)
  g_sum = gauss_2d(coord, a1, s1) + gauss_2d(coord, a2, s2);
  return g_sum;

def absorb(lum, psf, a = 0.28, d = 0.2, s = 3.1416, t = 0.68, e = 0.5):
  # lum is luminance distr (parameterized by x, y) in cd/m2
  # psf is point-spread function (param. by x, y) [a.u.]
  #a: cross-sectional area of the receptor in sq minutes
  #d: duration of stimulus in seconds
  #s: pupil area in sq mm
  #t: transmittance of ocular media
  #e: quantum efficiency of photoreceptor at lambda=555 nm

  # 2D convolution
	# mode/boundary/fillvalue are defaults, but we include them for clarity
	# why go with this? photorecepter, modeled by psf, is discrete in space
	# so, lum go beyond edge (mode=full), but nothing will be absorbed (fillvalue=0)
  l_p_conv = convolve2d(lum, psf, mode='full', boundary='fill', fillvalue=0);
  return a*d*s*t*e * 347.8 * np.sum(np.sum(l_p_conv)), l_p_conv;
