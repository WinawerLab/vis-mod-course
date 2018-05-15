Here I will be replicating some of the math in the following paper:

Watson, A. B. (2014). A formula for human retinal ganglion cell receptive field density as a 
function of visual field location. Journal of Vision, 14(7), 15-15.

The goal is to then use this paper to downsample/blur an image to match density of cones
and/or retinal ganglion cells across temporal/nasal/superior/inferior retina
---------
Done so far:
- replicated first 5 figures in watson paper
- used data from paper to blur an image given the different cone densities across all retinal locations
- create filters of retinal ganglion cells receptive field densities (using math in watson paper) and apply 
  to blurred cone image.
---------
I think I'm done!
---------
The code_blur folder has the imput/ouput images as well 
