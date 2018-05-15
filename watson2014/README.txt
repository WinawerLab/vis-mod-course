Here I will be replicating some of the math in the following paper:

Watson, A. B. (2014). A formula for human retinal ganglion cell receptive field density as a 
function of visual field location. Journal of Vision, 14(7), 15-15.

The goal is to then use this paper to downsample/blur an image to match density of cones
and/or retinal ganglion cells across temporal/nasal/superior/inferior retina
---------
Done so far:
- replicated first 5 figures in watson paper
- used data from paper to blur an image given the different cone densities across all retinal locations
- same thing but for retinal ganglion cell densities
  -add rgcd blur on top of cone blur
---------
I think I'm done!
---------
The code_blur folder has the imput/ouput images as well 
