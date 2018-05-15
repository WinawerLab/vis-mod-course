Made 4 filters for each retinal location. Summed them to 
make one filter and convolved with a grating of and high (8cpd) and low(1cpd) sf.

I did this for cones densities and then for retinal ganglion cell rf densities.

3 steps:
1 - blur image give cone densities using data provided by Watson 2014 from Curcio & Allen 1990
2 - replicate retinal ganglion cell receptive field density formula from Watson 2014 and create filters
3 - apply rgc rf density filter to already blurred cone image
