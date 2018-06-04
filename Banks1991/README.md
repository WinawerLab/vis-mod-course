#READ ME

This is the read me text file for the **Banks1991** folder, within the winawerlab/vis-mod-course github repository (https://github.com/WinawerLab/vis-mod-course/tree/master/Banks1991).

This repository is Matlab based and uses the ISETBio toolbox (https://github.com/isetbio/isetbio).

There is a ToolboxToolbox registration file in the forked Winawer lab ToolboxRegistry called 'VisualModelingGradCourse'. If you use the ToolboxToolbox (https://github.com/ToolboxHub/ToolboxToolbox), you can type tbUse('VisualModelingGradCourse') to add this repository to your paths.

Code in this folder makes an attempt to reproduce spatial frequency limits as a function of eccentricity, for an ideal observer model (Figure 3 and 5) from the Banks et al. (1991) paper.

Results can be obtained by running the script: s_mainAnalysis.m
Figure 3 and 5 can be obtained by running the script: s_visualizeResults.m


##Banks1991 folder structure:

- banksRootPath.m 		: function to set paths, needs to be at highest level Banks folder
- figures 			: folder to save figures, content will be ignored by .gitignore
- functions 			: folder with sub functions called by the s_mainAnalysis script
- results			: folder to save results (d prime, per eccentricity, per cone segment), 
				  content will be ignored by .gitignore
- scripts 			: folder containing the two main scripts (and html code for published scripts):
					```s_mainAnalysis
				  	s_visualizeResults```

##Next to do’s:
- Look into differences ISETBIO and Geisler computation of cone absorptions
- What target sizes are used for the stimuli at spatial frequencies and eccentricities not reported in table 1?
- Why is the ideal observer worse at the vary small spatial frequencies, compared to intermediate spatial frequencies?
- Implement retinal ganglion cell layer
- Implement option to simulate human observer results (computational observer model, add photon noise, human optics, etc)

###Full citation:
Peripheral spatial vision: limits imposed by optics, photoreceptors, and
receptor pooling. (1991) Banks, M.S., Sekuler, S.B., & Anderson, S.J.
Journal of the Optical Society of America Association (JOSAA). Vol 8, No 11, 1775-1787.