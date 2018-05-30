# Geisler 1984

This is an attempt to implement: Geisler, W. S. (1984). Physical limits
of acuity and hyperacuity. Journal of the Optical Society of America
A, 1(7), 775. http://dx.doi.org/10.1364/josaa.1.000775.

It currently reproduces figures 1 and (a modified version of) 2
without any problem. We can generate figure 4 and part of figure 5,
though the lines we end up creating are different, see the notebook
for an explanation.

## Install

This was created with python 3 and the required packages can be found
in `requirements.txt`.

If you are new to python, download python 3
from [Anaconda](https://www.anaconda.com/download/) and follow the
instructions. Then navigate to this directory on your command line and
do the following:

```
virtualenv --python=python3 geisler
source geisler/bin/activate
pip3 install -r requirements.txt
```

This will install all required packages, including Jupyter (see
the
[Jupyter lab docs](http://jupyterlab.readthedocs.io/en/stable/index.html) for
more info), which is necessary to view the notebook. Then, to open the
notebook, type `jupyter lab Geisler1984.ipynb`.

When you've finished, type `deactivate` to get your normal python path
back. All requirements are installed within the `geisler/` directory
created by the `virtualenv` call above and so you can remove that
directory when you're done to uninstall everything.

## Usage

Open the notebook (using `jupyter`). A variety of functions that we
call are stored within `helper_fcns.py`.
