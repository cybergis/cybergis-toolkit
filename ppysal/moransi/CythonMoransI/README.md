Warning:  This is unverified research code.

##Building##

The cython code needs to be compiled on the system that it will run on.  The contents of this directory have been compiled on an Ubuntu 12.04, 64bit box and may be portable to other architectures.

To build: `python setup.py build_ext --inplace`

We see two build warning that can be safely ignored.

##Running##

To run: `python lmorans.py`  

The code utilizes all available processing cores with an argument to set a defined number.

The test dataset has 3085 observations and 999 permutations are performed.  Performance is < 0.4 seconds on a dual core Ubuntu box.
