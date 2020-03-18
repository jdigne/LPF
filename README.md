# Local Probing Fields Analysis

Author: Julie Digne 2015-2019  
julie.digne@liris.cnrs.fr


##Description

This code is an implementation of the paper:

*Sparse Geometric Representation Through Local Shape Probing*, Julie Digne, Sébastien Valette, Raphaëlle Chaine, IEEE Transactions on Visualization and Computer Graphics, vol. 24, n. 7, pp2238-2250, July 2018.

(The paper was presented at Symposium on Geometry Processing - SGP - 2018)

* Project's webpage:   
https://perso.liris.cnrs.fr/julie.digne/articles/lpf.html  

* PDF paper:  
https://perso.liris.cnrs.fr/julie.digne/articles/lpf.pdf  

* IEEE final version:  
https://ieeexplore.ieee.org/document/7956272 

## License

This code is released under the GNU GPL v3 License.
See the attached [LICENSE](LICENSE) file.

## Build instructions

The build process relies on cmake which should be installed. Tested on Ubuntu 16.04, 18.04, 19.10.

1. Clone the project
+ cd LPF
+ mkdir build
+ cd build
+ cmake -DCMAKE_BUILD_TYPE=Release ..
+ make

The code is shipped with Eigen and nanoflann and has no other dependency.

##Usage
    lpf INPUT SEEDS OUTPUT -r <radius> -p -A <natoms>  
    lpf INPUT SEEDS OUTPUT -r <radius> -p -R <natoms>  
    lpf INPUT SEEDS OUTPUT -r <radius> -p -D <natoms>  

**INPUT**  (mandatory) ascii file containing the point set to analyze 

**SEEDS**  (mandatory) ascii file containing the seeds where lpf will be anchored

**OUTPUT** (mandatory) output file 
    
**-r**     `<radius>` (recommended) Radius that will be used to interpolate x,y positions. Otherwise r=1 will be used.  

**-p**     (recommended) will perform the computation in parallel.  

**-n**     `<nbins>` number of bins: will compute a unit-square grid with nxn points and keep only the points within radius 1. (All neighborhoods will be rescaled to 1). Alternatively use a specific generic pattern (option -g). If no file, an invalid file and no number of bins is specified a default number of bins (10) is used.  

**-l**     `<lambda>` lambda penalty for the sparse optimization. Default: 0.2.  

**-s**     Save immediately the descriptor (low memory consumption).  

**-A**     `<natoms>` analyse the descriptors.  

**-R**     `<natoms>` analyse the descriptors and resample.  

**-D**     `<natoms>` analyse the descriptors and denoise.  

**-v**     verify the reconstruction from the descriptors and outputs the point cloud. If -v and -s are both set, then -v will prevail.  

**-N**    `<n>` number of iterations for the dictionary computation. Default: 10.

**-g**     `<file>` sets the generic pattern to a precise file. See also '-n' option. If both -g and -n are used, and if the generic pattern file is valid, -g will prevail.  

**-d**     `<file>` input dictionary file  

**-c**     Compute parameterization instead of using a random initial seed orientation

### Example command line for denoising:

The seeds and the input file should be the same:

lpf noisy.txt noisy.txt output.txt -r 0.7 -n 16 -N 20 -D 32 -l 0.2


### Example command line for resampling

The seeds can a subset of the input points.

lpf noisy.txt seeds.txt output.txt -r 0.7 -n 16 -N 20 -R 32 -l 0.2


### Example command line for analysis

The seeds can a subset of the input points.

lpf noisy.txt seeds.txt output.txt -r 0.7 -n 16 -N 20 -A 32 -l 0.2


##Input data

This program takes as input a point cloud to analyze, a set of seeds, and optionally a generic pattern file and a dictionary file.

### Point cloud files

This includes the point cloud to analyse, the set of seeds, and the optional generic pattern files.
These files should be formatted as:  

x1	y1	z1  

x2	y2	z2  

.	.	.  

.	.	.  

xn	yn	zn


The coordinates for the generic pattern points are considered as scales.

###Dictionary file

The (optional) dictionary file should be formatted as a plain matrix. If the generic pattern contains M points and the dictionary contains D atoms, then this file should have *3M* lines and *D* columns.
