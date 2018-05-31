# Quantum-Free-Particle

## About 1D Quantum Free Particle

---

This program numerically integrates the Schrodinger equation on finite complex scalar fields for simulating interactions of quantum particles under varied observation.

This version implements a second-order in time finite difference method known as the "split-step" Crank-Nicolson method. By calculating energy states using the hamiltonian in both position and momentum space, this program is able to achieve numerically stable integration, which is necessary for finite difference methods.

Each time-iteration, the program evolves the wave function in the position basis. Then we apply a Fourier transform the wave function to evolve the non-linear term of the Hamiltonian in the momentum basis/phase-space. The waveform is then reverse Fourier transformed back into position space in order to repeat this evolution of the waveform.

  <img src="https://latex.codecogs.com/gif.latex?%5Cbg_black%20i%20%28%5Cfrac%7Bd%5Cpsi%7D%7Bdx%7D%29%20%3D%20-%5Cfrac%7B1%7D%7B2%7D%20%28%5Cfrac%7Bd%5Cpsi%7D%7Bdx%7D%29%5E%7B2%7D%20&plus;%20U%28x%29%5Cpsi%28x%29"/>
</p>



where:
<img src="https://latex.codecogs.com/gif.latex?%5Cpsi%28x%29%20%3D%20%5Cpsi_%7B%5Cmathbb%7BR%7D%7D%20%28x%29%20&plus;%20i%20%5Cast%20%5Cpsi_%7B%5Cmathbb%7BI%7D%7D%28x%29"/>
</p>


psi(x)= psireal(x) + i * psiimag(x)

<p align="center">
  <img src="https://github.com/mauckc/1D-Quantum-Free-Particle/blob/master/media/sample1.gif"/>
</p>
<p align="center">
    



## Pre-requisites:

fftw3 library
Mathematica for using "ListPlotter.nb" file for plotting C++ programs output
Folder named slices for storing output data

### Structure

*python/*  the code.

*models/*  contains

*media/*  contains images and video.

## Installing:
$ cd your/destination/
$ git clone https://...

<p align="center">
  <img src="https://github.com/mauckc/headpose/blob/master/media/sample2.gif"/>
</p>

## Compiling:

You must link this code with the fftw3 library. On Unix systems, link with -lfftw3 -lm.


### Dependencies
You need to have C++ and the following:

* [gcc GNU](https://gcc.gnu.org)
* [FFTw](http://fftw.org/)

<p align="center">
  <img src="https://github.com/mauckc/1D-Quantum-Free-Particle/blob/master/media/sample3.gif"/>
</p>


### License

All code in this project is provided as open source under the MIT license


---
-Ross Mauck
