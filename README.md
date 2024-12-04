
# **mestel2d** 

`mestel2d` is a simple yet effective C++ 2D particle-mesh N-body code designed to simulate the evolution of Mestel discs. It offers functionality to isolate and analyse specific harmonic contributions to the disc's dynamics.  

This tool employs a Cartesian grid for force computations and uses a cloud-in-cell scheme. For harmonic filtering, the code utilizes a finer polar grid for density estimation.  

It has notably been used in [Fouvry et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015A%26A...584A.129F/abstract), where its operation is extensively described, 
to study the collisional evolution of [Sellwood (2012)](https://ui.adsabs.harvard.edu/abs/2012ApJ...751...44S/abstract)'s disc [see also Roule et al. (in prep)].

---

## **Installation and Compilation**

Follow these steps to install and compile the code:  

1. Clone the repository:  
   ```bash
   git clone <repository_url>
   ```

2. **Dependencies:**  
   Ensure the libraries `fftw3` and `hdf5` are installed.  

3. Update the library paths:  
   Edit the first two lines of `mestel2d/Makefile` to reflect the locations of the installed libraries.  

4. Compile the code:  
   ```bash
   cd mestel2d
   make
   ```

---

## **Quick Test**

The repository includes a test simulation featuring an unstable [Zang N=4](https://dspace.mit.edu/handle/1721.1/27444) disc.  

### Prerequisites:  
- **Python** for generating initial conditions and analyzing/plotting results.  
- **FFmpeg** to create a simulation movie.  

If needed, modify the `test/runtest.sh` script to use your preferred Python environment.  

### Running the test:  
Execute the test with the following command:  
```bash
bash test/runtest.sh
```

---

## **Running Custom Simulations**

Simulation parameters are defined in the `mestel.cc` file. Follows the `test/runtest.sh` example script to run your own.

---

## **Authors**

John Magorrian - original version (main author)

Mathieu Roule - [@MathieuRoule](https://github.com/MathieuRoule) - softening kernel + marginal changes

