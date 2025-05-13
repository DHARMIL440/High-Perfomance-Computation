Here’s the updated `README.md` code with **"Charge deposition using Cloud-in-a-Cell (CIC)"** explicitly mentioned in the introduction and context. This version integrates your academic project purpose with scientific terminology and real-world relevance:

````markdown
# 2D Scattered Data Interpolation using Bilinear Scheme

This project implements bilinear interpolation of scattered 2D data points onto a structured mesh grid. The method used is a form of **Charge Deposition using Cloud-in-a-Cell (CIC)** — a widely used technique in particle-mesh simulations in physics and high performance computing (HPC).

It is modeled as a **"Similar Charge in Cloud"** problem, where values (e.g., charges or mass) from unstructured positions are spread over a structured computational domain. This concept is critical for efficient and accurate numerical simulation in various scientific and engineering applications.

## Problem Description

Given a set of scattered points \( (x_i, y_i, f_i) \), where \( f_i = 1 \), the objective is to interpolate the function values onto a structured grid of size \( M \times M \) using **bilinear interpolation**, based on the **CIC (Cloud-in-a-Cell)** method.

Each scattered point contributes to the values of the four surrounding grid points based on its relative distance to them.

## Files

- `PIC_interpolation.c`:  
  Serial implementation for bilinear interpolation. Can be extended to parallel code using OpenMP or MPI.

- `input_fileMaker.c`:  
  Utility to generate input binary file `input.bin` with random scattered points in a unit domain.

## Applications

This interpolation scheme has widespread applications:

1. **Computer Graphics and 3D Modeling**  
   • Used to generate smooth surfaces from scattered 3D points, such as those obtained from 3D scanning (e.g., LIDAR or photogrammetry).  
   • Example: Turning a point cloud of a scanned object into a mesh for rendering or animation.

2. **Scientific Visualization**  
   • In fields like meteorology or oceanography, scattered data points (e.g., temperature or pressure readings from sensors) are interpolated onto a mesh to visualize continuous fields over a region.

3. **Finite Element Analysis (FEA)**  
   • Engineers use this in simulations (e.g., structural analysis, fluid dynamics) to map scattered data (like material properties or experimental measurements) onto a structured mesh for numerical solving.

4. **Medical Imaging**  
   • Reconstructing surfaces or volumes from scattered data points, such as creating a 3D model of an organ from MRI or CT scan points.

5. **Machine Learning and Data Science**  
   • When dealing with sparse or unevenly sampled data, this scheme helps create a continuous representation for tasks like regression, visualization, or simulation.

## Cloud-in-a-Cell (CIC) – Charge Deposition

For each scattered point within a grid cell (with corner points \((X_i, Y_j)\), \((X_{i+1}, Y_j)\), etc.), the function value is deposited onto the grid using bilinear weights:

- \( w_{i,j} = (1 - dx)(1 - dy) \)  
- \( w_{i+1,j} = dx(1 - dy) \)  
- \( w_{i,j+1} = (1 - dx)dy \)  
- \( w_{i+1,j+1} = dxdy \)

Where:

- \( dx = \frac{x_i - X_i}{\Delta x} \)  
- \( dy = \frac{y_i - Y_j}{\Delta y} \)

These weights ensure smooth interpolation and conservation properties across the grid.

## Build and Run

### 1. Compile

```bash
gcc input_fileMaker.c -o input_gen
gcc PIC_interpolation.c -o interpolate -lm
````

### 2. Generate Input Data

```bash
./input_gen
```

This creates the file `input.bin` containing the scattered data.

### 3. Run Interpolation

```bash
./interpolate
```

### Output

The program outputs interpolated values on a structured grid. Output format depends on your implementation (e.g., binary or text).

## Test Configurations

Use the following setups for testing scalability and performance:

| Config | Nx   | Ny  | Points (Million) | Maxiter |
| ------ | ---- | --- | ---------------- | ------- |
| (a)    | 250  | 100 | 0.9              | 10      |
| (b)    | 250  | 100 | 5.0              | 10      |
| (c)    | 500  | 200 | 3.6              | 10      |
| (d)    | 500  | 200 | 20.0             | 10      |
| (e)    | 1000 | 400 | 14.0             | 10      |

## Performance and Optimization Goals

This project encourages not just accuracy but **high-performance parallel computation**:

* Implement and benchmark parallel versions (e.g., OpenMP, SIMD, MPI)
* Measure:

  * Execution time
  * Speedup vs thread count
  * Cache performance (valgrind, callgrind, Intel VTune)
* Evaluate:

  * Impact of hyperthreading
  * Thread scheduling strategies
  * Effect of data locality

## Notes

* Domain size is fixed to \[0, 1] × \[0, 1]
* Make sure interpolation results are validated against reference output
* Avoid race conditions in parallel implementations
