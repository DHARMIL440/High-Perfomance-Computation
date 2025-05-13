# 2D Scattered Data Interpolation using Bilinear Scheme

This project implements **bilinear interpolation** of scattered 2D data points onto a structured mesh grid using a method inspired by **Cloud-in-a-Cell (CIC)** — a widely adopted approach in particle-mesh simulations within physics and high-performance computing (HPC).

The problem is modeled as a **"Similar Charge in Cloud"** scenario, where values (e.g., charge, mass) from unstructured positions are deposited onto a structured computational domain. This concept is fundamental for efficient and accurate numerical simulations in scientific and engineering applications.

---

## 📌 Problem Description

Given a set of scattered points $(x_i, y_i, f_i)$, where $f_i = 1$, the goal is to interpolate function values onto a structured grid of size $M \times M$ using **bilinear interpolation** following the **CIC (Cloud-in-a-Cell)** method.

Each scattered point contributes to the values of the four nearest grid points based on its relative distance to them.

---

## 📁 Project Files

* **`PIC_interpolation.c`**
  Serial implementation of the bilinear interpolation algorithm. Can be extended for parallelism using OpenMP or MPI.

* **`input_fileMaker.c`**
  Utility to generate a binary input file `input.bin` with randomly distributed scattered points in the unit square domain.

---

## 🧠 Applications

This interpolation scheme finds broad application in the following domains:

1. **Computer Graphics & 3D Modeling**

   * Smooth surface generation from point clouds (e.g., from LIDAR or photogrammetry).
   * Example: Converting a scanned 3D object into a renderable mesh.

2. **Scientific Visualization**

   * Interpolating sensor data (e.g., temperature, pressure) for visualization in fields like meteorology or oceanography.

3. **Finite Element Analysis (FEA)**

   * Mapping scattered experimental or physical data onto a structured mesh for numerical simulations (e.g., stress analysis, CFD).

4. **Medical Imaging**

   * Reconstructing anatomical surfaces or volumes from scattered MRI/CT data.

5. **Machine Learning & Data Science**

   * Transforming sparse or uneven data into continuous grids for regression, analysis, or modeling tasks.

---

## 🧮 Cloud-in-a-Cell (CIC) – Charge Deposition

Each scattered point within a grid cell (bounded by corners $(X_i, Y_j)$, $(X_{i+1}, Y_j)$, etc.) contributes to surrounding grid points using bilinear weights:

$$
\begin{aligned}
w_{i,j} &= (1 - dx)(1 - dy) \\
w_{i+1,j} &= dx(1 - dy) \\
w_{i,j+1} &= (1 - dx)dy \\
w_{i+1,j+1} &= dxdy
\end{aligned}
$$

Where:

$$
dx = \frac{x_i - X_i}{\Delta x}, \quad dy = \frac{y_i - Y_j}{\Delta y}
$$

These weights ensure smooth, conservative deposition across the grid.

---

## ⚙️ Build and Run

### 1. Compile

```bash
gcc input_fileMaker.c -o input_gen
gcc PIC_interpolation.c -o interpolate -lm
```

### 2. Generate Input Data

```bash
./input_gen
```

Generates `input.bin` containing scattered point data.

### 3. Run Interpolation

```bash
./interpolate
```

### Output

The program outputs interpolated values over a structured grid. Output format may vary (binary/text), based on implementation.

---

## 🧪 Test Configurations

Suggested configurations to evaluate performance and scalability:

| Config | Nx   | Ny  | Points (Million) | Max Iterations |
| ------ | ---- | --- | ---------------- | -------------- |
| (a)    | 250  | 100 | 0.9              | 10             |
| (b)    | 250  | 100 | 5.0              | 10             |
| (c)    | 500  | 200 | 3.6              | 10             |
| (d)    | 500  | 200 | 20.0             | 10             |
| (e)    | 1000 | 400 | 14.0             | 10             |

---

## 🚀 Performance and Optimization Goals

This project not only emphasizes accuracy but also aims for **high-performance parallel computation**. Optimization directions include:

* Implement parallel versions (OpenMP, MPI, SIMD)
* Benchmark performance:

  * Execution time
  * Speedup vs. thread count
  * Cache behavior (Valgrind/Callgrind, Intel VTune)
* Evaluate:

  * Impact of hyperthreading
  * Thread scheduling strategies
  * Effects of memory locality

---

## 📝 Notes

* The simulation domain is fixed to $[0, 1] \times [0, 1]$.
* Always validate interpolation results against a reference or analytical benchmark.
* Take care to avoid **race conditions** when parallelizing the algorithm.
