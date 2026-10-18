# MSIPRM-OT

MATLAB implementation of the **multiscale primal–dual interior-point relaxation method (MSIPRM)** for large-scale discrete optimal transport with squared Euclidean costs.

The method combines a multiscale hierarchy, adaptive sparse active sets, and an interior-point relaxation solver. Coarse-grid solutions initialize finer-grid problems, while support identification reduces the number of transport variables handled by the solver.

This repository accompanies the manuscript:

> Shengyu Sun, Rui-Jin Zhang, Ruoyu Diao, and Yu-Hong Dai.  
> **A Multiscale Primal–Dual Interior-Point Relaxation Method for Large-Scale Optimal Transport Problems.**

## Implementation

### MATLAB version

The solver sources are located in the repository root.

- **Multiscale framework:** solves transport problems on progressively finer image grids.
- **Sparse active sets:** constructs and refines candidate transport supports using coarse-grid solutions and dual information.
- **Interior-point relaxation:** allows warm starts from solutions of preceding subproblems.
- **Reduced linear systems:** uses the sparse structure of the transport constraints when computing Newton directions.

| File or directory | Purpose |
| --- | --- |
| [`demo.m`](demo.m) | Multiscale experiment driver for DOTmark image pairs. |
| [`iprm4denseot.m`](iprm4denseot.m) | Dense IPRM solver used for the coarsest-grid initialization. |
| [`iprm4sparseot.m`](iprm4sparseot.m) | IPRM solver for sparse transport subproblems. |
| [`enlarge_Spt.m`](enlarge_Spt.m) | Constructs a finer-grid candidate support from a coarse-grid transport plan. |
| [`getSpts.m`](getSpts.m) | Searches occupied support blocks for negative reduced-cost entries. |
| [`extract_connected_components_bipartite.m`](extract_connected_components_bipartite.m) | Identifies connected components of the bipartite support graph. |
| [`data/`](data/) | Scripts generating the five synthetic test measures. |
| [`data_2048/`](data_2048/) | A large-scale experiment log and scripts plotting recorded iteration and timing data. |
| [`performanceprofile/`](performanceprofile/) | Stored solver timings and a performance-profile plotting script. |

## Requirements

- **MATLAB.** The experiments described in the manuscript used MATLAB R2024b on an Intel Core i7-14700 machine with 64 GB of RAM.
- **DOTmark data** for the image benchmark in `demo.m`. The dataset is distributed separately; see [DOTmark — A Benchmark for Discrete Optimal Transport](https://arxiv.org/abs/1610.03368).
- **Statistics and Machine Learning Toolbox** for the `randsample` call in the performance-profile script.

## Getting started

### Small solver example

Set MATLAB's current folder to the repository root. The following example solves transport between two measures on a 2-by-2 grid with unit grid spacing, using only the dense solver:

```matlab
a = ones(4, 1) / 4;
b = [1; 1; 3; 3] / 8;
massScale = norm([a; b]);

[x, s, lambda, z, status, iter] = iprm4denseot(4, 4, a, b, 0);

% Restore the mass scaling applied internally by the solver.
P = reshape(x, 4, 4) * massScale;

% Squared Euclidean costs in MATLAB's grid ordering.
index = (0:3)';
row = mod(index, 2);
col = floor(index / 2);
C = (row - row').^2 + (col - col').^2;

transportCost = sum(C(:) .* P(:));
marginalError = norm([sum(P, 1)' - a; sum(P, 2) - b]);
fprintf('Status: %d; iterations: %d\n', status, iter);
fprintf('Transport cost: %.8f; marginal error: %.3e\n', ...
        transportCost, marginalError);
```

The optimal quadratic transport cost is **0.25**. Here, columns of `P` correspond to `a` and rows correspond to `b`.

### Multiscale DOTmark experiments

The support-search routine `getSpts.m` is included in the repository root. It searches blocks occupied by the current support; it does not perform an exhaustive search over all transport pairs.

1. Download and extract the DOTmark data, preserving its category directories and CSV filenames.
2. In the local `loadData` function in `demo.m`, set `basePath` to the extracted `Data` directory. For example, this directory should contain `CauchyDensity/data32_1002.csv`.
3. Set `save_path` in `demo.m` to an existing output directory. The current script contains a machine-specific path and a fixed `1024` directory label; update both for the experiment being run.
4. Choose the resolution through `scale`, the image categories through the `for type = ...` loop, and the image pairs through `pairs`.
5. With the repository root as MATLAB's current folder, run:

```matlab
addpath(pwd);
run('demo.m');
```

The current defaults select **25 instances at 32-by-32 resolution**: five image pairs in each of the first five DOTmark categories. The multiscale hierarchy starts with a **16-by-16** dense initialization. Results are written to CSV files named `matlab-<type>-<scale>.csv` in the configured output directory.

## Experimental data

### Synthetic measures

The scripts in [`data/`](data/) generate normalized discrete measures and save the probability matrix as the variable `Z`.

| Script | Measure | Output file |
| --- | --- | --- |
| `SingleCauchy.m` | Single Cauchy-type component | `Singlecauchy.mat` |
| `SingleGauss.m` | Single Gaussian-type component | `Singlegauss.mat` |
| `dualCauchy.m` | Two Cauchy-type components | `dualcauchy.mat` |
| `dualGauss.m` | Two Gaussian-type components | `dualgauss.mat` |
| `CauchyGauss.m` | Cauchy–Gaussian mixture | `cauchygauss.mat` |

The default resolution is 1024-by-1024. Set `imageSize` in each script to generate another resolution. Use a separate output directory or filename for each resolution, since the scripts use fixed output filenames. The supplied `demo.m` reads DOTmark CSV files; a driver loading `Z` is needed to solve these synthetic instances.

### Timings and plots

- [`performanceprofile/solvingtime/`](performanceprofile/solvingtime/) contains timing arrays for several solvers and resolutions. These are recorded results; the competing solver implementations are separate dependencies.
- [`performanceprofile/plot_newcompare_lp.m`](performanceprofile/plot_newcompare_lp.m) processes timing arrays and plots performance profiles. Before use, update its input paths, sampling indices, shift parameter, time limit, and solver labels for the chosen experiment. The 128-resolution directory is currently named `sclale_128`.
- [`data_2048/Guass-Cauchy_Gauss.txt`](data_2048/Guass-Cauchy_Gauss.txt) records a large-scale synthetic run. [`Iter_time.m`](data_2048/Iter_time.m) and [`Fill_lineartime.m`](data_2048/Fill_lineartime.m) plot values recorded from that experiment.

Exact reproduction of aggregate timings requires the same test-instance indices and settings. The timing arrays do not provide the original instance indices for every sampled subset, and some manuscript experiments require additional drivers. Running times also depend on the hardware and MATLAB version.

## Citation

If this code is used in research, please cite the accompanying manuscript:

> Shengyu Sun, Rui-Jin Zhang, Ruoyu Diao, and Yu-Hong Dai. *A Multiscale Primal–Dual Interior-Point Relaxation Method for Large-Scale Optimal Transport Problems.* Manuscript.

```bibtex
@unpublished{sun_msiprm_ot,
  author = {Sun, Shengyu and Zhang, Rui-Jin and Diao, Ruoyu and Dai, Yu-Hong},
  title  = {A Multiscale Primal--Dual Interior-Point Relaxation Method
            for Large-Scale Optimal Transport Problems},
  note   = {Manuscript}
}
```

## Contact

**Shengyu Sun**  
E-mail: [sunshengyu@lsec.cc.ac.cn](mailto:sunshengyu@lsec.cc.ac.cn)

For questions about the implementation, please include the MATLAB version, input data or instance identifiers, parameter settings, and any error message.

## Copyright

Copyright © 2026 Shengyu Sun.

A license has not yet been specified for this repository. Please contact Shengyu Sun for licensing enquiries.
