# MSIPRM-OT

MATLAB implementation of the **multiscale primal–dual interior-point relaxation method (MSIPRM)** for large-scale optimal transport (OT) problems with the squared Euclidean cost.

See [Citation](#citation) for the accompanying paper.

## Implementation

- **Multiscale method:** constructs a hierarchy of OT problems at progressively finer levels and initializes the active set at each finer level from the support estimate at the previous level.
- **Active-set updates:** add indices corresponding to violated dual constraints to the active set at the current level.
- **Interior-point relaxation method:** solves each sparse subproblem using the solution of the previous subproblem as a warm start.
- **Schur complement system:** computes the Newton direction by solving a reduced Schur complement system derived from the normal equations.

| File or directory | Purpose |
| --- | --- |
| [`demo.m`](demo.m) | Numerical experiments on the DOTmark dataset. |
| [`iprm4denseot.m`](iprm4denseot.m) | Solves the dense OT problem at the coarsest level using IPRM. |
| [`iprm4sparseot.m`](iprm4sparseot.m) | Solves sparse subproblems using IPRM. |
| [`enlarge_Spt.m`](enlarge_Spt.m) | Constructs the initial active set at a finer level from the support estimate at the previous level. |
| [`getSpts.m`](getSpts.m) | Identifies violated dual constraints for active-set updates. |
| [`extract_connected_components_bipartite.m`](extract_connected_components_bipartite.m) | Counts the connected components of the bipartite graph and, when it is disconnected, adds edges from the original support set to repair the active set. |
| [`data/`](data/) | Constructs the five test measures described in Section 6.4 and Figure 4. |
| [`data_2048/`](data_2048/) | Computational details and plotting scripts for an OT problem on a $2048\times2048$ grid. |
| [`performanceprofile/`](performanceprofile/) | Running times and the script for plotting running-time performance profiles. |

For active-set updates, `getSpts.m` evaluates dual constraints in matrix blocks containing entries of the current active set. Dual constraints in the remaining blocks are not evaluated by this function.

## Experimental Setup

- **MATLAB.** The experiments in the paper were performed using MATLAB R2024b on a machine equipped with an Intel Core i7-14700 processor and 64 GB of RAM.
- **DOTmark dataset.** Download the dataset separately for the experiments in `demo.m`; see [DOTmark — A Benchmark for Discrete Optimal Transport](https://arxiv.org/abs/1610.03368).
- **Statistics and Machine Learning Toolbox.** Required for the `randsample` call in the performance-profile script.

### DOTmark experiments

1. Download and extract the DOTmark dataset, preserving the directories for the image classes and the CSV filenames.
2. Place `DOTmark_1.0` in the repository root so that, for example, `DOTmark_1.0/Data/CauchyDensity/data32_1002.csv` is available. If the dataset is stored elsewhere, update `basePath` in the local `loadData` function in `demo.m`.
3. Set `type` to the image class, `num1` and `num2` to the two image indices, and `scale` to the image resolution.
4. Set MATLAB's current folder to the repository root and run:

```matlab
addpath(pwd);
run('demo.m');
```

The current settings in `demo.m` solve one OT problem between images **1004 and 1005 in the ClassicImages class at resolution $512\times512$**. The OT problem at the coarsest level is solved on a **$16\times16$ grid**.

The batch loops and CSV output are currently commented out. For batch experiments, enable the `pairs` definition, the loops, and the corresponding result collection and CSV output code. Set `save_path` to an existing output directory and replace the fixed `1024` directory name with the actual resolution. Results are saved as `matlab-<type>-<scale>.csv` when CSV output is enabled.

## Numerical Results

### Test measures

The scripts in [`data/`](data/) discretize the two-dimensional probability distributions described in Section 6.4. Each script normalizes the discrete measure to have unit total mass and saves it as the variable `Z`.

| Script | Distribution in Figure 4 | Output file |
| --- | --- | --- |
| `SingleCauchy.m` | A single Cauchy distribution | `Singlecauchy.mat` |
| `SingleGauss.m` | A single Gaussian distribution | `Singlegauss.mat` |
| `dualCauchy.m` | Two Cauchy distributions | `dualcauchy.mat` |
| `dualGauss.m` | Two Gaussian distributions | `dualgauss.mat` |
| `CauchyGauss.m` | A mixture of a Cauchy distribution and a Gaussian distribution | `cauchygauss.mat` |

The current resolution in these scripts is $1024\times1024$. Set `imageSize` to generate the measures at another resolution. Use a separate output directory or filename for each resolution, since the scripts use fixed output filenames. The supplied `demo.m` reads DOTmark CSV files; experiments with these test measures require a script that loads `Z` from the MAT files.

### Running times and performance profiles

- [`performanceprofile/solvingtime/`](performanceprofile/solvingtime/) contains recorded running times for several solvers and resolutions. The implementations of the competing solvers are obtained separately.
- [`performanceprofile/plot_newcompare_lp.m`](performanceprofile/plot_newcompare_lp.m) computes the shifted geometric mean (SGM) and normalized SGM and plots running-time performance profiles. Before use, update the input paths, selected test instances, shift parameter, time limit, and solver names for the experiment. The directory for resolution $128\times128$ is currently named `sclale_128`.
- [`data_2048/Guass-Cauchy_Gauss.txt`](data_2048/Guass-Cauchy_Gauss.txt) records an experiment on a $2048\times2048$ grid. [`Iter_time.m`](data_2048/Iter_time.m) plots the number of IPRM iterations and running time. [`Fill_lineartime.m`](data_2048/Fill_lineartime.m) plots the average time per step for solving the linear system and the average filling rate of the active set.

For each image class, ten image pairs are selected at each of the $128\times128$ and $256\times256$ resolutions, and five image pairs are selected at the $512\times512$ resolution. The selected image pairs are listed in the `pairs` definitions in `demo.m`. Reproducing the reported SGM values requires the same test instances and parameter settings. Running times may vary with the machine and MATLAB version.

## Citation

If this code is used in research, please cite:

> Shengyu Sun, Rui-Jin Zhang, Ruoyu Diao, and Yu-Hong Dai.
> *A Multiscale Primal–Dual Interior-Point Relaxation Method for Large-Scale Optimal Transport Problems.*
> [arXiv:2608.12060](https://arxiv.org/abs/2608.12060), 2026.

```bibtex
@misc{sun_msiprm_ot,
  author        = {Sun, Shengyu and Zhang, Rui-Jin and Diao, Ruoyu and Dai, Yu-Hong},
  title         = {A Multiscale Primal--Dual Interior-Point Relaxation Method
                   for Large-Scale Optimal Transport Problems},
  year          = {2026},
  eprint        = {2608.12060},
  archivePrefix = {arXiv},
  primaryClass  = {math.OC},
  url           = {https://arxiv.org/abs/2608.12060}
}
```
