# MSIPRM-OT

MATLAB implementation of the **multiscale primal–dual interior-point relaxation method (MSIPRM)** for large-scale optimal transport (OT) problems with the squared Euclidean cost.

The multiscale framework constructs a hierarchy of OT problems at progressively finer levels. At each level, the OT problem is solved over a sequence of adaptively refined active sets, with the initial active set constructed from the support estimate at the previous level. The primal–dual interior-point relaxation method (IPRM) serves as the inner solver for each sparse subproblem.

This repository accompanies the manuscript:

> Shengyu Sun, Rui-Jin Zhang, Ruoyu Diao, and Yu-Hong Dai.  
> **A Multiscale Primal–Dual Interior-Point Relaxation Method for Large-Scale Optimal Transport Problems.**

## Implementation

MSIPRM is implemented in MATLAB. The code is located in the repository root.

- **Multiscale method:** constructs the initial active set for a finer-scale problem from the support estimate at a coarser scale.
- **Active-set updates:** identify violated dual constraints and refine the active set for the full OT problem.
- **Interior-point relaxation method:** uses the solution of the previous subproblem as a warm start.
- **Schur complement system:** computes the Newton direction by solving a reduced Schur complement system derived from the normal equations.

| File or directory | Purpose |
| --- | --- |
| [`demo.m`](demo.m) | Numerical experiments on the DOTmark dataset. |
| [`iprm4denseot.m`](iprm4denseot.m) | Solves the dense OT problem at the coarsest level using IPRM. |
| [`iprm4sparseot.m`](iprm4sparseot.m) | Solves sparse subproblems using IPRM. |
| [`enlarge_Spt.m`](enlarge_Spt.m) | Constructs the initial active set at a finer level from the support estimate at the previous level. |
| [`getSpts.m`](getSpts.m) | Identifies violated dual constraints for active-set updates. |
| [`extract_connected_components_bipartite.m`](extract_connected_components_bipartite.m) | Identifies the connected components of the bipartite graph associated with the active set. |
| [`data/`](data/) | Constructs the five test measures described in Section 6.4 and Figure 4. |
| [`data_2048/`](data_2048/) | Computational details and plotting scripts for an OT problem on a $2048\times2048$ grid. |
| [`performanceprofile/`](performanceprofile/) | Running times and the script for plotting running-time performance profiles. |

For active-set updates, `getSpts.m` evaluates dual constraints in matrix blocks containing entries of the current active set. Dual constraints in the remaining blocks are not evaluated by this function.

## Experimental Setup

- **MATLAB.** The experiments in the manuscript were performed using MATLAB R2024b on a machine equipped with an Intel Core i7-14700 processor and 64 GB of RAM.
- **DOTmark dataset.** Download the dataset separately for the experiments in `demo.m`; see [DOTmark — A Benchmark for Discrete Optimal Transport](https://arxiv.org/abs/1610.03368).
- **Statistics and Machine Learning Toolbox.** Required for the `randsample` call in the performance-profile script.

### DOTmark experiments

1. Download and extract the DOTmark dataset, preserving the directories for the image classes and the CSV filenames.
2. In the local `loadData` function in `demo.m`, set `basePath` to the extracted `Data` directory. For example, this directory should contain `CauchyDensity/data32_1002.csv`.
3. Set `save_path` in `demo.m` to an existing output directory. Update the path and the fixed `1024` directory name to match the experiment.
4. Set `scale` to the image resolution, choose the image classes in the `for type = ...` loop, and specify the image pairs in `pairs`.
5. Set MATLAB's current folder to the repository root and run:

```matlab
addpath(pwd);
run('demo.m');
```

The current settings in `demo.m` select **25 test instances at resolution $32\times32$**: five image pairs from each of the first five DOTmark image classes. The OT problem at the coarsest level is solved on a **$16\times16$ grid**. Results are saved as `matlab-<type>-<scale>.csv` in the specified output directory.

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

Reproducing the reported SGM values requires the same test instances and parameter settings. The recorded running times do not include the original image-pair indices for every set of randomly selected instances, and additional scripts are needed for some experiments in the manuscript. Running times also depend on the machine and MATLAB version.

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

For questions about the implementation, please include the MATLAB version, test instances, parameter settings, and any error message.

## Copyright

Copyright © 2026 Shengyu Sun.

A license has not yet been specified for this repository. Please contact Shengyu Sun for licensing enquiries.
