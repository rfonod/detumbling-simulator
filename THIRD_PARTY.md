# Third-Party Code

The MIT License in [`LICENSE`](LICENSE) covers the code written for this project. The files
listed below were written by others and are redistributed here for convenience; they remain
the work of their respective authors and are governed by their own terms, not by this
repository's license.

## Bundled MATLAB functions

| File | Author | Origin | Purpose here |
| :--- | :--- | :--- | :--- |
| [`putvar.m`](putvar.m) | John D'Errico | [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/27106-putvar-uigetvar) | Copies the simulation variables into the base workspace so that `plots.m` can read them |
| [`truncatedGaussian.m`](truncatedGaussian.m) | Bruno Luong | [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/23832-truncated-gaussian) | Draws the 3σ-truncated Gaussian dispersions used by the Monte Carlo campaigns |
| [`parfor_progress.m`](parfor_progress.m) | Jeremy Scheff | [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor-progress-bar-that-works-with-parfor) | Progress bar for the `parfor` loop in `run_mc.m` |

These three files carry an author attribution but no license text in their headers. They are
included as obtained from the File Exchange. If you redistribute this repository, check the
current terms on each submission's File Exchange page.

## NRLMSISE-00 atmosphere model

[`Future/NRLMSISE-00/`](Future/NRLMSISE-00/) contains the NRLMSISE-00 empirical atmosphere
model: the C implementation by **Dominik Brodowski**, with the MATLAB interface by
**Meysam Mahooti**. It is distributed under the 3-clause BSD license reproduced in
[`Future/NRLMSISE-00/license.txt`](Future/NRLMSISE-00/license.txt), and is used only by the
experimental `Future/` development line, not by the main simulator.

## Models and references

The simulator implements, rather than bundles, the following:

- **IGRF-12** geomagnetic field, evaluated through the Aerospace Toolbox function `igrfmagm`.
- The **B-dot gain derivation** follows G. Avanzini and F. Giulietti, *Magnetic Detumbling of a
  Rigid Spacecraft*, Journal of Guidance, Control, and Dynamics, 35(4), 2012.
