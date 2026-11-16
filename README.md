# MPACT

Metabolomics Peak Analysis Computational Tool.

> **The MPACT desktop app is deprecated.** MPACT now runs in the browser at
> **[mpactms.org](https://mpactms.org)**, with no installation. Its source is in
> [robertsamples/mpactweb](https://github.com/robertsamples/mpactweb).
>
> The desktop app is still available:
>
> - Download a release from the [Releases page](https://github.com/robertsamples/mpact/releases).
> - Its source code is on the [`legacy-desktop`](https://github.com/robertsamples/mpact/tree/legacy-desktop)
>   branch, along with its documentation.

## This repository

The `main` branch holds the MPACT analysis core: the modules that filter a
peak table, compute the statistics and search the databases. MPACT Web copies
these modules unmodified (`scripts/vendor_core.py` in the mpactweb repository)
and runs them on its server.

```
code/            the core modules
code/tests/      their tests
rawdata/PTY087I2 the example dataset the tests run on
benchmarks/      benchmarks of the filtering algorithms
```

To run the tests:

```
pip install -r requirements.txt pytest
python -m pytest code/tests
```

## Citation

If you use MPACT, please cite:

Samples, R. M.; Puckett, S. P.; Balunas, M. J. Metabolomics Peak Analysis
Computational Tool (MPACT): An Advanced Informatics Tool for Metabolomics and
Data Visualization of Molecules from Complex Biological Samples. *Anal. Chem.*
**2023**, *95* (23), 8770–8779.
[doi:10.1021/acs.analchem.2c04632](https://doi.org/10.1021/acs.analchem.2c04632)

## Credits

- **Robert M. Samples** (Center for Mass Spectrometry, Smith College): lead
  developer. Algorithms, UI design and benchmarking of the original MPACT;
  MPACT Web.
- **Yilu Chen** (Department of Food Science, University of Massachusetts
  Amherst): mass defect filter, adduct and neutral-loss grouping filter, and
  the revised isotope filter.
- **Sara P. Puckett**: UI design, testing and benchmarking of the original
  MPACT.
- **Marcy J. Balunas** (University of Michigan): supervision of the original
  project.

The original work was supported by National Science Foundation grant
IOS-1557914 and the Gordon and Betty Moore Foundation.

`code/pvclust.py` implements the multiscale bootstrap of Shimodaira (2004),
as in the R package pvclust (Suzuki and Shimodaira, 2006).

## History and license

MPACT was first released in 2022 as
[BalunasLab/mpact](https://github.com/BalunasLab/mpact). That release was
labelled GNU General Public License version 2, but it was built on PyQt5,
which is licensed under version 3, so version 3 was the license it needed.
This repository corrects the label and is released under the
[GNU General Public License version 3](LICENSE). It has been developed
independently since 2023. Files changed since the original release say so in
their headers.
