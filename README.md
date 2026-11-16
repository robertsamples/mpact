# MPACT (legacy desktop app)

> This is the deprecated desktop version of MPACT, kept for existing users.
> MPACT now runs in the browser at [mpactms.org](https://mpactms.org). Desktop
> builds are on the [Releases page](https://github.com/robertsamples/mpact/releases).

## Installation ##

1. Download and install Anaconda from https://www.anaconda.com/products/individual
2. Download and unzip repository (code>download ZIP) or clone with GitHub Desktop (code>open with GitHub Desktop)

## Startup ##

_Windows_

Double click MPACT shortcut in the main directory (note that nonfatal errors will cause program termination when run in this way and not through an IDE, if you wish to continue in the event of these errors launch through spyder using the Mac steps)

_Mac_

1. Launch spyder through anaconda navigator, or by typing 'spyder' into the anaconda prompt
2. Navigate to code directory in the  main project folder
3. Open main.py (not __main__.py with underscores) in spyder
4. click run button (green arrow in top toolbar)

Full instructions for data processing are included in the MPACT manual.

_Updates in Rev 23.05.15_

- Added data export tab
- Added GNPS peak table filtering functionality (experimental, only tested with FBMN export in MS-DIAL)
- Removed analysis info tab

_Updates in Rev 23.02.26_

- Added filtering of database hits by kingdom and genus
- Dependency check script spawns a window if any dependencies were installed informing user the kernal will restart.
- Bugfixes in group set interactions, and bruker data import

_Updates in Rev 23.02.19_

- Added support for Bruker Metaboscape peak lists
- Added support for MS-DIAL MSP files
- Added MSP file writer to export in-source fragmentation patterns
- Reduced lag when selecting features
- Added multithreading when generating feature abundance plots
- Fixed bug that resulted in feature highlights not being visible in data reanalysis
- Fixed issue preventing samples/injections from being selected and highlighted in NMDS plots
- Database hit images viewed are saved as images in a folder to eliminate rerendering

## History and license

MPACT was first released in 2022 as
[BalunasLab/mpact](https://github.com/BalunasLab/mpact). That release was
labelled GNU General Public License version 2, but it was built on PyQt5,
which is licensed under version 3, so version 3 was the license it needed.
This repository corrects the label and is released under the
[GNU General Public License version 3](LICENSE). Files changed since the
original release say so in their headers.

`code/pvclust.py` implements the multiscale bootstrap of Shimodaira (2004),
as in the R package pvclust (Suzuki and Shimodaira, 2006).
