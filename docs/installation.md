# Installation

MPACT runs on top of an Anaconda Python distribution and is launched as a
plain Python/PyQt5 script — there is no separate packaged installer.

## 1. Install Anaconda

Download and install Anaconda (or Miniconda) from
[anaconda.com/products/individual](https://www.anaconda.com/products/individual).
MPACT relies on the conda-compiled scientific stack (pandas, matplotlib,
scipy, NumPy) that ships with Anaconda's base environment.

## 2. Get the repository

Either:

- Download and unzip the repository (GitHub: **Code → Download ZIP**), or
- Clone it (`git clone https://github.com/robertsamples/mpact.git`, or with
  GitHub Desktop: **Code → Open with GitHub Desktop**).

It doesn't matter where you place the folder — MPACT's launcher script and
shortcut resolve all paths relative to wherever the repo actually lives.

## 3. First launch

The first time MPACT starts, it checks for a few dependencies that are not
part of a stock Anaconda install (`epam.indigo`, `UpSetPlot`, `squarify`)
and installs any that are missing. This can take a minute and may ask you
to restart the application/kernel once so the newly installed packages can
be imported. See [Getting Started](getting-started.md) for how to launch,
and [Troubleshooting](troubleshooting.md) if a dependency install ever
breaks your environment.

!!! warning "NumPy version hazard"
    If you ever see an error like `A module that was compiled using NumPy
    1.x cannot be run in NumPy 2.x`, your Anaconda environment's NumPy was
    upgraded to 2.x by some `pip install`, which breaks every
    conda-compiled module (pandas, matplotlib). Run:

    ```
    pip install "numpy<2"
    ```

    MPACT's own dependency installer is hardened against causing this (it
    always pins `numpy<2`), but a manual `pip install` of something else in
    the same environment can still trigger it. See
    [Troubleshooting](troubleshooting.md#numpy-2x-environment-break) for
    details.
