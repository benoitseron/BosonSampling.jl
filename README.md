# BosonSampling

This project implements standard and scattershot BosonSampling in Julia, including boson samplers and certification and optimization tools.

## Functionalities

A wide variety of tools are available:
* Boson-samplers, including partial distinguishability and loss
* Bunching tools and functions
* Various tools to validate experimental boson-samplers
* User-defined optical circuits built from optical elements
* Optimization functions over unitary matrices
* Photon counting tools for subsets and partitions of the output modes
* Tools to study permanent and generalized matrix function conjectures and counter-examples

## Installation

To install the package, launch a Julia REPL session and type

    julia> using Pkg; Pkg.add("BosonSampling")

Alternatively type on the `]` key. Then enter

    add BosonSampling

To use the package, write

    using BosonSampling

in your file.

## Authors

This package is written by Benoit Seron and Antoine Restivo. The original research presented in the package is done in collaboration with Dr. Leonardo Novo, Prof. Nicolas Cerf.
