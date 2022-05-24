# BosonSampling

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://AntoineRestivo.github.io/BosonSampling.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://AntoineRestivo.github.io/BosonSampling.jl/dev)

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

## Related package

  The present package takes advantage of efficient computation of matrix permanent from [Permanents.jl](https://github.com/benoitseron/Permanents.jl.git).  

## Authors & License

  - [Benoît Seron](mailto:benoitseron@gmail.com)
  - [Antoine Restivo](mailto:antoinerestivo@hotmail.fr)

Contact can be made by clicking on our names.

The original research presented in the package is done in collaboration with Dr. Leonardo Novo, Prof. Nicolas Cerf.

BosonSampling.jl is licensied under the [MIT license](https://github.com/benoitseron/BosonSampling.jl/blob/main/LICENSE).
