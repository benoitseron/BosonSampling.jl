---
title: 'BosonSampling.jl: a Julia package for linear quantum interferometry'
tags:
  - Julia
  - Boson sampling
  - Interferometry
  - Quantum optics
authors:
  - name: Benoit Seron.
    orcid: 0000-0002-1710-6587
    equal-contrib: true
    affiliation: "1"
    corresponding: true # (This is how to denote the corresponding author)
  - name: Antoine Restivo.
    orcid: xxx
    equal-contrib: true
    affiliation: "1"
    corresponding: false # (This is how to denote the corresponding author)
affiliations:
 - name: Université Libre de Bruxelles (ULB)
   index: 1
date: 7 June 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Quantum computing has grown extremely fast within the last decade. Although experimental realizations of practical general purpose quantum computers seems to be a distant wish, experimental realizations of quantum supremacy - a controlled task that cannot be realistically simulated by any classical super computer - were achieved in the last few years.

One of the main areas of research regards a primitive type of quantum computer, called a boson sampler. It is the quantum equivalent to Galton's board: photons are sent in a linear interferometer, and we observe their output patterns.

Only ten years elapsed between the theoretical proposal of Boson Sampling and going from manipulating only a few photons to hundreds today, thereby demonstrating a controlled task that no classical computer can realistically perform.

# Statement of need

`BosonSampling.jl` is a Julia package for simulating linear optics interferometry, including various models of Boson Samplers. Given the speed with which this research domain grows, Julia is an ideal language for researchers: it allows for easy to write code, while providing fast execution times out-of-the-box.

We laid out a basic type hierarchy that allows users to access a practical API while still being able to create new models. These new models can include anything from imperfect detectors to new types of Boson Sampling (nonlinear,...). Our package allows user to conveniently keep the same syntax even though the underlying models or calculations may be changed according to researcher's need.

A wide variety of practical functions are implemented in this package, including:
* Boson-samplers, including partial distinguishability and loss
* Bunching tools and functions
* Various tools to validate experimental boson-samplers
* User-defined optical circuits built from optical elements
* Optimization functions over unitary matrices
* Photon counting tools for subsets and partitions of the output modes
* Tools to study permanent and generalized matrix function conjectures and counter-examples

Speed of execution is a main interest of this package: as the battle for quantum supremacy rages, fast execution times are critical.

We also provide entierly novel tools regarding the validation of experiments: users can input experimental samples and obtain metrics on the quality of the experiment (photon indistinguishability,...). These too are crucial for asserting if an experiment lies in the quantum supremacy regime. To the knowledge of the authors, this is the first open-source package containing such tools.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Dr. Leonardo Novo, Prof. Nicolas Cerf, Carlos Fernandes and Dr. Fulvio Flamini during the genesis of this project.

# References
