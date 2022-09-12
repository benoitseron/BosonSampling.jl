# Embedding BosonSampling.jl in Python

[`Embedding Julia`](https://docs.julialang.org/en/v1/manual/embedding/) can be done in many ways and in several programming languages such as C/C++, C#, Fortran or Python. As an example we describe here how to integrate `BosonSampling.jl` to a python project by using [`PyJulia`](https://pyjulia.readthedocs.io/en/latest/installation.html) that provides a high-level interface to Julia through Python. In the following, we will assume that both Julia and Python are installed and added the `PATH`.

## Installation

The first step is to install the Julia module [`PyCall`](https://github.com/JuliaPy/PyCall.jl)

    julia> using Pkg
    julia> Pkg.add("PyCall")

Make sure that Julia uses your default Python distribution by using the Julia REPL in Shell mode (by entering `;`) through the command `which python`. If the answer is another Python distribution than your default one, reset the path before building `PyCall`:

    julia> ENV["PYTHON"] = "/path/to/python/binary"
    julia> Pkg.build("PyCall")

One can now install `PyJulia` via `pip` through

    $ python -m pip install --user julia

or directly via the Github repository: [`https://github.com/JuliaPy/pyjulia`]( https://github.com/JuliaPy/pyjulia).

Finally, install the remaining dependencies requiered by `PyJulia`:

    $ python
    >>> import julia
    >>> julia.install

It might still be possible that the Python interpreter used by `PyCall.jl` is not supported by `PyJulia`, throwing errors when importing modules from Julia. In this case, run the following lines

    >>> from julia.api import Julia
    >>> jl = Julia(compiled_modules=False)

## Example and workarounds

Once Julia imported in your Python project, any Julia module can be loaded via the standard `import`:

    >>> from julia import Main
    >>> from julia import Base
    >>> from julia import BosonSampling as bs

to finally use `BosonSampling.jl` as a Python module

    >>> r = bs.first_modes(4,4)
    >>> r
    <PyCall.jlwrap state = [1, 1, 1, 1]>

A key strength of `BosonSampling.jl` for efficient computing and modularity is its type architecture. As we have seen previously in the tutorial in [`Basic usage`](https://benoitseron.github.io/BosonSampling.jl/dev/tutorial/basic_usage.html), the definition of an [`Input`](@ref) requires to play with types and therefore use the delimiters `{}`. Those symbols are not valid identifiers for Python which therefore prevents to define any `Input` or take advantage of any type hierarchy in Julia. The easiest workaround is to define a function, e.g

    >>> curly = Main.eval("(a,b...) -> a{b...}")

This function defined once for all, one can now define [`Bosonic`](@ref) [`Input`](@ref) with the [`ModeOccupation`](@ref) defined here above

    >>> i = curly(bs.Input, bs.Bosonic)(r)

to finally choose an interferometer and run a boson sampling simulation via the [`cliffords_sampler`](@ref)

    >>> U = bs.RandHaar(4)
    >>> bs.cliffords_sampler(input=i, interf=U)
