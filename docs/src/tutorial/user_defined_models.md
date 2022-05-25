# Defining new models

One of the strengths of this package is the ease with which users may implement new models, for instance, new types of detectors (noisy, gaussian, threshold,...) without modifying the overall package structure: all calls will be similarly written, allowing for a very easy change of model with no extra coding (apart from the model itself).

For instance, one could implement a new type of [`OutputMeasurement`](@ref) that would be similar to the [`FockDetection`](@ref) but would account for random dark-counts in the detectors, or detector (in)efficiency.

Moreover, this is done without loss of performance given Julia's fast execution times. This would be much harder to do efficiently in a package linking a fast but lengthy to code programming language (C, Fortran,...) with a user interface language linking the fast routines that is easy to write but slow to execute (Python,...).

Note however that you can also use the above trick with our package if you wish to use our fast functions with a Python interface.

It is thus, in the authors' opinion, a good choice to use Julia for experimentalists who may want to account for a lot of subtleties not included in this package (or simply proper to their own experiment) as well as for theorists who may be interested in implementing new theoretical models, say nonlinear boson sampling.
