# Validation

Let's now see how we can use the tools of this package to validate boson samplers. Suppose we are given experimental results, how do we know that the boson sampler works?

Multiple methods exist to achieve this. We will focus here on a more general scheme as developed by the authors, which encompasses the majority of the methods used in the literature.

To do this, we look at photon counting in partition of the output modes of the interferometer. Generating this distribution is efficient and described in an other section of this tutorial. Through this, we can for instance try to see if the experimental dataset is more coherent with the hypothesis of indistinguishable particles, versus that of distinguishable ones.

We use a bayesian test to do this. Let's first see that this bayesian method can be used directly on the output results of the experiment (and not the photon counting in partitions )
