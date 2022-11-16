# Work in progress

# SpectralBoundaryIntegralMethod.jl: Spectral boundary integral method code

- This code was developed for part of [my dissertation](https://www.researchgate.net/publication/355033649_Simulations_of_Red_Blood_Cell_Flow_by_Boundary_Integral_Methods) to simulate red blood cell flow using boundary integral methods.
- This repository contains the code for the concepts and examples presented in Chapter 4 for analyzing the red blood cell motion and deformation in an unbounded domain.

## Citation

    @phdthesis{gurbuz2021Thesis,
    title={Simulations of Red Blood Cell Flow by Boundary Integral Methods},
    author={G\"urb\"uz, Ali},
    year={2021},
    school={State University of New York at Buffalo}
    }

### To (locally) reproduce this project, do the following:

1. Download this code base.
2. Open a Julia console and do:

```julia
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:

```julia
    using DrWatson
    @quickactivate "SpectralBoundaryIntegralMethod"
```

which auto-activates the project and enables local path handling from DrWatson.
