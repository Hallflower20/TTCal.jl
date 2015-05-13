# TTCal
![TTCal](ttcal.png)

[![Build Status](https://travis-ci.org/mweastwood/TTCal.jl.svg?branch=master)](https://travis-ci.org/mweastwood/TTCal.jl)
[![Coverage Status](https://coveralls.io/repos/mweastwood/TTCal.jl/badge.svg?branch=master)](https://coveralls.io/r/mweastwood/TTCal.jl?branch=master)

TTCal is a bandpass calibration routine developed for the OVRO LWA.

## Getting Started

TTCal requires the latest development version of the [Julia](http://julialang.org/) programming language.
```
$ git clone https://github.com/JuliaLang/julia.git
$ cd julia
$ make -j
```

To install TTCal, run:
```
$ julia
julia> Pkg.clone("https://github.com/mweastwood/CLI.jl.git")
julia> Pkg.clone("https://github.com/mweastwood/CasaCore.jl.git")
julia> Pkg.clone("https://github.com/mweastwood/TTCal.jl.git")
julia> Pkg.test("TTCal")
```
If all the tests pass, you are ready to begin using TTCal.
Simply add the `ttcal.jl` file to your `PATH` environment variable.
You can see the list of available commands by simpling running:

```
$ ttcal.jl 
 TTCal
=======
A calibration routine developed for the OVRO LWA.
Written by Michael Eastwood (mweastwood@astro.caltech.edu).

usage: ttcal.jl command options...

commands:
  bandpass        Solve for a bandpass calibration.
  polcal          Solve for a polarization calibration.
  applycal        Apply a calibration.
  diagnose        Diagnose a poor calibration.

Please provide one of the listed commands.
```

## Unpolarized Bandpass Calibration
## Polarized Bandpass Calibration
## Direction-Dependent Calibration
## Self-Calibration

## Sky Models

TTCal expects the sky model to be specified as a list of point sources in a JSON file.
Each source receives a name (`name`), J2000 position (`ra` and `dec`), and power-law spectrum.
The Stokes parameters (`I`, `Q`, `U`, and `V`) must be given at the reference frequency
(`freq`, measured in Hz). The spectral index is specified by `index`. Higher order terms
may also be specified. See the following example for Cas A and Cyg A using the published
spectra for Baars et al. 1977.

```
[
    {
        "ref": "Baars et al. 1977",
        "name": "Cas A",
        "ra": "23h23m24s",
        "dec": "58d48m54s",
        "I": 555904.26,
        "Q": 0.0,
        "U": 0.0,
        "V": 0.0,
        "freq": 1.0e6,
        "index": [-0.770]
    },
    {
        "ref": "Baars et al. 1977",
        "name": "Cyg A",
        "ra": "19h59m28.35663s",
        "dec": "+40d44m02.0970s",
        "I": 49545.02,
        "Q": 0.0,
        "U": 0.0,
        "V": 0.0,
        "freq": 1.0e6,
        "index": [+0.085,-0.178]
    }
]
```
