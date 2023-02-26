# BarcodedAmpliconDenoising

## Description

A Julia-based toolkit for consensus sequence generation from PacBio amplicon sequence data, compatible with flexible dual-barcode based multiplexing schemes. 

The core consensus generation method, `denoise_and_cluster()`, is based on the [RobustAmpliconDenoising.jl](https://github.com/MurrellGroup/RobustAmpliconDenoising.jl) package, but additionally clusters denoised sequence variants to collapse sequences likely generated via PCR errors during the amplification process. 

For ease of setting up new pipelines, we define a `SGA_pipeline()` function that preprocesses sequences, and calls this method. A standalone .ipynb running through the full analysis is included, along with several vizualization utilities at [BarcodedAmpliconDenoising.ipynb](notebooks/BarcodedAmpliconDenoising.ipynb).

## Setup

Activate the included environment from the Julia package manager. 

```julia
pkg> activate .
pkg> instantiate
pkg> precompile
```

Note: Installation depends on NextGenSeqUtils rev 1.5.3, which can be installed manually using

```juila
using Pkg
kg.add(PackageSpec(name="NextGenSeqUtils", rev="1.5.3", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))
```