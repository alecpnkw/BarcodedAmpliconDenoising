# BarcodedAmpliconDenoising

## Description

A Julia-based toolkit for consensus sequence generation from PacBio amplicon sequence data, compatible with flexible dual-barcode based multiplexing schemes. 

The core consensus generation method, `denoise_and_cluster()`, is based on the [RobustAmpliconDenoising.jl](https://github.com/MurrellGroup/RobustAmpliconDenoising.jl) package, but additionally clusters denoised sequence variants to collapse sequences likely generated via PCR errors during the amplification process. 

For ease of setting up new pipelines, we define a `SGA_pipeline()` function that preprocesses sequences, and calls this method. A standalone .ipynb running through the full analysis is included, along with several vizualization utilities at [BarcodedAmpliconDenoising.ipynb](notebooks/BarcodedAmpliconDenoising.ipynb).

## Installation

Installation depends on NextGenSeqUtils rev 1.5.3, which is not yet registered on Julia 1.0. To install, please run the following in the Julia REPL. 

```juila
using Pkg
Pkg.add(PackageSpec(name="NextGenSeqUtils", rev="1.5.3", url = "https://github.com/MurrellGroup/NextGenSeqUtils.jl.git"))
Pkg.add(PackageSpec(name="DPMeansClustering", rev="1.0", url = "https://github.com/MurrellGroup/DPMeansClustering.jl.git"))
Pkg.add(PackageSpec(name="RobustAmpliconDenoising", rev="1.0", url = "https://github.com/MurrellGroup/RobustAmpliconDenoising.jl.git"))
```