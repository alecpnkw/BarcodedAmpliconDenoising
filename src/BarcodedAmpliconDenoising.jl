module BarcodedAmpliconDenoising

using NextGenSeqUtils: double_primer_trim, read_fastq, write_fasta, kmer_count, corrected_kmer_dist, FAD
using RobustAmpliconDenoising: consensus_seq, denoise
using DPMeansClustering: dp_centers
using MultivariateStats: transform, MDS
using StatsBase
using PyPlot

include("denoise_and_cluster.jl")
include("SGA_pipeline.jl")
include("visualization.jl")

# denoise_and_cluster.jl...
export denoise_and_cluster,

# SGA_pipeline.jl...
SGA_pipeline,

# visualization.jl...
digital_gel_plot,
kmer_dist_mds_plot

end
