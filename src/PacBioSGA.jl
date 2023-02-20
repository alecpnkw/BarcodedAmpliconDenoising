module PacBioSGA

using NextGenSeqUtils: double_primer_trim, read_fastq, write_fasta, kmer_count, corrected_kmer_dist, FAD
using RobustAmpliconDenoising: consensus_seq, denoise
using DPMeansClustering: dp_centers
using MultivariateStats: transform, MDS
using StatsBase
using PyPlot

include("consensus.jl")
include("visualization.jl")

# consensus.jl...
export reconstruct_reads,

# visualization.jl...
digital_gel_plot,
kmer_dist_mds_plot

end
