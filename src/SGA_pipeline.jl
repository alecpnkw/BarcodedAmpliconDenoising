# generate reconstructed sequences for a FASTQ file
function SGA_pipeline(file::String, fwd_primer::String, rev_primer::String; max_depth = 1000, length_filter = (0.8, 1.2))
    seqID = replace(basename(file), ".fastq" => "")
    reads, _, _ = read_fastq(file)

    # min depth...
    if length(reads) < 3
        println("Insufficient read depth for $(seqID)!")
        return (String[],String[])
    end

    # max depth... 
    if !isnothing(max_depth) & (length(reads) > max_depth)
        reads = reads[1:max_depth]
    end

    # trim primers... 
    reads = double_primer_trim.(reads, uppercase(fwd_primer), uppercase(rev_primer))

    ## length filter...
    if !isnothing(length_filter)
        ls = length.(reads)
        l = StatsBase.median(ls)
        keeps = (ls .>= length_filter[1]*l) .& (ls .<= length_filter[2]*l)
        reads = reads[keeps]
    end

    # reconstruct...
    try 
        denoised_seqs, denoised_sizes = denoise_and_cluster(
            reads; 
            denoising_radius = 1.0, # radius for RAD denoising
            cluster_radius = 1e-3, # radius for RAD variant clustering
            cluster_proportion_thresh = 0.05 # proportion threshold to keep
        )
        denoised_names = ["$(seqID)_$(i)_count_$(s)" for (i,s) in enumerate(denoised_sizes)]
        return denoised_seqs, denoised_names
    catch
        @warn "$(seqID) failed to reconstruct!"
        return (String[],String[])
    end
end
