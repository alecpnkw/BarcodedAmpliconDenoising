# method to be one of :simple_consensus, :RAD, or :RobustChoice... 
function reconstruct_reads(seqs; method = :RobustChoice,
    fine_radius = 1.0, robust_cluster_radius = 0.01, RobustProportion = 0.05)

    out_seqs, out_names = String[], String[]

    # consensus
    if method == :simple_consensus
        global_consensus = consensus_seq(seqs);
        push!(out_seqs,global_consensus)
        push!(out_names,"CONSENSUS_length_$(length(global_consensus))_count_$(length(seqs))")
        return out_seqs, out_names
    elseif (method == :RAD) | (method == :RobustChoice)
        try
            # RAD
            rad_cons,rad_sizes = denoise(seqs; fine_radius = fine_radius, verbose=0);
            if method == :RAD
                for k in 1:length(rad_cons)
                    push!(out_seqs,rad_cons[k])
                    push!(out_names,"RAD_$(k)_length_$(length(rad_cons[k]))_count_$(rad_sizes[k])")
                end
            end
            # RobustChoice
            if method == :RobustChoice
                kmers = kmer_count.(rad_cons,6);
                _,_,clust_inds = dp_centers(kmers, robust_cluster_radius,
                    distfunc = corrected_kmer_dist,
                    center = mean,
                    cycle_lims = 30)
                for (i,inds) in enumerate(clust_inds)
                    best = argmax(rad_sizes[inds])
                    best_seq = rad_cons[inds[best]]
                    if rad_sizes[inds[best]]/sum(rad_sizes) >= RobustProportion
                        push!(out_seqs,best_seq)
                        push!(out_names,"RobustChoice_$(i)_length_$(length(best_seq))_count_$(rad_sizes[inds[best]])")
                    end
                end
            end
            return out_seqs, out_names
        catch
            @warn "Failed to reconstruct..."
            return String[], String[]
        end
    else
        @error "method not recognized!"
    end
end
