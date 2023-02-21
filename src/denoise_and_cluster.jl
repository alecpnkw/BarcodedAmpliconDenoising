# renamed method from 'RobustChoice'... 
function denoise_and_cluster(seqs; denoising_radius = 1.0, cluster_radius = 1e-3, cluster_proportion_thresh = 0.05)
    # RAD
    rad_cons,rad_sizes = denoise(seqs; fine_radius = denoising_radius, verbose=0);
    # RobustChoice
    out_seqs, out_sizes = String[], Int[]
    kmers = kmer_count.(rad_cons,6);
    _,_,clust_inds = dp_centers(kmers, cluster_radius,
        distfunc = corrected_kmer_dist,
        center = mean,
        cycle_lims = 30)
    for i in clust_inds
        best = argmax(rad_sizes[i])
        best_seq = rad_cons[i[best]]
        if rad_sizes[i[best]]/sum(rad_sizes) >= cluster_proportion_thresh
            push!(out_seqs,best_seq)
            push!(out_sizes,rad_sizes[i[best]])
        end
    end
    return out_seqs, out_sizes
end
