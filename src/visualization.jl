
function digital_gel_plot(files; labels = nothing, alpha = 0.025, height = 6, grid_lines_every = 1000)
    inds = Int[]
    x_pos = []
    points = []
    label_text = []
    count = length(files)
    for (i,file) in enumerate(files)
        seqs,_,_ = read_fastq(file)
        l = length(seqs)
        push!(inds,i)
        append!(points,length.(seqs))
        append!(x_pos,[(i - 0.85) + rand()*0.7 for s in seqs])
        if !isnothing(labels)
            push!(label_text,labels[i]*" $l")
        end
    end
    fig, ax = subplots(figsize=(count, height))
    for l in grid_lines_every:grid_lines_every:Int64(round(maximum(points)))
        plot([0,count],[l,l],alpha = 0.05, color = "red")
    end
    plot(x_pos,points,alpha = alpha,".",color = "grey")
    if !isnothing(labels)
        for i in 1:count
            text((i - 0.5) / count, 1.05, label_text[i] , rotation=90, transform = ax.transAxes)
        end
    end
    xlim(0,count);
    xticks([]);
    tight_layout(w_pad=0,h_pad=0,pad=0);
    return fig, ax 
end

function pairwise_dist(obs::Array, dist_func)
    l = length(obs)
    mat = zeros(l,l)
    for i in 1:l
        for j in i+1:l
            d = dist_func(obs[i], obs[j])
            #Assumes the distance is symmetric
            mat[i,j] = d
            mat[j,i] = d
        end
    end
    return mat
end

function kmer_dist_mds_plot(seqs::Array{String,1}; kwargs = Pair{Symbol,Any}[:alpha => 0.4])
    kmer_vecs = kmer_count.(vcat(seqs), 6);
    pdists = pairwise_dist(kmer_vecs, corrected_kmer_dist);
    fig, ax = subplots()
    # multidimensional scaling using MultivariateStats
    proj = transform(fit(MDS, pdists, maxoutdim=2, distances=true));
    ax.scatter(proj[1,:], proj[2,:]; alpha = 0.4, kwargs...)
    ax.set_xlabel("MDS 1")
    ax.set_ylabel("MDS 2")
    return fig, ax
end