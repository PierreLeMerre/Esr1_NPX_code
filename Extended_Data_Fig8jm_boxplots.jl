using DelimitedFiles
using Glob
using StatsPlots
using Plots
using HypothesisTests
using DataFrames

Genotype = "Esr"

col=["#1d1d1b" "#185c2a" "#F07F1A"]

# Need to change the input to the source data.....
src = "/Users/pielem/Library/CloudStorage/OneDrive-KI.SE/Pierre_Shared/LHA-LHb-PFC/Nature_Neuroscience_Revisions/tuning_score_"*Genotype*"NPY/"
filelist = glob("*.csv",src)

p = plot(layout=(4,2), size=(400,800))

for (f,file) in enumerate(filelist)
#    f = 1
 #   file = filelist[f]
    csv = readdlm(file,',')
    pupil = Vector{Float64}(undef, 1)
    id = Vector{Int64}(undef, 1)
    colors = Vector{String}(undef, 1)
    for b in 1:3
        append!(pupil,csv[2:end,b])
        append!(id,ones(length(csv[2:end,b])).*b)
        append!(colors,fill(col[b],length(csv[2:end,b])))
    end
    popfirst!(pupil)
    popfirst!(id)
    colors = colors[2:end]
    # remove nans
    deleteat!(colors, isnan.(pupil).==1)
    deleteat!(id, isnan.(pupil).==1)
    deleteat!(pupil, isnan.(pupil).==1)
    # Plot
    boxplot!(id, pupil , subplot=f, linewidth=2, xticks=(1:3,["All" "Ww" "Nw"]), xlabel="block", outliers=false, legend=false, fillcolor=["#1d1d1b" "#2babe2" "#9d9d9c"],
                 title="NPY "*Genotypes[f], ylabel="Tuning score", ylims=(-15, 15)) 
    #dotplot!(id, pupil, subplot=f, xticks=(1:3,["All" "Ww" "Nw"]), markersize=2,
           #color=colors, markerstrokewidth = 0, legend=false, title="Esr "*Genotypes[f], ylabel="Tuning score", ylims=(-15, 20))             
end
display(p)