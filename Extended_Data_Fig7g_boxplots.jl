using DelimitedFiles
using Glob
using StatsPlots
using Plots
using HypothesisTests

Genotypes = ["Esr"  "NPY" "VGlut2" "WT"]

src = "/Users/pierre/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Pierre_Shared/LHA-LHb-PFC/Nature_Neuroscience_Revisions/pupil/"
filelist = glob("Eye_closing_events_*.csv",src)

p = plot(layout=(2,2), size=(800,600))


for (f,file) in enumerate(filelist)
#    f = 1
 #   file = filelist[f]
    csv = readdlm(file,',')
    pupil = Vector{Float64}(undef, 1)
    id = Vector{Int64}(undef, 1)
    for b in 1:4
        append!(pupil,csv[2:end,b])
        append!(id,ones(length(csv[2:end,b])).*b)
    end
    popfirst!(pupil)
    popfirst!(id)
    # remove nans
    deleteat!(id, isnan.(pupil).==1)
    deleteat!(pupil, isnan.(pupil).==1)
    # Plot
    boxplot!(id, pupil , subplot=f, linewidth=2, xticks=(1:4,csv[1,:]), xlabel="block", outliers=false, legend=false, fillcolor=["#1d1d1b" "#2babe2" "#9d9d9c"],
                 title="Genotype "*Genotypes[f], ylabel="Eye closing events", ylims=(-10, 150)) 
    dotplot!(id, pupil, subplot=f, color=:white, markerstrokewidth = 1, legend=false)             
end
display(p)

# Statistical testing
csv = readdlm(filelist[1],',')
eye_b1 = csv[2:end,1] 
for (f,file) in enumerate(filelist[2:4])
    csv = readdlm(file,',')
    eye_b1 =hcat(eye_b1,[csv[2:end,1]; fill(NaN,5)])
end

csv = readdlm(filelist[1],',')
eye_b4 = csv[2:end,4] 
for (f,file) in enumerate(filelist[2:4])
    csv = readdlm(file,',')
    eye_b4 =hcat(eye_b4,[csv[2:end,4]; fill(NaN,5)])
end

#Genotypes = ["Esr"  "NPY" "VGlut2" "WT"]
# block 1 vs block 4
pvalue(SignedRankTest(float.(eye_b1[:,1]),float.(eye_b4[:,1])))
pvalue(SignedRankTest(float.(eye_b1[:,2]),float.(eye_b4[:,2])))
pvalue(SignedRankTest(float.(eye_b1[:,3]),float.(eye_b4[:,3])))
pvalue(SignedRankTest(float.(eye_b1[:,4]),float.(eye_b4[:,4])))

