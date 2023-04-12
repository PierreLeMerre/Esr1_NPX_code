using DelimitedFiles
using Glob
using StatsPlots
using Plots
using HypothesisTests

Genotypes = ["Esr"  "NPY" "VGlut2" "WT"]

# to adjust to run from the source data
src = "/Users/pierre/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Pierre_Shared/LHA-LHb-PFC/Nature_Neuroscience_Revisions/pupil/"
filelist = glob("Pupil_baseline_*.csv",src)

p = plot(layout=(2,2), size=(800,600))


for (f,file) in enumerate(filelist)
#    f = 1
 #   file = filelist[f]
    csv = readdlm(file,',')
    pupil = Vector{Float64}(undef, 1)
    id = Vector{Int64}(undef, 1)
    for b in 1:3
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
                 title="Genotype "*Genotypes[f], ylabel="Pupil area", ylims=(-1.5, 1.5)) 
                
end
display(p)

# Statistical testing
csv = readdlm(filelist[1],',')
pupil_b1 = csv[2:end,1] 
for (f,file) in enumerate(filelist[2:4])
    csv = readdlm(file,',')
    pupil_b1 =hcat(pupil_b1,csv[2:end,1])
end
# remove nans
[deleteat!(pupil_b1[:,i], isnan.(pupil_b1[:,i]).==1) for i in 1:4]

csv = readdlm(filelist[1],',')
pupil_b2 = csv[2:end,2] 
for (f,file) in enumerate(filelist[2:4])
    csv = readdlm(file,',')
    pupil_b2 =hcat(pupil_b2,csv[2:end,2])
end
# remove nans
[deleteat!(pupil_b2[:,i], isnan.(pupil_b2[:,i]).==1) for i in 1:4]

csv = readdlm(filelist[1],',')
pupil_b3 = csv[2:end,3] 
for (f,file) in enumerate(filelist[2:4])
    csv = readdlm(file,',')
    pupil_b3 =hcat(pupil_b3,csv[2:end,3])
end
# remove nans
[deleteat!(pupil_b3[:,i], isnan.(pupil_b3[:,i]).==1) for i in 1:4]

#Genotypes = ["Esr"  "NPY" "VGlut2" "WT"]
# block 1 vs block 2
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b1[:,1],isnan.(pupil_b1[:,1]))),float.(deleteat!(pupil_b2[:,1],isnan.(pupil_b2[:,1])))))
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b1[:,2],isnan.(pupil_b1[:,2]))),float.(deleteat!(pupil_b2[:,2],isnan.(pupil_b2[:,2])))))
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b1[:,3],isnan.(pupil_b1[:,3]))),float.(deleteat!(pupil_b2[:,3],isnan.(pupil_b2[:,3])))))
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b1[:,4],isnan.(pupil_b1[:,4]))),float.(deleteat!(pupil_b2[:,4],isnan.(pupil_b2[:,4])))))

# block 1 vs block 3
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b1[:,1],isnan.(pupil_b1[:,1]))),float.(deleteat!(pupil_b3[:,1],isnan.(pupil_b3[:,1])))))
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b1[:,2],isnan.(pupil_b1[:,2]))),float.(deleteat!(pupil_b3[:,2],isnan.(pupil_b3[:,2])))))
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b1[:,3],isnan.(pupil_b1[:,3]))),float.(deleteat!(pupil_b3[:,3],isnan.(pupil_b3[:,3])))))
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b1[:,4],isnan.(pupil_b1[:,4]))),float.(deleteat!(pupil_b3[:,4],isnan.(pupil_b3[:,4])))))

# block 2 vs block 3
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b2[:,1],isnan.(pupil_b2[:,1]))),float.(deleteat!(pupil_b3[:,1],isnan.(pupil_b3[:,1])))))
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b2[:,2],isnan.(pupil_b2[:,2]))),float.(deleteat!(pupil_b3[:,2],isnan.(pupil_b3[:,2])))))
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b2[:,3],isnan.(pupil_b2[:,3]))),float.(deleteat!(pupil_b3[:,3],isnan.(pupil_b3[:,3])))))
pvalue(MannWhitneyUTest(float.(deleteat!(pupil_b2[:,4],isnan.(pupil_b2[:,4]))),float.(deleteat!(pupil_b3[:,4],isnan.(pupil_b3[:,4])))))
