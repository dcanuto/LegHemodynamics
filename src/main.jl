importall LegModule

function main()

rstflag = "no" # restarting from scratch or previous simulation
hemoflag = "no" # 10% hemorrhage from left femoral artery
saveflag = "no" # save solution to .mat file
coupleflag = "no" # coupling to 3D liver tissue model
timeflag = "no" # solver timing
assimflag = "no" # patient data assimilation via EnKF

# build solution struct or generate ensemble
if assimflag == "no"
    if rstflag == "no"
        loadfile = "arterylist.txt"; # default artery data file for new sim
    elseif rstflag == "yes"
        loadfile = "lhoat2_$colnum.mat"; # restart file
    end
    system = CVModule.buildall(loadfile;numbeatstotal=10,restart=rstflag);
    savefile = "lhoat2_$colnum.mat" # filename for saving (only used if saveflag == "yes")
elseif assimflag == "yes"
    ensemblesize = 3;
    if rstflag == "no"
        loadfiles = ["arterytree.csv" for i=1:ensemblesize];
    elseif rstflag == "yes" loadfiles = ["test_1_$i.mat" for i=1:ensemblesize];
    end
    systems = pmap((a1)->CVModule.buildall(a1;numbeatstotal=1,restart=rstflag),loadfiles);
    savefiles = ["test_1_$i.mat" for i=1:ensemblesize];
end

if saveflag == "yes"
    if assimflag == "no"
        file = MAT.matopen(savefile, "w")
        write(file,"system",system)
        close(file)
    elseif assimflag == "yes"
        for i=1:ensemblesize
            file = MAT.matopen(savefiles[i],"w");
            write(file,"system",systems[i])
            close(file)
        end
    end
end

if assimflag == "no"
    if timeflag == "yes"
        return system, n, times
    elseif timeflag == "no"
        # return system, n
        return system
    end
elseif assimflag == "yes"
    return n,times,term_itr,split_itr
end

end
