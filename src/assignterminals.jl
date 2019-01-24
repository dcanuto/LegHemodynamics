function assignterminals!(system::LegSystem,R::Vector{Float64},C::Vector{Float64},
    old=Dict("a"=>0),restart="no",assim="no",sample="no")

    # construct and define terminals
    termctr = 0;
    for i in 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            termctr+=1;
            if assim == "no" && sample == "no"
                system.branches.term[i] = LegHemodynamics.ArterialTerminal();
            end
            # compliance, lower compartments
            if assim == "no" && sample == "no"
                if restart == "no"
                    push!(system.branches.term[i].C,C[termctr]);
                elseif restart == "yes"
                    temp = old[i]["C"];
                    push!(system.branches.term[i].C,temp[1]);
                end
            elseif assim == "yes" || sample == "yes"
                system.branches.term[i].C[1] = C[termctr];
            end
            # resistance, lower compartments
            if assim == "no" && sample == "no"
                if restart == "no"
                    push!(system.branches.term[i].R,system.solverparams.rho*
                        system.branches.c0[i][end]/system.branches.A0[i][end]);
                    push!(system.branches.term[i].R,R[termctr]-system.branches.term[i].R[1]);
                    println("Second resistance, artery $i: $(system.branches.term[i].R[end])")
                elseif restart == "yes"
                    temp = old[i]["R"];
                    push!(system.branches.term[i].R,temp[1]);
                    push!(system.branches.term[i].R,temp[2]);
                end
            elseif assim == "yes" || sample == "yes"
                system.branches.term[i].R[1] = system.solverparams.rho*
                    system.branches.c0[i][end]/system.branches.A0[i][end];
                system.branches.term[i].R[2] = R[termctr];
            end
        else
            if assim == "no" && sample == "no"
                system.branches.term[i] = LegHemodynamics.ArterialTerminal();
                if !isempty(system.branches.term[i].C)
                    system.branches.term[i].C[1] += NaN;
                else
                    push!(system.branches.term[i].C,NaN);
                end
                if !isempty(system.branches.term[i].R)
                    system.branches.term[i].R[1] += NaN;
                else
                    push!(system.branches.term[i].R,NaN);
                end
            end
        end
    end

end
