function assignterminals!(system::CVSystem,R::Vector{Float64},C::Vector{Float64},
    V::Vector{Float64},L::Float64,lowerflowfraction::Float64,venousfractionofR::Float64,
    venousfractionofL::Float64,venousfractionofC::Float64,venousfractionofV0::Float64,
    old=Dict("a"=>0),restart="no",assim="no",sample="no")
    R2Total = R[1];
    R3Total = R[2];
    R4Total = R[3];
    R5Total = R[4];
    C1Total = C[1];
    C2Total = C[2];
    C3Total = C[3];
    C4Total = C[4];
    C5Total = C[5];
    L5Total = L;
    V01 = V[1];
    V02 = V[2];
    V03 = V[3];
    V04 = V[4];
    V05 = V[5];

    numterminals = 0;
    numlower = 0;
    numupper = 0;

    for i = 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            numterminals += 1
            if system.branches.group[i] == "lower"
                numlower += 1
            else
                numupper += 1
            end
        end
    end

    # divide total peripheral impedance according to desired flow divisions
    R2Lower = numterminals*R2Total;
    R2Upper = numterminals*R2Total;
    R3Lower = numterminals*R3Total;
    R3Upper = numterminals*R3Total;
    R4Lower = numterminals*R4Total;
    R4Upper = numterminals*R4Total;
    C1Lower = C1Total/numterminals;
    C1Upper = C1Total/numterminals;
    C2Lower = C2Total/numterminals;
    C2Upper = C2Total/numterminals;
    C3Lower = C3Total/numterminals;
    C3Upper = C3Total/numterminals;
    C4Lower = C4Total/numterminals;
    C4Upper = C4Total/numterminals;
    V01Lower = V01/numterminals;
    V01Upper = V01/numterminals;
    V02Lower = V02/numterminals;
    V02Upper = V02/numterminals;
    V03Lower = V03/numterminals;
    V03Upper = V03/numterminals;
    V04Lower = V04/numterminals;
    V04Upper = V04/numterminals;

    C5Lower = venousfractionofC*C5Total/numterminals;
    C5Upper = venousfractionofC*C5Total/numterminals;
    V05Lower = venousfractionofV0*V05/numterminals;
    V05Upper = venousfractionofV0*V05/numterminals;

    system.svc.C = (1-lowerflowfraction)*(1-venousfractionofC)*C5Total;
    system.ivc.C = lowerflowfraction*(1-venousfractionofC)*C5Total;
    system.svc.V0 = 0.5*(1-venousfractionofV0)*V05;
    system.ivc.V0 = 0.5*(1-venousfractionofV0)*V05;

    RL = R5Total/lowerflowfraction;
    RU = R5Total/(1-lowerflowfraction);

    system.ivc.R = RL*(venousfractionofR/(1-venousfractionofR) + 1)^-1;
    system.svc.R = RU*(venousfractionofR/(1-venousfractionofR) + 1)^-1;

    Rl = RL - system.ivc.R;
    Ru = RU - system.svc.R;

    R5Lower = numlower*Rl;
    R5Upper = numupper*Ru;

    LL = L5Total/lowerflowfraction;
    LU = L5Total/(1-lowerflowfraction);

    system.ivc.L = LL*(venousfractionofL/(1-venousfractionofL) + 1)^-1;
    system.svc.L = LU*(venousfractionofL/(1-venousfractionofL) + 1)^-1;

    Ll = LL - system.ivc.L;
    Lu = LU - system.svc.L;

    L5Lower = numlower*Ll;
    L5Upper = numupper*Lu;

    # construct and define terminals
    for i in 1:length(system.branches.ID)
        if isempty(system.branches.children[i])
            if assim == "no" && sample == "no"
                system.branches.term[i] = CVModule.ArterialTerminal();
            end
            group = system.branches.group[i];
            if group == "lower"
                # compliance, lower compartments
                if assim == "no" && sample == "no"
                    if restart == "no"
                        push!(system.branches.term[i].C,C1Lower);
                        push!(system.branches.term[i].C,C2Lower);
                        push!(system.branches.term[i].C,C3Lower);
                        push!(system.branches.term[i].C,C4Lower);
                        push!(system.branches.term[i].C,C5Lower);
                    elseif restart == "yes"
                        temp = old[i]["C"];
                        push!(system.branches.term[i].C,temp[1]);
                        push!(system.branches.term[i].C,temp[2]);
                        push!(system.branches.term[i].C,temp[3]);
                        push!(system.branches.term[i].C,temp[4]);
                        push!(system.branches.term[i].C,temp[5]);
                    end
                elseif assim == "yes" || sample == "yes"
                    system.branches.term[i].C[1] = C1Lower;
                    system.branches.term[i].C[2] = C2Lower;
                    system.branches.term[i].C[3] = C3Lower;
                    system.branches.term[i].C[4] = C4Lower;
                    system.branches.term[i].C[5] = C5Lower;
                end
                # resistance, lower compartments
                if assim == "no" && sample == "no"
                    if restart == "no"
                        push!(system.branches.term[i].R,system.solverparams.rho*
                            system.branches.c0[i][end]/system.branches.A0[i][end]);
                        push!(system.branches.term[i].R,R2Lower);
                        push!(system.branches.term[i].R,R3Lower);
                        push!(system.branches.term[i].R,R4Lower);
                        push!(system.branches.term[i].R,R5Lower);
                    elseif restart == "yes"
                        temp = old[i]["R"];
                        push!(system.branches.term[i].R,temp[1]);
                        push!(system.branches.term[i].R,temp[2]);
                        push!(system.branches.term[i].R,temp[3]);
                        push!(system.branches.term[i].R,temp[4]);
                        push!(system.branches.term[i].R,temp[5]);
                    end
                elseif assim == "yes" || sample == "yes"
                    system.branches.term[i].R[1] = system.solverparams.rho*
                        system.branches.c0[i][end]/system.branches.A0[i][end];
                    system.branches.term[i].R[2] = R2Lower;
                    system.branches.term[i].R[3] = R3Lower;
                    system.branches.term[i].R[4] = R4Lower;
                    system.branches.term[i].R[5] = R5Lower;
                end
                if assim == "no" && sample == "no"
                    # unstressed volume, lower compartments
                    if restart == "no"
                        push!(system.branches.term[i].V0,V01Lower);
                        push!(system.branches.term[i].V0,V02Lower);
                        push!(system.branches.term[i].V0,V03Lower);
                        push!(system.branches.term[i].V0,V04Lower);
                        push!(system.branches.term[i].V0,V05Lower);
                    elseif restart == "yes"
                        temp = old[i]["V0"];
                        push!(system.branches.term[i].V0,temp[1]);
                        push!(system.branches.term[i].V0,temp[2]);
                        push!(system.branches.term[i].V0,temp[3]);
                        push!(system.branches.term[i].V0,temp[4]);
                        push!(system.branches.term[i].V0,temp[5]);
                    end
                    # inertance, lower veins
                    if restart == "no"
                        append!(system.branches.term[i].L,[zeros(4),L5Lower;])
                    elseif restart == "yes"
                        temp = old[i]["L"];
                        append!(system.branches.term[i].L,temp)
                    end
                elseif assim == "yes" || sample == "yes"
                    system.branches.term[i].V0[1] = V01Lower;
                    system.branches.term[i].V0[2] = V02Lower;
                    system.branches.term[i].V0[3] = V03Lower;
                    system.branches.term[i].V0[4] = V04Lower;
                    system.branches.term[i].V0[5] = V05Lower;
                end
            else
                # compliance, upper compartments
                if assim == "no" && sample == "no"
                    if restart == "no"
                        push!(system.branches.term[i].C,C1Upper);
                        push!(system.branches.term[i].C,C2Upper);
                        push!(system.branches.term[i].C,C3Upper);
                        push!(system.branches.term[i].C,C4Upper);
                        push!(system.branches.term[i].C,C5Upper);
                    elseif restart == "yes"
                        temp = old[i]["C"];
                        push!(system.branches.term[i].C,temp[1]);
                        push!(system.branches.term[i].C,temp[2]);
                        push!(system.branches.term[i].C,temp[3]);
                        push!(system.branches.term[i].C,temp[4]);
                        push!(system.branches.term[i].C,temp[5]);
                    end
                elseif assim == "yes" || sample == "yes"
                    system.branches.term[i].C[1] = C1Upper;
                    system.branches.term[i].C[2] = C2Upper;
                    system.branches.term[i].C[3] = C3Upper;
                    system.branches.term[i].C[4] = C4Upper;
                    system.branches.term[i].C[5] = C5Upper;
                end
                # resistance, upper compartments
                if assim == "no" && sample == "no"
                    if restart == "no"
                        push!(system.branches.term[i].R,system.solverparams.rho*
                            system.branches.c0[i][end]/system.branches.A0[i][end]);
                        push!(system.branches.term[i].R,R2Upper);
                        push!(system.branches.term[i].R,R3Upper);
                        push!(system.branches.term[i].R,R4Upper);
                        push!(system.branches.term[i].R,R5Upper);
                    elseif restart == "yes"
                        temp = old[i]["R"];
                        push!(system.branches.term[i].R,temp[1]);
                        push!(system.branches.term[i].R,temp[2]);
                        push!(system.branches.term[i].R,temp[3]);
                        push!(system.branches.term[i].R,temp[4]);
                        push!(system.branches.term[i].R,temp[5]);
                    end
                elseif assim == "yes" || sample == "yes"
                    system.branches.term[i].R[1] = system.solverparams.rho*
                        system.branches.c0[i][end]/system.branches.A0[i][end];
                    system.branches.term[i].R[2] = R2Upper;
                    system.branches.term[i].R[3] = R3Upper;
                    system.branches.term[i].R[4] = R4Upper;
                    system.branches.term[i].R[5] = R5Upper;
                end
                if assim == "no" && sample == "no"
                    # unstressed volume, upper compartments
                    if restart == "no"
                        push!(system.branches.term[i].V0,V01Upper);
                        push!(system.branches.term[i].V0,V02Upper);
                        push!(system.branches.term[i].V0,V03Upper);
                        push!(system.branches.term[i].V0,V04Upper);
                        push!(system.branches.term[i].V0,V05Upper);
                    elseif restart == "yes"
                        temp = old[i]["V0"];
                        push!(system.branches.term[i].V0,temp[1]);
                        push!(system.branches.term[i].V0,temp[2]);
                        push!(system.branches.term[i].V0,temp[3]);
                        push!(system.branches.term[i].V0,temp[4]);
                        push!(system.branches.term[i].V0,temp[5]);
                    end
                    # inertance, upper veins
                    if restart == "no"
                        append!(system.branches.term[i].L,[zeros(4),L5Upper;])
                    elseif restart == "yes"
                        temp = old[i]["L"];
                        append!(system.branches.term[i].L,temp)
                    end
                elseif assim == "yes" || sample == "yes"
                    system.branches.term[i].V0[1] = V01Upper;
                    system.branches.term[i].V0[2] = V02Upper;
                    system.branches.term[i].V0[3] = V03Upper;
                    system.branches.term[i].V0[4] = V04Upper;
                    system.branches.term[i].V0[5] = V05Upper;
                end
            end
        else
            if assim == "no" && sample == "no"
                system.branches.term[i] = CVModule.ArterialTerminal();
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
                if !isempty(system.branches.term[i].V0)
                    system.branches.term[i].V0[1] += NaN;
                else
                    push!(system.branches.term[i].V0,NaN);
                end
                push!(system.branches.term[i].L,NaN);
            end
        end
    end

end
