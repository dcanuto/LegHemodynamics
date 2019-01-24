function newtondist!(yout::Vector{Float64},iters::Vector{Int64},Ws::Float64,W2est::Float64,
    Pn::Float64,C::Float64,R0::Float64,R1::Float64,beta::Float64,c0::Float64,A0::Float64,
    W1::Float64,rho::Float64,f::Function,J::Function,maxiter::Int16,epsJ::Float64,
    epsN::Float64,maxval::Float64,h::Float64)

    # initial guess
    x0 = zeros(1);
    xx = zeros(1);
    x0[1] = W2est;

    # non-dimensionalize
    x0[1] /= Ws;

    # setup for iterations
    xx .= x0;
    x = zeros(1);
    N = 1;

    # allocation
    JJ = zeros(1);
    D = zeros(1);
    d = zeros(1);
    fvec = zeros(1);
    g = zeros(1);

    # max step size for line searches
    stpmax = 100*max(sqrt(norm(xx)),1);
    # println("Maximum step size: $stpmax")

    # useful parameter combinations
    a1 = 2*rho*(1+R1*C/h);
    a2 = (R0+R1+R0*R1*C/h)*(2*rho/beta)^2;

    # Newton iterations
    while N <= maxiter
        # determine Jacobian, check invertibility
        # println("Distal Jacobian:")
        J(JJ,xx,C,R0,R1,Ws,W1,rho,beta,c0,a1,a2,h);
        maximum!(d,abs.(JJ))
        D[1] = 1/d[1];
        # println(JJ)
        # println(D)
        # println(D*JJ)
        # println(cond(D*JJ))
        # println(det(D*JJ))
        if abs.((D.*JJ)[1]) < epsJ
            println(JJ)
            println(D.*JJ)
            error("Distal Newton Jacobian is singular.");
        end
        # compute gradient of line search objective function
        # println("Objective function:")
        f(fvec,xx,D,C,R0,R1,Ws,W1,rho,beta,c0,a1,a2,Pn,h);
        # println(f(xx,system,n,state))
        # println(fvec)
        fold = 0.5*dot(fvec,fvec);
        g = [(D.*JJ)'*fvec];
        # compute newton step
        # println(inv(D*JJ))
        # println("Newton step:")
        s = -inv((D.*JJ)[1])*fvec;
        # line search to update state vector
        # println("Line search:")
        fn, xn, check = LegHemodynamics.linedist(xx,fold,g,s,stpmax,f,J,
            maxiter,Ws,W1,rho,beta,C,R0,R1,c0,a1,a2,Pn,h,JJ,d,fvec,D);
        # println(xn[1]*vs)
        # check if sufficiently close to root
        # println("Root checking:")
        maximum!(d,abs.(JJ))
        D[1] = 1/d[1];
        f(fvec,xn,D,C,R0,R1,Ws,W1,rho,beta,c0,a1,a2,Pn,h);
        # println(fvec)
        if norm(fvec) <= 100*epsN
            x .= xn;
            x[1] = x[1]*Ws;
            # println("Converged estimate for W2: $(x[1])")
            # println(inv(D*JJ))
            # println(fvec)
            iters[1] += N;
            break
        end
        if check
            test = 0.
            den = max(fn,0.5*length(xn));
            for i = 1:length(xn)
                temp = abs(g[i])*max(abs(xn[i]),1.)/den;
                if temp > test
                    test = temp;
                end
            end
            if test < 1e-6
                check = true;
                println("Warning: gradient of objective function close to zero.")
            else
                check = false;
            end
        end
        # check for divergence
        if norm(fvec,Inf) >= maxval
            iters[1] += N;
            println(xn)
            println(fvec)
            error("Newton iteration diverged.");
        end
        N+=1;
        xx = xn;
        if N == maxiter
            println(JJ)
            # println(D)
            # println(D*JJ)
            println(xn)
            println(fvec)
            println(norm(fvec))
            error("Newton iteration failed to converge.");
        end
    end

    # update based on converged solution
    # println("Update:")
    yout[1] = (2*rho/beta)^2*(0.125*(W1-x[1])+c0)^4;
    yout[2] = yout[1]*0.5*(W1+x[1]);
    yout[3] = beta*(sqrt(yout[1]) - sqrt(A0)) - R0*yout[2];
    yout[4] = yout[3]/R1;
end
