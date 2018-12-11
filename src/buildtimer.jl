type CVTimer
    tfd::Float64
    tc::Float64
    ts::Float64
    tr::Float64
    tt::Float64
    function CVTimer()
        this = new();
        tfd = 0;
        tc = 0;
        ts = 0;
        tr = 0;
        tt = 0;
        return this
    end
end
