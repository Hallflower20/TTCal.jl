function xyz2uvw(x,y,z)
    Nant = length(x)
    Nbase = div(Nant*(Nant-1),2) + Nant
    u = Array{Float64}(undef, Nbase)
    v = Array{Float64}(undef, Nbase)
    w = Array{Float64}(undef, Nbase)
    α = 1
    for i = 1:Nant, j = i:Nant
        u[α] = x[j]-x[i]
        v[α] = y[j]-y[i]
        w[α] = z[j]-z[i]
        α += 1
    end
    u,v,w
end

function ant1ant2(Nant)
    Nbase = div(Nant*(Nant-1),2) + Nant
    ant1 = Array{Int32}(undef, Nbase)
    ant2 = Array{Int32}(undef, Nbase)
    α = 1
    for i = 1:Nant, j = i:Nant
        ant1[α] = i - 1
        ant2[α] = j - 1
        α += 1
    end
    ant1, ant2
end

function createms(Nant, Nfreq)
    Nbase = div(Nant*(Nant-1), 2) + Nant

    frame = ReferenceFrame()
    OVRO  = Position(pos"WGS84", 1207.969*u"m", -118.284441*u"°", 37.232271*u"°")
    pos   = measure(frame, OVRO, pos"ITRF")
    t = (2015. - 1858.) * 365. * 24. * 60. * 60. * u"s" # a rough current Julian date
    set!(frame, Epoch(epoch"UTC", t))
    set!(frame, pos)

    x = pos.x .+ 100randn(Nant)
    y = pos.y .+ 100randn(Nant)
    z = pos.z .+ 100randn(Nant)
    u, v, w = xyz2uvw(x, y, z)
    ν = range(40e6, stop=60e6, length=Nfreq) |> collect
    ant1, ant2 = ant1ant2(Nant)

    zenith = Direction(dir"AZEL", 0*u"°", 90*u"°")
    phase_dir = measure(frame, zenith, dir"J2000")

    name  = tempname()*".ms"
    table = Tables.create(name)

    subtable = Tables.create("$name/SPECTRAL_WINDOW")
    Tables.add_rows!(subtable, 1)
    subtable["CHAN_FREQ"] = reshape(ν, length(ν), 1)
    table[kw"SPECTRAL_WINDOW"] = subtable
    Tables.close(subtable)

    subtable = Tables.create("$name/ANTENNA")
    Tables.add_rows!(subtable, Nant)
    subtable["POSITION"] = [x y z]'
    table[kw"ANTENNA"] = subtable
    Tables.close(subtable)

    subtable = Tables.create("$name/FIELD")
    Tables.add_rows!(subtable,1)
    subtable["PHASE_DIR"] = reshape([ustrip(longitude(phase_dir));
                                     ustrip( latitude(phase_dir))], 2, 1)
    table[kw"FIELD"] = subtable
    Tables.close(subtable)

    Tables.add_rows!(table, Nbase)
    table["ANTENNA1"] = ant1
    table["ANTENNA2"] = ant2
    table["UVW"] = [u v w]'
    table["TIME"] = fill(ustrip(t), Nbase)
    table["FLAG_ROW"] = zeros(Bool, Nbase)
    table["FLAG"] = zeros(Bool, 4, Nfreq, Nbase)

    name, table
end

