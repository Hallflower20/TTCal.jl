function old_genvis_equivalent(meta, sources)
    model = zeros(JonesMatrix, TTCal.Nbase(meta), TTCal.Nfreq(meta))
    frame = TTCal.reference_frame(meta)
    phase_center = measure(frame, meta.phase_center, dir"ITRF")
    for source in sources
        azel = measure(frame, source.direction, dir"AZEL")
        az = longitude(azel)
        el =  latitude(azel)
        dir = measure(frame, source.direction, dir"ITRF")
        l = dir.x - phase_center.x
        m = dir.y - phase_center.y
        n = dir.z - phase_center.z
        for β = 1:TTCal.Nfreq(meta)
            ν = meta.channels[β]
            λ = TTCal.c / ν
            flux = source.spectrum(ν) |> TTCal.linear
            jones = meta.beam(ν, az, el)
            flux  = TTCal.congruence_transform(jones, flux)
            for α = 1:TTCal.Nbase(meta)
                antenna1 = meta.baselines[α].antenna1
                antenna2 = meta.baselines[α].antenna2
                r1 = meta.antennas[antenna1].position
                r2 = meta.antennas[antenna2].position
                u = (r1.x - r2.x) / λ
                v = (r1.y - r2.y) / λ
                w = (r1.z - r2.z) / λ
                model[α,β] += exp(2π*1im * (u*l + v*m + w*n)) * flux
            end
        end
    end
    model
end

@testset "genvis.jl" begin
    Nant  = 5
    Nfreq = 10
    name, ms = createms(Nant, Nfreq)
    meta = collect_metadata(ms, ConstantBeam())
    sources = readsources("sources.json")

    # test that genvis agrees with the above reference implementation
    visibilities = genvis(meta, sources)
    model = old_genvis_equivalent(meta, sources)
    @test vecnorm(visibilities.data - model) < sqrt(eps(Float64)) * vecnorm(visibilities.data)

    # a source at the phase center should have unity visibilities
    sources = [PointSource("Cristiano Ronaldo", meta.phase_center, PowerLaw(1,0,0,0,10e6,[0.0]))]
    visibilities = genvis(meta, sources)
    model = ones(JonesMatrix, TTCal.Nbase(meta), Nfreq)
    @test vecnorm(visibilities.data - model) < sqrt(eps(Float64)) * vecnorm(visibilities.data)

    # a point source and a Gaussian source with no width should have equal visibilities
    point = PointSource("point", Direction(dir"AZEL", 45degrees, 45degrees),
                                 PowerLaw(4, 3, 2, 1, 10e6, [1.0]))
    gaussian = GaussianSource("gaussian", Direction(dir"AZEL", 45degrees, 45degrees),
                                          PowerLaw(4, 3, 2, 1, 10e6, [1.0]), 0, 0, 0)
    point_visibilities = genvis(meta, point)
    gaussian_visibilities = genvis(meta, gaussian)
    @test vecnorm(point_visibilities.data - gaussian_visibilities.data) < eps(Float64) * vecnorm(point_visibilities.data)

    # a multi-component source should get the same visibilities as each of those sources together
    sources = [PointSource("point", Direction(dir"AZEL", 45degrees, 45degrees),
                                    PowerLaw(4, 3, 2, 1, 10e6, [1.0])),
               GaussianSource("gaussian", Direction(dir"AZEL", 45degrees, 45degrees),
                                          PowerLaw(4, 3, 2, 1, 10e6, [1.0]),
                                          deg2rad(500/3600), deg2rad(300/3600), deg2rad(50))]
    multi = MultiSource("multi", sources)
    vis1 = genvis(meta, sources)
    vis2 = genvis(meta, multi)
    @test vecnorm(vis1.data - vis2.data) < eps(Float64) * vecnorm(vis1.data)
end

