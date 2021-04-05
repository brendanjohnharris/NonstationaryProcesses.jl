p = 0.0:0.01:3.0

anim = @animate for ip ∈ p
    println(ip)
    S = simulate(waveDrivenHarmonicSim(X0 = [0.1, 0.0], parameter_profile_parameters = ((π,), (ip*π,), (π,), (π,))))
    #println(DN_HistogramMode_10(Array(timeseries(S, 1))))
    plot(S, downsample=3, size=(1000, 1000), title="Ω = $ip", border=:none, grid=:none)
end every 1
gif(anim, fps = 20)