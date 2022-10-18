function plot_fluc_results(solution::PowerGridSolution, fluc_node_idxs, ω_indices)
    plt_active_power = plot(solution, fluc_node_idxs, :p, lw = 3, ylabel = L"P[p.u.]", xlabel = L"t[s]", legend = false)

    plt_frequency = plot(solution, ω_indices, :x_1, legend = false, ylabel = L"ω[rad / s]", xlabel = L"t[s]")

    plt_voltage = plot(solution, :, :v, legend = false, ylabel = L"V [p.u.]", xlabel = L"t[s]")

    hist_voltage = histogram(vcat(solution(solution.dqsol.t, :, :v)...), legend = false, xlabel = L"V [p.u.]", color = colorant"darkgoldenrod1", linecolor = :match)
    hist_frequency = histogram(vcat(solution(solution.dqsol.t, ω_indices, :x_1)...), legend = false, xlabel =  L"ω[rad / s]", color = colorant"turquoise3", linecolor = :match)

    return plt_active_power, plt_frequency, plt_voltage, hist_voltage, hist_frequency
end