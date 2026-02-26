# Copyright (c) 2025–2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

integral_data_massive_FD = load(joinpath(integral_data_directory, "integrals_massive_FD.jld2"))

log_energy_density_vs_log_Tₘ_interpolation = LinearInterpolation(
    log.(integral_data_massive_FD["energy_density_FD_list"]),
    log.(integral_data_massive_FD["T_over_m_list"]);
    extrapolation=ExtrapolationType.Constant, cache_parameters=true
)
log_Tₘ_vs_log_energy_density_interpolation = LinearInterpolation(
    log.(integral_data_massive_FD["T_over_m_list"]),
    log.(integral_data_massive_FD["energy_density_FD_list"]);
    extrapolation=ExtrapolationType.Constant, cache_parameters=true
)
log_number_density_vs_log_Tₘ_interpolation = LinearInterpolation(
    log.(integral_data_massive_FD["number_density_FD_list"]),
    log.(integral_data_massive_FD["T_over_m_list"]);
    extrapolation=ExtrapolationType.Constant, cache_parameters=true
)
log_pressure_vs_log_Tₘ_interpolation = LinearInterpolation(
    log.(integral_data_massive_FD["pressure_FD_list"]),
    log.(integral_data_massive_FD["T_over_m_list"]);
    extrapolation=ExtrapolationType.Constant, cache_parameters=true
)

numeric_energy_density_FD(Tₘ) =
    if Tₘ < minimum(integral_data_massive_FD["T_over_m_list"])
        return 0
    elseif Tₘ > maximum(integral_data_massive_FD["T_over_m_list"])
        return (π^2 / 30) * (7 / 8) * Tₘ^4
    else
        (exp ∘ log_energy_density_vs_log_Tₘ_interpolation ∘ log)(Tₘ)
    end
numeric_Tₘ_FD(ρₘ) =
    if ρₘ < minimum(integral_data_massive_FD["energy_density_FD_list"])
        return 0
    elseif ρₘ > maximum(integral_data_massive_FD["energy_density_FD_list"])
        return (sqrt ∘ sqrt)(ρₘ / ((π^2 / 30) * (7 / 8)))
    else
        (exp ∘ log_Tₘ_vs_log_energy_density_interpolation ∘ log)(ρₘ)
    end

numeric_number_density_FD(Tₘ) =
    if Tₘ < minimum(integral_data_massive_FD["T_over_m_list"])
        return 0
    elseif Tₘ > maximum(integral_data_massive_FD["T_over_m_list"])
        return (ζ₃ / π^2) * (7 / 8) * Tₘ^3
    else
        (exp ∘ log_number_density_vs_log_Tₘ_interpolation ∘ log)(Tₘ)
    end

numeric_pressure_FD(Tₘ) =
    if Tₘ < minimum(integral_data_massive_FD["T_over_m_list"])
        return 0
    elseif Tₘ > maximum(integral_data_massive_FD["T_over_m_list"])
        return (π^2 / 90) * (7 / 8) * Tₘ^4
    else
        (exp ∘ log_pressure_vs_log_Tₘ_interpolation ∘ log)(Tₘ)
    end
