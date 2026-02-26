# Copyright (c) 2025–2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

g_star_data_Baumann = readdlm(joinpath(external_data_directory, "g_star_data_Baumann.csv"), ',')
log10_T_list = reverse!(log10.(g_star_data_Baumann[:, 1]))
g_star_list = reverse!(g_star_data_Baumann[:, 2])
g_star_interp_rslt = LinearInterpolation(g_star_list, log10_T_list; extrapolation=ExtrapolationType.Constant, cache_parameters=true)

N_ρ(T) = (g_star_interp_rslt ∘ log10 ∘ EUval)(GeV, T)

log10_T_in_GeV_list = -5:.01:3
log10_ρ_in_GeV_list = zero(log10_T_in_GeV_list)
for (ii, log10_T_in_GeV) ∈ enumerate(log10_T_in_GeV_list)
    T = (GeV ∘ exp10)(log10_T_in_GeV)
    ρ = (π^2 / 30) * N_ρ(T) * T^4
    log10_ρ_in_GeV_list[ii] = (log10 ∘ EUval)(GeV, ρ)
end
T_interp_rslt = LinearInterpolation(log10_T_in_GeV_list, log10_ρ_in_GeV_list; extrapolation=ExtrapolationType.Linear, cache_parameters=true)

function from_ρ_SM_to_T_SM(ρ_SM)
    ρ_SM_val = EUval(GeV, ρ_SM)
    ρ_SM_val < 0 && return EU(NaN)
    ρ_SM_val == 0 && return zero(EU)

    T_SM_val = (exp10 ∘ T_interp_rslt ∘ log10)(ρ_SM_val)

    return GeV(T_SM_val)
end
