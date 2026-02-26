# Copyright (c) 2025–2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

Br_ν_data = readdlm(joinpath(external_data_directory, "Br_neutrino.csv"), ',')
Br_pion_minus_data = readdlm(joinpath(external_data_directory, "Br_pion_minus.csv"), ',')
Br_kaon_minus_data = readdlm(joinpath(external_data_directory, "Br_kaon_minus.csv"), ',')
Br_kaon_zero_long_data = readdlm(joinpath(external_data_directory, "Br_kaon_zero_long.csv"), ',')

m_N_in_GeV_list = log10.(Br_ν_data[:, 1])
log10_Br_ν_list = log10.(Br_ν_data[:, 2])
Br_ν_interp_rslt = LinearInterpolation(log10_Br_ν_list, m_N_in_GeV_list; extrapolation=ExtrapolationType.Extension, cache_parameters=true)

log10_m_N_in_GeV_list = log10.(Br_pion_minus_data[:, 1])
log10_Br_pion_minus_list = log10.(Br_pion_minus_data[:, 2])
Br_pion_minus_interp_rslt = LinearInterpolation(log10_Br_pion_minus_list, log10_m_N_in_GeV_list; extrapolation=ExtrapolationType.Extension, cache_parameters=true)

log10_m_N_in_GeV_list = log10.(Br_kaon_minus_data[:, 1])
log10_Br_kaon_minus_list = log10.(Br_kaon_minus_data[:, 2])
Br_kaon_minus_interp_rslt = LinearInterpolation(log10_Br_kaon_minus_list, log10_m_N_in_GeV_list; extrapolation=ExtrapolationType.Extension, cache_parameters=true)

log10_m_N_in_GeV_list = log10.(Br_kaon_zero_long_data[:, 1])
log10_Br_kaon_zero_long_list = log10.(Br_kaon_zero_long_data[:, 2])
Br_kaon_zero_long_interp_rslt = LinearInterpolation(log10_Br_kaon_zero_long_list, log10_m_N_in_GeV_list; extrapolation=ExtrapolationType.Extension, cache_parameters=true)

Br_ν(m_N::EnergyUnit) = (exp10 ∘ Br_ν_interp_rslt ∘ log10 ∘ EUval)(GeV, m_N)
Br_pion_minus(m_N::EnergyUnit) = (exp10 ∘ Br_pion_minus_interp_rslt ∘ log10 ∘ EUval)(GeV, m_N)
Br_pion_plus(m_N::EnergyUnit) = Br_pion_minus(m_N)
Br_kaon_minus(m_N::EnergyUnit) = (exp10 ∘ Br_kaon_minus_interp_rslt ∘ log10 ∘ EUval)(GeV, m_N)
Br_kaon_zero_long(m_N::EnergyUnit) = (exp10 ∘ Br_kaon_zero_long_interp_rslt ∘ log10 ∘ EUval)(GeV, m_N)
