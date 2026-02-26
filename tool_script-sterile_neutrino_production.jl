# Copyright (c) 2025–2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

c_ν_e(T) = π * if zero(EU) ≤ T ≤ m_e
    7 / 27
elseif T ≤ m_μ
    91 / 216
elseif T ≤ Λ_QCD
    77 / 108
elseif T ≤ m_c
    203 / 108
elseif T ≤ m_τ
    244 / 81
elseif T ≤ m_b
    1981 / 648
elseif T ≤ M_W
    259 / 81
elseif T ≥ M_W
    259 / 81
else
    NaN
end

T_crit(m₄) = (m₄^2 / (2 * G_F^2))^(1/6)
function T_crit(m₄, c_ν::Function; relative_error=1e-6)
    this_T_crit = simple_T_crit = T_crit(m₄)
    new_T_crit = this_T_crit / c_ν(simple_T_crit)^(1/6)

    while abs((new_T_crit - this_T_crit) / this_T_crit) > relative_error
        this_T_crit = new_T_crit
        new_T_crit = simple_T_crit / c_ν(this_T_crit)^(1/6)
    end

    return new_T_crit
end

T_fo(m₄, U_a4_norm, c_ν::Function) = T_crit(m₄, c_ν) / cbrt(
    (1/4) * G_F * M_Pl * m₄ * U_a4_norm^2
)

BKG_coeff_peak(m₄, U_a4_norm) = G_F * M_Pl * m₄ * U_a4_norm^2 / 4
function BKG_coeff(T, m₄, U_a4_norm, c_ν::Function)
    this_T_crit = T_crit(m₄, c_ν)
    return BKG_coeff_peak(m₄, U_a4_norm) * (T ≥ this_T_crit ? (this_T_crit / T)^9 : (T / this_T_crit)^3)
end
