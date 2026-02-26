# Copyright (c) 2025–2026 Quan-feng WU <wuquanfeng@ihep.ac.cn>
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using Integrals
using JLD2

include("tool_script-directories.jl")
include("tool_script-geomspace.jl")

f_FD(pₘ, Tₘ) = inv(exp(sqrt(pₘ^2 + 1) / Tₘ) - 1)

integrand_energy_density_FD(pₘ, Tₘ) = sqrt(pₘ^2 + 1) * pₘ^2 * f_FD(pₘ, Tₘ) / (2 * π^2)
integrand_number_density_FD(pₘ, Tₘ) = pₘ^2 * f_FD(pₘ, Tₘ) / (2 * π^2)
integrand_pressure_FD(pₘ, Tₘ) = (pₘ^2 / (3 * sqrt(pₘ^2 + 1))) * pₘ^2 * f_FD(pₘ, Tₘ) / (2 * π^2)

Tₘ_list = geomspace(1.5e-3, 1e15, 10000)

energy_density_FD_list = zero(Tₘ_list)
number_density_FD_list = zero(Tₘ_list)
pressure_FD_list = zero(Tₘ_list)

for i_Tₘ ∈ eachindex(Tₘ_list)
    Tₘ = Tₘ_list[i_Tₘ]

    problem_energy_density_FD = IntegralProblem(integrand_energy_density_FD, (0., Inf), Tₘ)
    problem_number_density_FD = IntegralProblem(integrand_number_density_FD, (0., Inf), Tₘ)
    problem_pressure_FD = IntegralProblem(integrand_pressure_FD, (0., Inf), Tₘ)

    solution_energy_density_FD = solve(problem_energy_density_FD, QuadGKJL())
    solution_number_density_FD = solve(problem_number_density_FD, QuadGKJL())
    solution_pressure_FD = solve(problem_pressure_FD, QuadGKJL())

    @assert solution_energy_density_FD.retcode == SciMLBase.ReturnCode.Success
    @assert solution_number_density_FD.retcode == SciMLBase.ReturnCode.Success
    @assert solution_pressure_FD.retcode == SciMLBase.ReturnCode.Success

    energy_density_FD_list[i_Tₘ] = solution_energy_density_FD.u
    number_density_FD_list[i_Tₘ] = solution_number_density_FD.u
    pressure_FD_list[i_Tₘ] = solution_pressure_FD.u
end

jldopen(joinpath(integral_data_directory, "integrals_massive_FD.jld2"), "w+") do jld2_file
    jld2_file["T_over_m_list"] = Tₘ_list
    jld2_file["energy_density_FD_list"] = energy_density_FD_list
    jld2_file["number_density_FD_list"] = number_density_FD_list
    jld2_file["pressure_FD_list"] = pressure_FD_list
end
