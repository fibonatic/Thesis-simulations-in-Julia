using DifferentialEquations, Plots, LaTeXStrings
include("functions_list.jl")

time_span = (0.0,60.0)
parameters = param(Diagonal([1.0,2.0,3.0]))
omega_0 = [1;1.2;-1] * 33

figs = solve_and_plot(time_span, parameters, omega_0)
