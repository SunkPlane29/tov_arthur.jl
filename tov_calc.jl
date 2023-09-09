#
# This program will be used to 
# calculate the TOV and give possible other 
# quantities related to gravitational observables
#

# __author__: Arthur Efraim Bressan Pasqualotto
# __filiation__: PhD Student - UFSM - Santa Maria, RS, Brasil.
# __date__: 15:56 -=- 25/08/2023

# Objectives
# [ ] Solve the TOV equations from a table 
# [ ] Provide constraints for GW 
# TODO: Implementar a conversão G=c=1

# Importing the libraries
using DifferentialEquations

# Conditions to stop integration
# P = 0 
# m = M 
# φ = (1/2)*log(1-2GM/r)

# Get the data for the eos
# (firstly in MeV/fm³)
#-- Defining Constants --
begin
    G = 6.6708*1e-11                # Gravitational constant in SI
    c = 299792458.0                  # Speef of light in m/s
    nutosi = (1.60218*1e+32)          # MeV/fm^3 (Natural Units) -> J/m^3 (International Sysyem)
    nutogu = (1.60218*1e+32)*(G/(c^4)) # MeV/fm^3 (Natural Units) -> 1/m^2 (Geometrized Units)
    gutosi = (c^4.0)/G                    # 1/m^2 (Geometrized Units) -> J/m^3 (International System)
end

# Interpolate it to a function
# Defining the EOS

function isotropic(P,K,γ)
    # gravitational units
    return (P/K)^(1.0/γ)
end
# Defining the functions for dP, dM and dφ

function dφ(r,M,P)
    return (M + 4.0*π*r^3.0*P)/(r*(r-2.0*M))
end

function TOV!(du,u,p,t)
    P = u[1]
    M = u[2]
    ε = p(P)

    dphi = dφ(t,M,P)
    du[1] = -(ε + P)*dphi
    du[2] = 4π*t^2.0*ε
    du[3] = dphi
end

# Defining the stop condition 
function condition(u,t,integrator)
    # return u[3] - 0.5*log(1.0-2.0*u[2]/t)
    # return u[3] - 0.5*log(1.0-2*u[2]/t)
    return u[1]
end

affect!(integrator) = terminate!(integrator)

cb = ContinuousCallback(condition,affect!)

function solve_TOV(p0 = nutogu*20,tspan = (1.0e-15,5.0e4);EOS = P -> P/3.0)
    # Use it to solve the IVP (Initial Value Problem)
    # Initial Values:
    # P = P_0 
    # M = 0 
    # φ = 0 
    u0 = [p0,0.0,0.0]
    # tspan = (1.0e-10,5.0e4)
    p = EOS

    prob = ODEProblem(TOV!,u0,tspan,p,callback=cb)
    return solve(prob)
end

# This function returns the M x R diagram for 
# Compact objects 
function MR_diag(P=[1e-10,2e-10,3e-10,4e-10];EOS=eos)
    M = zeros(Float64,length(P))
    R = zeros(Float64,length(P))
    for i in 1:lastindex(P)
        sol = solve_TOV(P[i],EOS=eos)
        M[i] = sol[2,end]*(c^2/G)/1.989e30 # g.u. -> M_⊙
        R[i] = sol.t[end]
    end

    return M,R
end