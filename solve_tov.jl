begin
    include("tov_calc.jl")
    using Plots
    using CSV,DataFrames
end 

function lerp_1d(x::Array{Float64},y::Array{Float64},x0::Float64)::Float64
	idx = 1

	if x0 < x[1]

	if x[1] ≤ 0.
				return 0.
			else
	return x0*y[1]/x[1]
			end
	else

    for i in 1:length(x)
        if x[i] > x0
            idx = i - 1
            break
        end
    end

    Δx = x[idx+1]-x[idx]
    Δy = y[idx+1]-y[idx]

    a = Δy/Δx
    b = y[idx]-a*x[idx]

    return a*x0+b
	end
end

begin
    data = CSV.read("data/nov_07_24/data.csv",DataFrame,header=false) 
end
data[:,1]
function P(ε::Float64)::Float64
    x = data[:,2]
    y = data[:,1]
    return lerp_1d(x,y,ε)    
end

P_0 = data[1,2]
eos = x->(x/nutogu-4*57.5)*nutogu/3.0 # MeV/fm^3
4*57.5*nutogu
eos(nutogu)

solve_TOV(nutogu,EOS=eos)
plot(MR_diag(LinRange(nutogu,10*nutogu,10)))

eos(100)
nutogu