begin


# TOV diagram
# Esse programa resolve as equações de TOV e produz um diagrama M x Raio

using DelimitedFiles

G1 = 6.67430e-11                         # Constante Gravitacional em SI
c = 299792458.                            # Velocidade da Luz em SI m/2                    # Pi = 3.1415926...
nutosi = 1.60218e+32                 # MeV/fm^3 (Natural Units) -> J/m^3 (International Sysyem)
nutogu = 1.60218e+32*G1/(c^4)      # MeV/fm^3 (Natural Units) -> 1/m^2 (Geometrized Units)
gutosi = c^4/G1

# Variação da Massa Gravitacional em função do raio e da densidade de energia

function dmdr(r::Float64,e::Float64)::Float64
    return 4*π*(r^2)*e
end

# Variação da Pressão como função do raio e da densidade de energia

function dpdr(xp::Float64,xr::Float64,xm::Float64,xe::Float64)::Float64  # Variação da pressão por r
	if xr == 0.
		return 0.
	else
    a = (xp+xe)*(xm+4*π*(xr^3)*xp)/(xr*(xr-2*xm))
    return -a
	end
end

# interpolação linear de um ponto x0 de uma função definida com os vetores (x,y)

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

# Passo Runge-Kutta 4

function stepRK4TOV(p::Float64,m::Float64,r::Float64,dr::Float64,enp::Function)::Array{Float64}
    h = dr
    e = enp(p)

    k1p = dpdr(p,r,m,e)
    k1m = dmdr(r,e)

	#println(k1p,k1m)
    e = enp(p+h*k1p*0.5)

    k2p = dpdr(p+h*k1p*0.5,r+h*0.5,m+h*k1m*0.5,e)
    k2m = dmdr(r+h*0.5,e)

#println(k2p,k2m)

    e = enp(p+h*k2p*0.5)
    k3p = dpdr(p+h*k2p*0.5,r+h*0.5,m+h*k2m*0.5,e)
    k3m = dmdr(r+h*0.5,e)

#println(k3p,k3m)

    e = enp(p+h*k3p)
    k4p = dpdr(p+h*k3p,r+h,m+h*k3m,e)
    k4m = dmdr(r+h,e)

#println(k4p,k4m)

    dp = h*(k1p+2*k2p+2*k3p+k4p)/6
    dm = h*(k1m+2*k2m+2*k3m+k4m)/6

    return [m + dm ,p+dp,r+h]
end

function step_EULER(p::Float64,m::Float64,r::Float64,dr::Float64,enp::Function)::Array{Float64}

		h = dr
		e = enp(p)

		dm = dmdr(r,e)*h
		dp = dpdr(p,r,m,e)*h

		return [m+dm,p+dp,r+h]

	end

function TOV_solver(vP::Array{Float64},ve::Array{Float64},dr::Float64,P0::Float64)::Array{Float64}
    xP = P0
    xM = 0.0
    r = 0.0
	enp = x->lerp_1d(vP,ve,x)

	while xP ≥ 1e-12
        xM,xP,r = stepRK4TOV(xP,xM,r,dr,enp)
		#println(r," ",xP," ",xM)
    end

    return [xM*c^2/(G1*1.989e30),r]
end

	function MR_diagram(P::Array{Float64},e::Array{Float64},N::Int64)

		M = Array{Float64}(undef,N)
		R = Array{Float64}(undef,N)

		PP = LinRange(nutogu,P[length(P)]*nutogu,N)

		for i in 1:N

			M[i],R[i] = TOV_solver(P*nutogu,e*nutogu,1.0,PP[i])

		end

		#scatter(RR,MM)

		return R/1000,M

	end

	function File_eos(fname::String)

		ep = readdlm(fname,',',Float64);

		e = ep[:,1]
		P = ep[:,2]

		return e,P

	end


	function File_TOV(fname::String,N::Int)

		e,P = File_eos(fname)

		return MR_diagram(P,e,N)
	end

end

begin
	using Plots
	using DataFrames
	using CSV
end
begin 
begin
	data = File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data2.csv",1000)
	df = DataFrame((M=data[1],R=data[2]))
	CSV.write("MR_diag_TRS_14",df)
end

begin
	data = File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data1.csv",1000)
	df = DataFrame((M=data[1],R=data[2]))
	CSV.write("MR_diag_TRS_15",df)
end

begin
	data = File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data_MIT.csv",1000)
	df = DataFrame((M=data[1],R=data[2]))
	CSV.write("MR_diag_MIT",df)
end

begin
	data = File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data_mss13.csv",1000)
	df = DataFrame((M=data[1],R=data[2]))
	CSV.write("MR_diag_MSS13",df)
end

begin
	data = File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data_mss14.csv",1000)
	df = DataFrame((M=data[1],R=data[2]))
	CSV.write("MR_diag_MSS14",df)
end

begin
	data = File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data_mss15.csv",1000)
	df = DataFrame((M=data[1],R=data[2]))
	CSV.write("MR_diag_MSS15",df)
end
end 

begin 
	plot(File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data3.csv",1000),label="TRS 1.3")
	plot!(File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data2.csv",1000),label="TRS 1.4")
	plot!(File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data1.csv",1000),label="TRS 1.5")
	plot!(File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data_MIT.csv",1000),label="MIT")
	plot!(File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data_mss13.csv",1000),label="MSS 1.3")
	plot!(File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data_mss14.csv",1000),label="MSS 1.4")
	plot!(File_TOV("/home/arthuro/tov.jl/data/nov_07_24/data_mss15.csv",1000),label="MSS 1.5")
end