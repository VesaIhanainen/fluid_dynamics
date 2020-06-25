# using Pkg
# Pkg.add(DifferentialEquations)
# Pkg.Add(Plots)


using DifferentialEquations
using Plots

A = 2
N = 60
g = 9.81
gamma_0 = 1038

function problem(x,p,t)
    X = x[:,1]
    Y = x[:,2]
    G = g/gamma_0
    T = t/(2*π*A^2)/gamma_0
    x_prime = zeros(N)
    y_prime = zeros(N)
    for i in 1:N
        for j in 1:N
            if(i==j)
                continue
            end
            G_j = sqrt(1-(X[j]+1/N)^2)-sqrt(1-(X[j]-1/N)^2)
            x_prime[i] += G_j*(Y[i]-Y[j])/((X[i]-X[j])^2+(Y[i]-Y[j])^2)
            y_prime[i] -= G_j*(X[i]-X[j])/((X[i]-X[j])^2+(Y[i]-Y[j])^2)
        end
    end
    return_arr = zeros(N,N)
    return_arr[:,1] = x_prime
    return_arr[:,2] = y_prime
    return return_arr
end
function main()
   # Initial conditions
   u0 = zeros(N,N)
   u0[:,1] = range(-1,1,length=N)/2
   # Time span for solving
   tspan = (0.0,0.35)
   prob = ODEProblem(problem,u0,tspan)
   sol = solve(prob, DP5(),saveat=0.0001,reltol=1e-8,abstol=1e-8)
   print(length(sol[1,1,:]))
   GR.inline("png")
   anim = @animate for i=1:length(sol[1,1,:])
      plot(sol[:,1,i], sol[:,2,i], linestyle=:dash,marker=:dot)
      xlims!((0,0.6))
      ylims!((-0.6,0.1))
   end every 25
   gif(anim,"testi.gif",fps=30)
end

main()
