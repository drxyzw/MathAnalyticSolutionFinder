using Plots
f(x)=sin(x)/x
plot(f, xlims=(-10, 10), ylims=(-2, 2))
fa(x, n, s) = n == 0 ? 1.0 : (1-x^s/n^s/pi^s) * fa(x, n-1, s)
plot!(x -> fa(x, 1000, 2))
plot!(x -> fa(x, 1, 3))
plot!(x -> fa(x, 1000, 3))
