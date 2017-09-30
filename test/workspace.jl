workspace()

f1(x) = 1 / 2 * x
f2(x) = 1 / 2 * x + 1 / 2
funcs = [f1, f2]


@everywhere function worker(func, x)
    @parallel (hcat) for i = 1 : size(x)[2]
        func(x[:, i])
    end
end
# mapslices(funcs[k], x, 1)
initial = [1, 2]
max_iter = 2
num_funcs = length(funcs)
x = initial
for i = 1 : max_iter
    processes = [@spawn  worker(funcs[k], x) for k = 1 : num_funcs]
    x = reduce(hcat, fetch(processes[k]) for k = 1 : num_funcs)
end
