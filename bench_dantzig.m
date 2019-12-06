clear; clc;
%N = [5, 10, 20, 30, 40, 50, 75, 100, 125, 150, 200];
N = [5];
% Matriz para guardar:
% n | t | costo | iteraciones
bench = zeros(length(N), 4);
i = 1;
for num = N
    TSP;
    bench(i, :) = [num, t, costopt, iters];
    i  = i + 1;
end

writematrix(bench, 'benchmark_dantzig.txt', 'Delimiter', ' ');