%% ODE for double integrator
addpath("../src/");
clear all;
close all;
npos = 1;
t = 0:0.01:10;
x_0 = [1; 0];
u = 0;
v = [];
xs = [];
sj_c = [];
global ref; ref = [0.5; 0];
running_cost = @(x, u, ref)(quad_cost(x(1), ref(1)) + quad_cost(x(2), ref(1)) + quad_cost(u, ref(1)));
for i = 1:100
    x_0 = [rand; 0];
    ref = [rand; 0];
    [t, x] = ode45(@point_mass, t, x_0);
    xs = [xs, x];
    u1 = arrayfun(@(xi) -(xi(1)), x(:, 1));
    u2 = arrayfun(@(xi) -(sqrt(3) * x(2)), x(:, 2));
    u = u1 + u2 + ref(1);
    j = compute_cst2go(x, ref(1));
    temp = [x, j];
    temp = remove_mid_elems(temp, 0);
    sj_c = [sj_c, temp];
    v = [v, j];
    figure(1);
    hold on;
    plot3(x(:, 1), x(:,2), j);
end


v_sparse = remove_mid_elems(v, 3);
[v, idx] = sortrows(v, 1);
xs = xs(idx, :);

basis_funcs = {@(x)(x); @(x)(x.^2); @(x)(x(:,1:1) .* x(:, 1+1:end))};
A = build_basis_lib(xs, basis_funcs);

%% Computes weights with SR3
lam1 = 0.01; % good for l_1 regularizer
lam0 = 0.004; % good for l_0 regularizer
[x0, w0] = sr3(A, v, 'mode', '1', 'lam',lam0,'ptf',0);
w = [0, 0, sqrt(3), sqrt(3), 2];
err = norm((A * w0 - v).^2);
csvwrite("temp_v.csv", v);

figure();
plot(A * w0);
hold on;
plot(v);
title("SR3 Results");

figure();
plot(A * w') 
hold on;
plot(v);
title("Ground Truth Results");

%% Function Defs
function A = build_basis_lib(x, basis_f)
    A = [];
    [r, c] = size(basis_f);

    for i = 1:r
        f = basis_f(i); f = f{1};
        A = [A, f(x)];
    end
end

function xd = point_mass(t, x)
    global ref;
    u = ctrl_cb(x, ref(1));
    xd = [x(2); u];
end

function u = ctrl_cb(x, ref)
    u = -(x(1) + sqrt(3) * x(2) - ref);
end

function cst = quad_cost(x, ref)
    cst = (ref-x)^2;
end

function v = compute_values(x, u, ref)
    quad_cost = @(x)(x^2);
    pos_c = arrayfun(@quad_cost, x(:, 1), ref);
    vel_c = arrayfun(@quad_cost, x(:, 2), ref);
    ctr_c = arrayfun(@quad_cost, u(:, 1), ref);
    
    v= [];
    
    for i = 1:length(x)
        cst = sum(pos_c(i:end, 1)) + sum(vel_c(i:end, 1));
        v = [v; cst];
    end
end

function cst2go = compute_cst2go(x, ref)
    J = @(x)sqrt(3)*(ref - x(1))^2  + 2*x(2)*(ref - x(1)) + sqrt(3)*x(2)^2;
    cst2go = [];
    for i = 1:length(x)
        xi = x(i, :);
        cst2go = [cst2go; J(xi)]; 
    end 
end

function arr = remove_mid_elems(arr, times)
    for sparse = 1:times
        arr(1:2:end-1, :) = []; 
    end
end