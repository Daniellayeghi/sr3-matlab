%% Value learning with Mujoco results
% Use csv data from the initial and 
clear all;
addpath("../src/");
path     = "~/Repos/OptimisationBasedControl/data/";
x_files  = sort_by_date(dir(path + "*state_00.csv"));
v_files  = sort_by_date(dir(path + "*value_00.csv"));
sv_files = sort_by_date(dir(path + "*data_size*.csv"));
nstates  = 2;
nvalues  = 1;
ngoals   = 19;

%% Define library ensemble
basis_fs  = {@(x)(x); @(x)(x.^2); @(x)(x(:,1:1) .* x(:, 1+1:end)); @(x)ones(length(x), 1)};
candid_fs = [];

%% Computes weights with SR3
lam1 = 0.01; % good for l_1 regularizer
lam0 = 0.04; % good for l_0 regularizer

w_guess = [0, 0, 1.73, 1.73, 2, 0]';
for goal = 1:length(x_files)
    x = csvread(path + x_files(goal).name); x(1, :) = [];
    v = csvread(path + v_files(goal).name); v(1, :) = [];
    v = v/max(v); [v, idx] = sortrows(v, 1); x = x(idx, :);
    figure(1);
    hold on;
    plot3(x(:, 1), x(:, 2), v(:));
    A = build_basis_lib(x, basis_fs);
    [x0, w0] = sr3(A, v, 'mode', '1', 'lam',lam0,'ptf',0, 'w0', w_guess);
    w_guess = w0;
    candid_fs = [candid_fs; w0'];
end

%% Function Defs
function A = build_basis_lib(x, basis_f)
    A = [];
    [r, c] = size(basis_f);

    for i = 1:r
        f = basis_f(i); f = f{1};
        A = [A, f(x)];
    end
end

function s = sort_by_date(s)
    s(~[s.isdir]);
    [~,idx] = sort([s.datenum]);
    s = s(idx);
end

% 
% for goal = 1:length(sv_files)
%     s = csvread(path + sv_files(goal).name);
%     [r, c] = size(s); ninit = c/(nstates+nvalues);
%     w_guess = [0, 0, 1.73, 1.73, 2, 0]';
%     for init = 1:ninit
%         i = (nstates+nvalues) * init - nstates;
%         x = s(:, i:i+nstates-1);
%         v = s(:, i+nstates);
%         A = build_basis_lib(x, basis_fs);
%         [x0, w0] = sr3(A, v, 'mode', '1', 'lam',lam0,'ptf',0, 'w0', w_guess);
%         w_guess = w0;
%     end
%     candid_fs = [candid_fs; w0'];
% end