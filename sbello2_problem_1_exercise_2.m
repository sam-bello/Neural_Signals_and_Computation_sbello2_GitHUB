function [] = sbello2_problem_1_exercise_2(N)
%% Part A
% generate poisson pmf
g = gausswin(N,5) .* (cos(2*pi*(0:(N-1))/10).');
x = 2 * rand(N,1);
vec = exp(dot(g,x));
k = 1:100;
pmf = (vec.^k * exp(-vec)) ./ factorial(k);

% sample pmf to get histogram
r = [];
for i = 1:100*N
    r1 = 1;
    r2 = 1;
    while r2 > pmf(r1)
        r1 = randi([1 max(k)]);
        r2 = max(pmf)*rand;
    end
    r = [r; r1];
end
figure
histogram(r)
title("Histogram of r|x")

%% Part B
% M = N
M = N;
% generate input matrix
x_mat = 2 * rand(N,M);

% generate response vector
r = [];
for i = 1:M
    vec = exp(dot(g,x_mat(:,i)));
    k = 1:100;
    pmf = (vec.^k * exp(-vec)) ./ factorial(k);

    % sample pmf to get response
    r1 = 1;
    r2 = 1;
    while r2 > pmf(r1)
        r1 = randi([1 max(k)]);
        r2 = max(pmf)*rand;
    end
    r = [r; r1];
end

% least squares estimate of g
g_est1 = x_mat'\r;

% rmse of the estimate of g vs g
rmse = sqrt(sum((g_est1 - g) .^ 2) / length(g));

%% Part C
% minimixation cost function
part_c_cost_func = @(xx) cost_func(x_mat,r,xx);

% estimate g
g_est2 = fminunc(part_c_cost_func,rand(size(g)));
rmse2 = sqrt(sum((g_est2 - g) .^ 2) / length(g));

% repeat B and C for M = 2 * N
M = 2*N;
% generate input matrix
x_mat = 2 * rand(N,M);

% generate response vector
r = [];
for i = 1:M
    vec = exp(dot(g,x_mat(:,i)));
    k = 1:100;
    pmf = (vec.^k * exp(-vec)) ./ factorial(k);

    % sample pmf to get response
    r1 = 1;
    r2 = 1;
    while r2 > pmf(r1)
        r1 = randi([1 max(k)]);
        r2 = max(pmf)*rand;
    end
    r = [r; r1];
end

% least squares estimate of g
g_est3 = x_mat'\r;

% estimate g using minimization
part_c_cost_func = @(xx) cost_func(x_mat,r,xx);
g_est4 = fminunc(part_c_cost_func,rand(size(g)));

% repeat B and C for M = N / 2
M = round(N/2);
% generate input matrix
x_mat = 2 * rand(N,M);

% generate response vector
r = [];
for i = 1:M
    vec = exp(dot(g,x_mat(:,i)));
    k = 1:100;
    pmf = (vec.^k * exp(-vec)) ./ factorial(k);

    % sample pmf to get response
    r1 = 1;
    r2 = 1;
    while r2 > pmf(r1)
        r1 = randi([1 max(k)]);
        r2 = max(pmf)*rand;
    end
    r = [r; r1];
end

% least squares estimate of g
g_est5 = x_mat'\r;

% estimate g using minimization
part_c_cost_func = @(xx) cost_func(x_mat,r,xx);
g_est6 = fminunc(part_c_cost_func,rand(size(g)));

% plot curves
figure()
ax1 = subplot(1,3,1);
hold on
plot(g)
plot(g_est1)
plot(g_est2)
hold off
title("Estimation of G using M = N")
legend(["Actual g","Estimate of g using least-squares", "Estimate of g using minimization"])

ax2 = subplot(1,3,2);
hold on
plot(g)
plot(g_est3)
plot(g_est4)
hold off
title("Estimation of G using M = 2N")
legend(["Actual g","Estimate of g using least-squares", "Estimate of g using minimization"])

ax3 = subplot(1,3,3);
hold on
plot(g)
plot(g_est5)
plot(g_est6)
hold off
title("Estimation of G using M = N/2")
legend(["Actual g","Estimate of g using least-squares", "Estimate of g using minimization"])
linkaxes([ax1,ax2,ax3],'y')
%% Part D
% M = N
M = N;
% generate input matrix
x_mat = 2 * rand(N,M);

% generate response vector
r = [];
for i = 1:M
    vec = exp(dot(g,x_mat(:,i)));
    k = 1:100;
    pmf = (vec.^k * exp(-vec)) ./ factorial(k);

    % sample pmf to get response
    r1 = 1;
    r2 = 1;
    while r2 > pmf(r1)
        r1 = randi([1 max(k)]);
        r2 = max(pmf)*rand;
    end
    r = [r; r1];
end

% minimixation cost function with MAP
sigma = 1;
part_d_cost_func = @(xx) cost_func_map(x_mat,r,xx,sigma);

% estimate g
g_est7 = fminunc(part_d_cost_func,rand(size(g)));
rmse3 = sqrt(sum((g_est7 - g) .^ 2) / length(g));

%% Part E
A = 0.01;
M = 5000;
% generate input matrix
x_mat = 2 * rand(N,M);
x_mat = A * x_mat;

% generate response vector
r = [];
for i = 1:M
    vec = exp(dot(g,x_mat(:,i)));
    k = 1:100;
    pmf = (vec.^k * exp(-vec)) ./ factorial(k);

    % sample pmf to get response
    r1 = 1;
    r2 = 1;
    while r2 > pmf(r1)
        r1 = randi([1 max(k)]);
        r2 = max(pmf)*rand;
    end
    r = [r; r1];
end

% minimixation cost function with MAP
sigma = 1;
part_d_cost_func = @(xx) cost_func_map(x_mat,r,xx,sigma);

% estimate g
g_est8 = fminunc(part_d_cost_func,rand(size(g)));
% g_est8 = x_mat'\r;
rmse4 = sqrt(sum((g_est8 - g) .^ 2) / length(g));


figure()
hold on
plot(g)
plot(g_est8)
hold off
title("Estimation of G using at low signal (A = 0.01) with M = 100")
legend(["Actual g","Estimate of g using MAP estimation"])
end

function val = cost_func(A,y,x)
    val = 0;
    for i = 1:size(A,2)
        val = val + dot(A(:,1),x) - y(i)*log(dot(A(:,1),x));
    end
    % val = norm(A'*x - y).^2;
end

function val = cost_func_map(A,y,x,sigma)
    val = 0;
    for i = 1:size(A,2)
        val = val + dot(A(:,1),x) - y(i)*log(dot(A(:,1),x)) - log((1/(sigma*sqrt(2*pi)))*exp(-0.5*(dot(x,x)/sigma^2)));
    end
    % val = norm(A'*x - y).^2;
end