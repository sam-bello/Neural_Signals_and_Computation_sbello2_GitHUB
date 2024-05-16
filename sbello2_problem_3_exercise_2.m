function [] = sbello2_problem_3_exercise_2(data)
%% Part A
cond27 = data(27).A;
dt = 10e-3;

% plot heatmap
figure
imagesc(cond27')
colorbar
title("Condition 27")
xlabel("Time (10 ms)")
ylabel("Neuron")

%% Part B
% Estimate the A matrix using linear regression
dxdt = (cond27(2:end,:) - cond27(1:end - 1,:)) ./ dt;
A = cond27(1:end - 1,:) \ dxdt;

% Reconstruct data from A matrix
xt = zeros(size(cond27));
xt(2:end,:) = (inv(1 - A * dt) * (xt(1:end - 1,:))')';

% Plot error histogram
error = zeros(size(cond27,2),1);
for n = 1:size(cond27,2)
    error(n) = abs(norm((xt(:,n)) - norm(cond27(:,n))) / norm(cond27(:,n)));
end

figure()
histogram(error)
title('Histogram of Reconstruction Errors')
xlabel("Time (10 ms)")
ylabel("Error")

%% Part C
% Get data from all trials
data_all = [];
for i = 1:length(data)
    data_all = [data_all; data(i).A];
end

% Get the first 6 PCs
[coeff, score_all, ~, ~, ~, mu] = pca(data_all,'NumComponents',6);
score = reshape(score_all,[size(cond27,1),length(data),6]);

% get A matrix from cost function
pca_cost_func = @(A) cost_func(score,A); 
A = fminunc(pca_cost_func,zeros(6,6));

% reconstruct signal for condition 27 using A matrix and inverse PCA
x_pca = (squeeze(score(:,27,:)))';
x_pca(:,2:end) = dt * A * x_pca(:,1:end - 1) + x_pca(:,1:end - 1);
pca_recon = x_pca' * coeff' + repmat(mu,size(cond27,1),1);

% plot the original PCs vs reconstructed PCs
figure()
for i = 1:6
    subplot(2,3,i)
    hold on
    plot(score(:,27,i))
    plot(x_pca(i,:))
    hold off
    title(['PC ', int2str(i)])
    xlabel("PC 1")
    ylabel("PC 2")
end
sgtitle("Original PCs vs Reconstructed PCs")

% plot error histogram
error_pca = zeros(size(cond27,2),1);
for n = 1:size(cond27,2)
    error_pca(n) = abs(norm((pca_recon(:,n)) - norm(cond27(:,n))) / norm(cond27(:,n)));
end

figure()
histogram(error_pca)
title('Histogram of Reconstruction Errors using PCA')

%% Part D
% plot the first 2 dynamical principal dimensions
figure()
hold on
for i = 1:length(data)
    plot(score(:,i,1),score(:,i,2));
end
hold off
title("First 2 Dynamical Principal Dimensions")
xlabel("PC 1")
ylabel("PC 2")

%% Problem 3 Part E
% jPCA
jPCA_params.softenNorm = 5;
jPCA_params.suppressBWrosettes = true;
jPCA_params.suppressHistograms = true;

times = -50:10:150;
jPCA_params.numPCs = 6;
[Projection, Summary] = jPCA(data, times,jPCA_params);

phaseSpace(Projection, Summary);
printFigs(gcf,'.','-dpdf','Basic jPCA plot');

% eigenvalue spectra
A_jpca = Summary.jPCs;

A_jpca_eig = eig(A_jpca)
A_eig = eig(A)

end


%% cost function
function val = cost_func(xx,A)
    dt = 10e-3;
    sig = 0.1;
    val = 0;

    for i = 1:108
        x = squeeze(xx(:,i,:));
        val = val + sum(0.5*(x(2:end,:)- (x(1:end-1,:)*(A*dt) + x(1:end-1,:))) * ...
            ((sig^2*eye(6)*dt)^(-1)) * (x(2:end,:)-(x(1:end-1,:)*(A*dt) + x(1:end-1,:)))','all');
    end

end