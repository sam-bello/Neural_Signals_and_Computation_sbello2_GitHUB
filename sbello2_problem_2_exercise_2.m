function [] = sbello2_problem_2_exercise_2(data)
%% Part A
% compute PSTH by summing the spike counts for each neuron and time bin
% across all trials
PSTH_raw = zeros(size(data(1).spikes));
for i = 1:length(data)
    PSTH_raw = PSTH_raw + data(i).spikes;
end

% smoothen PSTH using GP
% i could get the prior distribution and MAP estimation in problem 1
% working quite right, so i used a guassian filter with sigma = cov(A) to
% smoothen the PSTH
PSTH_smooth = zeros(size(PSTH_raw));
for i = 1:size(PSTH_raw,1)
    k = cov(PSTH_raw(i,:));
    gauss = gausswin(10,k);
    gauss_filt = gauss/sum(gauss);
    PSTH_smooth(i,:) = conv(PSTH_raw(i,:),gauss_filt,'same');

end

% plot PSTH
figure()
for i = 1:size(PSTH_raw,1)
    ax1 = subplot(1,2,1);
    histogram(PSTH_raw(i,:))
    title(['Raw PSTH for Neuron ' int2str(i)])
    xlabel('T (1 ms bins)')
    ylabel('Spikes')

    ax2 = subplot(1,2,2);
    histogram(PSTH_smooth(i,:))
    title(['Smoothened PSTH for Neuron ' int2str(i)])
    xlabel('T (1 ms bins)')
    ylabel('Spikes')
    linkaxes([ax1,ax2],'xy')
end

%% Part B
% pca of raw and smooth PSTH
[coeff_raw, score_raw] = pca(PSTH_raw,'Centered',true);
[coeff_smooth, score_smooth] = pca(PSTH_smooth,'Centered',true);


% plot first 3 principal components
figure()
ax1 = subplot(1,2,1);
plot3(score_raw(:,1),score_raw(:,2),score_raw(:,3))
title('Raw PSTH in Space of First 3 PCs')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')

ax2 = subplot(1,2,2);
plot3(score_smooth(:,1),score_smooth(:,2),score_smooth(:,3))
title('Smoothed PSTH in Space of First 3 PCs')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
linkaxes([ax1,ax2],'xyz')

%% Part C
method = 'gpfa';
runIdx = 1;
xDim = 8;
kernSD = 30;

% Extract neural trajectories
result = neuralTraj(runIdx, data, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD, 'binWidth', 1);

% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);

% Plot neural trajectories in 3D space
plot3D(seqTrain, 'xorth', 'dimsToPlot', 1:3);
end