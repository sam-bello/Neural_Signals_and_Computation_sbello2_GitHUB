function [trace] = sbello2_problem_4_exercise_1(frames,mask,seed) % create a circle roi centered at the seed pixel 
n_frames = size(frames,4);
trace = zeros(n_frames,1);

% get trace as mean of each frame for the mask
for i = 1:n_frames
    frame = frames(:,:,:,i);
    trace(i) = mean(frame(mask),"all");
end

% plot trace
figure()
plot(1:n_frames,trace)
title(strcat("Time Trace of ROI with seed pixel: ",mat2str(seed)))
xlabel("Frame")
ylabel("Intensity")