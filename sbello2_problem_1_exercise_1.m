function [frames,corr14] = sbello2_problem_1_exercise_1(filename)
% get number of frames from file
fileinfo = imfinfo(filename);
n_frames = length(fileinfo);
% frames = zeros(fileinfo.Width,fileinfo.Height,1,n_frames);

% save each frame in the file in grayscale
for i = 1:n_frames
    frames(:,:,:,i) = mat2gray(imread(filename,i));
end

% play video using MATLAB Movie Player
implay(frames)

% compute correlation between frames 1 and 4
corr14 = normxcorr2(frames(:,:,:,1),frames(:,:,:,4));
padding = (length(frames(:,:,:,1)) - 1) / 2 + 1;
corr14 = corr14(padding:end - padding + 1,padding:end - padding + 1);