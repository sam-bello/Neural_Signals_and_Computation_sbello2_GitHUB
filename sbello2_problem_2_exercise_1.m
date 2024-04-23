function [frames,mean_im,median_im,var_im,range_im,std_im] = sbello2_problem_2_exercise_1(filename)
% get number of frames from file
fileinfo = imfinfo(filename);
n_frames = length(fileinfo);

% save each frame in the file in grayscale
for i = 1:n_frames
    frames(:,:,:,i) = mat2gray(imread(filename,i));
end

% get summary images
mean_im = mean(frames,4);
median_im = median(frames,4);
var_im = var(frames,0,4);

% plot summary images
figure()
imshow(mean_im)
title("Mean Image")

figure()
imshow(median_im)
title("Median Image")

figure()
imshow(10.*var_im)
title("Variance Image")

% test range and standard deviation summary images
range_im = range(frames,4);
std_im = std(frames,0,4);

figure()
imshow(range_im)
title("Range Image")

figure()
imshow(std_im)
title("Standard Deviation Image")