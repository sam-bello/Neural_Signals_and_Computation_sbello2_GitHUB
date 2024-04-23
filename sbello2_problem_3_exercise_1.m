function [mask] = sbello2_problem_3_exercise_1(sum_im,seed,min_rad) % create a circle roi centered at the seed pixel 
% initialize circular mask with minimum radius
[x,y] = ndgrid((1:size(sum_im,1)) - seed(1),(1:size(sum_im,2)) - seed(2));
mask = (x.^2 + y.^2) < (min_rad^2);
init_mask_mean = mean(sum_im(mask),"all"); % mean
init_mask_std = std(sum_im(mask),0,"all"); % standard deviations

% increas mask radius but 1 until outer ring being added has a mean that is
% more than one standard deviation away from the initial mask
fin = 0;
new_rad = min_rad + 1;
while ~fin
    % create new mask with new radius and get outer ring of new mask
    new_mask = (x.^2 + y.^2) < (new_rad^2);
    outer_ring = new_mask ~= mask;
    ring_mean = mean(sum_im(outer_ring),"all");
    
    % check if mean of outer ring is more than one standard deviation from
    % initial mask mean
    if abs(ring_mean - init_mask_mean) > init_mask_std
        fin = 1;
    end

    % update mask and radius
    mask = new_mask;
    new_rad = new_rad + 1;
end

% plot roi
roi = sum_im;
roi(mask ~= 1) = 0;
figure()
title(strcat("ROI with seed pixel: ",mat2str(seed)))
imshow(roi)