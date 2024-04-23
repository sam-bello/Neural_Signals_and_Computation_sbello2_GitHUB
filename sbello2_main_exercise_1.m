%% Exercise 1
clear
close all

%% Problem 1
movie1 = 'TEST_MOVIE_00001-small-motion.tif';
[frames1,corr14] = sbello2_problem_1_exercise_1(movie1);

%% Problem 2
movie2 = 'TEST_MOVIE_00001-small.tif';
[frames2,mean_im,median_im,var_im,range_im,std_im] = sbello2_problem_2_exercise_1(movie2);

%% Problem 3
seeds = [58 355;
        147 119;
        270 344
        398 415
        327 353];
min_rad = 5;
for i = 1:size(seeds,1)
    masks(:,:,i) = sbello2_problem_3_exercise_1(range_im,seeds(i,:),min_rad);
end

%% Problem 4
for i = 1:size(seeds,1)
    trace = sbello2_problem_4_exercise_1(frames2,masks(:,:,i),seeds(i,:));
end

%% Problem 5
sbello2_problem_5_exercise_1(frames2,range_im,5) % create a circle roi centered at the seed pixel 
