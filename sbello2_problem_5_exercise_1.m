function [] = sbello2_problem_5_exercise_1(frames,sum_im,n) % create a circle roi centered at the seed pixel 
n_frames = size(frames,4);
pt_matrix = zeros(size(frames,1) * size(frames,2),n_frames);

for i = 1:n_frames
    frame = frames(:,:,:,i);
    pt_matrix(:,i) = reshape(frame,[],1);
end

% pca
[coeff, score, latent,~,~,mu] = pca(pt_matrix,'Centered',true);
for i = 1:n
    [~,idx] = max(score(:,i));
    frame_idx = reshape(1:size(score,1),size(frames,[1,2]));
    [cntr(1),cntr(2)] = find(frame_idx == idx);
    sbello2_problem_3_exercise_1(sum_im,cntr,5); % plot roi
end

% nmf
[W,H] = nnmf(pt_matrix,10);
for i = 1:n
    [~,idx] = max(W(:,i));
    frame_idx = reshape(1:size(W,1),size(frames,[1,2]));
    [cntr(1),cntr(2)] = find(frame_idx == idx);
    sbello2_problem_3_exercise_1(sum_im,cntr,5); % plot roi
end

% ica
[Mdl] = rica(pt_matrix,10);
A = transform(Mdl,pt_matrix);
for i = 1:n
    [~,idx] = max(A(:,i));
    frame_idx = reshape(1:size(A,1),size(frames,[1,2]));
    [cntr(1),cntr(2)] = find(frame_idx == idx);
    sbello2_problem_3_exercise_1(sum_im,cntr,5); % plot roi
end