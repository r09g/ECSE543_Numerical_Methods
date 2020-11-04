clear; close all;

load r_evo.mat
l2norm = zeros(size(r,2),1);
linfnorm = zeros(size(r,2),1);

for i = 1:size(r,2)
    l2norm(i) = norm(r(:,i),2);
    linfnorm(i) = max(abs(r(:,i)));
end

figure,
hold on,
plot(0:size(r,2)-1, l2norm)
plot(0:size(r,2)-1, linfnorm)
xlim([0 size(r,2)-1])