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
plot(0:size(r,2)-1, l2norm, 'LineWidth', 1)
plot(0:size(r,2)-1, linfnorm, 'LineWidth', 1)
xlim([0 size(r,2)-1])
grid on
title("L2 Norm and L\infty Norm")
xlabel("Iterations")
ylabel("Value")
legend("L2 Norm","L\infty Norm")
