%% 任务1: 计算100维三对角矩阵的最大三个特征值
% 高级数值分析期中作业
% 矩阵A为100维三对角矩阵
% 主对角线：1, 2, 3, ..., 99, 100
% 次对角线：-1
% 右上角和左下角：-2

clear; close all; clc;

%% 第一步：构造矩阵A
n = 100;

% 方法1: 直接构造（用于较小规模矩阵）
A = zeros(n, n);

% 主对角线
for i = 1:n
    A(i, i) = i;
end

% 次对角线（上、下）
for i = 1:n-1
    A(i, i+1) = -1;
    A(i+1, i) = -1;
end

% 右上角和左下角
A(1, n) = -2;
A(n, 1) = -2;

fprintf('任务1：三对角矩阵的特征值计算\n');
fprintf('===========================================\n');
fprintf('矩阵大小: %d × %d\n', n, n);
fprintf('\n矩阵A的左上角5×5子矩阵:\n');
disp(A(1:5, 1:5));
fprintf('矩阵A的右下角5×5子矩阵:\n');
disp(A(n-4:n, n-4:n));

%% 第二步：计算特征值
% 使用eig函数计算所有特征值
eigenvalues = eig(A);

% 按降序排列
eigenvalues_sorted = sort(eigenvalues, 'descend');

% 提取最大的三个特征值
lambda_1 = eigenvalues_sorted(1);
lambda_2 = eigenvalues_sorted(2);
lambda_3 = eigenvalues_sorted(3);

fprintf('\n最大的三个特征值:\n');
fprintf('λ₁ = %.10f\n', lambda_1);
fprintf('λ₂ = %.10f\n', lambda_2);
fprintf('λ₃ = %.10f\n', lambda_3);

fprintf('\n精确到4位有效数字:\n');
fprintf('λ₁ ≈ %.4f\n', lambda_1);
fprintf('λ₂ ≈ %.4f\n', lambda_2);
fprintf('λ₃ ≈ %.4f\n', lambda_3);

%% 第三步：可视化
figure('Position', [100, 100, 1200, 500]);

% 子图1: 所有特征值分布
subplot(1, 2, 1);
plot(1:n, eigenvalues_sorted, 'bo', 'MarkerSize', 4, 'LineWidth', 1.5);
hold on;
plot(1:3, eigenvalues_sorted(1:3), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('特征值序号 (从大到小)');
ylabel('特征值大小');
title('100维三对角矩阵的所有特征值');
grid on;
legend('全部特征值', '最大的三个');

% 子图2: 矩阵结构可视化
subplot(1, 2, 2);
spy(A);
title('矩阵A的稀疏结构');

savefig('task1_eigenvalues.fig');
print('-dpng', 'task1_eigenvalues.png', '-r150');
close;

fprintf('\n图表已保存为 task1_eigenvalues.png\n');

%% 第四步：验证（可选）
fprintf('\n===========================================\n');
fprintf('验证: 检查计算结果的精度\n');
fprintf('===========================================\n');

% 验证特征值的正确性
residuals = zeros(3, 1);
[V, D] = eig(A);
[D_sorted, idx] = sort(diag(D), 'descend');

for i = 1:3
    v = V(:, idx(i));
    % 计算 ||Av - λv|| / ||λv||
    residuals(i) = norm(A*v - D_sorted(i)*v) / abs(D_sorted(i));
    fprintf('第%d个特征值的相对残差: %.2e\n', i, residuals(i));
end

% 所有结果至少达到4位精确有效数字
fprintf('\n===========================================\n');
fprintf('最终结果（至少4位精确有效数字）:\n');
fprintf('===========================================\n');
fprintf('λ₁ = %.6f\n', lambda_1);
fprintf('λ₂ = %.6f\n', lambda_2);
fprintf('λ₃ = %.6f\n', lambda_3);
