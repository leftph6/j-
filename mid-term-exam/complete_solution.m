%% 高级数值分析期中作业 - 完整解决方案
% 任务1: 计算100维三对角矩阵的最大三个特征值
% 任务2: 平衡状态下第2个小球的坐标

clear; close all; clc;

%% ======================== 任务1 开始 ========================
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║         高级数值分析期中作业 - MATLAB完整解决方案          ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

fprintf('【 任务1：三对角矩阵特征值计算 】\n');
fprintf('════════════════════════════════════════════════════════════\n');

% 第一步：构造矩阵A
n = 100;
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

fprintf('✓ 矩阵构造完成 (%d × %d)\n', n, n);
fprintf('  主对角线: 1, 2, 3, ..., 99, 100\n');
fprintf('  次对角线: -1\n');
fprintf('  角元素: -2\n\n');

% 第二步：计算特征值
eigenvalues = eig(A);
eigenvalues_sorted = sort(eigenvalues, 'descend');

% 提取最大的三个特征值
lambda1 = eigenvalues_sorted(1);
lambda2 = eigenvalues_sorted(2);
lambda3 = eigenvalues_sorted(3);

fprintf('✓ 特征值计算完成\n');
fprintf('  计算方法: 标准特征值分解 (eig function)\n\n');

fprintf('【 任务1 - 计算结果 】\n');
fprintf('────────────────────────────────────────────────────────────\n');
fprintf('最大的三个特征值:\n');
fprintf('  λ₁ = %.10f  ≈ %.6f\n', lambda1, lambda1);
fprintf('  λ₂ = %.10f  ≈ %.6f\n', lambda2, lambda2);
fprintf('  λ₃ = %.10f  ≈ %.6f\n', lambda3, lambda3);
fprintf('────────────────────────────────────────────────────────────\n\n');

% 第三步：绘制图表1
figure('Name', '任务1 - 特征值分析', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
plot(1:n, eigenvalues_sorted, 'bo', 'MarkerSize', 4, 'LineWidth', 1.5);
hold on;
plot(1:3, eigenvalues_sorted(1:3), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('特征值序号 (从大到小)', 'FontSize', 11);
ylabel('特征值大小', 'FontSize', 11);
title('100维三对角矩阵的所有特征值', 'FontSize', 12);
grid on;
legend('全部特征值', '最大的三个', 'FontSize', 10);

subplot(1, 2, 2);
spy(A);
title('矩阵A的稀疏结构', 'FontSize', 12);

sgtitle(['任务1结果: λ₁=', num2str(lambda1, '%.4f'), ...
         ', λ₂=', num2str(lambda2, '%.4f'), ...
         ', λ₃=', num2str(lambda3, '%.4f')], 'FontSize', 11);

print('-dpng', 'task1_eigenvalues.png', '-r150');
fprintf('✓ 任务1图表已保存: task1_eigenvalues.png\n\n');

% 验证
[V, D] = eig(A);
[D_sorted, idx] = sort(diag(D), 'descend');
fprintf('【 验证 - 特征值精度 】\n');
fprintf('────────────────────────────────────────────────────────────\n');
for i = 1:3
    v = V(:, idx(i));
    residual = norm(A*v - D_sorted(i)*v) / abs(D_sorted(i));
    fprintf('λ%d 相对残差: %.2e\n', i, residual);
end
fprintf('────────────────────────────────────────────────────────────\n\n');

%% ======================== 任务2 开始 ========================
fprintf('【 任务2：平衡状态下小球的坐标 】\n');
fprintf('════════════════════════════════════════════════════════════\n');

n_balls = 101;
n_unknowns = 99;

fprintf('✓ 系统参数\n');
fprintf('  小球总数: %d\n', n_balls);
fprintf('  未知变量: %d\n', n_unknowns);
fprintf('  弹簧原长: 0\n');
fprintf('  弹性系数: 100\n');
fprintf('  重力: (2+sin(i-1)) for i=2:101\n\n');

% ============ x方向方程组 ============
A_x = zeros(n_unknowns, n_unknowns);
b_x = zeros(n_unknowns, 1);

A_x(1, 1) = -2;
A_x(1, 2) = 1;
b_x(1) = 0;

for i = 2:n_unknowns-1
    A_x(i, i-1) = 1;
    A_x(i, i) = -2;
    A_x(i, i+1) = 1;
    b_x(i) = 0;
end

A_x(n_unknowns, n_unknowns-1) = 1;
A_x(n_unknowns, n_unknowns) = -2;
b_x(n_unknowns) = -100;

% ============ y方向方程组 ============
A_y = zeros(n_unknowns, n_unknowns);
b_y = zeros(n_unknowns, 1);

A_y(1, 1) = -2;
A_y(1, 2) = 1;
b_y(1) = (2 + sin(1)) / 100;

for i = 2:n_unknowns-1
    ball_index = i + 1;
    A_y(i, i-1) = 1;
    A_y(i, i) = -2;
    A_y(i, i+1) = 1;
    b_y(i) = (2 + sin(ball_index - 1)) / 100;
end

A_y(n_unknowns, n_unknowns-1) = 1;
A_y(n_unknowns, n_unknowns) = -2;
b_y(n_unknowns) = (2 + sin(100)) / 100;

fprintf('✓ 方程组构造完成\n');
fprintf('  系数矩阵大小: %d × %d (两个)\n', n_unknowns, n_unknowns);
fprintf('  方程类型: 线性三对角方程组\n\n');

% ============ 求解 ============
x_coords = A_x \ b_x;
y_coords = A_y \ b_y;

fprintf('✓ 方程组求解完成\n');
fprintf('  求解算法: LU分解\n\n');

% ============ 第2个小球的坐标 ============
ball_2_x = x_coords(1);
ball_2_y = y_coords(1);

fprintf('【 任务2 - 计算结果 】\n');
fprintf('────────────────────────────────────────────────────────────\n');
fprintf('第2个小球的坐标:\n');
fprintf('  x₂ = %.10f  ≈ %.6f\n', ball_2_x, ball_2_x);
fprintf('  y₂ = %.10f  ≈ %.6f\n', ball_2_y, ball_2_y);
fprintf('────────────────────────────────────────────────────────────\n\n');

% ============ 验证平衡条件 ============
all_x = [0; x_coords; 100];
all_y = [0; y_coords; 0];

fprintf('【 验证 - 平衡条件检查 】\n');
fprintf('────────────────────────────────────────────────────────────\n');
fprintf('选择前3个小球验证力的平衡:\n');

for i = 2:4
    force_x = 100*(all_x(i)-all_x(i-1)) + 100*(all_x(i+1)-all_x(i));
    force_y = 100*(all_y(i)-all_y(i-1)) + 100*(all_y(i+1)-all_y(i)) + (2 + sin(i-1));
    
    fprintf('\n小球%d:\n', i);
    fprintf('  x方向合力: %.2e (目标: 0)\n', force_x);
    fprintf('  y方向合力: %.2e (目标: 0)\n', force_y);
end
fprintf('\n────────────────────────────────────────────────────────────\n\n');

% ============ 绘制图表2 ============
figure('Name', '任务2 - 小球平衡分析', 'NumberTitle', 'off', ...
       'Position', [100, 100, 1400, 600]);

subplot(1, 2, 1);
plot(all_x, all_y, 'b-', 'LineWidth', 1);
hold on;
plot(all_x, all_y, 'ro', 'MarkerSize', 2);
plot(0, 0, 'g*', 'MarkerSize', 15, 'LineWidth', 2);
plot(100, 0, 'g*', 'MarkerSize', 15, 'LineWidth', 2);
plot(ball_2_x, ball_2_y, 'mo', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('x坐标', 'FontSize', 11);
ylabel('y坐标', 'FontSize', 11);
title('101个小球在平衡状态下的位置分布', 'FontSize', 12);
grid on;
legend('连接线', '小球', '固定点', '第2个小球', 'FontSize', 10);
axis equal;

subplot(1, 2, 2);
plot(0:100, all_y, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.5);
xlabel('小球编号', 'FontSize', 11);
ylabel('y坐标', 'FontSize', 11);
title('y坐标随小球编号的变化', 'FontSize', 12);
grid on;
hold on;
plot(2, ball_2_y, 'mo', 'MarkerSize', 10, 'LineWidth', 2);
text(2, ball_2_y-0.015, sprintf('第2个: %.4f', ball_2_y), ...
     'HorizontalAlignment', 'center', 'FontSize', 10);

sgtitle(['任务2结果: 第2个小球 (', num2str(ball_2_x, '%.4f'), ...
         ', ', num2str(ball_2_y, '%.4f'), ')'], 'FontSize', 11);

print('-dpng', 'task2_equilibrium.png', '-r150');
fprintf('✓ 任务2图表已保存: task2_equilibrium.png\n\n');

% ============ 统计信息 ============
fprintf('【 统计信息 】\n');
fprintf('────────────────────────────────────────────────────────────\n');
[y_min, min_idx] = min(all_y);
fprintf('最大下沉深度: %.6f (在第%d个小球)\n', abs(y_min), min_idx-1);
fprintf('x坐标范围: [%.6f, %.6f]\n', min(all_x), max(all_x));
fprintf('y坐标范围: [%.6f, %.6f]\n', min(all_y), max(all_y));
fprintf('────────────────────────────────────────────────────────────\n\n');

%% ======================== 最终总结 ========================
fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                    最终结果总结                             ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n');
fprintf('\n');

fprintf('【 任务1 】三对角矩阵最大特征值\n');
fprintf('  λ₁ = %.6f\n', lambda1);
fprintf('  λ₂ = %.6f\n', lambda2);
fprintf('  λ₃ = %.6f\n\n', lambda3);

fprintf('【 任务2 】第2个小球坐标\n');
fprintf('  (x₂, y₂) = (%.6f, %.6f)\n\n', ball_2_x, ball_2_y);

fprintf('【 精度说明 】\n');
fprintf('  所有结果均采用双精度浮点数 (精度 ~15-17位)\n');
fprintf('  满足题目要求: 至少4位精确有效数字\n\n');

fprintf('【 文件生成 】\n');
fprintf('  ✓ task1_eigenvalues.png\n');
fprintf('  ✓ task2_equilibrium.png\n');
fprintf('  ✓ 运行日志已显示在命令窗口\n\n');

fprintf('╔════════════════════════════════════════════════════════════╗\n');
fprintf('║                    程序运行完成！                          ║\n');
fprintf('╚════════════════════════════════════════════════════════════╝\n');
