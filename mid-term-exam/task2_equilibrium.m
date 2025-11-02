%% 任务2: 平衡状态下第2个小球的坐标
% 高级数值分析期中作业
% 101个小球由弹簧连接
% 第1个小球固定在(0,0)，第101个小球固定在(100,0)
% 弹簧原长为0，弹性系数为100
% 第i个小球受到大小为(2+sin(i-1))的向下重力

clear; close all; clc;

fprintf('任务2：平衡状态下小球的坐标计算\n');
fprintf('===========================================\n');

%% 第一步：建立方程组
% 对于第i个小球(i=2到100)的平衡条件：
% x方向：x_{i-1} - 2x_i + x_{i+1} = 0
% y方向：y_{i-1} - 2y_i + y_{i+1} = (2+sin(i-1))/100
% 边界条件：(x_1, y_1) = (0, 0), (x_101, y_101) = (100, 0)

n_balls = 101;  % 总共101个小球
n_unknowns = 99;  % 99个未知数（第2到100个小球）

% ============ x方向方程组 ============
A_x = zeros(n_unknowns, n_unknowns);
b_x = zeros(n_unknowns, 1);

% 第2个小球(i=2)：x_1 - 2x_2 + x_3 = 0
% -2x_2 + x_3 = -x_1 = 0
A_x(1, 1) = -2;
A_x(1, 2) = 1;
b_x(1) = 0;

% 中间的小球(i=3到100)
for i = 2:n_unknowns-1
    A_x(i, i-1) = 1;
    A_x(i, i) = -2;
    A_x(i, i+1) = 1;
    b_x(i) = 0;
end

% 第100个小球(i=100)：x_99 - 2x_100 + x_101 = 0
% x_99 - 2x_100 = -x_101 = -100
A_x(n_unknowns, n_unknowns-1) = 1;
A_x(n_unknowns, n_unknowns) = -2;
b_x(n_unknowns) = -100;

% 求解x方向
x_coords = A_x \ b_x;

% ============ y方向方程组 ============
A_y = zeros(n_unknowns, n_unknowns);
b_y = zeros(n_unknowns, 1);

% 第2个小球(i=2)：y_1 - 2y_2 + y_3 = (2+sin(1))/100
% -2y_2 + y_3 = (2+sin(1))/100 - y_1 = (2+sin(1))/100
A_y(1, 1) = -2;
A_y(1, 2) = 1;
b_y(1) = (2 + sin(1)) / 100;

% 中间的小球(i=3到100)
for i = 2:n_unknowns-1
    ball_index = i + 1;  % 小球编号
    A_y(i, i-1) = 1;
    A_y(i, i) = -2;
    A_y(i, i+1) = 1;
    b_y(i) = (2 + sin(ball_index - 1)) / 100;
end

% 第100个小球(i=100)：y_99 - 2y_100 + y_101 = (2+sin(100))/100
% y_99 - 2y_100 = (2+sin(100))/100 - y_101 = (2+sin(100))/100
A_y(n_unknowns, n_unknowns-1) = 1;
A_y(n_unknowns, n_unknowns) = -2;
b_y(n_unknowns) = (2 + sin(100)) / 100;

% 求解y方向
y_coords = A_y \ b_y;

% ============ 第2个小球的坐标 ============
ball_2_x = x_coords(1);
ball_2_y = y_coords(1);

fprintf('\n第2个小球的坐标:\n');
fprintf('x₂ = %.10f\n', ball_2_x);
fprintf('y₂ = %.10f\n', ball_2_y);

fprintf('\n精确到4位有效数字:\n');
fprintf('x₂ ≈ %.4f\n', ball_2_x);
fprintf('y₂ ≈ %.4f\n', ball_2_y);

%% 第二步：验证和分析
fprintf('\n===========================================\n');
fprintf('验证: 检查平衡条件\n');
fprintf('===========================================\n');

% 构造所有小球的坐标
all_x = [0; x_coords; 100];
all_y = [0; y_coords; 0];

% 验证前几个小球的平衡条件
fprintf('\n前3个小球的平衡条件验证:\n');
for i = 2:4
    % 力的平衡方程
    force_x = 100*(all_x(i)-all_x(i-1)) + 100*(all_x(i+1)-all_x(i));
    force_y = 100*(all_y(i)-all_y(i-1)) + 100*(all_y(i+1)-all_y(i)) + (2 + sin(i-1));
    
    fprintf('小球%d:\n', i);
    fprintf('  x方向合力 = %.2e (应接近0)\n', force_x);
    fprintf('  y方向合力 = %.2e (应接近0)\n', force_y);
end

%% 第三步：可视化
figure('Position', [100, 100, 1400, 600]);

% 子图1: 所有小球的位置
subplot(1, 2, 1);
plot(all_x, all_y, 'b-', 'LineWidth', 1);
hold on;
plot(all_x, all_y, 'ro', 'MarkerSize', 2);
plot(0, 0, 'g*', 'MarkerSize', 15, 'LineWidth', 2);
plot(100, 0, 'g*', 'MarkerSize', 15, 'LineWidth', 2);
plot(ball_2_x, ball_2_y, 'mo', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('x坐标');
ylabel('y坐标');
title('101个小球在平衡状态下的位置分布');
grid on;
legend('连接线', '小球', '固定点', '第2个小球');
axis equal;

% 子图2: y坐标的变化趋势
subplot(1, 2, 2);
plot(0:100, all_y, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.5);
xlabel('小球编号');
ylabel('y坐标');
title('y坐标随小球编号的变化');
grid on;
hold on;
plot(2, ball_2_y, 'mo', 'MarkerSize', 10, 'LineWidth', 2);
text(2, ball_2_y-0.02, sprintf('第2个: %.4f', ball_2_y), 'HorizontalAlignment', 'center');

savefig('task2_equilibrium.fig');
print('-dpng', 'task2_equilibrium.png', '-r150');
close;

fprintf('\n\n图表已保存为 task2_equilibrium.png\n');

%% 第四步：统计信息
fprintf('\n===========================================\n');
fprintf('统计信息:\n');
fprintf('===========================================\n');
fprintf('最大下沉深度(|y_min|) = %.6f (在第%d个小球)\n', min(all_y), find(all_y == min(all_y))-1);
fprintf('小球x坐标范围: [%.4f, %.4f]\n', min(all_x), max(all_x));
fprintf('小球y坐标范围: [%.6f, %.6f]\n', min(all_y), max(all_y));

fprintf('\n===========================================\n');
fprintf('最终结果（至少4位精确有效数字）:\n');
fprintf('===========================================\n');
fprintf('第2个小球坐标: (%.6f, %.6f)\n', ball_2_x, ball_2_y);
fprintf('===========================================\n');
