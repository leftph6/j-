% 快速参考指南 - 期中作业MATLAB解决方案
% 
% 本文件包含所有关键代码片段，可直接复制使用
% ===========================================================

%% 【快速开始】运行完整解决方案
% 在MATLAB命令窗口输入：
% >> run('complete_solution.m')
% 
% 或者：
% >> cd到文件目录
% >> complete_solution

%% ==================== 任务1 快速代码 ====================

% 【最小化代码版本】
clear; clc;

% 构造矩阵
A = diag(1:100) + diag(-ones(99,1),1) + diag(-ones(99,1),-1);
A(1,100) = -2;
A(100,1) = -2;

% 计算特征值
evals = sort(eig(A), 'descend');

% 显示结果
fprintf('最大三个特征值:\n');
fprintf('λ1 = %.10f\n', evals(1));
fprintf('λ2 = %.10f\n', evals(2));
fprintf('λ3 = %.10f\n', evals(3));

%% ==================== 任务2 快速代码 ====================

% 【最小化代码版本】
clear; clc;

% 构造系数矩阵（x方向）
n = 99;
A_x = -2*eye(n) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
b_x = zeros(n,1);
b_x(n) = -100;

% 构造系数矩阵（y方向）
A_y = -2*eye(n) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
b_y = zeros(n,1);
for i = 1:n
    idx = i + 1;
    b_y(i) = (2 + sin(idx-1)) / 100;
end
b_y(n) = (2 + sin(100)) / 100;

% 求解
x = A_x \ b_x;
y = A_y \ b_y;

% 显示结果
fprintf('第2个小球的坐标:\n');
fprintf('x2 = %.10f\n', x(1));
fprintf('y2 = %.10f\n', y(1));

%% ==================== 调试技巧 ====================

% 【验证矩阵构造】
% spy(A)  % 显示矩阵稀疏结构
% A(1:5, 1:5)  % 显示左上角子矩阵

% 【验证特征值计算】
% [V, D] = eig(A);
% residual = norm(A*V(:,1) - D(1,1)*V(:,1)) / abs(D(1,1));
% % residual应该很小（< 1e-10）

% 【验证平衡条件】
% force = A*x;  % 应该接近b
% norm(A*x - b) / norm(b)  % 相对误差应该< 1e-12

%% ==================== 常用MATLAB命令 ====================

% 矩阵构造
% A = zeros(n, n)        % 零矩阵
% A = eye(n)             % 单位矩阵
% A = diag(v)            % 对角矩阵
% A = diag(v, k)         % k为偏移（1=上,−1=下）
% A = [1,2; 3,4]         % 直接输入

% 矩阵操作
% A'                     % 转置
% A + B                  % 加法
% A * B                  % 矩阵乘法
% A .* B                 % 元素乘法
% A \ b                  % 求解 Ax=b

% 特征值特征向量
% eig(A)                 % 所有特征值
% [V, D] = eig(A)        % V=特征向量，D=特征值

% 排序
% sort(x)                % 升序排列
% sort(x, 'descend')     % 降序排列
% [v, idx] = sort(x)     % 返回值和索引

% 向量操作
% norm(x)                % 向量范数
% norm(A, 'fro')         % Frobenius范数
% dot(a, b)              % 内积

%% ==================== 性能对比 ====================

% 测试不同大小的矩阵
sizes = [50, 100, 200, 500];
times = zeros(size(sizes));

for i = 1:length(sizes)
    n = sizes(i);
    A = diag(1:n) + diag(-ones(n-1,1),1) + diag(-ones(n-1,1),-1);
    A(1,n) = -2;
    A(n,1) = -2;
    
    tic;
    eig(A);
    times(i) = toc;
end

% 绘制性能曲线
% figure;
% loglog(sizes, times, 'bo-');
% hold on;
% loglog(sizes, (sizes/100).^3 * times(2), 'r--');
% xlabel('矩阵维数');
% ylabel('计算时间(秒)');
% legend('实际时间', 'O(n^3)参考');

%% ==================== 输出格式控制 ====================

% 改变数字显示格式
% format short        % 5位有效数字（默认）
% format long         % 15位有效数字
% format scientific   % 科学计数法
% format rat          % 有理数

% 使用fprintf精确控制
% fprintf('%d\n', 5)           % 整数
% fprintf('%.4f\n', 3.14159)   % 4位小数
% fprintf('%.4e\n', 0.0001)    % 科学计数法

%% ==================== 错误处理 ====================

% try-catch语句
try
    x = A \ b;
catch ME
    fprintf('错误: %s\n', ME.message);
end

% 检查矩阵奇异性
if rcond(A) < 1e-10
    warning('矩阵可能为奇异矩阵');
end

% 验证解
residual = norm(A*x - b) / norm(b);
if residual > 1e-10
    warning('解的精度可能不够好');
end

%% ==================== 导出结果 ====================

% 保存为文本文件
% save results.txt lambda1 lambda2 lambda3 -ascii

% 保存为MATLAB文件
% save results.mat A evals x y

% 导出到Excel
% writematrix([evals(1:3)], 'results.xlsx', 'Sheet', 1)

%% ==================== 常见问题排查 ====================

% Q: 为什么计算结果有差异?
% A: 浮点数精度、舍入误差、算法差异等。使用长格式对比：
%    format long
%    result

% Q: 矩阵太大怎么办？
% A: 使用稀疏矩阵：
%    A_sparse = sparse(A);
%    eig(A_sparse)

% Q: 特征值有复数？
% A: 非对称矩阵可能有复数特征值。题目矩阵对称，不会出现此情况。

% Q: 线性方程组无解？
% A: 检查矩阵是否奇异：
%    det(A)      % 接近0则奇异
%    rank(A)     % 应等于矩阵维数

%% ==================== 进阶优化 ====================

% 使用稀疏矩阵加快计算
A_sparse = spdiags([[-2*ones(100,1)], [-ones(100,1)], [1:100]', ...
                    [-ones(100,1)]], -1:1, 100, 100);
A_sparse(1,100) = -2;
A_sparse(100,1) = -2;

% 效果：对于大矩阵可节省内存和时间

% 使用并行计算
% parfor i = 1:1000
%     result(i) = compute_something(data(i));
% end

%% ==================== 可视化技巧 ====================

% 绘制特征值分布
% histogram(evals, 20)

% 绘制热力图
% imagesc(A)
% colorbar

% 绘制3D曲线
% plot3(x, y, z)

% 绘制等高线
% contour(X, Y, Z, 20)

%% ==================== 参数敏感性分析 ====================

% 改变矩阵维数的影响
for n = [50, 100, 200]
    A = diag(1:n) + diag(-ones(n-1,1),1) + diag(-ones(n-1,1),-1);
    evals = sort(eig(A), 'descend');
    fprintf('n=%d: λ1=%.4f\n', n, evals(1));
end

% 改变重力的影响
for gravity_scale = [0.5, 1.0, 2.0]
    b_y_scaled = gravity_scale * b_y;
    y_scaled = A_y \ b_y_scaled;
    fprintf('gravity_scale=%.1f: y2=%.4f\n', gravity_scale, y_scaled(1));
end

%% ==================== 最终检查清单 ====================

% 提交前检查：
% ✓ 矩阵维数正确 (100×100)
% ✓ 特征值至少4位有效数字
% ✓ 坐标至少4位有效数字
% ✓ 代码注释完整
% ✓ 图表已生成并保存
% ✓ 结果已验证
% ✓ 没有警告或错误信息

%% ==================== 参考资源 ====================

% MATLAB官方文档
% https://www.mathworks.com/help/matlab/

% 常用函数快速查询
% help eig
% help linsolve
% help sparse

% 更多例子见完整代码文件：
% - complete_solution.m        (完整解决方案)
% - task1_eigenvalues.m        (任务1详细版)
% - task2_equilibrium.m        (任务2详细版)

% END OF QUICK REFERENCE
