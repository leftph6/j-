
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.linalg import eigh
from scipy.optimize import fsolve

# ==================== 任务1: 计算三对角矩阵的特征值 ====================
# 构造100维三对角矩阵
n = 100

# 主对角线
main_diag = np.arange(1, n+1, dtype=float)

# 次对角线（上下均为-1）
sub_diag = -np.ones(n-1)

# 构造矩阵
A = np.diag(main_diag) + np.diag(sub_diag, 1) + np.diag(sub_diag, -1)

# 右上角和左下角为-2
A[0, n-1] = -2
A[n-1, 0] = -2

print("任务1：三对角矩阵的特征值计算")
print("="*50)
print(f"矩阵大小: {n}×{n}")
print(f"\n矩阵A的左上角5×5子矩阵:")
print(A[:5, :5])
print(f"\n矩阵A的右下角5×5子矩阵:")
print(A[-5:, -5:])

# 计算所有特征值
eigenvalues = np.linalg.eigvalsh(A)
# 按降序排列
eigenvalues_sorted = np.sort(eigenvalues)[::-1]

# 最大的三个特征值
top_3_eigenvalues = eigenvalues_sorted[:3]

print(f"\n最大的三个特征值:")
print(f"λ₁ = {top_3_eigenvalues[0]:.10f}")
print(f"λ₂ = {top_3_eigenvalues[1]:.10f}")
print(f"λ₃ = {top_3_eigenvalues[2]:.10f}")

# ==================== 任务2: 平衡状态下小球的坐标 ====================
print("\n" + "="*50)
print("任务2：平衡状态下小球的坐标计算")
print("="*50)

# 101个小球，位置编号为1到101
# 第1个小球固定在(0,0)，第101个小球固定在(100,0)
# 相邻小球由弹簧连接，原长为0，弹性系数为100

# 对于第i个小球(i=2到100)，设其坐标为(x_i, y_i)
# 受力平衡条件：
# x方向：受到左右弹簧的水平拉力
# y方向：受到左右弹簧的竖直拉力+重力

# 设第i个小球坐标为(x_i, y_i)
# 第1个：(0, 0)  [固定]
# 第2个到第100个：待求
# 第101个：(100, 0)  [固定]

# 弹簧力向量：从i到i+1的位移向量乘以100（弹性系数）
# 对于第i个小球：
#   从i-1到i的弹簧力：100*[(x_i-x_{i-1}, y_i-y_{i-1})]
#   从i到i+1的弹簧力：100*[(x_{i+1}-x_i, y_{i+1}-y_i)]
# 重力：(0, -(2+sin(i-1)))，其中i从2到101

# 整理方程：对第i个小球(i=2到100)
# x方向：100(x_i - x_{i-1}) + 100(x_{i+1} - x_i) = 0
#        x_{i-1} - 2x_i + x_{i+1} = 0
# y方向：100(y_i - y_{i-1}) + 100(y_{i+1} - y_i) = 2+sin(i-1)
#        y_{i-1} - 2y_i + y_{i+1} = (2+sin(i-1))/100

# 边界条件：
# x_1 = 0, y_1 = 0
# x_{101} = 100, y_{101} = 0

# 对于x方向，方程为：x_{i-1} - 2x_i + x_{i+1} = 0，i=2到100
# 这是一个线性方程组，可以直接求解

# 构造x方向的方程组
# 有99个未知数（x_2到x_100），有99个方程
A_x = np.zeros((99, 99))
b_x = np.zeros(99)

for i in range(99):
    # 对于第(i+2)个小球（下标从2开始）
    if i == 0:
        # 第2个小球：x_1 - 2x_2 + x_3 = 0
        # -2x_2 + x_3 = -x_1 = 0
        A_x[0, 0] = -2
        A_x[0, 1] = 1
        b_x[0] = 0
    elif i == 98:
        # 第100个小球：x_99 - 2x_100 + x_101 = 0
        # x_99 - 2x_100 = -x_101 = -100
        A_x[98, 97] = 1
        A_x[98, 98] = -2
        b_x[98] = -100
    else:
        # 中间的小球
        A_x[i, i-1] = 1
        A_x[i, i] = -2
        A_x[i, i+1] = 1
        b_x[i] = 0

x_coords = np.linalg.solve(A_x, b_x)

# 对于y方向，方程为：y_{i-1} - 2y_i + y_{i+1} = (2+sin(i-1))/100，i=2到100
A_y = np.zeros((99, 99))
b_y = np.zeros(99)

for i in range(99):
    # 对于第(i+2)个小球（下标从2开始）
    ball_index = i + 2  # 小球编号
    if i == 0:
        # 第2个小球：y_1 - 2y_2 + y_3 = (2+sin(1))/100
        # -2y_2 + y_3 = (2+sin(1))/100 - y_1 = (2+sin(1))/100
        A_y[0, 0] = -2
        A_y[0, 1] = 1
        b_y[0] = (2 + np.sin(1)) / 100
    elif i == 98:
        # 第100个小球：y_99 - 2y_100 + y_101 = (2+sin(100))/100
        # y_99 - 2y_100 = (2+sin(100))/100
        A_y[98, 97] = 1
        A_y[98, 98] = -2
        b_y[98] = (2 + np.sin(100)) / 100
    else:
        # 中间的小球
        A_y[i, i-1] = 1
        A_y[i, i] = -2
        A_y[i, i+1] = 1
        b_y[i] = (2 + np.sin(ball_index - 1)) / 100

y_coords = np.linalg.solve(A_y, b_y)

# 第2个小球的坐标（对应x_coords[0]和y_coords[0]）
ball_2_x = x_coords[0]
ball_2_y = y_coords[0]

print(f"\n第2个小球的坐标:")
print(f"x₂ = {ball_2_x:.10f}")
print(f"y₂ = {ball_2_y:.10f}")

# 验证结果
print(f"\n验证结果（精确到4位有效数字）:")
print(f"x₂ ≈ {ball_2_x:.4f}")
print(f"y₂ ≈ {ball_2_y:.4f}")

# 绘制所有小球的位置
all_x = np.concatenate([[0], x_coords, [100]])
all_y = np.concatenate([[0], y_coords, [0]])

plt.figure(figsize=(14, 6))
plt.plot(all_x, all_y, 'b-', linewidth=1, label='连接线')
plt.plot(all_x, all_y, 'ro', markersize=3, label='小球位置')
plt.plot([0, 100], [0, 0], 'g*', markersize=15, label='固定小球')
plt.plot(ball_2_x, ball_2_y, 'go', markersize=8, label=f'第2个小球 ({ball_2_x:.4f}, {ball_2_y:.4f})')
plt.xlabel('x坐标')
plt.ylabel('y坐标')
plt.title('101个小球在平衡状态下的位置分布')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig('balls_equilibrium.png', dpi=150)
plt.close()

print("\n图表已保存为 'balls_equilibrium.png'")

# 保存结果供MATLAB代码参考
results = {
    'eigenvalue_1': top_3_eigenvalues[0],
    'eigenvalue_2': top_3_eigenvalues[1],
    'eigenvalue_3': top_3_eigenvalues[2],
    'ball_2_x': ball_2_x,
    'ball_2_y': ball_2_y,
    'all_x': all_x,
    'all_y': all_y
}

print("\n" + "="*50)
print("计算完成！")
