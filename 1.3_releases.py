import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

#发动机参数
n_min = 600 #发动机最低转速
n_max = 4000 #发动机最高转速
n_max_power = 3800 # 最大功率转速，假设为3800转，需要从发动机特性确定。
n_max_torque = 2500  # 最大扭矩转速，假设为2500转，需要从发动机特性确定。

#车辆质量参数
m_load = 2000 #装载质量 (kg)
m_empty = 1800 #整车整备质量 (kg)
m_total = 3880 #总质量 (kg)

# 速度参数
v_max = 120  # 最大车速 (km/h)
v_min = 0  # 最小车速 (km/h)
v_min_2=5.14229 # 二档最低车速
# 车轮参数
r = 0.367  # 车轮半径 (m)

# 传动系统参数
eta_t = 0.85  # 传动系统机械效率
f = 0.013  # 滚动阻力系数
CdA = 2.77  # 空气阻力系数 × 迎风面积 (m^2)
i0 = 5.83  # 主减速器传动比

# 转动惯量
I_f = 0.218  # 飞轮转动惯量 (kg*m^2)
I_w1 = 1.798  # 前轮转动惯量 (kg*m^2)
I_w2 = 3.598  # 后轮转动惯量 (kg*m^2)

#变速箱传动比(4档)
ig_1_4 = 6.09
ig_2_4 = 3.09
ig_3_4 = 1.71
ig_4_4 = 1.00

#变速箱传动比(5档)
ig_1_5 = 5.56
ig_2_5 = 2.769
ig_3_5 = 1.644
ig_4_5 = 1.00
ig_5_5 = 0.793

#轴距 & 质心位置
L = 3.2 #轴距
a = 1.947 #质心至前轴距离（满载）
h_g = 0.9 #质心高

#空气密度
rho = 1.225  # 空气密度 (kg/m^3) 在标准大气条件下

g = 9.81 #重力加速度
G = m_total * g #重力

# 考虑转动惯量的牛二
def cal_delta_m(gear):
    return m_total + (2*I_w1 + 2*I_w2 + I_f *i0 **2 *gear_ratios[gear-1] **2 )/(r**2)

# 定义函数：计算发动机转矩
def engine_torque(n):
    """n: 发动机转速 (r/min)"""
    n_scaled = n / 1000
    Tq = -19.313 + 295.27 * n_scaled - 165.44 * n_scaled**2 + 40.874 * n_scaled**3 - 3.8445 * n_scaled**4
    return Tq

def traction_force(n, ig):
    """
    n: 发动机转速 (r/min)
    ig: 变速器传动比
    """
    Tq = engine_torque(n)
    Tw = Tq * ig * i0 * eta_t  # 车轮转矩
    F = Tw / r  # 驱动力
    return F

def resistance_force(v):
    """v: 车速 (km/h)"""
    v_mps = v / 3.6  # 将 km/h 转换为 m/s
    Fr = G * f  # 滚动阻力
    Fa = 0.5 * rho * CdA * v_mps**2  # 空气阻力
    F_total = Fr + Fa  # 总阻力
    return F_total

# 选择变速箱类型 (4 或 5)
gearbox_type = 5  # 选择 5 档变速箱

# 定义传动比列表
if gearbox_type == 4:
    gear_ratios = [ig_1_4, ig_2_4, ig_3_4, ig_4_4]
    num_gears = 4
elif gearbox_type == 5:
    gear_ratios = [ig_1_5, ig_2_5, ig_3_5, ig_4_5, ig_5_5]
    num_gears = 5
else:
    raise ValueError("Invalid gearbox type. Choose 4 or 5.")

def calculate_max_speed_iterative(gear_ratios, n_max, r, i0, eta_t, f, CdA, rho, G, v_min, v_max, resistance_force_values): # 计算最高车速
    max_speed = 0
    v_range = np.linspace(v_min, v_max, 100)
    for ig in gear_ratios:
        gear_max_speed = 0  # 当前档位的最高速度
        for v in v_range:
            # 计算发动机转速
            n = (v / 3.6) * 60 * i0 * ig / (2 * np.pi * r)

            # 检查转速是否超过最大转速
            if n > n_max:
                gear_max_speed = (n_max / 60) * (2 * np.pi * r) * (3.6 / (i0 * ig))
                break  # 超过最大转速，退出循环
            
            # 计算驱动力和阻力
            index = np.argmin(np.abs(v_range - v))
            F_traction  = traction_force(n, ig) 
            F_resistance = resistance_force_values[index]  # 使用之前定义的函数

            # 检查驱动力是否小于等于阻力
            if F_traction <= F_resistance and v>80:
                gear_max_speed = v
                break  # 找到平衡点，退出循环
        max_speed = max(max_speed, gear_max_speed)

    return max_speed


def calculate_max_grade(gear_ratios, n_max_torque, r, i0, eta_t, f, G):
    """计算最大爬坡度"""
    ig_1 = gear_ratios[0]
    Tq_max = engine_torque(n_max_torque)
    F_traction_max = (Tq_max * ig_1 * i0 * eta_t) / r
    theta_max=0
    for theta in np.linspace(0, np.pi/2, 100):
        if G * np.sin(theta)+G * f * np.cos(theta) > F_traction_max:
            break
        else:
            theta_max=theta
    max_grade = np.tan(theta_max) * 100 
    return max_grade

def calculate_adhesion_rate(m_total, g, a, h_g, L, theta, F_traction_max):
    """计算附着率"""
    F_z = (m_total * g * (a * np.cos(theta) + h_g * np.sin(theta))) / L
    mu = F_traction_max / F_z
    return mu

# 创建图表
plt.figure(figsize=(18, 12))  # 调整图表大小, 更大的图表，容纳更多子图

# 绘制驱动力-车速曲线
plt.subplot(3, 2, 1)  # 创建一个子图，调整为 3 行 2 列
v_range = np.linspace(v_min, v_max, 100)  # 使用统一的车速范围
for i in range(num_gears):
    ig = gear_ratios[i]
    # 生成发动机转速范围
    n_range = np.linspace(n_min, n_max, 100)
    # 计算车速 (km/h)
    v_range_temp = n_range * r * 2 * np.pi / 60 / i0 / ig * 3.6
    # 计算驱动力
    traction_force_values = [traction_force(n, ig) for n in n_range]
    plt.plot(v_range_temp, traction_force_values, label=f'Gear {i+1}')

plt.xlabel('Vehicle Speed (km/h)')
plt.ylabel('Traction Force (N)')
plt.title('Traction Force vs. Speed ({} speed)'.format(gearbox_type))
plt.grid(True)
plt.legend()  # 显示图例
plt.xlim(0, v_max)  # 限制 x 轴范围
plt.ylim(0, 15000) #限制y轴范围

# 绘制扭矩-转速曲线
plt.subplot(3, 2, 2)  # 创建一个子图
n_range = np.linspace(n_min, n_max, 100)  # 在最小转速和最大转速之间生成 100 个点
torque_values = [engine_torque(n) for n in n_range]
plt.plot(n_range, torque_values)
plt.xlabel('Engine Speed (r/min)')
plt.ylabel('Torque (N*m)')
plt.title('Engine Torque vs. Speed')
plt.grid(True)

# 绘制阻力-速度曲线
plt.subplot(3, 2, 3)  # 创建一个子图
v_range = np.linspace(v_min, v_max, 100)  # 使用统一的车速范围
resistance_force_values = [resistance_force(v) for v in v_range]
plt.plot(v_range, resistance_force_values)
plt.xlabel('Vehicle Speed (km/h)')
plt.ylabel('Resistance Force (N)')  # 修改 y 轴标签
plt.title('Resistance Force vs. Speed') # 修改标题
plt.grid(True)
plt.xlim(0, v_max)  # 限制 x 轴范围

# 绘制驱动力与阻力平衡图
plt.subplot(3, 2, 4)  # 创建一个子图

v_range = np.linspace(v_min, v_max, 100)  # 使用统一的车速范围
resistance_force_values = [resistance_force(v) for v in v_range]  # 计算阻力

# 绘制每一档的驱动力曲线和阻力曲线
for i in range(num_gears):
    ig = gear_ratios[i]
    # 计算驱动力
    n_range = np.linspace(n_min, n_max, 100)
    v_range_temp = n_range * r * 2 * np.pi / 60 / i0 / ig * 3.6
    traction_force_values = [traction_force(n, ig) for n in n_range]

    # 截断驱动力曲线，只保留与阻力曲线车速范围重叠的部分
    valid_indices = (v_range_temp >= v_min) & (v_range_temp <= v_max)
    v_range_truncated = v_range_temp[valid_indices]
    traction_force_values_truncated = np.array(traction_force_values)[valid_indices]

    plt.plot(v_range_truncated, traction_force_values_truncated, label=f'Gear {i+1}')  # 绘制驱动力曲线

plt.plot(v_range, resistance_force_values, label='Resistance', color='black', linestyle='--')  # 绘制阻力曲线，黑色虚线
plt.xlabel('Vehicle Speed (km/h)')
plt.ylabel('Force (N)')
plt.title('Traction & Resistance Force Balance')
plt.grid(True)
plt.xlim(0, v_max)  # 限制 x 轴范围
plt.ylim(0, 15000)
plt.legend()

# 添加加速度倒数曲线图
plt.subplot(3, 2, 5)
v_range = np.linspace(v_min, 70, 200) # 速度范围从 v_min 到 70 km/h
# ig = gear_ratios[1] # 第二个传动比对应2档
def inverse_acceleration(v): #定义加速度倒数的函数
        gear=2
        ig = gear_ratios[gear-1]
        n = (v / 3.6) * 60 * i0 * ig / (2 * np.pi * r)
        while(n>n_max):
            gear+=1
            ig = gear_ratios[gear-1]
            n = (v / 3.6) * 60 * i0 * ig / (2 * np.pi * r)
        F_traction = traction_force(n, ig)
        F_resistance = resistance_force(v)
        F_net = F_traction - F_resistance
        delta_m=cal_delta_m(gear)
        a = F_net / delta_m
        if a <= 0:
           return np.inf
        return 1/a
            
acceleration_inverse_values = [inverse_acceleration(v) for v in v_range]
plt.plot(v_range, acceleration_inverse_values)
plt.xlabel('Vehicle Speed (km/h)')
plt.ylabel('Inverse of Acceleration (s^2/m)')
plt.title('Inverse of Acceleration vs. Speed (Gear 2&3&4)')
plt.grid(True)
plt.xlim(v_min, 70)


# 计算最大速度
v_range = np.linspace(v_min, v_max, 100)
resistance_force_values = [resistance_force(v) for v in v_range]
speed_interval = 0.1  # 速度间隔 (km/h)
max_speed_iterative = calculate_max_speed_iterative(gear_ratios, n_max, r, i0, eta_t, f, CdA, rho, G, v_min, v_max, resistance_force_values)
print(f"\n迭代法计算最高车速: {max_speed_iterative:.2f} km/h")

# 计算最大爬坡度
max_grade = calculate_max_grade(gear_ratios, n_max_torque, r, i0, eta_t, f, G)
print(f"理论最大爬坡度: {max_grade:.2f} %")

# 爬坡度对应附着率检查
theta = np.arctan(max_grade / 100)  # 坡度角
ig_1 = gear_ratios[0]
Tq_max = engine_torque(n_max_torque)
F_traction_max = (Tq_max * ig_1 * i0 * eta_t) / r
adhesion_rate = calculate_adhesion_rate(m_total, g, a, h_g, L, theta, F_traction_max)
print(f"最大爬坡度对应附着率: {adhesion_rate:.2f}")

# 使用积分计算加速时间
time, error = integrate.quad(inverse_acceleration, v_min_2, 70)
print(f"\n用2档起步加速行驶至70km/h的加速时间 (积分): {time:.2f} s")
plt.tight_layout()
plt.show()
