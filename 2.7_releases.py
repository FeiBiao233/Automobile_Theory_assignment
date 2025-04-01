import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import pandas as pd

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

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
gear_ratios = [5.56, 2.769, 1.644, 1.00, 0.793]  # 五档变速器传动比
num_gears=5 #档位数
gear_top = gear_ratios[-1]      # 最高档传动比
gear_second = gear_ratios[-2]   # 次高档传动比

# 转动惯量
I_f = 0.218  # 飞轮转动惯量 (kg*m^2)
I_w1 = 1.798  # 前轮转动惯量 (kg*m^2)
I_w2 = 3.598  # 后轮转动惯量 (kg*m^2)

#几何参数
L = 3.2 #轴距
a = 1.947 #质心至前轴距离（满载）
h_g = 0.9 #质心高

#其他参数
rho = 1.225  # 空气密度 (kg/m^3) 在标准大气条件下
Qid = 0.299         # 怠速油耗 (mL/s)
fuel_density = 0.74 # 燃油密度 (g/mL)
g = 9.81 #重力加速度
G = m_total * g #重力

# ============================== 发动机特性参数 ==============================
engine_coeffs = pd.DataFrame({
    'n':    [815,   1207,   1614,   2012,   2603,   3006,   3403,   3804],
    'B0': [1326.8, 1354.7, 1284.4, 1222.9, 1141.0, 1051.2, 1233.9, 1129.7],
    'B1': [-416.46,-303.98,-189.75,-121.59,-98.893,-73.714,-84.478,-45.291],
    'B2': [72.379, 36.657, 14.524, 7.0035,4.4763, 2.8593, 2.9788, 0.71113],
    'B3': [-5.8629,-2.0553,-0.51184,-0.18517,-0.091077,-0.05138,-0.047449,-0.00075215],
    'B4': [0.17768,0.043072,0.0068164,0.0018555,0.00068906,0.00035032,0.00028230,-0.000038568]
})

# ============================== 核心计算函数 ==============================
def interpolate_coefficients(n):
    """根据转速n插值获取B0-B4系数"""
    n_values = engine_coeffs['n'].values
    if n < n_values[0] or n > n_values[-1]:
        return None
    idx = np.searchsorted(n_values, n)
    idx = max(1, min(idx, len(n_values)-1))
    n_low, n_high = n_values[idx-1], n_values[idx]
    w_low = (n_high - n) / (n_high - n_low)
    w_high = (n - n_low) / (n_high - n_low)
    coeffs = w_low * engine_coeffs.iloc[idx-1, 1:] + w_high * engine_coeffs.iloc[idx, 1:]
    return coeffs.values

# 考虑转动惯量的牛二
def cal_delta_m(gear):
    return m_total + (2*I_w1 + 2*I_w2 + I_f *i0 **2 *gear_ratios[gear-1] **2 )/(r**2)

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

#速度(km/h)计算发动机转速 (r/min)
def cal_rpm(v_kmh, gear):
    v = v_kmh /3.6  # km/h → m/s
    wheel_rps = v / (2 * np.pi * r)  # 车轮转/秒
    return wheel_rps * gear_ratios[gear-1] * i0 * 60  # 发动机转速

def cal_vkmh(n, gear): #计算车速
    n_scaled = n / (i0 * gear_ratios[gear-1])
    v = (2 * np.pi * r) * n_scaled /60
    v_kmh=v*3.6 
    return v_kmh

def engine_torque(n):
    """n: 发动机转速 (r/min)"""
    n_scaled = n / 1000
    Tq = -19.313 + 295.27 * n_scaled - 165.44 * n_scaled**2 + 40.874 * n_scaled**3 - 3.8445 * n_scaled**4
    return Tq

def cal_resist_power(v_kmh, acceleration, slope_angle,gear): #计算发动机功率需求 (kW)
    v = v_kmh /3.6  # m/s
    F_roll = G * f * np.cos(slope_angle)
    F_air = 0.5 * rho * CdA * v**2  # 修改后的空气阻力计算
    F_accel = cal_delta_m(gear) * acceleration
    F_gradient = G * np.sin(slope_angle)
    F_total = F_roll + F_air + F_accel + F_gradient
    return (F_total * v) / (eta_t * 1000)

def calculate_fuel_rate(n, Pe):
    # 计算燃油消耗率 (mL/s)
    if Pe <= 0 or n < n_min:
        return Qid
    coeffs = interpolate_coefficients(n)
    if coeffs is None:
        return Qid
    B0, B1, B2, B3, B4 = coeffs
    b = B0 + B1*Pe + B2*Pe**2 + B3*Pe**3 + B4*Pe**4  # g/(kW·h)
    return (b * Pe) / (3600 * fuel_density)

def cal_fuel_time(t,v_start,accel): #建立在6种工况中油耗(mL/s)与时间的函数
    v=v_start+accel*t
    n=cal_rpm(v,gear=3)
    Pe=cal_resist_power(v,accel,0,gear=3)
    fuel_rate=calculate_fuel_rate(n,Pe) #油耗(mL/s)
    return fuel_rate/1000 # 转换为L/s

# ============================== 六工况油耗计算 ==============================
six_conditions = [
    {'type': '匀速', 'duration':7.2,   'v_start':25, 'v_end':25, 'accel':0},
    {'type': '匀加速', 'duration':16.7,  'v_start':25, 'v_end':40, 'accel':0.25},
    {'type': '匀速', 'duration':22.5,  'v_start':40, 'v_end':40, 'accel':0},
    {'type': '匀加速', 'duration':14.0,  'v_start':40, 'v_end':50, 'accel':0.20},
    {'type': '匀速', 'duration':18.0,  'v_start':50, 'v_end':50, 'accel':0},
    {'type': '匀减速', 'duration':19.3,  'v_start':50, 'v_end':25, 'accel':-0.36},
]

def simulate_six_cycles():
    total_fuel=0  # 总油耗 (mL)
    for cond in six_conditions:
        duration = cond['duration']
        accel = cond['accel']
        v_start = cond['v_start']
        v_end = cond['v_end']
        fuel=integrate.quad(cal_fuel_time, 0, duration, args=(v_start,accel))[0]
        total_fuel+=fuel
    return total_fuel

# ============================== 功率平衡图 ==============================
def plot_power_balance():
    plt.figure(figsize=(12, 8))
    v_range = np.linspace(0, 120, 1000)  # 调整车速范围
    n_range = np.linspace(n_min, n_max, 1000)  # 发动机转速范围
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # 各档位驱动力功率
    for idx in range(num_gears):
        ig=gear_ratios[idx]  # 当前档位传动比
        v_range_temp = cal_vkmh(n_range, idx+1)  # 计算车速范围
        Pe = [traction_force(n,gear_ratios[idx])*cal_vkmh(n,idx+1)/3.6/1000 for n in n_range]  # 计算功率
        plt.plot(v_range_temp, Pe, label=f'{idx+1}档 (ig={ig:.2f})', linewidth=2)
    
    # 行驶阻力功率
    resist_power = [(0.5*rho*CdA*(v/3.6)**3 + G*f*(v/3.6))/1000 for v in v_range]
    plt.plot(v_range, resist_power, label='行驶阻力', linewidth=3, linestyle='--', color='#2c3e50')
    plt.title('货车功率平衡图', fontsize=16)
    plt.xlabel('车速 (km/h)', fontsize=12)
    plt.ylabel('功率 (kW)', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.xlim(0, 120)
    plt.ylim(0, 80)
    plt.show()

# ============================== 等速油耗曲线 ==============================
def plot_constant_speed_consumption():
    plt.figure(figsize=(10, 6))  # 首先创建 figure

    # 计算四档
    v_range_4 = np.linspace(20, 90, 1000)
    fuel_consumption_4 = []
    for v in v_range_4:
        n = cal_rpm(v, 4)
        Pe = cal_resist_power(v, 0, 0, 4)
        fuel_rate = calculate_fuel_rate(n, Pe)
        #  正确的百公里油耗计算方法
        fuel_ml_per_km = (fuel_rate * 3600) *100 /v
        fuel_consumption_4.append(fuel_ml_per_km / 1000)  # L/100km
    plt.plot(v_range_4, fuel_consumption_4, label=f'4档 (ig={gear_ratios[4-1]:.2f})')

    # 计算五档
    v_range_5 = np.linspace(20, 120, 1000)
    fuel_consumption_5 = []
    for v in v_range_5:
        n = cal_rpm(v, 5)
        Pe = cal_resist_power(v, 0, 0, 5)
        fuel_rate = calculate_fuel_rate(n, Pe)
        #  正确的百公里油耗计算方法
        fuel_ml_per_km = (fuel_rate * 3600) *100 /v
        fuel_consumption_5.append(fuel_ml_per_km / 1000)  # L/100km
    plt.plot(v_range_5, fuel_consumption_5, label=f'5档 (ig={gear_ratios[5-1]:.2f})')

    plt.xlabel('车速 (km/h)')
    plt.ylabel('油耗 (L/100km)')
    plt.title('等速百公里油耗曲线')
    plt.legend()
    plt.grid()
    plt.show()


# ============================== 主函数 ==============================
if __name__ == '__main__':
    #计算六工况油耗
    total_fuel_100km = simulate_six_cycles()/1.075*100 # 将总油耗转换为百公里油耗
    print(f"六工况循环百公里油耗: {total_fuel_100km:.2f} L/100km")

    # 绘制功率平衡图
    plot_power_balance()

    # 绘制等速油耗曲线
    plot_constant_speed_consumption()