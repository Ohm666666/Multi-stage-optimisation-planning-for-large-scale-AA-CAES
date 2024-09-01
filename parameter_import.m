function parameter_import
%% excel文件写入数据
filename = 'D:\文章-大规模CAES多阶段优化规划-程序\0.公用参数\文章所用数据V20240822.xlsx';
par_G = readmatrix(filename,'Sheet','火电机组参数设置','Range','B2:P15');               %从EXCEL导入数据--火电机组参数
par_W = readmatrix(filename,'Sheet','风电场参数设置','Range','B2:B7');                  %从EXCEL导入数据--风电场参数
par_P = readmatrix(filename,'Sheet','光伏电站参数设置','Range','B12:B17');              %从EXCEL导入数据--光伏电站参数
par_CAES = readmatrix(filename,'Sheet','多容量等级CAES参数设置','Range','C3:F44');      %从EXCEL导入数据--CAES参数
par_sys = readmatrix(filename,'Sheet','系统参数设置','Range','B2:B15');                 %从EXCEL导入数据--系统参数
par_grow = readmatrix(filename,'Sheet','系统参数设置','Range','F2:G6');                 %从EXCEL导入数据--源荷增长比例

%% 火电机组数据处理
P_G_max = par_G(1,:)';          %火电机组--额定出力
P_G_min = par_G(2,:)';          %火电机组--最小技术出力
s_G = par_G(3,:)';              %火电机组--启停成本
sst_G = par_G(4,:)';            %火电机组--最小启停机时间--单位h
ramp_G = par_G(5,:)';           %火电机组--爬坡率
a_G = par_G(6,:)';              %火电机组--购电成本系数a
b_G = par_G(7,:)';              %火电机组--购电成本系数b
c_G = par_G(8,:)';              %火电机组--正备用成本系数c
d_G = par_G(9,:)';              %火电机组--负备用成本系数d
adj_G = par_G(10,:)';           %火电机组--调差系数
kp_G = par_G(11,:)';            %火电机组--调速器比例增益系数
ki_G = par_G(12,:)';            %火电机组--调速器积分增益系数
H_G = par_G(13,:)';             %火电机组--惯性时间常数
r_G_fm = par_G(14,:)';          %火电机组--限幅比例

%% 风电场--光伏电站数据处理
P_WS_N = par_W(1);               %风电场--装机容量
P_PV_N = par_P(1);               %光伏电站--装机容量
c_wind_waste_pen = par_W(2);     %弃风惩罚成本系数
K_WS_1 = par_W(4)/50;            %风电场--一次调频系数
K_WS_2 = par_W(3)/50;            %风电场--惯量响应系数
r_WS_fm_up = par_W(5);           %风电场--正向一次调频限幅系数
r_WS_fm_down = par_W(6);         %风电场--负向一次调频限幅系数

%% CAES电站数据处理
P_CAES_cN = par_CAES(1,:)';          %CAES电站--压缩机--额定功率MW
P_CAES_c_min = par_CAES(2,:)';       %CAES电站--压缩机--最小出力比%
yita_CAES_c = par_CAES(5,:)';        %CAES电站--压缩机--额定效率%
P_CAES_dN = par_CAES(7,:)';          %CAES电站--膨胀机--额定功率MW
P_CAES_d_min = par_CAES(8,:)';       %CAES电站--膨胀机--最小出力比%
yita_CAES_d = par_CAES(11,:)';       %CAES电站--膨胀机--额定效率%
V_CAES_N = par_CAES(13,:)';          %CAES电站--储气室--体积(万m^3)
p_CAES_max = par_CAES(14,:)';        %CAES电站--储气室--最大气压
p_CAES_min = par_CAES(15,:)';        %CAES电站--储气室--最小气压
p_CAES_0 = par_CAES(16,:)';          %CAES电站--储气室--初始气压
kc_1 = par_CAES(25,:)';              %CAES电站--压缩系数--1高效率段
kc_2 = par_CAES(26,:)';              %CAES电站--压缩系数--2中效率段
kc_3 = par_CAES(27,:)';              %CAES电站--压缩系数--3低效率段
kd_1 = par_CAES(28,:)';              %CAES电站--发电系数--1高效率段
kd_2 = par_CAES(29,:)';              %CAES电站--发电系数--2中效率段
kd_3 = par_CAES(30,:)';              %CAES电站--发电系数--3低效率段
adj_CAES = par_CAES(32,:)';          %CAES电站--调差系数
kp_CAES = par_CAES(33,:)';           %CAES电站--调速器比例增益系数
ki_CAES = par_CAES(34,:)';           %CAES电站--调速器积分增益系数
H_CAES = par_CAES(35,:)';            %CAES电站--惯性时间常数
r_CAES_fm = par_CAES(36,:)';         %CAES电站--限幅比例

c_CAES_c = par_CAES(37,:)';          %CAES电站--压缩机--建设成本系数 ($/MW)
c_CAES_d = par_CAES(38,:)';          %CAES电站--膨胀机--建设成本系数 ($/MW)
c_CAES_h = par_CAES(39,:)';          %CAES电站--储气室--建设成本系数 ($/万m^3)
c_CAES_OM_fix = par_CAES(40,:)';     %CAES电站--单位固定运维成本$/(MW·a)
c_CAES_OM_var = par_CAES(41,:)';     %CAES电站--单位可变运维成本$/(MW)
T_CAES_span = par_CAES(42,:)';       %CAES电站--期望寿命

%%%%%待定%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_CAES = 21.5;                       %CAES电站正备用成本系数c
% d_CAES = 20;                         %CAES电站负备用成本系数d
%%%%%待定%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 系统数据处理
gama = par_sys(1);                   %折现率
delta_sal = par_sys(2);              %净残值率
r_fe_wind = par_sys(3);              %风电预测误差比例
r_fe_load = par_sys(4);              %负荷预测误差比例
df_steady_lim = par_sys(5)/50;       %稳态频差上限标幺值(0.2Hz)
df_max_lim = par_sys(6)/50;          %最大频差上限标幺值(0.5Hz)
vf_lim = par_sys(7)/50;              %频率变化率上限标幺值(0.35Hz/s)
df_deadzone = par_sys(8)/50;         %系统调频死区标幺值(0.03Hz)
k_load = par_sys(9);                 %负荷频率响应系数
kai = par_sys(10:13);                %典型日一年中天数占比
T_all = par_sys(14);                 %规划工程期限

%% 源荷增长数据处理
r_wind_grow = par_grow(:,1);         %新能源增长比例
r_load_grow = par_grow(:,2);         %负荷增长比例

%% 数据保存
save parameter.mat
end