function [Cost_1,P_CAES_c_sum] = Underlevel_CAES(x,st,ss)    
%下层日前优化调度程序：场景A
% V20240412创建;V20240418开始编写;V20240421暂不考虑详细版CAES备用约束，暂不考虑储热约束(默认储热能力不受限);
% 输入参数x：配置方案(其中每一个元素表示：某一时段某一地理点接入某一类型CAES的个数),st：规划阶段,ss：典型日标志
% 输出参数Cost：该配置方案下+该阶段+该典型日条件下：计算得到的日前调度运行成本,P_CAES：CAES充电容量
yalmip('clear')
%% 参数导入
load parameter.mat                                     %导入系统参数
load typical_day.mat                                   %导入典型日参数

%% 数据选择与处理
P_TD = typical_day(ss).d;                                           %选择典型日
P_wind_f = r_wind_grow(st).*P_TD(:,1:3)';                           %计及源荷增长比率--风电场出力预测数据
P_load_f = r_load_grow(st).*P_TD(:,4)';                             %计及源荷增长比率--负荷预测数据
%生成源荷预测误差
P_wind_real = P_wind_f;P_wind_real(P_wind_real~=0)=1;
P_wind_up = P_wind_real.*(40-st*5 + 0.5.*r_fe_wind.*P_wind_f);      %风电预测上波动[N_wind*Horizon]
P_wind_up = min(P_wind_up,r_wind_grow(st).*[P_WS_N;P_WS_N;P_PV_N]-P_wind_f);         %不越限
P_wind_down = P_wind_real.*(40-st*5 + 0.5.*r_fe_wind.*P_wind_f);    %风电预测下波动[N_wind*Horizon]
P_wind_down = min(P_wind_down,P_wind_f);                            %不越限
P_load_up = 90-st*5 + 0.5.*r_fe_load.*P_load_f;                     %日前负荷预测数据-上波动
P_load_down = 100-st*5 + 0.5.*r_fe_load.*P_load_f;                  %日前负荷预测数据-下波动
%生成净负荷及其预测误差
P_nl = P_load_f - sum(P_wind_f,1);                                  %净负荷--功率值    
P_nl_up = P_load_down + sum(P_wind_up,1);                           %净负荷--上预测误差
P_nl_dn = P_load_up + sum(P_wind_down,1);                           %净负荷--下预测误差
%生成系统灵活爬坡调节需求
R_nl_up = max(0,(min((P_nl(2:end)-P_nl(1:end-1)),0) + (P_nl_up(2:end) + P_nl_dn(1:end-1))));
R_nl_dn = max(0,(min(P_nl(1:end-1)-(P_nl(2:end)),0) + (P_nl_up(1:end-1) + P_nl_dn(2:end))));
R_nl_up = [R_nl_up,P_nl_up(end)];
R_nl_dn = [R_nl_dn,P_nl_dn(end)];

%生成预想功率扰动
temp_0 = [(P_load_up + sum(P_wind_down,1));( ...
    P_load_down + sum(P_wind_up,1))];
P_LW_error = max(temp_0);                                           %系统可能发生的预想功率扰动

%根据规划阶段选择系统内各地理点已建设的CAES电站个数x_st=[4*9],V20240419仅允许每一节点只能配置一台机组
if st == 1
    x_st = x(1:4,:);                           %1阶段
elseif  st == 2
    x_st = x(1:4,:) + x(5:8,:);                %2阶段
elseif  st == 3
    x_st = x(1:4,:) + x(5:8,:) + x(9:12,:);    %3阶段
else
    x_st = x(1:4,:) + x(5:8,:) + x(9:12,:) + x(13:16,:);    %其他阶段
end
x_real = double((1 - sum(x_st)==0)');
%% 参数设置
M = 1000000000;                                             %大M
Horizon = 24;                                               %调度点数
t_60 = 60;                                                  %明确调度步长(60min)
delta_t = 1;                                                %明确调度步长(60min)
t_15 = 15;                                                  %明确事故备用时长(15min)
t_5 = 5;                                                    %明确灵活爬坡时长(5min)

N_G = 15;                                                   %火电机组--数量
N_WS = 3;                                                   %风电场--（10聚合1）数量
N_CAES = nnz(x_real);                                       %CAES电站--数量 V20240419仅允许每一节点只能配置一台机组,否则需要！增加！节点数
%CAES参数处理--去除置零点
P_CAES_c_N = sum(x_st.*P_CAES_cN)';                         %CAES电站--压缩机--额定功率MW[N_CAES*1]
P_CAES_c_N = Rz(P_CAES_c_N);
P_CAES_d_N = sum(x_st.*P_CAES_dN)';                         %CAES电站--膨胀机--额定功率MW[N_CAES*1]
P_CAES_d_N = Rz(P_CAES_d_N);
CAES_VN = sum(x_st.*V_CAES_N)';                             %CAES电站--储气室体积[N_CAES*1]
CAES_VN = Rz(CAES_VN);
yita_CAES_c = sum(x_st.*yita_CAES_c);                       %CAES电站--压缩机--额定效率%
yita_CAES_c = Rz(yita_CAES_c);
yita_CAES_d = sum(x_st.*yita_CAES_d);                       %CAES电站--膨胀机--额定效率%
yita_CAES_d = Rz(yita_CAES_d);
P_CAES_c_min = sum(x_st.*P_CAES_c_min.*P_CAES_cN)';         %CAES电站--压缩机--最小出力MW[N_CAES*1]
P_CAES_c_min = Rz(P_CAES_c_min);
P_CAES_d_min = sum(x_st.*P_CAES_d_min.*P_CAES_dN)';         %CAES电站--膨胀机--最小出力MW[N_CAES*1]
P_CAES_d_min = Rz(P_CAES_d_min);

p_CAES_max = sum(x_st.*p_CAES_max)';                        %CAES电站--储气室--最大气压
p_CAES_max = Rz(p_CAES_max);
p_CAES_min = sum(x_st.*p_CAES_min)';                        %CAES电站--储气室--最小气压
p_CAES_min = Rz(p_CAES_min);
p_CAES_0 = sum(x_st.*p_CAES_0)';                            %CAES电站--储气室--初始气压
p_CAES_0 = Rz(p_CAES_0);
kc_1 = sum(x_st.*kc_1)';                                    %CAES电站--压缩系数--1高效率段
kc_1 = Rz(kc_1);
kc_2 = sum(x_st.*kc_2)';                                    %CAES电站--压缩系数--2中效率段
kc_2 = Rz(kc_2);
kc_3 = sum(x_st.*kc_3)';                                    %CAES电站--压缩系数--3低效率段
kc_3 = Rz(kc_3);
kd_1 = sum(x_st.*kd_1)';                                    %CAES电站--发电系数--1高效率段
kd_1 = Rz(kd_1);
kd_2 = sum(x_st.*kd_2)';                                    %CAES电站--发电系数--2中效率段
kd_2 = Rz(kd_2);
kd_3 = sum(x_st.*kd_3)';                                    %CAES电站--发电系数--3低效率段
kd_3 = Rz(kd_3);

adj_CAES = sum(x_st.*adj_CAES)';                            %CAES电站--调差系数
adj_CAES = Rz(adj_CAES);
kp_CAES = sum(x_st.*kp_CAES)';                              %CAES电站--调速器比例增益系数
kp_CAES = Rz(kp_CAES);
ki_CAES = sum(x_st.*ki_CAES)';                              %CAES电站--调速器积分增益系数
ki_CAES = Rz(ki_CAES);
H_CAES = sum(x_st.*H_CAES)';                                %CAES电站--惯性时间常数
H_CAES = Rz(H_CAES);
r_CAES_fm = sum(x_st.*r_CAES_fm)';                          %CAES电站--限幅比例
r_CAES_fm = Rz(r_CAES_fm);

%CAES变工况参数处理
P_CAES_c_max = P_CAES_c_N;
P_CAES_c_max_1 = P_CAES_c_max;                              %CAES电站--压缩工况下--1高效率段--功率上限
P_CAES_c_min_1 = (80.01*10^-2).*P_CAES_c_N;                 %CAES电站--压缩工况下--1高效率段--功率下限
P_CAES_c_max_2 = (80.00*10^-2).*P_CAES_c_N;                 %CAES电站--压缩工况下--2中效率段--功率上限
P_CAES_c_min_2 = (60.01*10^-2).*P_CAES_c_N;                 %CAES电站--压缩工况下--2中效率段--功率下限
P_CAES_c_max_3 = (60.00*10^-2).*P_CAES_c_N;                 %CAES电站--压缩工况下--3低效率段--功率上限
P_CAES_c_min_3 = P_CAES_c_min;                              %CAES电站--压缩工况下--3低效率段--功率下限

P_CAES_d_max = P_CAES_d_N;
P_CAES_d_max_1 = P_CAES_d_max;                              %CAES电站--发电工况下--1高效率段--功率上限
P_CAES_d_min_1 = (80.01*10^-2).*P_CAES_d_N;                 %CAES电站--发电工况下--1高效率段--功率下限
P_CAES_d_max_2 = (80.00*10^-2).*P_CAES_d_N;                 %CAES电站--发电工况下--2中效率段--功率上限
P_CAES_d_min_2 = (60.01*10^-2).*P_CAES_d_N;                 %CAES电站--发电工况下--2中效率段--功率下限
P_CAES_d_max_3 = (60.00*10^-2).*P_CAES_d_N;                 %CAES电站--发电工况下--3低效率段--功率上限
P_CAES_d_min_3 = P_CAES_d_min;                              %CAES电站--发电工况下--3低效率段--功率下限

c_CAES_run = repmat(x_real.*3.*[1.50;1.45;1.40;1.35;1.30;1.25;1.20;1.15;1.10],1,Horizon);         %计算CAES组运行成本
c_CAES_run(all(c_CAES_run == 0,2),:)=[];
c_CAES_res = repmat(x_real.*2.*[1.50;1.45;1.40;1.35;1.30;1.25;1.20;1.15;1.10],1,Horizon);
c_CAES_res(all(c_CAES_res == 0,2),:)=[];

lambda_1 = 1*10^-4;                                         %系统--缺电率限值
lambda_2 = 5*10^-2;                                         %系统--消纳率限值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%待修改%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
miu_p = 0.05;                               %储气室始末气压最大允许偏差
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%待修改%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mar_G = min(repmat(P_G_max.*r_G_fm,1,Horizon), ...
        repmat(P_G_max./adj_G*df_steady_lim,1,Horizon));     %火电机组极限功率扰动情况下一次调频变化量[N_G * Horizon]——确定值

%% 变量定义
% 常规机组变量
u_G = binvar(N_G,Horizon,'full');                           %火电机组--运行标志变量，取值为1则为发电状态，取值为0则为停机状态
P_G = sdpvar(N_G,Horizon,'full');                           %火电机组--有功出力功率
R_G_up = sdpvar(N_G,Horizon,'full');                        %火电机组--正爬坡
R_G_down = sdpvar(N_G,Horizon,'full');                      %火电机组--负爬坡
R_G_em = sdpvar(N_G,Horizon,'full');                        %火电机组--事故备用

% CAES电站变量
P_CAES_c = sdpvar(N_CAES,Horizon,'full');                   %CAES电站--压缩工况--有功用电功率
P_CAES_d = sdpvar(N_CAES,Horizon,'full');                   %CAES电站--发电工况--有功出力功率
u_CAES_c = binvar(N_CAES,Horizon,'full');                   %CAES电站--压缩工况--标志变量,取值为1则为压缩工况
u_CAES_d = binvar(N_CAES,Horizon,'full');                   %CAES电站--发电工况--标志变量,取值为1则为发电工况
u_CAES_s = binvar(N_CAES,Horizon,'full');                   %CAES电站--待机工况--标志变量,取值为1则为待机工况

P_CAES_c_1 = sdpvar(N_CAES,Horizon,'full');                 %CAES电站--压缩工况--有功用电功率--1高效率段
P_CAES_c_2 = sdpvar(N_CAES,Horizon,'full');                 %CAES电站--压缩工况--有功用电功率--2中效率段
P_CAES_c_3 = sdpvar(N_CAES,Horizon,'full');                 %CAES电站--压缩工况--有功用电功率--3低效率段
P_CAES_d_1 = sdpvar(N_CAES,Horizon,'full');                 %CAES电站--发电工况--有功出力功率--1高效率段
P_CAES_d_2 = sdpvar(N_CAES,Horizon,'full');                 %CAES电站--发电工况--有功用电功率--2中效率段
P_CAES_d_3 = sdpvar(N_CAES,Horizon,'full');                 %CAES电站--发电工况--有功用电功率--3低效率段

u_CAES_c_1 = binvar(N_CAES,Horizon,'full');                 %CAES电站--压缩工况--标志变量--1高效率段
u_CAES_c_2 = binvar(N_CAES,Horizon,'full');                 %CAES电站--压缩工况--标志变量--2中效率段
u_CAES_c_3 = binvar(N_CAES,Horizon,'full');                 %CAES电站--压缩工况--标志变量--3低效率段
u_CAES_d_1 = binvar(N_CAES,Horizon,'full');                 %CAES电站--发电工况--标志变量--1高效率段
u_CAES_d_2 = binvar(N_CAES,Horizon,'full');                 %CAES电站--发电工况--标志变量--2中效率段
u_CAES_d_3 = binvar(N_CAES,Horizon,'full');                 %CAES电站--发电工况--标志变量--3低效率段

CAES_pst = sdpvar(N_CAES,Horizon,'full');                   %CAES电站--储气室--气压变量,类似SOC
% H_CAES_tes = sdpvar(N_CAES,Horizon,'full');                 %CAES电站--储热室--热量变量,类似SOC  V20240421不考虑储热

R_CAES_up_c = sdpvar(N_CAES,Horizon,'full');                %CAES--压缩工况--正爬坡
R_CAES_down_c = sdpvar(N_CAES,Horizon,'full');              %CAES--压缩工况--负爬坡
R_CAES_up_d = sdpvar(N_CAES,Horizon,'full');                %CAES--发电工况--正爬坡
R_CAES_down_d = sdpvar(N_CAES,Horizon,'full');              %CAES--发电工况--负爬坡

R_CAES_em_c = sdpvar(N_CAES,Horizon,'full');                %CAES--压缩工况--事故备用
R_CAES_em_s = sdpvar(N_CAES,Horizon,'full');                %CAES--待机工况--事故备用
R_CAES_em_d = sdpvar(N_CAES,Horizon,'full');                %CAES--发电工况--事故备用

% 风电变量
P_wind = sdpvar(N_WS,Horizon,'full');                       %风电场--实际出力
% Mar_W = sdpvar(N_WS,Horizon,'full');                      %风电参与一次调频的变化量限值
P_load = sdpvar(1,Horizon,'full');                          %负荷--实际供电

%% 目标函数 
C_G_run = sum(sum(repmat(a_G,1,Horizon).*P_G*delta_t+repmat(b_G,1,Horizon).*u_G));                                                %计算火电运行成本
C_G_start = sum(sum(repmat(s_G,1,Horizon-1).*(u_G(:,2:Horizon)-u_G(:,1:Horizon-1)+abs(u_G(:,2:Horizon)-u_G(:,1:Horizon-1)))/2));  %计算火电开机成本
C_G_stop = sum(sum(repmat(s_G,1,Horizon-1).*(u_G(:,1:Horizon-1)-u_G(:,2:Horizon)+abs(u_G(:,1:Horizon-1)-u_G(:,2:Horizon)))/2));   %计算火电关机成本
C_G_re = sum(sum(repmat(c_G,1,Horizon).*R_G_em))*delta_t;                                                                         %计算火电事故备用成本
C_G_rl = sum(sum(repmat(c_G,1,Horizon).*R_G_up+repmat(d_G,1,Horizon).*R_G_down))*delta_t;                                         %计算火电灵活爬坡成本
C_wind = sum(c_wind_waste_pen*sum(P_wind_f - P_wind))*delta_t;                                                                    %计算弃风成本 
C_load = sum(2*c_wind_waste_pen*sum(P_load_f - P_load))*delta_t;                                                                  %计算缺电惩罚成本 
C_CAES_run = sum(sum(c_CAES_run.*(P_CAES_d+P_CAES_c)));                                                                           %计算CAES组运行成本
C_CAES_rl = sum(sum(c_CAES_res.*(R_CAES_em_c+R_CAES_em_d+R_CAES_em_s)));                                                          %计算CAES组事故备用成本
C_CAES_re = sum(sum(c_CAES_res.*(R_CAES_up_c+R_CAES_up_d+R_CAES_down_c+R_CAES_down_d)));                                          %计算CAES组灵活爬坡成本

% 1
Cst=[];
%% 约束条件--常规机组
%常规机组出力约束
Cst=[Cst,u_G.*repmat(P_G_min,1,Horizon) <= P_G];
Cst=[Cst,P_G <= u_G.*repmat(P_G_max,1,Horizon)];
%常规机组爬坡率约束
Cst=[Cst,-t_60*repmat(ramp_G,1,Horizon-1) - (1-u_G(:,2:Horizon)+1-u_G(:,1:Horizon-1))*M <= P_G(:,2:Horizon)-P_G(:,1:Horizon-1)];
Cst=[Cst,P_G(:,2:Horizon)-P_G(:,1:Horizon-1) <= t_60*repmat(ramp_G,1,Horizon-1) + (1-u_G(:,2:Horizon)+1-u_G(:,1:Horizon-1))*M];
%常规机组事故备用约束
Cst=[Cst,R_G_em <= t_15*repmat(ramp_G,1,Horizon).*u_G];           %事故备用-爬坡率约束
Cst=[Cst,P_G + R_G_em <= repmat(P_G_max,1,Horizon).*u_G];         %事故备用-功率约束

%常规机组灵活爬坡约束
Cst=[Cst,R_G_up <= t_5*repmat(ramp_G,1,Horizon).*u_G];            %正灵活爬坡-爬坡率约束
Cst=[Cst,R_G_down <= t_5*repmat(ramp_G,1,Horizon).*u_G];          %负灵活爬坡-爬坡率约束
Cst=[Cst,P_G + R_G_up <= repmat(P_G_max,1,Horizon).*u_G];         %正灵活爬坡-功率约束
Cst=[Cst,P_G - R_G_down >= repmat(P_G_min,1,Horizon).*u_G];       %负灵活爬坡-功率约束
%常规机组一次调频约束
Cst=[Cst,P_G + u_G.*Mar_G <= repmat(P_G_max,1,Horizon).*u_G];
Cst=[Cst,P_G - u_G.*Mar_G >= repmat(P_G_min,1,Horizon).*u_G];
%常规机组启停机时间约束
for n = 1: N_G
    for t = 1: Horizon-sst_G(n)
        Cst=[Cst,sum(u_G(n , t+1:t+sst_G(n)))>=sst_G(n)*(u_G(n,t+1)-u_G(n,t))];
        Cst=[Cst,sum(1-u_G(n , t+1:t+sst_G(n)))>=sst_G(n)*(u_G(n,t)-u_G(n,t+1))];
    end
    for t = Horizon-sst_G(n)+1: Horizon-1
        Cst=[Cst,sum(u_G(n , t+1:Horizon))+sum(u_G(n , 1:t+sst_G(n)-Horizon))>=sst_G(n)*(u_G(n,t+1)-u_G(n,t))];
        Cst=[Cst,sum(1-u_G(n , t+1:Horizon))+sum(1-u_G(n , 1:t+sst_G(n)-Horizon))>=sst_G(n)*(u_G(n,t)-u_G(n,t+1))];
    end
    Cst=[Cst,sum(u_G(n,1:sst_G(n)))>=sst_G(n)*(u_G(n,1)-u_G(n,Horizon))];
    Cst=[Cst,sum(1-u_G(n,1:sst_G(n)))>=sst_G(n)*(u_G(n,Horizon)-u_G(n,1))];
end
Cst=[Cst,[R_G_em;R_G_up;R_G_down] >=0];
%% 约束条件--风电场
Cst=[Cst,P_wind <= P_wind_f];           %风电出力约束
Cst=[Cst,P_wind >= 0];                  %风电出力约束
% 2
%% 约束条件--CAES电站
%CAES电站--考虑变工况效率--总出力/总状态与分出力/分状态间的关系约束
Cst=[Cst,P_CAES_c == (P_CAES_c_1 + P_CAES_c_2 + P_CAES_c_3)];
Cst=[Cst,P_CAES_d == (P_CAES_d_1 + P_CAES_d_2 + P_CAES_d_3)];
Cst=[Cst,u_CAES_c == (u_CAES_c_1 + u_CAES_c_2 + u_CAES_c_3)];
Cst=[Cst,u_CAES_d == (u_CAES_d_1 + u_CAES_d_2 + u_CAES_d_3)];
Cst=[Cst,u_CAES_c_1 + u_CAES_c_2 + u_CAES_c_3 <= 1];
Cst=[Cst,u_CAES_d_1 + u_CAES_d_2 + u_CAES_d_3 <= 1];
Cst=[Cst,u_CAES_c + u_CAES_d <= 1];
%CAES电站--待机工况
Cst=[Cst,u_CAES_s <= (1 - u_CAES_c - u_CAES_d)];
%CAES电站--出力上下限约束
Cst=[Cst,u_CAES_c_1.*repmat(P_CAES_c_min_1,1,Horizon) <= P_CAES_c_1];
Cst=[Cst,P_CAES_c_1 <= u_CAES_c_1.*repmat(P_CAES_c_max_1,1,Horizon)];
Cst=[Cst,u_CAES_c_2.*repmat(P_CAES_c_min_2,1,Horizon) <= P_CAES_c_2];
Cst=[Cst,P_CAES_c_2 <= u_CAES_c_2.*repmat(P_CAES_c_max_2,1,Horizon)];
Cst=[Cst,u_CAES_c_3.*repmat(P_CAES_c_min_3,1,Horizon) <= P_CAES_c_3];
Cst=[Cst,P_CAES_c_3 <= u_CAES_c_3.*repmat(P_CAES_c_max_3,1,Horizon)];
Cst=[Cst,u_CAES_d_1.*repmat(P_CAES_d_min_1,1,Horizon) <= P_CAES_d_1];
Cst=[Cst,P_CAES_d_1 <= u_CAES_d_1.*repmat(P_CAES_d_max_1,1,Horizon)];
Cst=[Cst,u_CAES_d_2.*repmat(P_CAES_d_min_2,1,Horizon) <= P_CAES_d_2];
Cst=[Cst,P_CAES_d_2 <= u_CAES_d_2.*repmat(P_CAES_d_max_2,1,Horizon)];
Cst=[Cst,u_CAES_d_3.*repmat(P_CAES_d_min_3,1,Horizon) <= P_CAES_d_3];
Cst=[Cst,P_CAES_d_3 <= u_CAES_d_3.*repmat(P_CAES_d_max_3,1,Horizon)];
%CAES电站--气压变化等式约束
% 1时段CAES电站--气压变化等式约束
Cst=[Cst,CAES_pst(:,1) == (p_CAES_0 + kc_1.*P_CAES_c_1(:,1) + kc_2.*P_CAES_c_2(:,1) + ...
    kc_3.*P_CAES_c_3(:,1) - kd_1.*P_CAES_d_1(:,1) - kd_2.*P_CAES_d_2(:,1) - kd_3.*P_CAES_d_3(:,1))];
% 2~24时段CAES电站--气压变化等式约束
for t=2:Horizon
    Cst=[Cst,CAES_pst(:,t) == (CAES_pst(:,t-1)+ kc_1.*P_CAES_c_1(:,t) + kc_2.*P_CAES_c_2(:,t) + ...
        kc_3.*P_CAES_c_3(:,t) - kd_1.*P_CAES_d_1(:,t) - kd_2.*P_CAES_d_2(:,t) - kd_3.*P_CAES_d_3(:,t))];
end
%CAES电站--气压上下限不等式约束
Cst=[Cst,repmat(p_CAES_min,1,Horizon) <= CAES_pst];
Cst=[Cst,CAES_pst <= repmat(p_CAES_max,1,Horizon)];
%CAES电站--气压始末保持不变约束
Cst=[Cst,abs(CAES_pst(:,end) - p_CAES_0) <= miu_p];
%CAES电站--事故备用约束
%V20240521要求压缩工况提供极限事故备用
Cst=[Cst,R_CAES_em_c <= u_CAES_c.*repmat((P_CAES_c_max + P_CAES_d_max),1,Horizon)];      %压缩工况--事故备用
Cst=[Cst,R_CAES_em_c <= P_CAES_c + repmat(P_CAES_d_max,1,Horizon).*u_CAES_c];            %压缩工况--事故备用
Cst=[Cst,R_CAES_em_c <= (CAES_pst - repmat(p_CAES_min,1,Horizon))./(repmat(kd_1,1,Horizon).*t_15/t_60) + ...
    P_CAES_c + (1 - u_CAES_s).*repmat(P_CAES_d_max,1,Horizon)];                          %压缩工况--事故备用

Cst=[Cst,R_CAES_em_s >= u_CAES_s.*repmat((P_CAES_d_min),1,Horizon)];                     %待机工况--事故备用
Cst=[Cst,R_CAES_em_s <= u_CAES_s.*repmat((P_CAES_d_max),1,Horizon)];                     %待机工况--事故备用
Cst=[Cst,R_CAES_em_s <= (CAES_pst - repmat(p_CAES_min,1,Horizon))./(repmat(kd_1,1,Horizon).*t_15/t_60) - ...
    0 + (1 - u_CAES_s).*repmat(P_CAES_d_max,1,Horizon)];                                 %待机工况--事故备用

Cst=[Cst,R_CAES_em_d <= u_CAES_d.*repmat((P_CAES_d_max - P_CAES_d_min),1,Horizon)];      %发电工况--事故备用
Cst=[Cst,R_CAES_em_d <= repmat(P_CAES_d_max,1,Horizon) - P_CAES_d];                      %发电工况--事故备用
Cst=[Cst,R_CAES_em_d <= (CAES_pst - repmat(p_CAES_min,1,Horizon))./(repmat(kd_1,1,Horizon).*t_15/t_60) - ...
    P_CAES_d + (1 - u_CAES_d).*repmat(P_CAES_d_max,1,Horizon)];                          %发电工况--事故备用

%CAES电站--灵活爬坡调节约束
%简化版不允许跨工况提供爬坡调节，也就是说CAES停机没有爬坡调节能力，发电及压缩工况只能继续保持当前工况,略去气压对爬坡调节的影响
%压缩工况
Cst=[Cst,R_CAES_up_c <= u_CAES_c.*repmat((P_CAES_c_max - P_CAES_c_min),1,Horizon)];      %压缩工况--正爬坡
Cst=[Cst,R_CAES_up_c <= P_CAES_c - repmat(P_CAES_c_min,1,Horizon).*u_CAES_c];            %压缩工况--正爬坡
Cst=[Cst,R_CAES_down_c <= repmat(P_CAES_c_max,1,Horizon) - P_CAES_c];                    %压缩工况--负爬坡
Cst=[Cst,R_CAES_down_c <= u_CAES_c.*repmat((P_CAES_c_max - P_CAES_c_min),1,Horizon)];    %压缩工况--负爬坡
Cst=[Cst,R_CAES_down_c <= (repmat(p_CAES_max,1,Horizon) - CAES_pst)./(repmat(kc_1,1,Horizon).*t_5/t_60) - ...
    P_CAES_c + (1 - u_CAES_c).*repmat(P_CAES_c_max,1,Horizon)];                          %压缩工况--负爬坡
%发电工况
Cst=[Cst,R_CAES_up_d <= u_CAES_d.*repmat((P_CAES_d_max - P_CAES_d_min),1,Horizon)];      %发电工况--正爬坡
Cst=[Cst,R_CAES_up_d <= repmat(P_CAES_d_max,1,Horizon) - P_CAES_d];                      %发电工况--正爬坡
Cst=[Cst,R_CAES_up_d <= (CAES_pst - repmat(p_CAES_min,1,Horizon))./(repmat(kd_1,1,Horizon).*t_5/t_60) - ...
    P_CAES_d + (1 - u_CAES_d).*repmat(P_CAES_d_max,1,Horizon)];                          %发电工况--正爬坡
Cst=[Cst,R_CAES_down_d <= u_CAES_d.*repmat((P_CAES_d_max - P_CAES_d_min),1,Horizon)];    %发电工况--负爬坡
Cst=[Cst,R_CAES_down_d <= P_CAES_d - u_CAES_d.*repmat(P_CAES_d_min,1,Horizon)];          %发电工况--负爬坡
%停机工况不配拥有爬坡
Cst=[Cst,[R_CAES_up_c;R_CAES_down_c;R_CAES_up_d;R_CAES_down_d;R_CAES_em_c;R_CAES_em_s;R_CAES_em_d] >= 0];

% 3
%% 约束条件--系统
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%待修改%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%系统--保供促消约束
Cst=[Cst,P_load <= P_load_f];           %负荷供电约束
Cst=[Cst,P_load >= 0];                  %负荷供电约束
%保供电约束
Cst=[Cst,sum(P_load_f-P_load)/sum(P_load_f) <= lambda_1];
%保消纳约束
Cst=[Cst,sum(sum(P_wind_f-P_wind))/sum(sum(P_wind_f)) <= lambda_2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%待补充%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%系统--功率平衡约束
Cst=[Cst,sum(P_G,1) + sum(P_wind,1) + sum(P_CAES_d,1) == sum(P_load,1) + sum(P_CAES_c,1)];
%系统--事故备用约束
Cst=[Cst,sum(R_G_em,1) + sum(R_CAES_em_c,1) + sum(R_CAES_em_s,1) + sum(R_CAES_em_d,1) >= 0.05.*sum(P_load,1)];      %事故备用(V20240521设置为5%负荷) 
%系统--灵活爬坡约束
Cst=[Cst,sum(R_G_up,1) + sum(R_CAES_up_d,1) + sum(R_CAES_up_c,1) >= R_nl_up];         %正爬坡(净负荷增+) 
Cst=[Cst,sum(R_G_down,1) + sum(R_CAES_down_d,1) + sum(R_CAES_down_c,1) >= R_nl_dn];   %负爬坡(净负荷降-)
%系统--频率安全--稳态频差约束V20240421计及CAES的影响
Cst=[Cst,P_LW_error./(df_steady_lim-df_deadzone) - k_load*sum(P_load,1) - k_load*sum(P_CAES_c,1) - ...
    K_WS_1.*sum(P_wind,1) <= sum(u_G.*repmat(P_G_max./adj_G,1,Horizon),1) + sum(u_CAES_d.*repmat(P_CAES_d_max./adj_CAES,1,Horizon),1)];
%系统--频率安全--频率变化率约束V20240421计及CAES的影响
Cst=[Cst,P_LW_error./vf_lim - K_WS_2.*sum(P_wind,1) <= 2*sum(u_G.*repmat(P_G_max.*H_G,1,Horizon),1) + ...
    sum(u_CAES_d.*repmat(P_CAES_d_max.*H_CAES,1,Horizon),1)];

% Cst=[Cst,V_node(1,:) == 1];                                                                                          %首节点电压恒定为1
% Cst=[Cst,phase_node(1,:) == 0];                                                                                      %首节点为参考节点，相角为0
% Cst=[Cst,V_min <= V_node <= V_max];                                                                                  %节点电压约束   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%待修改%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4
%% 求解
obj = C_G_run+C_G_start+C_G_stop+C_G_re+C_G_rl+C_wind + C_load+C_CAES_run+C_CAES_rl+C_CAES_re; 
options = sdpsettings('solver','gurobi','verbose',2,'gurobi.TimeLimit',1200);
options.gurobi.MIPGap = 0.0001;
optimize(Cst,obj,options);

C_G_run = value(C_G_run);
C_G_start = value(C_G_start);
C_G_stop = value(C_G_stop);
C_G_rl = value(C_G_rl);
C_G_re = value(C_G_re);
C_wind = value(C_wind);
obj = value(obj);

u_G = value(u_G);
P_G = value(P_G);
R_G_em= value(R_G_em);
R_G_up = value(R_G_up);
R_G_down = value(R_G_down);

u_CAES_c = value(u_CAES_c);
u_CAES_d = value(u_CAES_d);
P_CAES_c = value(P_CAES_c);
P_CAES_d = value(P_CAES_d);
CAES_pst = value(CAES_pst);
R_CAES_up_c = value(R_CAES_up_c);
R_CAES_up_d = value(R_CAES_up_d);
R_CAES_down_c = value(R_CAES_down_c);
R_CAES_down_d = value(R_CAES_down_d);
%% 传参（跑配置程序用）
Cost_1 = C_G_run+C_G_start+C_G_stop+C_G_rl+C_G_re+C_wind;
P_CAES_c_sum = sum(sum(abs(P_CAES_c)));
P_wind = value(P_wind);
P_load = value(P_load);

%% 传参（跑plot程序用）
c_sum = C_G_run+C_G_start+C_G_stop+C_G_rl+C_wind;
Cost_2 = [C_G_run;C_G_start;C_G_stop;C_G_rl;C_wind;c_sum];

% 绘图数据整理（跑配置程序时隐去）
P_wind_consumption = sum(P_wind)/sum(P_wind_f);
%系统稳态频差计算
df_steady = P_LW_error./(k_load*sum(P_load,1) + K_WS_1.*sum(P_wind,1) + ...
    sum(u_G.*repmat(P_G_max./adj_G,1,Horizon),1) + sum(u_CAES_d.*repmat(P_CAES_d_max./adj_CAES,1,Horizon),1));
df_steady = 50.*value(df_steady);
%系统频率变化率计算
vf = P_LW_error./(2*sum(u_G.*repmat(P_G_max.*H_G,1,Horizon),1) + ...
    2*sum(u_CAES_d.*repmat(P_CAES_d_max.*H_CAES,1,Horizon),1) + K_WS_2.*sum(P_wind,1));
vf = 50.*value(vf);
%数据整理
Fre_data_da = [vf;df_steady];                                                               %用于展示频率指标情况
Gen_data_da = [u_G;P_G;R_G_up;R_G_down];                                                     %用于展示火电机组启停计划、出力情况、备用情况
CAES_data_da = [u_CAES_c;u_CAES_d;P_CAES_c;P_CAES_d;...
    CAES_pst;R_CAES_up_c+R_CAES_up_d;R_CAES_down_c+R_CAES_down_d];                           %用于展示CAES启停计划、出力情况、备用情况

%% 绘图
% figure(1)
% heatmap(Gen_data_da(1:15,:));
% title('火电开机状态');
% 
% figure(2)
% ribbon(Gen_data_da(16:30,:)')
% title('火电出力情况');
% 
% figure(3)
% heatmap(CAES_data_da(1:9,:)-CAES_data_da(10:18,:));
% colormap("gray")
% title('CAES运行工况(压缩为正/发电为负)');
% 
% figure(4)
% ribbon(CAES_data_da(19:27,:)'-CAES_data_da(28:36,:)');
% title('CAES出力情况(压缩为正/发电为负)');
% 
% figure(5)
% ribbon(CAES_data_da(37:45,:)');
% title('CAES储气室气压');
% zlim([40 55])
% 
% figure(6)
% bar3(sum(CAES_data_da(19:27,:)-CAES_data_da(28:36,:)));
% title('CAES总出力情况(压缩为正/发电为负)');
% 
% figure(7)
net_load_1 = P_load - sum(P_wind,1);
% bar3(net_load_1);
% title('净负荷1(负荷-新能源)');
% 
% figure(8)
net_load_2 = P_load - sum(P_wind,1) + ...
    sum(CAES_data_da(19:27,:)-CAES_data_da(28:36,:));
% bar3(net_load_2);
% title('净负荷2(负荷-新能源-储能)');

%% 数据分析
% v2p_1 = max(net_load_1) - min(net_load_1);
% v2p_2 = max(net_load_2) - min(net_load_2);
% ramp_max_1 = max(diff(net_load_1));
% ramp_max_2 = max(diff(net_load_2));
if C_G_run == 0
    Cost_1 = M*M;
end
end