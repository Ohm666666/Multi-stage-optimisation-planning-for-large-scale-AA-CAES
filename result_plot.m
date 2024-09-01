clc,clear
%% 参数导入
x = [47;44;52;28]*2;
% x = [0;0;0;0];
% parameter_import  %每当更新设备参数时使用该程序，调用"参数导入"函数更新parameter.mat
load parameter.mat

ESS_EN = x .* E_Bat_N;
ESS_PN = x .* P_Bat_N;

%% 计算目标函数——总成本年值
F_1 = 0;
F_ESS = zeros(1,1);                             %储能投资成本
F_Run = zeros(1,1);                             %负荷用电成本
r_cal = r*((1+r)^T_pro)/((1+r)^T_pro-1);        %折算系数

for flag=1:4
    [cost,P_ess] = Underlevel_FM(x,flag);  %调用下层函数计算典型日的运行成本
    Cost(:,flag) = cost;
    P_ESS(flag) = P_ess;
end

P_ESS_sum = sum(P_ESS .* [91,92,91,92]);                             %计算储能年充放电功率MW
C_ins = (c_ESS_E*sum(ESS_EN) + c_ESS_P*sum(ESS_PN));                 %计算储能投资成本                                 
C_rep = (1/(1+r)^T_Bat) * C_ins;                                     %计算储能更新成本                               
C_om_y = c_ESS_OM_fix*sum(ESS_PN) + c_ESS_OM_var*P_ESS_sum ;         %计算储能年运维成本
F_Run = sum((Cost .* [91,92,91,92]),2);                              %计算年运行成本                              
F_ESS = r_cal*(C_ins + C_rep) + C_om_y;                              %计算储能总成本
F_1 = F_ESS + F_Run(6);                                              %计算总成本年值，即适应度函数    

Cost_data = [F_1;F_Run;F_ESS;C_ins*r_cal;C_rep*r_cal;C_om_y];