function Fitness = Costfun_CAES(x)
% function Cost_data = Costfun_CAES(x)
% 中层计算适应度函数：场景A
% V20240412创建,V20240416开始编写,V20240417/18调试修改--1代,局部跑通
% V20240522/23调试修改--2代,全局系统调试
% V20240604增加展示成本项部分Cost_data
% 输入参数x：配置方案(其中每一个元素表示：某一时段某一地理点接入某一类型CAES的个数)
% 输出参数Fitness：该配置方案下计算得到的适应度函数，即总成本
%% 参数导入
% parameter_import          %每当更新设备参数时使用该程序，调用"参数导入"函数更新parameter.mat
load('parameter.mat')       %导入设备参数
x_time = sum(x,2);          %仅考虑时间属性，忽略地理属性，用于计算成本[4*1;4*1;4*1;4*1]分为四阶段

%% 参数设置/处理
T_all = T_all;                                  %总规划周期40y
T_es = 10;                                      %单阶段规划周期10y
N_State = 4;                                    %规划阶段数
c_CAES_con = P_CAES_cN.*c_CAES_c + ...
    P_CAES_dN.*c_CAES_d + V_CAES_N.*c_CAES_h;           %单机建设成本
M = 10^10;                                             %大M

%% 成本项定义
Fitness = 0;                                    %定义适应度函数--总成本项
C_con = 0;                                      %定义成本项1：储能建设成本
C_om = 0;                                       %定义成本项2：储能维护成本
C_op = 0;                                       %定义成本项3：电网运行成本
B_sal = 0;                                      %定义成本项4：残值回收成本
C_op_ys = 0;                                    %定义辅助成本项：年电网运行成本

config_lim = [450,900,1250,1650];
for i = 1: 4
    judge_1(i) = (sum(x_time(1:4*i,:).*repmat(P_CAES_dN,i,1)) >= config_lim(i));
    judge_2(i) = sum(x_time(4*i-3:4*i,:)) > 0;
end

if ~all(any([judge_1,judge_2],1))
    Fitness = M^2.5;                          %计算总成本，即适应度函数
else
    %% 调用下层函数
    %调用函数计算日前调度结果
    parfor state = 1: N_State+1
        for flag = 1: 4
            [cost_d,P_CAES_d] = Underlevel_CAES(x,state,flag);    %调用下层函数计算典型日的运行成本
            Cost(state,flag) = cost_d;                            %存储！日！调度运行成本
            P_CAES(state,flag) = P_CAES_d;                        %存储！日！CAES充放电容量
        end
    end

    %% 数据处理计算
    %计算成本项1：储能建设成本
    for s = 1: N_State
        xs(s,:) = x_time(4*s-3:4*s)';                       %该阶段下需要配置的CAES台数[1*4]
        gr = (1+gama)^(-(s-1)*T_es);                        %该阶段下折现系数计算
        C_con = C_con + gr*sum(xs(s,:).* c_CAES_con');      %累加该阶段储能建设成本
                        pin_con(s,:) = gr*sum(xs(s,:).* c_CAES_con');  %各阶段储能建设成本
    end

    %计算成本项2：储能维护成本--固定部分
    xi = 1:1:T_all;                                         %辅助变量：各年[1*40]
    for s = 1: N_State
        C_om_fix_y((s-1)*T_es+1:s*T_es) = repmat(sum(c_CAES_OM_fix'.*P_CAES_dN'.*sum(xs(1:s,:))),1,T_es);    %计算各年的！年！固定维护成本 [1*40]
    end
    C_om_fix = sum(((1+gama).^(-xi)).*C_om_fix_y);          %计算！总！固定维护成本

    %计算成本项2：储能维护成本--可变部分
    x = [1,T_es+1,2*T_es+1,3*T_es+1,4*T_es];                %辅助变量：各规划阶段初年+终年[1*5]
    C_om_var_ys = sum(P_CAES.* ...
        repmat(c_CAES_OM_var',5,1).*kai',2)';               %计算各规划阶段初年的！年！可变维护成本[1*5]
    xi = 1:1:T_all;                                         %辅助变量：各年[1*40]
    C_om_var_y = interp1(x,C_om_var_ys,xi,"spline");        %计算各年的！年！可变维护成本 [1*40]
    C_om_var = sum(((1+gama).^(-xi)).*C_om_var_y);          %计算！总！可变维护成本

    %计算成本项2：储能维护成本--总
    C_om = C_om_fix + C_om_var;                                     %计算！总！储能维护成本

    %计算成本项3：电网运行成本
    x = [1,T_es+1,2*T_es+1,3*T_es+1,4*T_es];                        %辅助变量：各规划阶段初年+终年[1*5]
    C_op_ys = sum(Cost.*kai',2)';                                   %计算各规划阶段初年的！年！调度运行成本[1*5]
    xi = 1:1:T_all;                                                 %辅助变量：各年[1*40]
    C_op_y = interp1(x,C_op_ys,xi,"spline");                        %计算各年的！年！调度运行成本 [1*40]
    C_op = sum(((1+gama).^(-xi)).*C_op_y);                          %计算！总！调度运行成本

    %计算成本项4：残值回收成本
    T_use = repmat(T_es.*[4;3;2;1],1,4);                                %各阶段CAES电站使用时间
    k_dep = (1-T_use.*(1-delta_sal)./repmat(T_CAES_span',N_State,1));   %计算折旧系数[4*4]
    B_sal = (1+gama)^(-T_all)*sum(sum(xs.*c_CAES_con'.*k_dep));         %计算！总！残值回收成本

    %% 输出结果
    F_ESS = C_con + C_om - B_sal;                                   %计算储能总成本
    Fitness = C_con + C_om + C_op - B_sal;                          %计算总成本，即适应度函数

    for s = 1: N_State
        gr = (1+gama)^(-(s-1)*T_es);
        pin_con(s,:) = gr*sum(xs(s,:).* c_CAES_con');                          %各阶段储能建设成本
        pin_om_fix(s,:) = sum(((1+gama).^ ...
            (-xi((s-1)*T_es+1:s*T_es))).*C_om_fix_y((s-1)*T_es+1:s*T_es));     %各阶段固定维护成本
        pin_om_var(s,:) = sum(((1+gama).^ ...
        (-xi((s-1)*T_es+1:s*T_es))).*C_om_var_y((s-1)*T_es+1:s*T_es));         %各阶段可变维护成本
        pin_op(s,:) = sum(((1+gama).^ ...
            (-xi((s-1)*T_es+1:s*T_es))).*C_op_y((s-1)*T_es+1:s*T_es));         %各阶段调度运行成本
    end

    Cost_data = [C_op;C_con;C_om;B_sal;Fitness];

% %     filename = 'D:\文章-大规模CAES多阶段优化规划-程序\算例EP2.xlsx';
% %     writematrix(Cost_data,filename,'Sheet','规划层成本分析','Range','B2:B6')
    xlswrite("D:\文章-大规模CAES多阶段优化规划-程序\算例EP1.xlsx",Cost_data,'规划层成本分析','B2:B6');
    xlswrite("D:\文章-大规模CAES多阶段优化规划-程序\算例EP1.xlsx",pin_op,'电网运行成本分析','B2:B5');
    xlswrite("D:\文章-大规模CAES多阶段优化规划-程序\算例EP1.xlsx",pin_con,'储能建设运维成本分析','B2:B5');
    xlswrite("D:\文章-大规模CAES多阶段优化规划-程序\算例EP1.xlsx",pin_om_fix,'储能建设运维成本分析','G2:G5');
    xlswrite("D:\文章-大规模CAES多阶段优化规划-程序\算例EP1.xlsx",pin_om_var,'储能建设运维成本分析','L2:L5');
end
end