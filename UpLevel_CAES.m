%上层优化规划程序，用于生成配置方案
% V20240412创建：上层配置模型：场景A——主程序
% V20240522修改调试全系统跑通（规模15*10,个数*代数）
% V20240528修改调试更大规模寻优（规模30*20,个数*代数）
clear
clc
%% 参数设置
%模型参数
% % % par_CAES = xlsread('G:\王佳旭\文章-大规模CAES多阶段优化规划-程序\0.公用参数\文章数据收资V20240228.xlsx','多容量等级CAES参数设置','B3:E20');%从EXCEL导入数据，CAES参数
State = 4;                                          %规划时段N个
Geography = 9;                                      %规划地理选点M个
Type = 4;                                           %规划容量等级K个,容量等级按从小到大的方式排列
variable = State*Geography*Type;                    %粒子变量维度：N*M*K
var_max = repmat([1;1;1;1],State,Geography);        %某一时段某一地理点接入某一类型CAES的最大个数
var_min = repmat([0;0;0;0],State,Geography);        %某一时段某一地理点接入某一类型CAES的最小个数

%算法参数
global partysize;
partysize = 50;                                     %PSO算法--粒子数目
Dmax = 30;                                          %PSO算法--最大迭代次数
Matersmax = var_max;                                %PSO算法--粒子上限
Matersmin = var_min;                                %PSO算法--粒子下限
step = repmat([1;1;1;1],State,Geography);           %PSO算法--粒子步长
len = Matersmax - Matersmin;                        %PSO算法--粒子搜索区间
v_max = repmat([1;1;1;1],State,Geography);          %PSO算法--粒子搜索速度上限
tolerance = 15;                                     %PSO算法--收敛误差  
tol_D = 8;                                          %PSO算法--连续tol_D次全局最优值变化小于tolerance则认为收敛
C1 = 2;                                             %PSO算法--自我学习因子
C2 = 1;                                             %PSO算法--群体学习因子
wholetime = zeros(Dmax,1);                          %PSO算法--存储每一次迭代所需时间
k = 0;                                              %PSO算法--收敛代数计次用
global_f = zeros(1,Dmax);                           %PSO算法--全局最优解

%%  种群初始化
var = struct('Matersmin',Matersmin,'step',step, ...
    'Matersmax',Matersmax,'len',len,'v_max',v_max); %PSO算法--传参用
for i = 1: partysize
%     if (i==1||i==2)%种子1
%         mater = zeros(16,9);
%         mater(3,4) = ones(1,1);mater(4,8) = ones(1,1);
%         mater(8,5) = ones(1,1);mater(8,7) = ones(1,1);
%         mater(10,2) = ones(1,1);mater(15,9) = ones(1,1);
%         mater(16,1) = ones(1,1);
%         Maters(i) = struct('p',mater);                  %Maters会随着迭代更改
%     elseif (i==3||i==4)%种子2 
%         mater = zeros(16,9);
%         mater(3,4) = ones(1,1);mater(4,8) = ones(1,1);
%         mater(8,5) = ones(1,1);mater(8,7) = ones(1,1);
%         mater(10,2) = ones(1,1);mater(15,9) = ones(1,1);
%         mater(16,1) = ones(1,1);
%         Maters(i) = struct('p',mater);                  %Maters会随着迭代更改
%     elseif (i==5||i==6)%种子3
%         mater = zeros(16,9);
%         mater(4,4) = ones(1,1);mater(3,8) = ones(1,1);
%         mater(8,5) = ones(1,1);mater(8,7) = ones(1,1);
%         mater(11,2) = ones(1,1);mater(15,9) = ones(1,1);
%         mater(16,1) = ones(1,1);
%         Maters(i) = struct('p',mater);                  %Maters会随着迭代更改 
%     elseif (i==7||i==8)%种子3
%         mater = zeros(16,9);
%         mater(4,4) = ones(1,1);mater(3,8) = ones(1,1);
%         mater(8,5) = ones(1,1);mater(8,7) = ones(1,1);
%         mater(10,2) = ones(1,1);mater(15,9) = ones(1,1);
%         mater(16,1) = ones(1,1);
%         Maters(i) = struct('p',mater);                  %Maters会随着迭代更改

%     else
        mater = create_x(var);                          %随机生成粒子，初始化位置
        Maters(i) = struct('p',mater);                  %Maters会随着迭代更改
%     end
end

%% 粒子群寻优
for d = 1: Dmax
    tstart = tic;
    %% 计算适应度函数
    parfor j = 1: partysize
        Copp(j) = Costfun_CAES(Maters(j).p);  %调用中下层函数--计算总成本年值
        disp(strcat('第',num2str(d),'代第',num2str(j),'个粒子求解完成,适应度函数为',num2str(Copp(j))));
    end

    %% 个体最优位置/值的更新
    for j = 1: partysize
        fitness(j) = Copp(j);               %适应度函数
        if d == 1
            personal_M(j).p = Maters(j).p;  %第一次迭代：个体最优位置为初始位置
            personal_f(j) = fitness(j);     %第一次迭代：个体最优值为初始值
        elseif fitness(j) < personal_f(j)   %当某次迭代适应度函数比个体最优值小时，更新个体最优值
            personal_M(j).p = Maters(j).p;  %2+次迭代：个体最优位置--更新位置
            personal_f(j) = fitness(j);     %2+次迭代：个体最优值--更新值
        end
    end

    %% 全局最优位置/值的更新
    [global_f(d),jj] = min(personal_f);     %全局最优解--值：把所有粒子（也就是个体最优解们）里面的（适应度函数）最小值找到
    global_M(d).p = personal_M(jj).p;       %全局最优解的位置(即粒子的取值)

    %% 收敛判据
    if d>=2 && global_f(d-1)-global_f(d)<tolerance
        k = k+1;
        if k >= tol_D                       %条件：如果连续tol_D次全局最优值变化不超过收敛误差，则结束寻优
            break;
        end
    else
        k = 0;
    end

    %% 速度更新
    w = 0.8-(0.8-0.6)/Dmax*d;                               %迭代次数d越大，惯性权重越小，增大局部寻优能力
%     for j = 1: partysize
%         v_n = v(j).p;                                       %定义速度更新辅助变量 v_n [12*9]
%         v_n = round(w.*v_n + ...                            %速度项1：原始值
%             C1.*rand.*(personal_M(j).p-Maters(j).p) + ...   %速度项2：自我学习
%             C2.*rand.*(global_M(d).p-Maters(j).p));         %速度项3：群体学习
%         v_n = min(v_n,v_max);                               %速度区间约束
%         v_n = max(v_n,-v_max);                              %更新后的速度不超过最大值
%         v(j).p = v_n;                                       %粒子速度更新
%     end

    %% 位置更新
    for j = 1: partysize
        Maters(j).p = update_x(Maters(j).p,personal_M(j).p,global_M(d).p,w,C1,C2);
        Maters(j).p = min(Maters(j).p,Matersmax);   %位置区间约束
        Maters(j).p = max(Maters(j).p,Matersmin);
    end
    wholetime(d,1) = toc(tstart);
end
F_ULtotal = global_f(d);            %最终优化结果：最优适应度函数    
P_ULall = global_M;                 %存储每次迭代结果：历次全局最优位置