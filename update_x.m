%函数用于：更新粒子位置
%V20240522创建，需要修改满足上层约束
function [var_new] = update_x(var_0,var_person,var_global,w,C1,C2)
n1 = size(var_0, 1);                  %规划时段数
n2 = size(var_0, 2);                  %地理点个数
%三粒子地理点配置上限
Geo_up_0 = nnz(var_0);                %自我--地理点配置上限
Geo_up_1 = nnz(var_person);           %个体--地理点配置上限
Geo_up_2 = nnz(var_global);           %全局--地理点配置上限

%三粒子数据处理--非零元素取值（1~2）--零元素取值（0~1.5）
g_0 = (1-Geo_up_0).*1.5*rand(1, n2) + Geo_up_0.*(1+rand(1, n2));
g_1 = (1-Geo_up_1).*1.5*rand(1, n2) + Geo_up_1.*(1+rand(1, n2));
g_2 = (1-Geo_up_2).*1.5*rand(1, n2) + Geo_up_2.*(1+rand(1, n2));

x_0 = (1-var_0).*1.5.*rand(n1, n2) + var_0.*(1+rand(n1, n2));
x_1 = (1-var_person).*1.5.*rand(n1, n2) + var_person.*(1+rand(n1, n2));
x_2 = (1-var_global).*1.5.*rand(n1, n2) + var_global.*(1+rand(n1, n2));

%三粒子合并处理
Geo_up = round((w*Geo_up_0+C1*Geo_up_1+C2*Geo_up_2)/(w+C1+C2));%加权平均处理
g_temp = w.*g_0 + C1.*g_1 + C2.*g_2;
x_temp = w.*x_0 + C1.*x_1 + C2.*x_2;
if Geo_up > n2
    Geo_up = n2;%判断是否越限
elseif Geo_up <= 7
    Geo_up = 7;%判断是否越限
end

for t = 1:20
    %生成新粒子
    for i = 1 : n2
        xa(:,i) = double(x_temp(:,i)==max(x_temp(:,i)));
    end
    xb = rand(1, n2);
    temp_a = sort(xb,"descend");
    temp_b = temp_a(Geo_up);
    xb = double(xb>=temp_b);

    var_new = xa.*xb;

    ft = feasibility_test(var_new);
    if ft == 1
        break;
    else
        Geo_up = Geo_up+1;
        Geo_up = min(Geo_up,n2);
    end
end
end