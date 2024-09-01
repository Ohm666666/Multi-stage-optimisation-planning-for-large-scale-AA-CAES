%函数用于：随机生成粒子，初始化位置
%V20240415创建，需要修改满足上层约束
function [x] = create_x(var)
n1 = size(var.Matersmin, 1);
n2 = size(var.Matersmin, 2);        %地理点个数
Geo_up = randi([7, n2], 1);         %生成地理点配置上限
for t = 1:20
    %粒子处理第一步--保证各地理点只建一台机组
    x0 = rand(n1, n2);                  %生成随机矩阵【n1,n2】
    for i = 1 : n2
        x1(:,i) = double(x0(:,i)==max(x0(:,i)));
    end
    %粒子处理第二步--生成应该配置的地理点
    x2 = rand(1, n2);
    temp_a = sort(x2,"descend");
    if Geo_up == 0
        temp_b = inf;
    else
        temp_b = temp_a(Geo_up);
    end
    x2 = double(x2>=temp_b);
    %粒子处理第三步--产生粒子
    x = x1.*x2;
    ft = feasibility_test(x);
    if ft == 1
        break;
    else
        Geo_up = Geo_up+1;
        Geo_up = min(Geo_up,n2);
    end
end
end