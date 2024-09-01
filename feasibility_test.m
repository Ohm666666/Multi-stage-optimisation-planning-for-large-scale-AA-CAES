function y = feasibility_test(x)
x_time = sum(x,2);          %仅考虑时间属性，忽略地理属性，用于计算成本[4*1;4*1;4*1;4*1]分为四阶段
config_lim = [450,900,1250,1650];
load('parameter.mat')       %导入设备参数
for i = 1: 4
    judge_1(i) = (sum(x_time(1:4*i,:).*repmat(P_CAES_dN,i,1)) >= config_lim(i));
    judge_2(i) = sum(x_time(4*i-3:4*i,:)) > 0;
end
if ~all(any([judge_1,judge_2],1))
    y = 0;
else
    y = 1;
end