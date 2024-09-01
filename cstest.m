%用于计算经济性成本
clc,clear
% load('result_202407020044.mat')
% Costfun_CAES(global_M(27).p)

mater = zeros(16,9);
mater(4,8) = ones(1,1);mater(4,9) = ones(1,1);
mater(8,6) = ones(1,1);mater(8,7) = ones(1,1);
mater(12,5) = ones(1,1);
mater(14,4) = ones(1,1);mater(15,2) = ones(1,1);
state = 4;
flag = 2;
x = mater;
Fitness = Costfun_CAES(x);