%用于计算日前调度计划
clc,clear
% load('result_202407020044.mat')
% mater = global_M(27).p;

mater = zeros(16,9);
mater(4,8) = ones(1,1);mater(4,9) = ones(1,1);
mater(8,6) = ones(1,1);mater(8,7) = ones(1,1);
mater(12,5) = ones(1,1);
mater(14,4) = ones(1,1);mater(15,2) = ones(1,1);

state = 4;          %阶段
flag = 1;           %典型日

x = mater;
[cost_d,P_CAES_d] = Underlevel_CAES_EP3(x,state,flag);    %调用下层函数计算典型日的运行成本
