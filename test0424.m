clear
clc
        mater = zeros(16,9);
% mater(1,7) = ones(1,1);
%         mater(2,8:8) = ones(1,1);
% mater(2,4) = ones(1,1);
%         mater(3,8) = ones(1,1);
%         mater(2,1:1) = ones(1,1);
        mater(4,1:3) = ones(1,3);
% %         mater(8,7) = ones(1,1);mater(8,9) = ones(1,1);
% %         mater(12,6) = ones(1,1);
% %         mater(14,8) = ones(1,1);
state = 3;          %阶段
flag = 1;          %典型日

x = mater;
[cost_d,P_CAES_d] = Underlevel_CAES(x,state,flag);    %调用下层函数计算典型日的运行成本
