function [Cost_2,P_ess] = Underlevel_test(x,st,ss)    

Cost_2 = sum(sum(x+st+ss));
P_ess = Cost_2;

end