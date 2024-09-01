function typical_day_import
%导入典型日源荷数据,V20240423使用24点数据
%% excel文件写入数据
filename = 'F:\王佳旭\文章-大规模CAES多阶段优化规划-程序\0.公用参数\新能源负荷预测数据V20240601.xlsx';
par_tday_1 = readmatrix(filename,'Sheet','典型日1源荷数据','Range','B2:E25');               %从EXCEL导入数据--典型日1源荷数据
par_tday_2 = readmatrix(filename,'Sheet','典型日2源荷数据','Range','B2:E25');               %从EXCEL导入数据--典型日2源荷数据
par_tday_3 = readmatrix(filename,'Sheet','典型日3源荷数据','Range','B2:E25');               %从EXCEL导入数据--典型日3源荷数据
par_tday_4 = readmatrix(filename,'Sheet','典型日4源荷数据','Range','B2:E25');               %从EXCEL导入数据--典型日4源荷数据

%% 典型日参数以结构体形式存储
typical_day(1) = struct('d',par_tday_1);
typical_day(2) = struct('d',par_tday_2);
typical_day(3) = struct('d',par_tday_3);
typical_day(4) = struct('d',par_tday_4);

%% 数据保存
save ('typical_day.mat','typical_day')
end