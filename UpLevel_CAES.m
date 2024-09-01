%�ϲ��Ż��滮���������������÷���
% V20240412�������ϲ�����ģ�ͣ�����A����������
% V20240522�޸ĵ���ȫϵͳ��ͨ����ģ15*10,����*������
% V20240528�޸ĵ��Ը����ģѰ�ţ���ģ30*20,����*������
clear
clc
%% ��������
%ģ�Ͳ���
% % % par_CAES = xlsread('G:\������\����-���ģCAES��׶��Ż��滮-����\0.���ò���\������������V20240228.xlsx','�������ȼ�CAES��������','B3:E20');%��EXCEL�������ݣ�CAES����
State = 4;                                          %�滮ʱ��N��
Geography = 9;                                      %�滮����ѡ��M��
Type = 4;                                           %�滮�����ȼ�K��,�����ȼ�����С����ķ�ʽ����
variable = State*Geography*Type;                    %���ӱ���ά�ȣ�N*M*K
var_max = repmat([1;1;1;1],State,Geography);        %ĳһʱ��ĳһ��������ĳһ����CAES��������
var_min = repmat([0;0;0;0],State,Geography);        %ĳһʱ��ĳһ��������ĳһ����CAES����С����

%�㷨����
global partysize;
partysize = 50;                                     %PSO�㷨--������Ŀ
Dmax = 30;                                          %PSO�㷨--����������
Matersmax = var_max;                                %PSO�㷨--��������
Matersmin = var_min;                                %PSO�㷨--��������
step = repmat([1;1;1;1],State,Geography);           %PSO�㷨--���Ӳ���
len = Matersmax - Matersmin;                        %PSO�㷨--������������
v_max = repmat([1;1;1;1],State,Geography);          %PSO�㷨--���������ٶ�����
tolerance = 15;                                     %PSO�㷨--�������  
tol_D = 8;                                          %PSO�㷨--����tol_D��ȫ������ֵ�仯С��tolerance����Ϊ����
C1 = 2;                                             %PSO�㷨--����ѧϰ����
C2 = 1;                                             %PSO�㷨--Ⱥ��ѧϰ����
wholetime = zeros(Dmax,1);                          %PSO�㷨--�洢ÿһ�ε�������ʱ��
k = 0;                                              %PSO�㷨--���������ƴ���
global_f = zeros(1,Dmax);                           %PSO�㷨--ȫ�����Ž�

%%  ��Ⱥ��ʼ��
var = struct('Matersmin',Matersmin,'step',step, ...
    'Matersmax',Matersmax,'len',len,'v_max',v_max); %PSO�㷨--������
for i = 1: partysize
%     if (i==1||i==2)%����1
%         mater = zeros(16,9);
%         mater(3,4) = ones(1,1);mater(4,8) = ones(1,1);
%         mater(8,5) = ones(1,1);mater(8,7) = ones(1,1);
%         mater(10,2) = ones(1,1);mater(15,9) = ones(1,1);
%         mater(16,1) = ones(1,1);
%         Maters(i) = struct('p',mater);                  %Maters�����ŵ�������
%     elseif (i==3||i==4)%����2 
%         mater = zeros(16,9);
%         mater(3,4) = ones(1,1);mater(4,8) = ones(1,1);
%         mater(8,5) = ones(1,1);mater(8,7) = ones(1,1);
%         mater(10,2) = ones(1,1);mater(15,9) = ones(1,1);
%         mater(16,1) = ones(1,1);
%         Maters(i) = struct('p',mater);                  %Maters�����ŵ�������
%     elseif (i==5||i==6)%����3
%         mater = zeros(16,9);
%         mater(4,4) = ones(1,1);mater(3,8) = ones(1,1);
%         mater(8,5) = ones(1,1);mater(8,7) = ones(1,1);
%         mater(11,2) = ones(1,1);mater(15,9) = ones(1,1);
%         mater(16,1) = ones(1,1);
%         Maters(i) = struct('p',mater);                  %Maters�����ŵ������� 
%     elseif (i==7||i==8)%����3
%         mater = zeros(16,9);
%         mater(4,4) = ones(1,1);mater(3,8) = ones(1,1);
%         mater(8,5) = ones(1,1);mater(8,7) = ones(1,1);
%         mater(10,2) = ones(1,1);mater(15,9) = ones(1,1);
%         mater(16,1) = ones(1,1);
%         Maters(i) = struct('p',mater);                  %Maters�����ŵ�������

%     else
        mater = create_x(var);                          %����������ӣ���ʼ��λ��
        Maters(i) = struct('p',mater);                  %Maters�����ŵ�������
%     end
end

%% ����ȺѰ��
for d = 1: Dmax
    tstart = tic;
    %% ������Ӧ�Ⱥ���
    parfor j = 1: partysize
        Copp(j) = Costfun_CAES(Maters(j).p);  %�������²㺯��--�����ܳɱ���ֵ
        disp(strcat('��',num2str(d),'����',num2str(j),'������������,��Ӧ�Ⱥ���Ϊ',num2str(Copp(j))));
    end

    %% ��������λ��/ֵ�ĸ���
    for j = 1: partysize
        fitness(j) = Copp(j);               %��Ӧ�Ⱥ���
        if d == 1
            personal_M(j).p = Maters(j).p;  %��һ�ε�������������λ��Ϊ��ʼλ��
            personal_f(j) = fitness(j);     %��һ�ε�������������ֵΪ��ʼֵ
        elseif fitness(j) < personal_f(j)   %��ĳ�ε�����Ӧ�Ⱥ����ȸ�������ֵСʱ�����¸�������ֵ
            personal_M(j).p = Maters(j).p;  %2+�ε�������������λ��--����λ��
            personal_f(j) = fitness(j);     %2+�ε�������������ֵ--����ֵ
        end
    end

    %% ȫ������λ��/ֵ�ĸ���
    [global_f(d),jj] = min(personal_f);     %ȫ�����Ž�--ֵ�����������ӣ�Ҳ���Ǹ������Ž��ǣ�����ģ���Ӧ�Ⱥ�������Сֵ�ҵ�
    global_M(d).p = personal_M(jj).p;       %ȫ�����Ž��λ��(�����ӵ�ȡֵ)

    %% �����о�
    if d>=2 && global_f(d-1)-global_f(d)<tolerance
        k = k+1;
        if k >= tol_D                       %�������������tol_D��ȫ������ֵ�仯�����������������Ѱ��
            break;
        end
    else
        k = 0;
    end

    %% �ٶȸ���
    w = 0.8-(0.8-0.6)/Dmax*d;                               %��������dԽ�󣬹���Ȩ��ԽС������ֲ�Ѱ������
%     for j = 1: partysize
%         v_n = v(j).p;                                       %�����ٶȸ��¸������� v_n [12*9]
%         v_n = round(w.*v_n + ...                            %�ٶ���1��ԭʼֵ
%             C1.*rand.*(personal_M(j).p-Maters(j).p) + ...   %�ٶ���2������ѧϰ
%             C2.*rand.*(global_M(d).p-Maters(j).p));         %�ٶ���3��Ⱥ��ѧϰ
%         v_n = min(v_n,v_max);                               %�ٶ�����Լ��
%         v_n = max(v_n,-v_max);                              %���º���ٶȲ��������ֵ
%         v(j).p = v_n;                                       %�����ٶȸ���
%     end

    %% λ�ø���
    for j = 1: partysize
        Maters(j).p = update_x(Maters(j).p,personal_M(j).p,global_M(d).p,w,C1,C2);
        Maters(j).p = min(Maters(j).p,Matersmax);   %λ������Լ��
        Maters(j).p = max(Maters(j).p,Matersmin);
    end
    wholetime(d,1) = toc(tstart);
end
F_ULtotal = global_f(d);            %�����Ż������������Ӧ�Ⱥ���    
P_ULall = global_M;                 %�洢ÿ�ε������������ȫ������λ��