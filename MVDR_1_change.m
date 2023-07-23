% Title:    最小方差无畸变响应波束形成器，MVDR
% Author:   Xu Zhe
% Date:     2020-07-25
% Brief:    一维均匀线阵波束形成，MVDR波束形成器

% 系统初始化 Initialization
close all
clear
clc

% 参数 Parameters
M   = 12;               % 阵元数目
c   = 3e8;              % 光速
f   = 10e9;             % 频率
l   = c/f;              % 波长
d   = l/2;              % 阵元间隔 --关注
phi = [10,30,0]*pi/180;   % 信号和干扰角度
Kt  = 1;                % 目标信号数目
snr = 20;               % 信噪比
inr = 30;               % 干噪比
Ac  = 3601;             % 角度采样数

% 计算阵元位置
R = d * (0:M-1)';

% 计算方向向量
K = length(phi);        % 信号数目
A = zeros(M,K);         % 初始化方向矩阵 Steering Matrix
for k = 1:K
    for m = 1:M
        A(m,k) = exp(1j * 2 * pi * R(m,1)* sin(phi(k))/l);
    end
end

% 构造信号
Ts = 1024;              % 快拍数
Kj = K - Kt;            % 干扰信号数目
% 一个干扰信号，一个目标
s = zeros(K,Ts);        % 信号
sd = zeros(1,Ts);       % 期望信号
for k = 1:Kt
    for i = 0:Ts-1
        s(k,i+1) = sqrt(10^(snr/10)) * exp(1j * 2 * pi * f * i/(2*Ts)); % 3*Ts
    end
    sd = sd + s(k,:);
end
sd = sd/Kt;
if Kj ~= 0
    for k = Kt+1:K
        for i = 0:Ts-1
            s(k,i+1) = sqrt(10^(inr/10)/2) * (randn(1) + 1j*randn(1)) * exp(1j * 2 * pi * f * i/(2*Ts));
        end
    end
end

% 噪声
n = (randn(M,Ts)+1j*randn(M,Ts))/sqrt(2);
% 接收数据向量
x = A * s+n;

% MVDR  
J = A * s; % 不再只考虑干扰
% Rnj代替Rx
Rnj = 1/Ts * (J*J'+n*n');
% inv,pinv求逆
Rpriv = pinv(Rnj);
w = Rpriv*A(:,1)/(A(:,1)'*Rpriv*A(:,1));

% 方向图
theta = linspace(-pi,pi,Ac);
y = zeros(1,Ac);
for i = 1:Ac
    a = exp(1j * 2 * pi * R(:,1)* sin(theta(i))/l);  % a [M,1]为对应theta角的导向矢量
    y(i) = w' * a;
end

ydB = 20*log10(abs(y)/max(abs(y)));

% 绘图
figure('name','波束方向图','color',[1,1,1]);
set(gcf,'position',[458,342,290,200])
plot(theta * 180/pi,ydB,'color','b','linestyle','-','linewidth',1);
grid on
xlim([-90,90]);
ylim([-60,0]);
set(gca,'xtick',-90:30:90);
set(gca,'FontName','Times New Roman','FontSize',10.5);
xlabel('Direction of Arrival of the Signal (Degree)','FontSize',10.5);
ylabel('Magnitude Response (dB)','FontSize',10.5);

% 保存
% print('-depsc','E:\07 Latex\Paper Note\figures\MVDR_1.eps');
