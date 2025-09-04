close all;
clear;
clc;
% load('Pnvskappanew.mat')
figure
width_cm = 10;
height_cm = 10;
set(gcf,'unit','centimeters','position',[20 5 width_cm height_cm])
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [1.6/height_cm 1.45/width_cm],[1.2/height_cm 0.8/height_cm], [1.5/width_cm 0.3/width_cm]);
if ~make_it_tight,  clear subplot;  end       %%%%%%%%%% gap [h, w], h [down, up], w [left right]

tic;% 开始计时

%参数
delta = 2*pi*0; %MHz  delta = Delta_1 + Delta_2
Delta = 2*pi*0.1; %MHz  Delta1
Omega1 = 2*pi*1; %MHz
Omega2 = 2*pi*1; %MHz
thetacL = pi/2;
thetacR = -pi/2;
theta = -pi/2;
%phi1 = -pi/3;
%phi2 = pi/9;

%定义微分方程组
odefunL = @(tL, yL) [
    -1.0j*0.5*Omega1*(yL(2)+yL(3));  % dy1/dt
    -1.0j*Delta*yL(2)-1.0j*0.5*Omega1*yL(1)-1.0j*0.5*Omega2*exp(-1.0j*theta)*yL(4);  % dy2/dt
    -1.0j*Delta*yL(3)-1.0j*0.5*Omega1*yL(1)+1.0j*0.5*Omega2*exp(-1.0j*thetacL)*yL(4); % dy3/dt
    -1.0j*delta*yL(4)-1.0j*0.5*Omega2*exp(1.0j*theta)*yL(2)+1.0j*0.5*Omega2*exp(1.0j*thetacL)*yL(3) % dy4/dt
];

odefunR = @(tR, yR) [
    -1.0j*0.5*Omega1*(yR(2)+yR(3));  % dy1/dt
    -1.0j*Delta*yR(2)-1.0j*0.5*Omega1*yR(1)-1.0j*0.5*Omega2*exp(-1.0j*theta)*yR(4);  % dy2/dt
    -1.0j*Delta*yR(3)-1.0j*0.5*Omega1*yR(1)+1.0j*0.5*Omega2*exp(-1.0j*thetacR)*yR(4); % dy3/dt
    -1.0j*delta*yR(4)-1.0j*0.5*Omega2*exp(1.0j*theta)*yR(2)+1.0j*0.5*Omega2*exp(1.0j*thetacR)*yR(3) % dy4/dt
];

% odefunL = @(tL, yL) [
%     -1.0j*0.5*Omega1*exp(-1.0j*phi1)*(yL(2)+yL(3));  % dy1/dt
%     -1.0j*Delta*yL(2)-1.0j*0.5*Omega1*exp(1.0j*phi1)*yL(1)-1.0j*0.5*Omega2*exp(-1.0j*phi2)*exp(-1.0j*theta)*yL(4);  % dy2/dt
%     -1.0j*Delta*yL(3)-1.0j*0.5*Omega1*exp(1.0j*phi1)*yL(1)+1.0j*0.5*Omega2*exp(-1.0j*phi2)*exp(-1.0j*thetacL)*yL(4); % dy3/dt
%     -1.0j*delta*yL(4)-1.0j*0.5*Omega2*exp(1.0j*phi2)*exp(1.0j*theta)*yL(2)+1.0j*0.5*Omega2*exp(1.0j*phi2)*exp(1.0j*thetacL)*yL(3) % dy4/dt
% ];
% 
% odefunR = @(tR, yR) [
%     -1.0j*0.5*Omega1*exp(-1.0j*phi1)*(yR(2)+yR(3));  % dy1/dt
%     -1.0j*Delta*yR(2)-1.0j*0.5*Omega1*exp(1.0j*phi1)*yR(1)-1.0j*0.5*Omega2*exp(-1.0j*phi2)*exp(-1.0j*theta)*yR(4);  % dy2/dt
%     -1.0j*Delta*yR(3)-1.0j*0.5*Omega1*exp(1.0j*phi1)*yR(1)+1.0j*0.5*Omega2*exp(-1.0j*phi2)*exp(-1.0j*thetacR)*yR(4); % dy3/dt
%     -1.0j*delta*yR(4)-1.0j*0.5*Omega2*exp(1.0j*phi2)*exp(1.0j*theta)*yR(2)+1.0j*0.5*Omega2*exp(1.0j*phi2)*exp(1.0j*thetacR)*yR(3) % dy4/dt
% ];

% 初始条件
tspan = [0, 20];         % 时间范围
yL0 = [1; 0; 0; 0]; 
yR0 = [1; 0; 0; 0];

% 求解微分方程组
[tL, yL] = ode45(odefunL, tspan, yL0);
[tR, yR] = ode45(odefunR, tspan, yR0);

% 归一化系数
CL= abs(yL(:, 1)).^2+abs(yL(:, 2)).^2+abs(yL(:, 3)).^2+abs(yL(:, 4)).^2;
CR= abs(yR(:, 1)).^2+abs(yR(:, 2)).^2+abs(yR(:, 3)).^2+abs(yR(:, 4)).^2;

% % 计算模平方
P4L = (abs(yL(:, 4)).^2)./CL; %PL_gamma
P4R = (abs(yR(:, 4)).^2)./CR; %PR_gamma

% 绘图
subplot(2,1,1)
h1 = plot(tL, P4L, 'r-', 'DisplayName', 'P4L');
axis([0 20 0 1]);%x轴y轴范围
xlabel('$t\,(\mu s)$','interpreter', 'latex','FontSize',10);
ylabel('$P_{\gamma}^{(L)}$','interpreter', 'latex','FontSize',10);
%legend('show');
grid on;

subplot(2,1,2)
h2 = plot(tR, P4R, 'b-', 'DisplayName', 'P4R');
axis([0 20 0 1]);%x轴y轴范围
xlabel('$t\,(\mu s)$','interpreter', 'latex','FontSize',10);
ylabel('$P_{\gamma}^{(R)}$','interpreter', 'latex','FontSize',10);
%legend('show');
grid on;

toc; % 结束计时并显示经过的时间
set(gcf, 'PaperSize', [width_cm height_cm]);
%saveas(gcf, '/Users/zoufen/Desktop/Fig3a01.pdf')

% 提取数据
xData1 = get(h1, 'XData');
yData1 = get(h1, 'YData');
xData2 = get(h2, 'XData');
yData2 = get(h2, 'YData');

% 保存到 CSV 文件
dataL = [xData1' yData1'];
dataR = [xData2' yData2'];
writematrix(dataL, '/Users/zoufen/Desktop/plot_data_thetaR_L.csv');
writematrix(dataR, '/Users/zoufen/Desktop/plot_data_thetaR_R.csv');



function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
if (nargin<4) || isempty(gap),    gap=0.01;  end
if (nargin<5) || isempty(marg_h),  marg_h=0.05;  end
if (nargin<5) || isempty(marg_w),  marg_w=marg_h;  end
if isscalar(gap),   gap(2)=gap;  end
if isscalar(marg_h),  marg_h(2)=marg_h;  end
if isscalar(marg_w),  marg_w(2)=marg_w;  end
gap_vert   = gap(1);
gap_horz   = gap(2);
marg_lower = marg_h(1);
marg_upper = marg_h(2);
marg_left  = marg_w(1);
marg_right = marg_w(2);
[subplot_col,subplot_row]=ind2sub([n,m],p);  
subplot_cols=1+max(subplot_col)-min(subplot_col);
subplot_rows=1+max(subplot_row)-min(subplot_row);
height=(1-(marg_lower+marg_upper)-(m-1)*gap_vert)/m;
width =(1-(marg_left+marg_right)-(n-1)*gap_horz)/n;
merged_height=subplot_rows*( height+gap_vert )- gap_vert;
merged_width= subplot_cols*( width +gap_horz )- gap_horz;
merged_bottom=(m-max(subplot_row))*(height+gap_vert) +marg_lower;
merged_left=(min(subplot_col)-1)*(width+gap_horz) +marg_left;
pos_vec=[merged_left merged_bottom merged_width merged_height];
h=subplot('Position',pos_vec,varargin{:});
if (nargout < 1),  clear h;  end
end
