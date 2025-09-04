close all;
clear;
clc;
figure
width_cm = 10;
height_cm = 10;
set(gcf,'unit','centimeters','position',[20 5 width_cm height_cm])
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [1.6/height_cm 1.45/width_cm],[1.2/height_cm 0.8/height_cm], [1.5/width_cm 0.3/width_cm]);
if ~make_it_tight,  clear subplot;  end       %%%%%%%%%% gap [h, w], h [down, up], w [left right]

tic;% 开始计时
%参数
Delta1 = 2*pi*0.1; %MHz  Delta1
Delta2 = 2*pi*0.4; %MHz  Delta1
Omega1 = 2*pi*1; %MHz
Omega2 = 2*pi*1; %MHz
% Delta1 = 0; %MHz  Delta1
% Delta2 = -(2*pi*50)/sqrt(2); %MHz  Delta1
% Omega1 = 2*pi*50; %MHz
% Omega2 = 2*pi*0.5; %MHz
thetacL = pi/2;
thetacR = -pi/2;

% dtheta的范围
dtheta_vals = linspace(-pi/2, pi/2, 501);

%预分配数组
PL = zeros(length(dtheta_vals), 1);
PR = zeros(length(dtheta_vals), 1);
ee = zeros(length(dtheta_vals), 1);

for i = 1:length(dtheta_vals)
    theta  = thetacL + dtheta_vals(i);

    % 定义微分方程组
    odefunL = @(tL, yL) [
        -1.0j*0.5*Omega1*(yL(2)+yL(3));  % dy1/dt
        -1.0j*Delta1*yL(2)-1.0j*0.5*Omega1*yL(1)-1.0j*0.5*Omega2*exp(-1.0j*theta)*yL(4);  % dy2/dt
        -1.0j*Delta1*yL(3)-1.0j*0.5*Omega1*yL(1)+1.0j*0.5*Omega2*exp(-1.0j*thetacL)*yL(4); % dy3/dt
        -1.0j*(Delta1+Delta2)*yL(4)-1.0j*0.5*Omega2*exp(1.0j*theta)*yL(2)+1.0j*0.5*Omega2*exp(1.0j*thetacL)*yL(3) % dy4/dt
    ];

    odefunR = @(tR, yR) [
        -1.0j*0.5*Omega1*(yR(2)+yR(3));  % dy1/dt
        -1.0j*Delta1*yR(2)-1.0j*0.5*Omega1*yR(1)-1.0j*0.5*Omega2*exp(-1.0j*theta)*yR(4);  % dy2/dt
        -1.0j*Delta1*yR(3)-1.0j*0.5*Omega1*yR(1)+1.0j*0.5*Omega2*exp(-1.0j*thetacR)*yR(4); % dy3/dt
        -1.0j*(Delta1+Delta2)*yR(4)-1.0j*0.5*Omega2*exp(1.0j*theta)*yR(2)+1.0j*0.5*Omega1*exp(1.0j*thetacR)*yR(3) % dy4/dt
    ];

    % 初始条件
    tspan = [0,320]; %T 应该足够大以近似无穷大
    yL0 = [1; 0; 0; 0]; 
    yR0 = [1; 0; 0; 0];

    % 求解微分方程组
    [tL, yL] = ode45(odefunL, tspan, yL0);
    [tR, yR] = ode45(odefunR, tspan, yR0);

    % 归一化系数
    CL = sum(abs(yL).^2, 2);
    CR = sum(abs(yR).^2, 2);

    % 计算模平方
    P4L = (abs(yL(:, 4)).^2)./CL; %PL_gamma
    P4R = (abs(yR(:, 4)).^2)./CR; %PR_gamma
    
    % 使用 trapz 进行数值积分并计算时间平均
    PL = trapz(tL, P4L) / tL(end); % 使用 tL 的最后一个元素作为 T
    PR = trapz(tR, P4R) / tR(end); % 使用 tR 的最后一个元素作为 T
    ee(i) = abs(PL-PR)/(PL+PR);
end


% 绘图
subplot(1,1,1)
h1 = plot(dtheta_vals, ee(:), 'r-', 'DisplayName', 'P4L');
axis([-pi/2 pi/2 0 1]);%x轴y轴范围
xlabel('$\delta\theta$','interpreter', 'latex','FontSize',10);
ylabel('$v_{\mathrm{EPT}}$','interpreter', 'latex','FontSize',10);
set(gca, 'xtick', [-pi/2 -pi/4 0 pi/4 pi/2], 'xticklabels', {'-\pi/2' '-\pi/4' '0' '\pi/4' '\pi/2'},'FontName','Times New Roman','fontsize',10)
set(gca, 'ytick', [0 0.5 1], 'yticklabels', {'0' '0.5' '1'},'FontName','Times New Roman','fontsize',10)
%legend('show');
grid on;

toc; % 结束计时并显示经过的时间
%set(gcf, 'PaperSize', [width_cm height_cm]);
%saveas(gcf, '/Users/zoufen/Desktop/vEPTvsdtheta.pdf')

% 提取数据
xData1 = get(h1, 'XData');
yData1 = get(h1, 'YData');

% 保存到 CSV 文件
%data = [xData1' yData1'];
%writematrix(data, '/Users/zoufen/Desktop/vEPTvsdtheta.csv');

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