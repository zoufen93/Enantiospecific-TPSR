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
hbar = 1.055*10^(-34);
n = 3.22*10^(18);%m^-3
d0 = 3*3.33564*10^(-30); %C*m
ep0 = 8.854*10^(-12); %C^2/(N*m)
omega = 2*pi*6*10^(14); %Hz
c = 3*10^(8);%m/s
Delta2 = 2*pi*0.4*10^6/sqrt(2); %Hz  delta = Delta_1 + Delta_2
Delta1 = 2*pi*0.1*10^6; %Hz     Delta1
Omega1 = 2*pi*1*10^6; %Hz
Omega2 = 2*pi*0.01*10^6; %Hz
kappa1 = 2*pi*10*10^6;
kappa2 = 2*pi*10*10^6;
thetacL = pi/2;
thetacR = -pi/2;
npoint = 201;
step = 0.01;
%nth = 0;

% 初始化数据表
thetaLdata = zeros(npoint,2);
thetaRdata = zeros(npoint,2);

% 定义态矢量
p = 4;
Id = eye(p);
ket = @(i) Id(:, i);
bra = @(i) Id(i, :);
aket = ket(1); abra = bra(1);
bpket = ket(2); bpbra = bra(2);
bmket = ket(3); bmbra = bra(3);
gket = ket(4); gbra = bra(4);

% 符号矩阵
rhoLss = sym('rhoLss', [p p], 'real');
rhoRss = sym('rhoRss', [p p], 'real');
    
    
% 把符号矩阵改成变量形式
for i = 1:npoint
    theta = pi*(i*step-1.01);

    % 定义 Hamiltonians
    HL = Delta1*(bpket*bpbra+bmket*bmbra)+(Delta1+Delta2)*gket*gbra ...
         +0.5*Omega1*(bpket*abra+aket*bpbra+bmket*abra+aket*bmbra) ...
         +0.5*Omega2*(exp(1.0j*theta)*gket*bpbra+bpket*gbra*exp(-1.0j*theta)) ...
         -0.5*Omega2*(exp(1.0j*thetacL)*gket*bmbra+bmket*gbra*exp(-1.0j*thetacL));

    HR = Delta1*(bpket*bpbra+bmket*bmbra)+(Delta1+Delta2)*gket*gbra ...
         +0.5*Omega1*(bpket*abra+aket*bpbra+bmket*abra+aket*bmbra) ...
         +0.5*Omega2*(exp(1.0j*theta)*gket*bpbra+bpket*gbra*exp(-1.0j*theta)) ...
         -0.5*Omega2*(exp(1.0j*thetacR)*gket*bmbra+bmket*gbra*exp(-1.0j*thetacR));

    % Liouvillian
    L0L = -1.0j*kron(eye(p),HL) + 1.0j*kron(HL.',eye(p));
    L0R = -1.0j*kron(eye(p),HR) + 1.0j*kron(HR.',eye(p));

    % 添加耗散
    C_ops = {bpket*gbra, bmket*gbra, aket*gbra, aket*bpbra, aket*bmbra};
    kappa_vals = [kappa1, kappa1, kappa1, kappa2, kappa2]; % 每个C的耗散速率
    for cc = 1:length(C_ops)
        C = C_ops{cc};
        kappa_c = kappa_vals(cc); % 当前C对应的kappa
        CL = kappa_c/2 * (2*kron(conj(C), C) ...
             - kron(eye(p), C'*C) - kron((C'*C).', eye(p)));
        L0L = L0L + CL;
        L0R = L0R + CL;
    end
    
    % ---- 稳态方程 L*rho = 0  + trace(rho)=1 ----
    % 归一化行向量：1×p^2
    norm_row = reshape(eye(p), 1, []);   % 等价于 vec(I)' 

    % 左体系
    A_L = [L0L; norm_row];               % (p^2+1) × p^2
    b_L = [zeros(p^2,1); 1];             % (p^2+1) × 1
    rhoL_vec = A_L\b_L;                  % p^2 × 1
    rhoLss_sol = reshape(rhoL_vec, p, p);

    % 右体系
    A_R = [L0R; norm_row];
    b_R = [zeros(p^2,1); 1];
    rhoR_vec = A_R\b_R;
    rhoRss_sol = reshape(rhoR_vec, p, p);

    % ---- 数值清理与兜底，保证不出现 NaN/非厄米/trace误差 ----
    % 逼近厄米
    rhoLss_sol = (rhoLss_sol + rhoLss_sol')/2;
    rhoRss_sol = (rhoRss_sol + rhoRss_sol')/2;

    % 归一化 trace=1（防止最小二乘误差）
    trL = trace(rhoLss_sol);
    trR = trace(rhoRss_sol);
    if abs(trL) < 1e-12 || ~isfinite(trL), rhoLss_sol = eye(p)/p; else, rhoLss_sol = rhoLss_sol/real(trL); end
    if abs(trR) < 1e-12 || ~isfinite(trR), rhoRss_sol = eye(p)/p; else, rhoRss_sol = rhoRss_sol/real(trR); end

    % 若仍有非有限元素，则用均匀混合态兜底
    if any(~isfinite(rhoLss_sol(:))), rhoLss_sol = eye(p)/p; end
    if any(~isfinite(rhoRss_sol(:))), rhoRss_sol = eye(p)/p; end

    % 提取元
    rhoLgbp = rhoLss_sol(4,2);
    rhoLgbm = rhoLss_sol(4,3);
    rhoRgbp = rhoRss_sol(4,2);
    rhoRgbm = rhoRss_sol(4,3);

    % 存储数据
    thetaLdata(i, :) = [theta/pi, ((omega/c)*(n*d0^(2)/ep0)*imag(2*exp(-1.0j*thetacL)*rhoLgbm - 2*exp(-1.0j*theta)*rhoLgbp))/(hbar*Omega2)];
    thetaRdata(i, :) = [theta/pi, ((omega/c)*(n*d0^(2)/ep0)*imag(2*exp(-1.0j*thetacR)*rhoRgbm - 2*exp(-1.0j*theta)*rhoRgbp))/(hbar*Omega2)];
end

subplot(1,1,1)
h1= plot(thetaLdata(:, 1), thetaLdata(:, 2), '-');
hold on;
h2= plot(thetaRdata(:, 1), thetaRdata(:, 2), '--');
hold off;
legd=legend('$L$','$R$','location','southeast');%图例
set(legd,'Box','off','interpreter', 'latex','FontName','Times New Roman','FontSize',10);%设置图例
legd.ItemTokenSize = [15,10];
%title('$|a\rangle=|1,0\rangle$, $|b\rangle=|4,-1\rangle$, $|c\rangle=|2,0\rangle$, $|d\rangle=|4,1\rangle$', 'interpreter', 'latex', 'FontSize', 10);
%axis([-pi pi 0 1.1]);%x轴y轴范围
xlabel('$\theta/\pi$','interpreter', 'latex','FontSize',10);%设置x轴标签
ylabel('$A(\theta)\,[m^{-1}]$','interpreter', 'latex','FontSize',10);%设置y轴标签
set(gca,'linewidth',0.8);%设置坐标轴线宽
%set(gca, 'xtick', [0 40 80], 'xticklabels', {'0' '40' '80'},'FontName','Times New Roman','fontsize',10)
%set(gca, 'ytick', [-1 -0.5 0 0.5 1], 'yticklabels', {'-1' '-0.5' '0' '0.5' '1'},'FontName','Times New Roman','fontsize',10)
%set(gca,'yminortick','on');%设置小刻度标签
grid on;

toc; % 结束计时并显示经过的时间

%%% Save figure
%set(gcf, 'PaperSize', [width_cm height_cm]);
%saveas(gcf, '/Users/zoufen/Desktop/Fig3d.pdf')

% 提取数据
xData1 = get(h1, 'XData');
yData1 = get(h1, 'YData');
yData2 = get(h2, 'YData');

% 保存到 CSV 文件
data = [xData1' yData1' yData2'];
writematrix(data, '/Users/zoufen/Desktop/Avstheta.csv');



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