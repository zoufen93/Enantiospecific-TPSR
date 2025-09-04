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
Delta2 = 2*pi*0.4; %MHz  Delta2
Delta1 = 2*pi*0.1; %MHz  Delta1
Omega1 = 2*pi*1; %MHz
Omega2 = 2*pi*1; %MHz
thetacL = pi/2;
thetacR = -pi/2;
theta = thetacL;
tlist = 0:0.01:400;

Gamma2 = 0; %MHz

% 定义态矢量
p = 4;
Id = eye(p);
ket = @(i) Id(:, i);
bra = @(i) Id(i, :);
aket = ket(1); abra = bra(1);
bpket = ket(2); bpbra = bra(2);
bmket = ket(3); bmbra = bra(3);
gket = ket(4); gbra = bra(4);

% 定义 Hamiltonians
HL = Delta1*(bpket*bpbra+bmket*bmbra)+(Delta1+Delta2)*gket*gbra...
     +0.5*Omega1*(bpket*abra+aket*bpbra+bmket*abra+aket*bmbra)...
     +0.5*Omega2*(exp(1.0j*theta)*gket*bpbra+bpket*gbra*exp(-1.0j*theta))...
     -0.5*Omega2*(exp(1.0j*thetacL)*gket*bmbra+bmket*gbra*exp(-1.0j*thetacL));

HR = Delta1*(bpket*bpbra+bmket*bmbra)+(Delta1+Delta2)*gket*gbra...
     +0.5*Omega1*(bpket*abra+aket*bpbra+bmket*abra+aket*bmbra)...
     +0.5*Omega2*(exp(1.0j*theta)*gket*bpbra+bpket*gbra*exp(-1.0j*theta))...
     -0.5*Omega2*(exp(1.0j*thetacR)*gket*bmbra+bmket*gbra*exp(-1.0j*thetacR));
 
H1p = bpket*bpbra+bmket*bmbra+gket*gbra;
H2p = gket*gbra;

% Liouvillian
L0_commL = -1.0j*kron(eye(p),HL) + 1.0j*kron(HL.',eye(p));
L0_commR = -1.0j*kron(eye(p),HR) + 1.0j*kron(HR.',eye(p));
L1p_comm = -1.0j*kron(eye(p),H1p) + 1.0j*kron(H1p.',eye(p));
L2p_comm = -1.0j*kron(eye(p),H2p) + 1.0j*kron(H2p.',eye(p));

% 初始条件 (interaction picture)
%rho0_tilde = rho0;  % 因为 t=0 时 e^{-L0*0} = I
rho0 = aket*abra;
rho0_vec = rho0(:); %列向量化

% 二阶微扰项
L1p_tilde_L_rho2 = @(t1) expm(-L0_commL * t1)* L1p_comm * L1p_comm * expm(L0_commL * t1) * rho0_vec; 
L1p_tilde_R_rho2 = @(t1) expm(-L0_commR * t1)* L1p_comm * L1p_comm * expm(L0_commR * t1) * rho0_vec;
L2p_tilde_L_rho2 = @(t1) expm(-L0_commL * t1)* L2p_comm * L2p_comm * expm(L0_commL * t1) * rho0_vec;
L2p_tilde_R_rho2 = @(t1) expm(-L0_commR * t1)* L2p_comm * L2p_comm * expm(L0_commR * t1) * rho0_vec;

%beta0的范围
Gamma1_vals = linspace(0, 0.1, 51);

% 初始化数据表
PL = zeros(length(Gamma1_vals), 1);
PR = zeros(length(Gamma1_vals), 1);
ee = zeros(length(Gamma1_vals), 1);

% 开启并行
parpool('local');

parfor i = 1:length(Gamma1_vals)
    Gamma1  = 2*pi*Gamma1_vals(i);
    
    
    rhoL_vec_2 = 2*Gamma1*integral(L1p_tilde_L_rho2, 0, 320, ...
        'ArrayValued', true, 'AbsTol', 1e-12, 'RelTol', 1e-9) ...
        + 2*Gamma2*integral(L2p_tilde_L_rho2, 0, 320, ...
        'ArrayValued', true, 'AbsTol', 1e-12, 'RelTol', 1e-9);

    rhoR_vec_2 = 2*Gamma1*integral(L1p_tilde_R_rho2, 0, 320, ...
        'ArrayValued', true, 'AbsTol', 1e-12, 'RelTol', 1e-9) ...
        + 2*Gamma2*integral(L2p_tilde_R_rho2, 0, 320, ...
        'ArrayValued', true, 'AbsTol', 1e-12, 'RelTol', 1e-9);
   
    rhoL_vec_2 = real(rhoL_vec_2);
    rhoR_vec_2 = real(rhoR_vec_2);

    % 预计算 rho^{t} 的矩阵指数方法
    rhotp_cache_L = cell(length(tlist), 1);
    rhotp_cache_R = cell(length(tlist), 1); 
    for k = 1:length(tlist)
        t = tlist(k);
        rhoL_withNoise = expm(L0_commL * t)*(rho0_vec + rhoL_vec_2);
        rhoR_withNoise = expm(L0_commR * t)*(rho0_vec + rhoR_vec_2);
        
        % === 归一化 + 防负值修正 ===
        rhotp_cache_L{k} = enforcePhysicalState(rhoL_withNoise, p);
        rhotp_cache_R{k} = enforcePhysicalState(rhoR_withNoise, p);
    end
    
    % 计算 P_L(t), P_R(t)
    PL_t = zeros(size(tlist));
    PR_t = zeros(size(tlist));

    for k = 1:length(tlist)
        % 还原矩阵形式
        rhoL_t = reshape(rhotp_cache_L{k}, p, p);
        rhoR_t = reshape(rhotp_cache_R{k}, p, p);

        % 假设 P_L(t) / P_R(t) 定义为在 |g> 态的概率
        PL_t(k) = real(gbra * rhoL_t * gket);
        PR_t(k) = real(gbra * rhoR_t * gket);
    end
    
    % 计算时间平均
    Tmax = tlist(end);
    PL_avg = trapz(tlist, PL_t) / Tmax;
    PR_avg = trapz(tlist, PR_t) / Tmax;
    PL(i) = PL_avg;
    PR(i) = PR_avg;
    ee(i) = abs(PL_avg-PR_avg)/(PL_avg+PR_avg);
end

delete(gcp('nocreate'));

subplot(2,1,1)
h1 =plot(Gamma1_vals, PL(:), 'r-', 'DisplayName', 'PgammaL');
hold on;
h2 =plot(Gamma1_vals, PR(:), 'b--', 'DisplayName', 'PgammaL');
hold off;
legd=legend('L','R','location','northeast');%图例
set(legd,'Box','off','interpreter', 'latex','FontName','Times New Roman','FontSize',10);%设置图例
legd.ItemTokenSize = [15,10];
%axis([0 4 0 1]);%x轴y轴范围
title('$\Gamma_{2}/2\pi=0.01\,\mathrm{MHz}$', 'interpreter', 'latex', 'FontSize', 10);
xlabel('$\Gamma_{1}/2\pi$','interpreter', 'latex','FontSize',10);
ylabel('$\bar{P}_{\gamma}^{(L,R)}$','interpreter', 'latex','FontSize',10);
%legend('show');
grid on;

subplot(2,1,2)
h3 = plot(Gamma1_vals, ee(:), 'r-', 'DisplayName', 'PgammaL');
hold on;
% legd=legend('Without phase noise','With phase noise','location','northeast');%图例
% set(legd,'Box','off','interpreter', 'latex','FontName','Times New Roman','FontSize',10);%设置图例
% legd.ItemTokenSize = [15,10];
axis([0 0.1 0 1]);%x轴y轴范围
%title('$\Gamma_{2}/2\pi=0.01\,\mathrm{MHz}$', 'interpreter', 'latex', 'FontSize', 10);
xlabel('$\Gamma_{1}/2\pi$','interpreter', 'latex','FontSize',10);
ylabel('$v_{\mathrm{EPT}}$','interpreter', 'latex','FontSize',10);
%legend('show');
grid on;

toc; % 结束计时并显示经过的时间

%%% Save figure
set(gcf, 'PaperSize', [width_cm height_cm]);
saveas(gcf, '/Users/zoufen/Desktop/vEPTvsGamma1.pdf')

% 提取数据
xData1 = get(h1, 'XData');
yData1 = get(h1, 'YData');
yData2 = get(h2, 'YData');
yData3 = get(h3, 'YData');

% 保存到 CSV 文件
data = [xData1' yData1' yData2'];
writematrix(data, '/Users/zoufen/Desktop/vEPTvsGamma1.csv');

% ============ 工具函数 ============
function rho_vec_fixed = enforcePhysicalState(rho_vec, p)
    rho_mat = reshape(rho_vec, p, p);

    % 保证厄米性
    rho_mat = (rho_mat + rho_mat')/2;

    % 去掉虚部残差
    rho_mat = real(rho_mat);

    % 负特征值截断
    [V,D] = eig(rho_mat);
    D = diag(D);
    D(D < 0) = 0;
    rho_mat = V*diag(D)*V';

    % 归一化
    tr_val = trace(rho_mat);
    if tr_val ~= 0
        rho_mat = rho_mat / tr_val;
    else
        rho_mat = eye(p)/p;
    end

    rho_vec_fixed = reshape(rho_mat, [], 1);
end

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