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
Delta2 = -2*pi*10.2; %MHz  Delta2
Delta1 = 2*pi*10; %MHz  Delta1
Omega1 = 2*pi*1; %MHz
Omega2 = 2*pi*1; %MHz
thetacL = pi/2;
thetacR = -pi/2;
theta = thetacL;
tlist = 0:0.01:320;

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

L0_commL = liouvillian_commutator(HL);
L0_commR = liouvillian_commutator(HR);

%beta0的范围
beta0_vals = linspace(0, 0.5, 51);

% 初始化数据表
% PL = zeros(length(beta0_vals), 1);
% PR = zeros(length(beta0_vals), 1);
ee = zeros(length(beta0_vals), 1);

% 开启并行
parpool('local');

parfor i = 1:length(beta0_vals)
    beta0  = beta0_vals(i);
    
    % 初始条件 (interaction picture)
    rho0 = (1-beta0)*aket*abra + 0.5*beta0*bpket*bpbra + 0.5*beta0*bmket*bmbra;
    rho0_vec = rho0(:); %列向量化

    % 预计算 rho^{t} 的矩阵指数方法
    rhot_cache_L = cell(length(tlist), 1); 
    rhot_cache_R = cell(length(tlist), 1); 
    for k = 1:length(tlist)
        t = tlist(k);
        rhot_cache_L{k} = expm(L0_commL*t)*rho0_vec;
        rhot_cache_R{k} = expm(L0_commR*t)*rho0_vec;
    end 

    % 计算 P_L(t), P_R(t) 
    PL_t = zeros(size(tlist));
    PR_t = zeros(size(tlist));

    for k = 1:length(tlist)
        % 还原矩阵形式
        rhoL_t = reshape(rhot_cache_L{k}, p, p);
        rhoR_t = reshape(rhot_cache_R{k}, p, p);

        % 假设 P_L(t) / P_R(t) 定义为在 |g> 态的概率
        PL_t(k) = real(gbra * rhoL_t * gket);
        PR_t(k) = real(gbra * rhoR_t * gket);
    end
    % 计算时间平均
    Tmax = tlist(end);
    PL_avg = trapz(tlist, PL_t) / Tmax;
    PR_avg = trapz(tlist, PR_t) / Tmax;
    ee(i) = abs(PL_avg-PR_avg)/(PL_avg+PR_avg);
end

delete(gcp('nocreate'));

% 绘图
subplot(1,1,1)
h1 = plot(beta0_vals, ee(:), 'r-', 'DisplayName', 'PgammaL');
axis([0 0.5 0 1]);%x轴y轴范围
%title('$\beta_0$', 'interpreter', 'latex', 'FontSize', 10);
xlabel('$\beta_0$','interpreter', 'latex','FontSize',10);
ylabel('$v_{\mathrm{EPT}}$','interpreter', 'latex','FontSize',10);
%legend('show');
grid on;


toc; % 结束计时并显示经过的时间
%set(gcf, 'PaperSize', [width_cm height_cm]);
%saveas(gcf, '/Users/zoufen/Desktop/vEPTvsbeta0.pdf')

% 提取数据
xData1 = get(h1, 'XData');
yData1 = get(h1, 'YData');

% 保存到 CSV 文件
dataL = [xData1' yData1'];
writematrix(dataL, '/Users/zoufen/Desktop/vEPTvsbeta0.csv');

function L = liouvillian_commutator(H)
% LIOUVILLIAN_COMMUTATOR  返回 -i[H, ·] 的超算符矩阵（列向量化）
%   H : d x d Hamiltonian (Hermitian)
%   L : d^2 x d^2 matrix such that vec(-i[H,rho]) = L * vec(rho)

d = size(H,1);
L = -1i * (kron(eye(d), H) - kron(H.', eye(d)));
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
