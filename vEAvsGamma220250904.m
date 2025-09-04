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
hbar = 1.055*10^(-34);
n = 3.22*10^(18);%m^-3
d0 = 3*3.33564*10^(-30); %C*m
ep0 = 8.854*10^(-12); %C^2/(N*m)
omega = 2*pi*6*10^(14); %Hz
c = 3*10^(8);%m/s
Delta1 = 2*pi*0.1*10^(6); %MHz  
Delta2 = 2*pi*0.4*10^(6); %MHz  
Omega1 = 2*pi*1*10^(6); %MHz
Omega2 = 2*pi*0.01*10^(6); %MHz
kappa1 = 2*pi*10*10^6;
kappa2 = 2*pi*10*10^6;
thetacL = pi/2;
thetacR = -pi/2;
theta = thetacL;
tlist = 0:0.01:320; % us（单位只影响横轴刻度，不影响算子）

% 扰动强度扫描 (MHz -> Hz)
Gamma1        = 0;           % Hz
Gamma2_vals_M = linspace(0, 0.1, 21);   % MHz
Gamma2_vals   = 2*pi*Gamma2_vals_M*10^(6);

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

% 添加耗散
C_ops = {bpket*gbra, bmket*gbra, aket*gbra, aket*bpbra, aket*bmbra};
kappa_vals = [kappa1, kappa1, kappa1, kappa2, kappa2]; % 每个C的耗散速率
for cc = 1:length(C_ops)
    C = C_ops{cc};
    kappa_c = kappa_vals(cc); % 当前C对应的kappa
    CL = kappa_c/2 * (2*kron(conj(C), C) ...
         - kron(eye(p), C'*C) - kron((C'*C).', eye(p)));
    L0_commL = L0_commL + CL;
    L0_commR = L0_commR + CL;
end
    
% 初始条件 (interaction picture)
%rho0_tilde = rho0;  % 因为 t=0 时 e^{-L0*0} = I
rho0 = aket*abra;
rho0_vec = rho0(:); %列向量化 p^2 × 1

% 初始化数据表
AL = zeros(length(Gamma2_vals), 1);
AR = zeros(length(Gamma2_vals), 1);
ee = zeros(length(Gamma2_vals), 1);

% ===== 并行池（不覆盖变量名）=====
if isempty(gcp('nocreate'))
    parpool('local');
end

parfor ii = 1:numel(Gamma2_vals)

    Gamma2 = Gamma2_vals(ii);   % Hz

    % ---- 安全积分，得到二阶修正向量 ----
    L1p_tilde_L_rho2 = @(t1) safeExpmTimes(-L0_commL,t1, L1p_comm*L1p_comm,  L0_commL, t1, rho0_vec);
    L1p_tilde_R_rho2 = @(t1) safeExpmTimes(-L0_commR,t1, L1p_comm*L1p_comm,  L0_commR, t1, rho0_vec);
    L2p_tilde_L_rho2 = @(t1) safeExpmTimes(-L0_commL,t1, L2p_comm*L2p_comm,  L0_commL, t1, rho0_vec);
    L2p_tilde_R_rho2 = @(t1) safeExpmTimes(-L0_commR,t1, L2p_comm*L2p_comm,  L0_commR, t1, rho0_vec);

    rhoL_vec_2 = 2*Gamma1*safeIntegral(L1p_tilde_L_rho2, 0, 320) ...
               + 2*Gamma2*safeIntegral(L2p_tilde_L_rho2, 0, 320);
    rhoR_vec_2 = 2*Gamma1*safeIntegral(L1p_tilde_R_rho2, 0, 320) ...
               + 2*Gamma2*safeIntegral(L2p_tilde_R_rho2, 0, 320);

    % 数值清理
    rhoL_vec_2 = real(rhoL_vec_2);
    rhoR_vec_2 = real(rhoR_vec_2);

    % ---- 随时间演化并修正到物理态 ----
    rhotp_cache_L_end = zeros(p^2,1);
    rhotp_cache_R_end = zeros(p^2,1);

    for kk = 1:numel(tlist)
        t = tlist(kk);

        vL = safeExpmTimes(L0_commL, t, [], [], 0, rho0_vec + rhoL_vec_2);
        vR = safeExpmTimes(L0_commR, t, [], [], 0, rho0_vec + rhoR_vec_2);

        % 防 NaN/Inf + Hermitian + 正定 + 归一化
        vL = enforcePhysicalState(vL, p);
        vR = enforcePhysicalState(vR, p);

        % 只需最终时刻
        if kk == numel(tlist)
            rhotp_cache_L_end = vL;
            rhotp_cache_R_end = vR;
        end
    end

    % ---- 恢复矩阵、取元并计算 A^(L/R) ----
    rhoLss = reshape(rhotp_cache_L_end, p, p);
    rhoRss = reshape(rhotp_cache_R_end, p, p);

    rhoLgbp = rhoLss(4,2);  % <g|rho|bp>
    rhoLgbm = rhoLss(4,3);  % <g|rho|bm>
    rhoRgbp = rhoRss(4,2);
    rhoRgbm = rhoRss(4,3);

    pref = (omega/c)*(n*d0^2/ep0)/(hbar*Omega2);
    aL   = pref * imag( 2*exp(-1i*thetacL)*rhoLgbm - 2*exp(-1i*theta)*rhoLgbp );
    aR   = pref * imag( 2*exp(-1i*thetacR)*rhoRgbm - 2*exp(-1i*theta)*rhoRgbp );

    % 防 NaN/Inf（极端情况下退回 0）
    if ~isfinite(aL), aL = 0; end
    if ~isfinite(aR), aR = 0; end

    AL(ii) = aL;
    AR(ii) = aR;

    ee(ii) = abs(aL-aR)/(aL+aR);
end

subplot(2,1,1)
h1 = plot(Gamma2_vals_M, AL(:), 'r-', 'DisplayName', 'PgammaL');
hold on;
h2 = plot(Gamma2_vals_M, AR(:), 'b--', 'DisplayName', 'PgammaL');
hold off;
legd=legend('L','R','location','northeast');%图例
set(legd,'Box','off','interpreter', 'latex','FontName','Times New Roman','FontSize',10);%设置图例
legd.ItemTokenSize = [15,10];
%axis([0 0.1 0 2000]);%x轴y轴范围
%title('$\Gamma_{2}/2\pi=0.01\,\mathrm{MHz}$', 'interpreter', 'latex', 'FontSize', 10);
xlabel('$\Gamma_{1}/2\pi$','interpreter', 'latex','FontSize',10);
ylabel('$A^{(L,R)}$','interpreter', 'latex','FontSize',10);
%legend('show');
grid on;

subplot(2,1,2)
h3 = plot(Gamma2_vals_M, ee(:), 'r-', 'DisplayName', 'PgammaL');
hold on;
% legd=legend('Without phase noise','With phase noise','location','northeast');%图例
% set(legd,'Box','off','interpreter', 'latex','FontName','Times New Roman','FontSize',10);%设置图例
% legd.ItemTokenSize = [15,10];
axis([0 0.1 0 1]);%x轴y轴范围
%title('$\Gamma_{2}/2\pi=0.01\,\mathrm{MHz}$', 'interpreter', 'latex', 'FontSize', 10);
xlabel('$\Gamma_{1}/2\pi$','interpreter', 'latex','FontSize',10);
ylabel('$v_{\mathrm{EA}}$','interpreter', 'latex','FontSize',10);
%legend('show');
grid on;

toc; % 结束计时并显示经过的时间

%%% Save figure
%set(gcf, 'PaperSize', [width_cm height_cm]);
%saveas(gcf, '/Users/zoufen/Desktop/vEAvsGamma2.pdf')

% 提取数据
xData1 = get(h1, 'XData');
yData1 = get(h1, 'YData');
yData2 = get(h2, 'YData');
yData3 = get(h3, 'YData');

% 保存到 CSV 文件
data = [xData1' yData1' yData2'];
writematrix(data, '/Users/zoufen/Desktop/vEAvsGamma2.csv');

% =========================
% ===== 本地工具函数 =====
% =========================

function rho_vec_fixed = enforcePhysicalState(rho_vec, p)
% 将长度为 p^2 的列向量修正为密度矩阵：
% 1) NaN/Inf 拦截 -> I/p
% 2) Hermitian 化
% 3) 特征值截断为 >=0（防负）
% 4) 归一化 trace=1

    if any(~isfinite(rho_vec))
        rho_vec_fixed = (eye(p)/p);
        rho_vec_fixed = rho_vec_fixed(:);
        return;
    end

    rho = reshape(rho_vec, p, p);

    % Hermitian
    rho = (rho + rho')/2;

    % Eig 修正到 PSD
    [V,D] = eig((rho+rho')/2);        % 再次确保 Hermitian
    d = real(diag(D));
    d(d < 0) = 0;                     % 去掉小负值
    rho_psd = V*diag(d)*V';

    % trace 归一化，避免 0/0
    trv = real(trace(rho_psd));
    if ~isfinite(trv) || trv <= 0
        rho_psd = eye(p)/p;
    else
        rho_psd = rho_psd/trv;
    end

    rho_vec_fixed = rho_psd(:);
end

function val = safeIntegral(fhandle, a, b)
% 对 vector-valued 被积函数做稳健积分；失败时返回 0 向量（保持长度与 f(a) 一致）

    try
        val = integral(fhandle, a, b, 'ArrayValued', true, ...
            'AbsTol', 1e-12, 'RelTol', 1e-9);
        if any(~isfinite(val))
            % 触发兜底
            x0  = fhandle(a);
            val = zeros(size(x0));
        end
    catch
        x0  = fhandle(a);
        val = zeros(size(x0));
    end
end

function out = safeExpmTimes(L, t1, M, R, t2, v)
% 计算  expm(L*t1) * ( M * expm(R*t2) * v )
% 若 M、R、t2 为空，则退化为 expm(L*t1) * v
    try
        if isempty(M)
            y = expm(L*t1) * v;
        else
            y = expm(R*t2) * v;
            y = M * y;
            y = expm(L*t1) * y;
        end
        if any(~isfinite(y)), error('inf_or_nan'); end
        out = y;
    catch
        out = zeros(size(v));
    end
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