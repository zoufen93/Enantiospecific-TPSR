close all;
clear;
clc;
figure
width_cm = 14;
height_cm = 5;
set(gcf,'unit','centimeters','position',[20 5 width_cm height_cm])
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [1.6/height_cm 1.45/width_cm],[1.0/height_cm 0.5/height_cm], [1.4/width_cm 0.25/width_cm]);
if ~make_it_tight,  clear subplot;  end       %%%%%%%%%% gap [h, w], h [down, up], w [left right]

tic;% 开始计时
%选取约化普朗克常数hbar=1
a = 2*pi*8572.05; % MHz
b = 2*pi*3640.10; % MHz
c = 2*pi*2790.96; % MHz
debye = 2*pi*503.3; % 1 Debye = h*0.5033 GHz/(kV/cm)
da = 1.201*debye; %\mu_a
db = 1.916*debye;
dc = 0.365*debye;
d0 = da;
dp = -(db+1.0j*dc)/sqrt(2);
dm = (db-1.0j*dc)/sqrt(2);
Jc = 5;
nc = (Jc + 1)^2 - 1;

% 定义函数 H1, H2, H3
H1 = @(J, k) 0.5 * (b + c) * (J * (J + 1) - k^2) + a * k^2;
H2 = @(J, k) 0.25 * (b - c) * sqrt((J * (J + 1) - k * (k + 1)) * (J * (J + 1) - (k + 1) * (k + 2)));
H3 = @(J, k) 0.25 * (b - c) * sqrt((J * (J + 1) - k * (k - 1)) * (J * (J + 1) - (k - 1) * (k - 2)));

% 定义函数 HR
HR = @(J) computeHR(J, H1, H2, H3);

% 初始化 Hr 矩阵
Hr1 = cell(Jc + 1, 1);
Hr0 = cell(Jc + 1, 1);

% 填充 Hr 矩阵
for J = 1:Jc
    Hr1{J + 1} = HR(J);
end

for J = 0:Jc
    Hr0{J + 1} = HR(J);
end

% 将Hr转换为块对角矩阵Hrd
Hrd1 = blkdiag(Hr1{:});
Hrd0 = blkdiag(Hr0{:});

% 定义 Wigner3j 函数
wigner_value = @(Ja, Ma, Jb, Mb, s) Wigner3j([Ja, 1, Jb], [Ma, -s, -Mb]);

% 定义 G[Ja, ka, Jb, kb] 函数
G_value = @(Ja, ka, Jb, kb) sqrt((2 * Ja + 1) * (2 * Jb + 1)) *(d0 * (-1)^(-kb) * wigner_value(Ja, ka, Jb, kb, 0)+ ...
                 dp * (-1)^(-kb+1) * wigner_value(Ja, ka, Jb, kb, 1) + dm * (-1)^(-kb-1) * wigner_value(Ja, ka, Jb, kb, -1));

% 定义 hi[Ja, ka, Jb, kb] 函数
hi_value_p = @(Ja, ka, Jb, kb) (-1)^(1+1) * wigner_value(Ja, 1, Jb, 1, 0) * G_value(Ja, ka, Jb, kb); %M=1
hi_value_m = @(Ja, ka, Jb, kb) (-1)^(-1+1) * wigner_value(Ja, -1, Jb, -1, 0) * G_value(Ja, ka, Jb, kb); %M=-1
hi_value_0 = @(Ja, ka, Jb, kb) (-1)^(0+1) * wigner_value(Ja, 0, Jb, 0, 0) * G_value(Ja, ka, Jb, kb); %M=0

%计算 Hi 数组
Him_p = zeros(nc, nc);
Him_m = zeros(nc, nc);
ii=1;
for Ja = 1:Jc
    for ka = -Ja:Ja
        jj=1;
        for Jb = 1:Jc
            for kb = -Jb:Jb
                Him_p(ii,jj) = hi_value_p(Ja, ka, Jb, kb);
                Him_m(ii,jj) = hi_value_m(Ja, ka, Jb, kb);
                jj=jj+1;
            end
        end
      ii=ii+1;
    end
end

Him_0 = zeros(nc+1, nc+1);
ii=1;
for Ja = 0:Jc
    for ka = -Ja:Ja
        jj=1;
        for Jb = 0:Jc
            for kb = -Jb:Jb
                Him_0(ii,jj) = hi_value_0(Ja, ka, Jb, kb);
                jj=jj+1;
            end
        end
      ii=ii+1;
    end
end

% e0的范围
e0 = 20;
HpL = Hrd1 - e0 * Him_p; %M=1
HmL = Hrd1 - e0 * Him_m; %M=-1
H0L = Hrd0 - e0 * Him_0; %M=0

HpR = Hrd1 + e0 * Him_p; %M=1
HmR = Hrd1 + e0 * Him_m; %M=-1
H0R = Hrd0 + e0 * Him_0; %M=0

HpL = (HpL + HpL') / 2;  %确保矩阵对称
HmL = (HmL + HmL') / 2;
H0L = (H0L + H0L') / 2;
HpR = (HpR + HpR') / 2; 
HmR = (HmR + HmR') / 2;
H0R = (H0R + H0R') / 2;

% P_{J,K,M;L}的布居
% P000L = abs(A(1, 1, H0L))^2;
% P1m10L = abs(A(1, 2, H0L))^2;
% P100L = abs(A(1, 3, H0L))^2;
% P1p10L = abs(A(1, 4, H0L))^2;
% P2m20L = abs(A(1, 5, H0L))^2;
% P2m10L = abs(A(1, 6, H0L))^2;
% P200L = abs(A(1, 7, H0L))^2;
% P2p10L = abs(A(1, 8, H0L))^2;
% P2p20L = abs(A(1, 9, H0L))^2;

% P_{J,M;L}的布居
% P00L = abs(A(1, 1, H0L))^2;
% P10L = abs(A(1, 2, H0L))^2 + abs(A(1, 3, H0L))^2 + abs(A(1, 4, H0L))^2;
% P20L = abs(A(1, 5, H0L))^2 + abs(A(1, 6, H0L))^2 + abs(A(1, 7, H0L))^2 + abs(A(1, 8, H0L))^2 + abs(A(1, 9, H0L))^2;
% P30L = abs(A(1, 10, H0L))^2 + abs(A(1, 11, H0L))^2 + abs(A(1, 12, H0L))^2 + abs(A(1, 13, H0L))^2 + abs(A(1, 14, H0L))^2 + abs(A(1, 15, H0L))^2 + abs(A(1, 16, H0L))^2;


% 定义每个 J 态的起止索引（手动或由其他函数生成）
J_vals_0 = 0:Jc;
J_indices_0 = {1, 2:4, 5:9, 10:16 17:25 26:36};  % 可扩展为自动生成
P10_L = zeros(1, Jc+1);  % 存放每个 J 的布居
P40_L = zeros(1, Jc+1);  % 存放每个 J 的布居
P10_R = zeros(1, Jc+1);  % 存放每个 J 的布居
P40_R = zeros(1, Jc+1);  % 存放每个 J 的布居
for j = 1:length(J_vals_0)
    idx = J_indices_0{j};  % 获取第 j 个 J 态对应的索引
    P10_L(j) = sum(abs(A(1, idx, H0L)).^2);  % 概率总和
    P40_L(j) = sum(abs(A(4, idx, H0L)).^2);  % 概率总和
    P10_R(j) = sum(abs(A(1, idx, H0R)).^2);  % 概率总和
    P40_R(j) = sum(abs(A(4, idx, H0R)).^2);  % 概率总和
end

J_vals_1 = 1:Jc;
J_indices_1 = {1:3, 4:8, 9:15 16:24 25:35};  % 可扩展为自动生成
P1p_L = zeros(1, Jc);  % 存放每个 J 的布居
P1m_L = zeros(1, Jc);  % 存放每个 J 的布居
P3p_L = zeros(1, Jc);  % 存放每个 J 的布居
P3m_L = zeros(1, Jc);  % 存放每个 J 的布居
for j = 1:length(J_vals_1)
    idx = J_indices_1{j};  % 获取第 j 个 J 态对应的索引
    P1p_L(j) = sum(abs(A(1, idx, HpL)).^2);  % 概率总和
    P1m_L(j) = sum(abs(A(1, idx, HmL)).^2);  % 概率总和
    P3p_L(j) = sum(abs(A(3, idx, HpL)).^2);  % 概率总和
    P3m_L(j) = sum(abs(A(3, idx, HmL)).^2);  % 概率总和
end

disp(['J=5, P10 = ', num2str(P10_L(6))]);
disp(['J=5, P40 = ', num2str(P40_L(6))]);
disp(['J=5, P31 = ', num2str(P3p_L(5))]);
disp(['J=5, P11 = ', num2str(P1p_L(5))]);

% % 绘图
% subplot(2,2,1)
% bar(J_vals_0, P10_L);
% xlabel('$J$', 'interpreter', 'latex', 'FontSize', 10);
% ylabel('$P_{|J,0;L\rangle}^{\alpha}$', 'interpreter', 'latex', 'FontSize', 12);
% title('$\alpha = 1$ and $M=0$', 'interpreter', 'latex', 'FontSize', 10);
% 
% subplot(2,2,2)
% bar(J_vals_0, P40_L);
% xlabel('$J$', 'interpreter', 'latex', 'FontSize', 10);
% ylabel('$P_{|J,0;L\rangle}^{\gamma}$', 'interpreter', 'latex', 'FontSize', 12);
% title('$\gamma = 4$ and $M=0$' , 'interpreter', 'latex', 'FontSize', 10);
% 
% subplot(2,2,3)
% bar(J_vals_1, P3p_L);
% xlabel('$J$', 'interpreter', 'latex', 'FontSize', 10);
% ylabel('$P_{|J,1;L\rangle}^{\beta}$', 'interpreter', 'latex', 'FontSize', 12);
% title('$\beta = 3$ and $M=1$', 'interpreter', 'latex', 'FontSize', 10);
% 
% subplot(2,2,4)
% bar(J_vals_1, P3m_L);
% xlabel('$J$', 'interpreter', 'latex', 'FontSize', 10);
% ylabel('$P_{|J,-1;L\rangle}^{\beta}$', 'interpreter', 'latex', 'FontSize', 12);
% title('$\beta = 3$ and $M=-1$' , 'interpreter', 'latex', 'FontSize', 10);

% 绘图 - 四组数据合并成一个散点图
subplot(1,2,1)
% α = 1, M=0
semilogy(J_vals_0, P10_L, 'o', 'MarkerSize', 5, 'DisplayName', '$|\alpha,0;L\rangle$', 'LineWidth', 1.0);
hold on;

%  β = 1, M=-1
semilogy(J_vals_1, P1m_L, '^', 'MarkerSize', 5, 'DisplayName', '$|\beta,\pm1;L\rangle$', 'LineWidth', 1.0);
hold on;

% β = 1, M=1
% semilogy(J_vals_1, P3p_L, 'd', 'MarkerSize', 5, 'DisplayName', '$\beta = 3, M=1$', 'LineWidth', 1.0);
% hold on;

% γ = 4, M=0
semilogy(J_vals_0, P40_L, 's', 'MarkerSize', 5, 'DisplayName', '$|\gamma,0;L\rangle$', 'LineWidth', 1.0);
hold off;
% 统一设置
xlabel('$J$', 'interpreter', 'latex', 'FontSize', 10);
ylabel('$P_{J}$', 'interpreter', 'latex', 'FontSize', 10);
legend('Box','off','Interpreter', 'latex', 'FontSize', 9, 'Location', 'best');
axis([-0.2 5.2 10^(-7) 3]);%x轴y轴范围
set(gca,'linewidth',0.8);%设置坐标轴线宽
set(gca, 'xtick', [0 1 2 3 4 5], 'xticklabels', {'0' '1' '2' '3' '4' '5'},'FontName','Times New Roman','fontsize',10)
set(gca, 'ytick', [10^(-6) 10^(-4) 10^(-2) 10^(0)], 'yticklabels', {'10^{-6}' '10^{-4}' '10^{-2}' '1'},'FontName','Times New Roman','fontsize',10)
grid off;


% 绘图 - 四组数据合并成一个散点图
subplot(1,2,2)
% α = 1, M=0
semilogy(J_vals_0, P10_L, 'o', 'MarkerSize', 5, 'DisplayName', '$|\alpha,0;L\rangle$', 'LineWidth', 1.0);
hold on;

%  β = 3, M=-1
semilogy(J_vals_1, P3m_L, '^', 'MarkerSize', 5, 'DisplayName', '$|\beta,\pm1;L\rangle$', 'LineWidth', 1.0);
hold on;

% β = 3, M=1
% semilogy(J_vals_1, P3p_L, 'd', 'MarkerSize', 5, 'DisplayName', '$\beta = 3, M=1$', 'LineWidth', 1.0);
% hold on;

% γ = 4, M=0
semilogy(J_vals_0, P40_L, 's', 'MarkerSize', 5, 'DisplayName', '$|\gamma,0;L\rangle$', 'LineWidth', 1.0);
hold off;
% 统一设置
xlabel('$J$', 'interpreter', 'latex', 'FontSize', 10);
ylabel('$P_{J}$', 'interpreter', 'latex', 'FontSize', 10);
legend('Box','off','Interpreter', 'latex', 'FontSize', 9, 'Location', 'best');
axis([-0.2 5.2 10^(-7) 3]);%x轴y轴范围
set(gca,'linewidth',0.8);%设置坐标轴线宽
set(gca, 'xtick', [0 1 2 3 4 5], 'xticklabels', {'0' '1' '2' '3' '4' '5'},'FontName','Times New Roman','fontsize',10)
set(gca, 'ytick', [10^(-6) 10^(-4) 10^(-2) 10^(0)], 'yticklabels', {'10^{-6}' '10^{-4}' '10^{-2}' '1'},'FontName','Times New Roman','fontsize',10)
grid off;

%%% Save figure
set(gcf, 'PaperSize', [width_cm height_cm]);
saveas(gcf, '/Users/zoufen/Desktop/FigS0.pdf')

% Define HR function
function HR_matrix = computeHR(J, H1, H2, H3)
    n = 2 * J + 1;
    % Initialize matrix with size n x n
    HR_matrix = zeros(n,n);
    
    % Fill the main diagonal
    for i = 1:n
        k = i - (J + 1);
        HR_matrix(i, i) = H1(J, k);
    end
    
    % Fill the superdiagonal (2 positions above the main diagonal)
    for i = 1:n-2
        k = i - (J + 1);
        HR_matrix(i, i + 2) = H2(J, k);
    end
    
    % Fill the subdiagonal (2 positions below the main diagonal)
    for i = 3:n
        k = i - (J + 1);
        HR_matrix(i, i - 2) = H3(J, k);
    end
end

function [eigenvalues, eigenvectors] = computeEig(Hi)
    [V, D] = eig(Hi); % 计算特征向量和特征值
    [eigenvalues, idx] = sort(diag(D), 'ascend'); % 提取特征值并从小到大排序
    eigenvectors = V(:, idx); % 提取排序后对应的特征向量
end

% 定义 A[i, t, Hi] 函数
function value = A(i, t, Hi)
    [~, eigenvectors] = computeEig(Hi);
    vector = eigenvectors(:, i); % 提取特定的特征向量
    normalized_vector = vector / norm(vector); % 归一化特征向量
    value = normalized_vector(t); % 提取特定索引的值
end

function w = Wigner3j(j123, m123)

j1 = j123(1); j2 = j123(2); j3 = j123(3);
m1 = m123(1); m2 = m123(2); m3 = m123(3);

% Input error checking
if any( j123 < 0 )
    error( 'The j must be non-negative' )
elseif any( rem( [j123, m123], 0.5 ) )
    error( 'All arguments must be integers or half-integers' )
elseif any( rem( (j123 - m123), 1 ) )
    error( 'j123 and m123 do not match' );
end

% Selection rules
if ( j3 > (j1 + j2) ) || ( j3 < abs(j1 - j2) ) ... % j3 out of interval
   || ( m1 + m2 + m3 ~= 0 ) ... % non-conserving angular momentum
   || any( abs( m123 ) > j123 ) % m is larger than j
    w = 0;
    return
end

% Simple common case
if ~any( m123 ) && rem( sum( j123 ), 2 ) % m1 = m2 = m3 = 0 & j1 + j2 + j3 is odd
    w = 0;
    return
end

% Evaluation
t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;

tmin = max( 0,  max( t1, t2 ) );
tmax = min( t3, min( t4, t5 ) );

t = tmin : tmax;
w = sum( (-1).^t .* exp( -ones(1,6)  ...
* gammaln( [t; t-t1; t-t2; t3-t; t4-t; t5-t] +1 ) + ...
                         gammaln( [j1+j2+j3+1, j1+j2-j3, ...
                         j1-j2+j3, -j1+j2+j3, j1+m1, ...
                          j1-m1, j2+m2, j2-m2, j3+m3, ...
                           j3-m3] +1 )* ...
                           [-1; ones(9,1)] * 0.5 ) ) ...
                            * (-1)^( j1-j2-m3 );

% Warnings
if isnan( w )
    warning( 'MATLAB:Wigner3j:NaN', 'Wigner3J is NaN!' )
elseif isinf( w )
    warning( 'MATLAB:Wigner3j:Inf', 'Wigner3J is Inf!' )
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