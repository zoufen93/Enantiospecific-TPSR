close all;
clear;
clc;
figure
width_cm = 20;
height_cm = 15;
set(gcf,'unit','centimeters','position',[20 5 width_cm height_cm])
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [1.6/height_cm 1.45/width_cm],[1.2/height_cm 0.8/height_cm], [1.5/width_cm 0.25/width_cm]);
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

%Raman电场强度
Es = 1;

% e0的范围
e0_vals = linspace(0.001, 20, 600);

%动态分配 eigvals_plot 数组大小
thetacL = zeros(length(e0_vals), 1);
thetacR = zeros(length(e0_vals), 1);
Dc = zeros(length(e0_vals), 1);

% 开启并行
parpool('local');

parfor i = 1:length(e0_vals)
    e0 = e0_vals(i);
    HpL = Hrd1 - e0 * Him_p; %M=1
    HmL = Hrd1 - e0 * Him_m; %M=-1
    H0L = Hrd0 - e0 * Him_0;%M=0

    HpR = Hrd1 + e0 * Him_p; %M=1
    HmR = Hrd1 + e0 * Him_m; %M=-1
    H0R = Hrd0 + e0 * Him_0; %M=0

    HpL = (HpL + HpL') / 2; % 确保矩阵对称
    HmL = (HmL + HmL') / 2;
    H0L = (H0L + H0L') / 2;
    HpR = (HpR + HpR') / 2; 
    HmR = (HmR + HmR') / 2;
    H0R = (H0R + H0R') / 2;
    
    %computeOmegaP(Jc,i,j,d0,dp,dm,-Es,Hi,Hj)
    OmegaL_da = computeOmegaP(Jc,4,1,d0,dp,dm,-Es,HpL,H0L); %a-->alpha; b,d-->beta; c-->gamma
    OmegaR_da = computeOmegaP(Jc,4,1,d0,dp,dm,Es,HpR,H0R);

    OmegaL_dc = computeOmegaP(Jc,4,2,d0,dp,dm,-Es,HpL,H0L);
    OmegaR_dc = computeOmegaP(Jc,4,2,d0,dp,dm,Es,HpR,H0R);

    OmegaL_ba = computeOmegam(Jc,4,1,d0,dp,dm,-Es,HmL,H0L);
    OmegaR_ba = computeOmegam(Jc,4,1,d0,dp,dm,Es,HmR,H0R);

    OmegaL_bc = computeOmegam(Jc,4,2,d0,dp,dm,-Es,HmL,H0L);
    OmegaR_bc = computeOmegam(Jc,4,2,d0,dp,dm,Es,HmR,H0R);

    OmegaL_cb = conj(OmegaL_bc);
    OmegaL_ad = conj(OmegaL_da);
    OmegaR_cb = conj(OmegaR_bc);
    OmegaR_ad = conj(OmegaR_da);
    OmegaL = angle(-OmegaL_ad*OmegaL_dc*OmegaL_cb*OmegaL_ba)/2;
    OmegaR = angle(-OmegaR_ad*OmegaR_dc*OmegaR_cb*OmegaR_ba)/2;
    thetacL(i) = 2*OmegaL;
    thetacR(i) = 2*OmegaR;
    %Dc(i) = 1 - (abs(sin(OmegaL)*sin(OmegaR)+cos(OmegaL)*cos(OmegaR)))^2;  %点乘
    Dc(i) = (abs(-sin(OmegaL)*cos(OmegaR)+cos(OmegaL)*sin(OmegaR)))^2; %叉乘
end

delete(gcp('nocreate'));

% 绘制前 9 个排序的特征值作为 e0 的函数
subplot(2,1,1)
h1 = plot(e0_vals, thetacL(:), '-', 'LineWidth', 1);
hold on;
h2 = plot(e0_vals, thetacR(:), '--', 'LineWidth', 1);
hold off;
legd=legend('$\theta_{c}^L$','$\theta_{c}^R$','location','southeast');%图例
set(legd,'Box','off','interpreter', 'latex','FontName','Times New Roman','FontSize',10);%设置图例
legd.ItemTokenSize = [15,10];
title('$|a\rangle=|1,0\rangle$, $|b\rangle=|3,-1\rangle$, $|d\rangle=|3,1\rangle$, $|c\rangle=|4,0\rangle$', 'interpreter', 'latex', 'FontSize', 10);
%axis([0 80 -1.1 1.1]);%x轴y轴范围
xlabel('$E_{0}$ (kV/cm)','interpreter', 'latex','FontSize',10);%设置x轴标签
ylabel('$\theta_{c}^{L,R}$','interpreter', 'latex','FontSize',10);%设置y轴标签
set(gca,'linewidth',0.8);%设置坐标轴线宽
%set(gca, 'xtick', [0 40 80], 'xticklabels', {'0' '40' '80'},'FontName','Times New Roman','fontsize',10)
%set(gca, 'ytick', [-1 -0.5 0 0.5 1], 'yticklabels', {'-1' '-0.5' '0' '0.5' '1'},'FontName','Times New Roman','fontsize',10)
%set(gca,'yminortick','on');%设置小刻度标签
grid on;

subplot(2,1,2)
h3 = plot(e0_vals, Dc(:), '-', 'LineWidth', 1);
hold off;
%title('$|a\rangle=|1,0\rangle$, $|b\rangle=|1,-1\rangle$, $|c\rangle=|3,0\rangle$, $|d\rangle=|1,1\rangle$', 'interpreter', 'latex', 'FontSize', 10);
%axis([0 80 -1.1 1.1]);%x轴y轴范围
xlabel('$E_{0}$ (kV/cm)','interpreter', 'latex','FontSize',10);%设置x轴标签
ylabel('$D$','interpreter', 'latex','FontSize',10);%设置y轴标签
set(gca,'linewidth',0.8);%设置坐标轴线宽
%set(gca, 'xtick', [0 40 80], 'xticklabels', {'0' '40' '80'},'FontName','Times New Roman','fontsize',10)
%set(gca, 'ytick', [-1 -0.5 0 0.5 1], 'yticklabels', {'-1' '-0.5' '0' '0.5' '1'},'FontName','Times New Roman','fontsize',10)
%set(gca,'yminortick','on');%设置小刻度标签
grid on;

toc; % 结束计时并显示经过的时间
 
% 提取数据
xData1 = get(h1, 'XData');
yData1 = get(h1, 'YData');
xData2 = get(h2, 'XData');
yData2 = get(h2, 'YData');
xData3 = get(h3, 'XData');
yData3 = get(h3, 'YData');

%保存到 CSV 文件
data = [xData1' yData1' yData2' yData3'];
writematrix(data, '/Users/zoufen/Desktop/plot_data142.csv');

% %%% Save figure
% set(gcf, 'PaperSize', [width_cm height_cm]);
% saveas(gcf, '/Users/zoufen/Desktop/Fig1.pdf')

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

% 定义 OmegaP函数<M=1|dp|M=0>
function s = computeOmegaP(Jc,i,j,d0,dp,dm,Es,Hp,H0)
s1 = 0;
na = (Jc+1)^2-1;
nb = (Jc+1)^2;
for ta = 1:na
      Ja = ceil(sqrt(ta+1) - 1); %向上取整
      ka = ta - (Ja^2+Ja);
      for tb = 1:nb
            Jb = ceil(sqrt(tb) - 1);
            kb = tb - (Jb^2+Jb+1);
            s1= s1 + sqrt(2*Ja+1)*sqrt(2*Jb+1)*Es*conj(A(i,ta,Hp))*A(j,tb,H0)*...
                    (-1)^(0)*Wigner3j([Ja,1,Jb],[1,-1,0])*...
                    (d0*(-1)^(-kb+0)*Wigner3j([Ja,1,Jb],[ka,0,-kb])...
                    +dp*(-1)^(-kb+1)*Wigner3j([Ja,1,Jb],[ka,-1,-kb])...
                    +dm*(-1)^(-kb-1)*Wigner3j([Ja,1,Jb],[ka,1,-kb]));
      end
end
s = s1;
end

% 定义 Omegam函数<M=-1|dm|M=0>
function s = computeOmegam(Jc,i,j,d0,dp,dm,Es,Hm,H0)
s1=0;
na = (Jc+1)^2-1;
nb = (Jc+1)^2;
for ta = 1:na
      Ja = ceil(sqrt(ta+1) - 1); %向上取整
      ka = ta - (Ja^2+Ja);
      for tb = 1:nb
            Jb = ceil(sqrt(tb) - 1);
            kb = tb - (Jb^2+Jb+1);
            s1= s1+ sqrt(2*Ja+1)*sqrt(2*Jb+1)*Es*conj(A(i,ta,Hm))*A(j,tb,H0)*...
                    (-1)^(0)*Wigner3j([Ja,1,Jb],[-1,1,0])*...
                    (d0*(-1)^(-kb+0)*Wigner3j([Ja,1,Jb],[ka,0,-kb])...
                    +dp*(-1)^(-kb+1)*Wigner3j([Ja,1,Jb],[ka,-1,-kb])...
                    +dm*(-1)^(-kb-1)*Wigner3j([Ja,1,Jb],[ka,1,-kb]));
      end
end
s = s1;
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


function w = Wigner3j( j123, m123 )

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