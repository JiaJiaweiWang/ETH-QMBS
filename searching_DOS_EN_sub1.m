load('E','E')
load('N_exp.mat','N_exp')
load('validp.mat','validp')
load('V.mat','V')
load('positiveSym.mat','nump')
% load('HSnorm.mat','HSnorm')
load('sn.mat','sn')
load('tn.mat','tn')


d=3;
N = 10; % 格点数

dE=E-E';
dN=N_exp-N_exp';
E_bar=(E+E')/2;
N_bar=(N_exp+N_exp')/2;


%%

E_std=std(E(:));
N_std=std(N_exp(:));

Omega=zeros(validp,1);

for i=1:validp
    Omega(i)=density_sts(E(i),N_exp(i),E_std,N_std,1);
end

[Eg,Ng]=meshgrid(linspace(min(E), max(E), 30), linspace(min(N_exp), max(N_exp), 30));

omeg_mesh=zeros(size(Eg,1),size(Eg,2));
for i=1:size(Eg,1)
    for k=1:size(Eg,2)
        omeg_mesh(i,k)=density_sts(Eg(i,k),Ng(i,k),E_std,N_std,1);
    end
end
%%

scatter3(E(sn),N_exp(sn),Omega(sn),50, 'MarkerEdgeColor', 'r','LineWidth', 1.5)
hold on
scatter3(E,N_exp,Omega,30, Omega,'filled')
hold on
surf(Eg, Ng, omeg_mesh, 'FaceAlpha', 0.1, 'EdgeColor', 'b','FaceColor','black','EdgeAlpha', 0.5)
xlim([min(E), max(E)])
ylim([min(N_exp)-0.02, max(N_exp)])
zlim([0,max(Omega)+0.02])
xticks(-5:5:5)
yticks(-0.5:1:2.5)
zticks(0:0.05:0.1)
text(1,1,0.15, {sprintf('$h_1$=std($E_i$)=%.3f ', ...
    E_std); ...
    sprintf('$h_2$=std($N_i$)=%.3f ', ...
    N_std)}, ...
    'FontName', 'Times New Roman', 'FontSize', 14, ...
    'BackgroundColor', 'yellow', 'EdgeColor', 'black', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'Interpreter', 'latex');
box on
set(gca,'linewidth',1.1)
set(gca,'FontName','Times New Roman','FontSize',18)
set(gcf, 'Position', [500, 300, 550, 400]);
%% 子系统的划分

A=1;
Ab=N-A;

vresh=zeros(d^A,d^Ab,validp);
for i=1:validp
    vresh(:,:,i)=reshape(V(:,i),d^A,d^Ab);
end

vsub=zeros(d^A,d^A,validp);
for i=1:validp
    vsub(:,:,i)=vresh(:,:,i)*vresh(:,:,i)';
end

%% Hilbert-Schmidt norm

HSnorm=zeros(validp,validp);
for i=1:validp
    for k=1:i
        rhom=vresh(:,:,i)*vresh(:,:,k)';
        HSnorm(i,k)=(trace(rhom'*rhom))^(1/2);
        HSnorm(k,i)=HSnorm(i,k);
    end
end

save('HSnorm1','HSnorm')

%%


numE=50;
numN=50;
degrid=linspace(-16.4,16.4 ...
    ,numE);
dngrid=linspace(-2.6,2.6,numN);
dE_grid=zeros((numE-1)*(numN-1),1);
dN_grid=zeros((numE-1)*(numN-1),1);
dcount_grid=zeros((numE-1)*(numN-1),1);

validcountt=0;
for i=1:numE-1
    for k=1:numN-1
        condition=dE<degrid(1,i+1) & dE>degrid(1,i) & dN<dngrid(1,k+1) & dN>dngrid(1,k) & abs(dE)>0.0000001;
        validcountt=validcountt+1;
        if any(condition(:))
            dE_grid(validcountt)=mean(dE(condition));
            dN_grid(validcountt)=mean(dN(condition));
            dcount_grid(validcountt)=sum(condition(:));
        else
            dE_grid(validcountt)=(degrid(1,i+1)+degrid(1,i))/2;
            dN_grid(validcountt)=(dngrid(1,k+1)+dngrid(1,k))/2;
            dcount_grid(validcountt)=0;
        end
    end
end

dcount_grid=dcount_grid/max(dcount_grid(:));

[fitresult, gof] = dEdN_Gaussian(dE_grid, dN_grid, dcount_grid);

b=fitresult.b;
c=fitresult.c;
d=fitresult.d;
renorm_factor=c*d/pi;

densidEdN= (renorm_factor/b) * fitresult(dE,dN);


%%


[E_mesh, N_mesh] = meshgrid(linspace(-16.4, 16.4, 100), linspace(-2.6, 2.6, 100));
density_mesh = (renorm_factor/b) *fitresult(E_mesh, N_mesh);

p=0.1;

% 方法1：直接指定具体的等高线
contour_levels = [renorm_factor*(1-p),renorm_factor*(1-2*p),renorm_factor*(1-3*p),renorm_factor*(1-4*p),renorm_factor*(1-5*p),renorm_factor*(1-6*p)]; % 您想要的具体高度值
scatter(dE(:),dN(:),5,densidEdN(:),'filled')
hold on
contour(E_mesh, N_mesh, density_mesh, contour_levels,'LineWidth',1.5,'EdgeColor','k');
colorbar;
xlabel('dE');
ylabel('dN');
title('指定高度的等高线图');
% hold on
% scatter(dE(ranran),dN(ranran),5,'filled')

%% 初始值之赋予
E_std=std(E(:));
N_std=std(N_exp(:));

initial_guess = [E_std, N_std, 120000]; % [h1_initial, h2_initial, a_initial]
lb = [0.1, 0.1, 50];   % 参数下界
ub = [10, 10, 10000000];      % 参数上界

%% 设置优化选项
options = optimset('Display', 'iter', 'MaxIter', 1000, 'MaxFunEvals', 3000);

%% 循环参数设置
numofCycle=9;

p=0.1;
%%
H1=zeros(numofCycle,1);
H2=zeros(numofCycle,1);
H1_std=zeros(numofCycle,1);
H2_std=zeros(numofCycle,1);

R2=zeros(numofCycle,1);

for k=1:numofCycle
    densi_stt=renorm_factor;
    densi_end=renorm_factor*(1-k*p);
    rang=densidEdN > densi_end & densidEdN <densi_stt & abs(dE)>0.000001;

    R2D2=0;
    H11=0;
    H22=0;
    parastd=zeros(3,1);
    for xx=0:3
        X = linspace(-8.6, 8.6, 17-xx);
        Y = linspace(0.15, 2.8, 14-xx);
        Ebargrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
        Nbargrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
        normgrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
        gridcount=0;
        
        for ii=1:size(X,2)-1
            for kk=1:size(Y,2)-1
                rxy=E_bar>X(ii) & E_bar<=X(ii+1) & N_bar>Y(kk) & N_bar<Y(kk+1) & rang;
                if any(rxy(:))
                    gridcount=gridcount+1;
                    Ebargrid(gridcount)=mean(E_bar(rxy));
                    Nbargrid(gridcount)=mean(N_bar(rxy));
                    normgrid(gridcount)=sqrt(mean(HSnorm(rxy) .^ 2));
                end
            end
        end
        
        x_data=Ebargrid(1:gridcount,1);
        y_data=Nbargrid(1:gridcount,1);
        normgrid=normgrid(1:gridcount,1);
        den_data=normgrid .^ (-2);

        fprintf('开始优化...\n');

        % 定义目标函数（残差平方和）
        objective_function = @(params) calculate_residuals(params, x_data, y_data, den_data);
        
        % 执行非线性最小二乘优化
        [optimal_params, ~, residual, ~, ~] = lsqnonlin(objective_function, initial_guess, lb, ub, options);
        
        fprintf('优化完成！\n');

        R_squared = 1 - sum(residual.^2) / sum((den_data - mean(den_data)).^2);
        if R_squared>R2D2
            R2D2=R_squared;
            H11=optimal_params(1);
            H22=optimal_params(2);
            RMSE = sqrt(mean(residual.^2));
            jacobian_matrix = calculate_jacobian(optimal_params, x_data, y_data);
            parastd = sqrt(diag(inv(jacobian_matrix' * jacobian_matrix))) * RMSE;
        end
        if R2D2>0.93
            break
        end
    end
    
    H1(k)=H11;
    H2(k)=H22;
    R2(k)=R2D2;
    H1_std(k)=parastd(1);
    H2_std(k)=parastd(2);
end

%%

H1mean = mean(H1(:));
H2mean = mean(H2(:));


errorbar(1:length(H1(:)), H1(:), H1_std(:), 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'h1', 'LineWidth', 1);
hold on;
errorbar(1:length(H2(:)), H2(:), H2_std(:), 's', 'Color', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'h2', 'LineWidth', 1);
hold on
line([0,length(H1(:))+1], [H1mean, H1mean], 'LineWidth', 1.5, 'LineStyle','-.', 'Color', 'b')
hold on
line([0,length(H2(:))+1], [H2mean, H2mean], 'LineWidth', 1.5, 'LineStyle','-.', 'Color', 'r')
hold on
scatter(1:length(H1(:)),R2(:),40 ,'Marker','square','MarkerEdgeColor','k','LineWidth',1.5)
legend({'$h_1$', '$h_2$', ...
    ['$\overline{h}_1$ = ', num2str(H1mean, '%.2f')], ...
    ['$\overline{h}_2$ = ', num2str(H2mean, '%.2f')], ...
    '$R^2$'}, ...
    'Location', 'northeast','Interpreter','latex','NumColumns',3);
xlim([0,length(H1(:))+1])
ylim([0,1.5])
% text(6, 1.4, {sprintf('mean$(h_1) = %.3f$', ...
%     H1mean); ...
%     sprintf('mean$(h_2) = %.3f$', ...
%     H2mean)}, ...
%     'FontName', 'Times New Roman', 'FontSize', 14, ...
%     'BackgroundColor', 'yellow', 'EdgeColor', 'black', ...
%     'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
%     'Interpreter', 'latex');
set(gca,'linewidth',1.1)
set(gca,'FontName','Times New Roman','FontSize',18)
set(gcf, 'Position', [500, 300, 550, 400]);
% ylim([0,1])

% save('H1mean','H1mean')
% save('H2mean','H2mean')

%% 计算区域展示

[E_mesh, N_mesh] = meshgrid(linspace(-16.4, 16.4, 100), linspace(-2.6, 2.6, 100));
density_mesh = (renorm_factor/b) *fitresult(E_mesh, N_mesh);

p=0.1;

% 方法1：直接指定具体的等高线
contour_levels = zeros(1, 9);
for i = 1:9
    contour_levels(i) = renorm_factor * (1 - i * p);
end % 您想要的具体高度值
scatter(dE(:),dN(:),5,densidEdN(:),'filled')
hold on
contour(E_mesh, N_mesh, density_mesh, contour_levels,'LineWidth',0.8,'EdgeColor','k');
colorbar off;
box on
xlim([min(dE(:)),max(dE(:))])
ylim([min(dN(:)),max(dN(:))])
xticks(-16:8:16)
yticks(-2:2:2)
set(gca,'linewidth',1.1)
set(gca,'FontName','Times New Roman','FontSize',18)
set(gcf, 'Position', [500, 300, 550, 400]);% hold on
% scatter(dE(ranran),dN(ranran),5,'filled')



%%

rang=densidEdN > renorm_factor*(1-0.5) & abs(dE)>0.000001;
numofsample=size(E_bar(rang),1);
X = linspace(-8.6, 8.6, 15);
Y = linspace(0.15, 2.8, 12);

Ebargrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
Nbargrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
normgrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
gridcount=0;
minsample=10000;
gridnum=zeros((size(X,2)-1)*(size(Y,2)-1),1);

for i=1:size(X,2)-1
    for k=1:size(Y,2)-1
        rxy=E_bar>X(i) & E_bar<=X(i+1) & N_bar>Y(k) & N_bar<Y(k+1) & rang;
        if any(rxy(:))
            gridcount=gridcount+1;
            Ebargrid(gridcount)=mean(E_bar(rxy));
            Nbargrid(gridcount)=mean(N_bar(rxy));
            normgrid(gridcount)=sqrt(mean(HSnorm(rxy) .^ 2));
            gridnum(gridcount)=sum(rxy(:));
            if gridnum(gridcount)<minsample
                minsample=gridnum(gridcount);
            end
        end
    end
end

Ebargrid=Ebargrid(1:gridcount,1);
Nbargrid=Nbargrid(1:gridcount,1);
normgrid=normgrid(1:gridcount,1);
fnorm=normgrid .^ (-2);
gridnum=gridnum(1:gridcount,1);


x_data=Ebargrid(:);
y_data=Nbargrid(:);
den_data=fnorm(:);

figure;

subplot(1, 2, 1);
scatter3(x_data,y_data,den_data)
% zlim([0,10000])

subplot(1, 2, 2);
scatter3(x_data,y_data,log(gridnum))
title(sprintf('单个格子内最小样本数目 %d 总样本数目 %d', minsample, numofsample))
set(gcf, 'Position', [100, 100, 1600, 600]);

%% 执行优化
fprintf('开始优化...\n');

% 定义目标函数（残差平方和）
objective_function = @(params) calculate_residuals(params, x_data, y_data, den_data);

% 设置优化选项
options = optimset('Display', 'iter', 'MaxIter', 1000, 'MaxFunEvals', 3000);

% 执行非线性最小二乘优化
[optimal_params, resnorm, residual, exitflag, output] = lsqnonlin(objective_function, initial_guess, lb, ub, options);

fprintf('优化完成！\n');

%% 计算拟合评估指标
fprintf('计算拟合评估指标...\n');

% 计算预测值
predicted_den = zeros(size(den_data));
for i = 1:length(x_data)
    predicted_den(i) = density_sts(x_data(i), y_data(i), optimal_params(1), optimal_params(2), optimal_params(3));
end

% 计算各种评估指标
RMSE = sqrt(mean(residual.^2));
MAE = mean(abs(residual));
R_squared = 1 - sum(residual.^2) / sum((den_data - mean(den_data)).^2);
Adjusted_R_squared = 1 - (1 - R_squared) * (length(den_data) - 1) / (length(den_data) - length(optimal_params) - 1);
MSE = mean(residual.^2);

% 计算参数标准误差和置信区间
jacobian_matrix = calculate_jacobian(optimal_params, x_data, y_data);
parameter_std = sqrt(diag(inv(jacobian_matrix' * jacobian_matrix))) * RMSE;
confidence_intervals = [optimal_params' - 1.96*parameter_std, optimal_params' + 1.96*parameter_std];

%% 绘制结果
fprintf('绘制拟合结果...\n');

figure;

% 原始数据散点图
subplot(2, 3, 1);
scatter3(x_data, y_data, den_data, 20, den_data, 'filled');
colorbar;
xlabel('x');
ylabel('y');
zlabel('密度');
title('原始数据');
grid on;

% 拟合曲面
subplot(2, 3, 2);
% 创建网格用于绘制曲面
[x_grid, y_grid] = meshgrid(linspace(min(x_data), max(x_data), 30), ...
                           linspace(min(y_data), max(y_data), 30));
z_fit = zeros(size(x_grid));

for i = 1:numel(x_grid)
    z_fit(i) = density_sts(x_grid(i), y_grid(i), optimal_params(1), optimal_params(2), optimal_params(3));
end

surf(x_grid, y_grid, z_fit, 'EdgeColor', 'none');
hold on;
scatter3(x_data, y_data, den_data, 20, 'r', 'filled');
xlabel('x');
ylabel('y');
zlabel('密度');
title('拟合曲面 vs 原始数据');
legend('拟合曲面', '原始数据', 'Location', 'best');
grid on;

% 残差图
subplot(2, 3, 3);
plot(predicted_den, residual, 'o');
hold on;
plot([min(predicted_den), max(predicted_den)], [0, 0], 'r--', 'LineWidth', 2);
xlabel('预测值');
ylabel('残差');
title('残差图');
grid on;

% 预测 vs 实际值
subplot(2, 3, 4);
plot(den_data, predicted_den, 'o');
hold on;
plot([min(den_data), max(den_data)], [min(den_data), max(den_data)], 'r--', 'LineWidth', 2);
xlabel('实际值');
ylabel('预测值');
title(sprintf('预测 vs 实际 (R² = %.4f)', R_squared));
grid on;

% 第一幅新图：参数和拟合信息展示
subplot(2, 3, 5);
axis off; % 关闭坐标轴

% 创建信息文本
info_text = {
    sprintf('拟合参数结果:');
    sprintf('');
    sprintf('h1 = %.6f ± %.6f', optimal_params(1), parameter_std(1));
    sprintf('h2 = %.6f ± %.6f', optimal_params(2), parameter_std(2)); 
    sprintf('a  = %.6f ± %.6f', optimal_params(3), parameter_std(3));
    sprintf('');
    sprintf('拟合评估指标:');
    sprintf('R^2 = %.4f', R_squared);
    sprintf('调整R^2 = %.4f', Adjusted_R_squared);
    sprintf('RMSE = %.6f', RMSE);
    sprintf('MAE = %.6f', MAE);
    sprintf('MSE = %.6f', MSE);
    sprintf('');
    sprintf('数据点数: %d', length(x_data))
};

% 显示文本
text(0.1, 0.9, info_text, 'VerticalAlignment', 'top', 'FontSize', 10, ...
     'FontName', 'FixedWidth', 'Interpreter', 'none');
title('参数估计和拟合信息');

% 第二幅新图：log(density) vs x_data
subplot(2, 3, 6);

% 计算每个数据点的对数密度
log_density = zeros(validp,1);
for i = 1:validp
    density_val = density_sts(E(i), N_exp(i), optimal_params(1), optimal_params(2), 1);
    % 避免对非正数取对数
    log_density(i) = log(density_val);
end

% 绘制散点图
scatter(E, log_density, 30, 'filled', 'MarkerFaceColor', [0.2 0.6 0.8]);
xlabel('x');
ylabel('log(密度)');
xlim([-8,8]);
title('对数密度 vs x');
grid on;

% 添加趋势线（如果数据足够）
valid_idx = ~isnan(log_density);
if sum(valid_idx) > 1
    hold on;
    [E_sorted, sort_idx] = sort(E(valid_idx));
    log_density_sorted = log_density(valid_idx);
    log_density_sorted = log_density_sorted(sort_idx);

    % 使用移动平均平滑曲线
    window_size = min(5, floor(length(E_sorted)/4));
    if window_size > 1
        smooth_log_density = movmean(log_density_sorted, window_size);
        plot(E_sorted, smooth_log_density, 'r-', 'LineWidth', 2, ...
             'DisplayName', sprintf('%d点移动平均', window_size));
        legend('Location', 'best');
    end
end

% 调整图形大小以适应新增的子图
set(gcf, 'Position', [100, 100, 1600, 900]);




%%

scatter3(x_data, y_data, den_data, 50, den_data,'filled','MarkerEdgeColor', 'k','LineWidth', 1)
colormap("cool")
% clim([min(abs(predicted_den-den_data)),max(abs(predicted_den-den_data))])
hold on
surf(x_grid, y_grid, z_fit, 'FaceAlpha', 0.3, 'EdgeColor', 'b','FaceColor','interp','EdgeAlpha', 0.5,'FaceLighting', 'gouraud')
% 主光源 - 从斜上方照射
light('Position', [max(x_data) max(y_data) max(den_data(:))*2], 'Style', 'infinite');
% 补光 - 从另一角度照射
light('Position', [min(x_data) max(y_data) max(den_data(:))*1.5], 'Style', 'infinite');
% 底部补光 - 减少阴影
light('Position', [0 0 min(z_fit(:))-1], 'Style', 'local');
material shiny;
xlim([min(x_data), max(x_data)])
ylim([min(y_data)-0.02, max(y_data)])
zlim([0,max(den_data)+500])
xticks(-5:5:5)
yticks(0.5:1:2.5)
zticks(0:1000:2000)
legend('$\mathcal{T}_{grid}^{-1}$', '$a\cdot\Omega(\mathcal{E},\mathcal{N})$', 'Location', 'northwest')
set(legend, 'Interpreter', 'latex')
box on
set(gca,'linewidth',1.1)
set(gca,'FontName','Times New Roman','FontSize',18)
set(gcf, 'Position', [500, 300, 550, 400]);
% grid off

%%
scatter(den_data, predicted_den, 50,"yellow" ,'filled','MarkerEdgeColor', 'k','LineWidth', 0.5);
hold on;
plot([0, 2000], [0, 2000], 'r--', 'LineWidth', 2);
box on
xlim([0, 1750])
ylim([0, 1750])
xticks(0:1000:2000)
yticks(0:1000:2000)
set(gca,'linewidth',1.1)
set(gca,'FontName','Times New Roman','FontSize',18)
set(gcf, 'Position', [500, 300, 550, 400])
text(150, 1650, {sprintf('$R^2 = %.3f$', ...
    R_squared); ...
    sprintf('$h_1 = %.3f \\pm %.3f$', ...
    optimal_params(1), parameter_std(1)); ...
    sprintf('$h_2 = %.3f \\pm %.3f$', ...
    optimal_params(2), parameter_std(2))}, ...
    'FontName', 'Times New Roman', 'FontSize', 14, ...
    'BackgroundColor', 'yellow', 'EdgeColor', 'black', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'Interpreter', 'latex');
% xlabel('实际值')
% ylabel('预测值')
% title(sprintf('预测 vs 实际 (R² = %.4f)', R_squared))
% grid on;

%% 计算残差

function residuals = calculate_residuals(params, x_data, y_data, den_data)
% 计算残差
    h1 = params(1);
    h2 = params(2);
    a = params(3);
    
    predicted_den = zeros(size(den_data));
    for i = 1:length(x_data)
        predicted_den(i) = density_sts(x_data(i), y_data(i), h1, h2, a);
    end
    
    residuals = predicted_den - den_data;
end



%% 态密度计算

function [den] = density_sts(x, y, h1, h2, a)
    % 添加持久变量以避免重复加载
    persistent validp Ee Nn
    
    if isempty(validp)
        % 加载您的数据文件
        load('N_exp.mat', 'N_exp');
        load('E.mat', 'E');  % 假设您的E数据保存在E.mat中
        
        Ee=E;
        Nn=N_exp;
        validp = length(E);
        
        fprintf('数据加载完成，validp = %d\n', validp);
    end
    
    den = 0;
    for i = 1:validp
        distance = (((x - Ee(i)) .^ 2) / (h1^2)) + (((y - Nn(i)) .^ 2) / (h2^2));
        den = den + exp(-distance / 2);
    end
    den = a * den / (2 * pi * validp * h1 * h2);
end


%% Jacobian矩阵

function J = calculate_jacobian(params, x_data, y_data)
% 计算雅可比矩阵（数值微分）
    epsilon = 1e-6;
    n_params = length(params);
    n_data = length(x_data);
    J = zeros(n_data, n_params);
    
    for i = 1:n_params
        params_plus = params;
        params_plus(i) = params_plus(i) + epsilon;
        
        f_plus = zeros(n_data, 1);
        f_minus = zeros(n_data, 1);
        
        for j = 1:n_data
            f_plus(j) = density_sts(x_data(j), y_data(j), params_plus(1), params_plus(2), params_plus(3));
            f_minus(j) = density_sts(x_data(j), y_data(j), params(1), params(2), params(3));
        end
        
        J(:, i) = (f_plus - f_minus) / epsilon;
    end
end

