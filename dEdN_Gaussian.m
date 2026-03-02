function [fitresult, gof] = dEdN_Gaussian(dE_grid, dN_grid, dcount_grid)
%CREATEFIT1(DE_GRID,DN_GRID,DCOUNT_GRID)
%  创建一个拟合。
%
%  要进行 '无标题拟合 1' 拟合的数据:
%      X 输入: dE_grid
%      Y 输入: dN_grid
%      Z 输出: dcount_grid
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 29-Nov-2025 13:30:25 自动生成


%% 拟合: '无标题拟合 1'。
[xData, yData, zData] = prepareSurfaceData( dE_grid, dN_grid, dcount_grid );

% 设置 fittype 和选项。
ft = fittype( 'b*exp(-(c*x)^2-(d*y)^2)', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.743132468124916 0.392227019534168 0.706046088019609];

% 对数据进行模型拟合。
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );


% 绘制数据拟合图。
figure( 'Name', '无标题拟合 1' );
h = plot( fitresult, [xData, yData], zData );
legend( h, '无标题拟合 1', 'dcount_grid vs. dE_grid, dN_grid', 'Location', 'NorthEast', 'Interpreter', 'none' );
% 为坐标区加标签
xlabel( 'dE_grid', 'Interpreter', 'none' );
ylabel( 'dN_grid', 'Interpreter', 'none' );
zlabel( 'dcount_grid', 'Interpreter', 'none' );
grid on
view( 6.1, 90.0 );


