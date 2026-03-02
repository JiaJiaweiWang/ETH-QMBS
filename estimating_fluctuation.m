load('E.mat','E')
load('validp.mat','validp')
load('N_exp.mat','N_exp')
load('V.mat','V')
load('sn.mat','sn')
load('tn.mat','tn')
load('positiveSym.mat','nump')
load('proj_perms.mat','proj_perms')
d=3;
N=10;

N_bar=(N_exp+N_exp')/2;
E_bar=(E+E')/2;
dE=E-E';
dN=N_exp-N_exp';


%%

psi0=product_states(proj_perms(2,:),N);

c=V'*psi0;
c=c/sqrt(abs(c'*c));
c2=abs(c) .^2;

%%
colorde=[0.0863    0.9804    0.8627];

scatter(E,log(c2),40,'filled','MarkerFaceColor',colorde,'MarkerEdgeColor','b','LineWidth',0.1)
hold on
scatter(E(sn),log(c2(sn)),50,'filled','MarkerFaceColor',colorde,'LineWidth',2,'MarkerEdgeColor','r')
xlim([-9,9])
ylim([-27 ,0])
xticks(-8:4:8)
yticks(-20:10:0)
box on
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 300, 550, 400]);



%%
omega=zeros(validp,validp);
for k=1:validp-1
    for i=k+1:validp
        omega(k,i)=density_sts(E_bar(k,i),N_bar(k,i));
        omega(i,k)=omega(k,i);
    end
end

omega_reverse=omega .^ (-1);
omega_reverse(1:validp+1:end)=0;

save('Omega&Omega_reverse','omega','omega_reverse')

load('Omega&Omega_reverse.mat')


%% 
rdd2=speye(d^(N-2));

psi_R1=kron(spin(0),spin(-1));
psi_R2=kron(spin(1),spin(0));
Rk=psi_R1*psi_R1'-psi_R2*psi_R2';
Rk=kron(Rk,rdd2);

%% R_k矩阵元的平方

R_matrix=V'*Rk*V;
R_matrix_2=R_matrix .^ 2;

%%
Insym=zeros(validp,1);
for i=1:nump
    Insym(i)=1;
end

%%

rp=zeros(validp,1);
rm=zeros(validp,1);
rp_deno=zeros(validp,1);
rm_deno=zeros(validp,1);
r_nume=zeros(validp,1);

for i=1:validp
    r_up=0;
    rm_down=0;
    rp_down=0;
    for k=1:validp
        if Insym(i)~=Insym(k)
        % if i~=k
            r_up=r_up+R_matrix_2(i,k);
            rp_down=rp_down+(R_matrix_2(i,k)/((dE(k,i)-1)^2));
            rm_down=rm_down+(R_matrix_2(i,k)/((dE(k,i)+1)^2));
        end
    end
    r_nume(i)=r_up;
    rp_deno(i)=rp_down;
    rm_deno(i)=rm_down;
    rp(i)=r_up/rp_down;
    rm(i)=r_up/rm_down;
end

%%
scatter(E,log(rp),40,rp_deno,'filled','MarkerEdgeColor','b','LineWidth',0.1)
colormap("cool")
colorbar
clim([0, max(rp_deno)])
set(colorbar, 'Ticks', 0:0.4:0.8);
hold on
scatter(E(sn),log(rp(sn)),50,'LineWidth',2,'MarkerEdgeColor','r')
xlim([-9,9])
ylim([-4.1 ,1.8])
xticks(-8:4:8)
yticks(-4:2:2)
box on
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 300, 550, 400]);

%%
uisetcolor
%%
colorde=[0.0863    0.9804    0.8627];
scatter(E,r_nume,40,'filled','MarkerFaceColor',colorde,'MarkerEdgeColor','b','LineWidth',0.1)
% colormap("cool")
% colorbar
% clim([0, max(rp_deno)])
% set(colorbar, 'Ticks', 0:0.4:0.8);
hold on
scatter(E(sn),r_nume(sn),50,'filled','MarkerFaceColor',colorde,'LineWidth',2,'MarkerEdgeColor','r')
xlim([-9,9])
ylim([0,0.061])
xticks(-8:4:8)
yticks(0:0.03:0.06)
box on
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 300, 550, 400]);
%%
colorde=[0.0863    0.9804    0.8627];
scatter(E,rp_deno,40,log(rp),'filled','MarkerEdgeColor','b','LineWidth',0.1)
colormap("cool")
colorbar
clim([min(log(rp)), max(log(rp))])
set(colorbar, 'Ticks', -4:4:0);
hold on
scatter(E(sn),rp_deno(sn),50,'filled','MarkerFaceColor',colorde,'LineWidth',2,'MarkerEdgeColor','r')
xlim([-9,9])
ylim([0 ,1])
xticks(-8:4:8)
yticks(0:0.5:1)
box on
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 300, 550, 400]);
%%
scatter(E(sn),r_nume(sn))
hold on
scatter(E(tn),r_nume(tn))

%%
InsymM=Insym+Insym';
%%
range_dE=dE<1.6 & dE>0.4 & InsymM==1 & R_matrix_2>10^(-30);

% scatter3(E_bar(range_dE),N_bar(range_dE),R_matrix_2(range_dE) .^ (-1));


numofsample=size(E_bar(range_dE),1);
X = linspace(-8.6, 8.6, 16);
Y = linspace(0.15, 2.8, 9);

Ebargrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
Nbargrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
Rm_2grid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
gridcount=0;

for i=1:size(X,2)-1
    for k=1:size(Y,2)-1
        rxy=E_bar>X(i) & E_bar<=X(i+1) & N_bar>Y(k) & N_bar<Y(k+1) & range_dE;
        if any(rxy(:))
            gridcount=gridcount+1;
            Ebargrid(gridcount)=mean(E_bar(rxy));
            Nbargrid(gridcount)=mean(N_bar(rxy));
            Rm_2grid(gridcount)=mean(R_matrix_2(rxy));
        end
    end
end

Ebargrid=Ebargrid(1:gridcount,1);
Nbargrid=Nbargrid(1:gridcount,1);
Rm_2grid=Rm_2grid(1:gridcount,1);
Rm_nega2=Rm_2grid .^ (-1);


x_data=Ebargrid(:);
y_data=Nbargrid(:);
den_data=Rm_nega2(:);

scatter3(x_data,y_data,den_data)


%%

dden=density_sts(x_data,y_data);

a = dden \ den_data;

[x_grid, y_grid] = meshgrid(linspace(min(x_data), max(x_data), 30), ...
                           linspace(min(y_data), max(y_data), 30));

z_grid=a* density_sts(x_grid,y_grid);
%%
scatter3(x_data, y_data, den_data, 50, den_data,'filled','MarkerEdgeColor', 'k','LineWidth', 1)
colormap("cool")
% clim([min(abs(predicted_den-den_data)),max(abs(predicted_den-den_data))])
hold on
surf(x_grid, y_grid, z_grid, 'FaceAlpha', 0.3, 'EdgeColor', 'b','FaceColor','interp','EdgeAlpha', 0.5,'FaceLighting', 'gouraud')
% 主光源 - 从斜上方照射
light('Position', [max(x_data) max(y_data) max(den_data(:))*2], 'Style', 'infinite');
% 补光 - 从另一角度照射
light('Position', [min(x_data) max(y_data) max(den_data(:))*1.5], 'Style', 'infinite');
% 底部补光 - 减少阴影
light('Position', [0 0 min(z_grid(:))-1], 'Style', 'local');
material shiny;
xlim([min(x_data), max(x_data)])
ylim([min(y_data)-0.02, max(y_data)])
zlim([0,max(den_data)+3000])
xticks(-5:5:5)
yticks(0.5:1:2.5)
zticks(0:10000:20000)
legend('$\mathcal{R}_{grid}^{-1}$', '$a\cdot\Omega(\mathcal{E},\mathcal{N})$', 'Location', 'northwest')
set(legend, 'Interpreter', 'latex')
box on
set(gca,'linewidth',1.1)
set(gca,'FontName','Times New Roman','FontSize',18)
set(gcf, 'Position', [500, 300, 550, 400]);
%%

f2=R_matrix_2 ./ omega_reverse;
f2(1:validp+1:end)=0;

f2_w=f2 ./ ((dE-1) .^2 );

%%
range_dE=abs(dE)>0.0000000001 & f2_w>10^(-10);
% scatter3(dE(range_dE),dN(range_dE),log(f2_w(range_dE)))

scatter(dE(range_dE),f2_w(range_dE))
%%


X = linspace(min(dE(range_dE)), max(dE(range_dE)), 500);
% Y = linspace(min(dN(range_dE)), max(dN(range_dE)), 20);

dEgrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
% dNgrid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
f2grid=zeros((size(X,2)-1)*(size(Y,2)-1),1);
gridcount=0;

for i=1:size(X,2)-1
    % for k=1:size(Y,2)-1
    
    % rxy=dE>X(i) & dE<=X(i+1) & dN>Y(k) & dN<Y(k+1) & range_dE;
    rxy=dE>X(i) & dE<=X(i+1) & range_dE;
    if any(rxy(:))
        gridcount=gridcount+1;
        dEgrid(gridcount)=mean(dE(rxy));
        % dNgrid(gridcount)=mean(dN(rxy));
        f2grid(gridcount)=mean(f2_w(rxy));
    end

    % end
end

dEgrid=dEgrid(1:gridcount,1);
% dNgrid=dNgrid(1:gridcount,1);
f2grid=f2grid(1:gridcount,1);

%%

scatter(dEgrid,f2grid,'filled')

%% fit f_coarse with E and N
options = optimoptions('lsqcurvefit', ...
    'OptimalityTolerance', 1e-10, ...  % 默认1e-6，改为更小的值
    'MaxFunctionEvaluations', 5000, ...  % 提高函数评估上限（默认1000）
    'FunctionTolerance', 1e-16, ...  % 降低函数容差（默认1e-6）
    'StepTolerance', 1e-18, ...      % 降低步长容差（默认1e-6）
    'MaxIterations', 5000, ...       % 增加最大迭代次数（默认400）
    'Display', 'iter');              % 显示迭代过程

% p0=[10^(-4),0.1,0.7,1*(10^(-5)),0.1,2.8,0.5,0.5];
% p0=[10^(-4),0.1,0.7,0.5];
p0=[10^(-4),0.1,0.7];
% X=zeros(size(dEgrid,1)*size(dEgrid,2),2);
% X(:,1)=dEgrid(:);
% X(:,2)=dNgrid(:);
fvalue=f2grid(:);
p_fit = lsqcurvefit(@f_est, p0, dEgrid(:), fvalue, [], [], options);
% p_fit = lsqcurvefit(@f_est, p0, dE(range_dE), f2_w(range_dE), [], [], options);


f_est_value=f_est(p_fit,X(:));


%%

scatter(dEgrid(:),f2grid(:),30,'filled','MarkerFaceColor',colorde,'MarkerEdgeColor','b')
hold on
plot(X(:),f_est_value,'LineStyle','-.','LineWidth',1.5)
xlim([-4,4])
ylim([0 ,0.0015])
xticks(-4:4:4)
% yticks(0:0.5:1)
box on
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 300, 550, 400]);
%%

est_f2=zeros(validp,validp);
for i=1:validp
    for k=1:validp
        xx=[dE(i,k),dN(i,k)];
        est_f2(i,k)=f_est(p_fit,xx);
    end
end
scatter3(dE(:),dN(:),est_f2(:))
%%
est_R2=omega_reverse .* est_f2;

%%
rp_d_est=zeros(validp,1);

for i=1:validp
    rp_down=0;
    for k=1:validp
        if Insym(i)~=Insym(k)
        % if i~=k
            rp_down=rp_down+est_R2(k,i);
        end
    end
    rp_d_est(i)=rp_down;
end

%%
scatter(E,log(rp_d_est),40,log(rp),'filled','MarkerEdgeColor','b','LineWidth',0.1)
colormap('cool')
hold on
scatter(E(sn),log(rp_d_est(sn)),50,'filled','MarkerFaceColor',colorde,'LineWidth',2,'MarkerEdgeColor','r')
xlim([-8.6,8.6])
ylim([-3,1])
xticks(-8:4:8)
yticks(-3:2:1)
box on
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 300, 550, 400]);

%%
% function [fvalue]=f_est(p,X)
% de=X(:,1);
% f_de1=p(1)*(exp(-(((de-p(3)) ./ p(2)) .^ 2))+ exp(-(((de+p(3)) ./ p(2)) .^ 2)));
% f_de2=p(4)*(exp(-(((de-p(6)) ./ p(5)) .^ 2))+ exp(-(((de+p(6)) ./ p(5)) .^ 2)));
% dn=X(:,2);
% f_dn1=exp(-((dn ./ p(7)) .^ 2));
% f_dn2=exp(-((dn ./ p(8)) .^ 2));
% fvalue=f_de1 .* f_dn1 + f_de2 .* f_dn2;
% end

%%
% function [fvalue]=f_est(p,X)
% de=X(:,1);
% f_de1=p(1)*exp(-(((de-p(3)) ./ p(2)) .^ 2));
% dn=X(:,2);
% f_dn1=exp(-((dn ./ p(4)) .^ 2));
% fvalue=f_de1 .* f_dn1;
% end

%%

function [fvalue]=f_est(p,X)
de=X(:,1);
f_de1=p(1)*exp(-(((de-p(3)) ./ p(2)) .^ 2));
fvalue=f_de1;
end