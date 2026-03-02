load('V.mat','V')
load('E.mat','E')
load('proj_perms.mat','proj_perms')
load('validp.mat','validp')
load('N_exp.mat','N_exp')
%% 可观测量
d=3;
N=10;

s1=sparse([0,1/(2^(1/2)),0;1/(2^(1/2)),0,0;0,0,0]);
s2=sparse([0,0,0;0,0,1/(2^(1/2));0,1/(2^(1/2)),0]);
l=sparse([0,0,0;0,0,0;0,0,1]);
r=sparse([1,0,0;0,0,0;0,0,0]);
sx=s1+s2;
sz=sparse([1,0,0;0,0,0;0,0,-1]);
s2z=sparse([1,0,0;0,1,0;0,0,0]);
id3=speye(d);
idd=speye(d^N);
rdd1=speye(d^(N-1));
rdd2=speye(d^(N-2));

O1=kron(s2z,rdd1);
O2=kron(kron(sz,sz),rdd2);

%% 时间

T=linspace(1000,4000,30000);
T=T';

%% 演化以及时间平均，时间涨落

O_mean1=zeros(validp,1);
O_fluct1=zeros(validp,1);

for k=1:validp
    psi0=product_states(proj_perms(k,:),N);
    O_exp=evolve(E,V,psi0,O1,T);
    O_mean1(k)=mean(O_exp);
    O_fluct1(k)=std(O_exp);
end

%%

save('O_fluct_s2z','O_fluct1')
save('O_mean_s2z','O_mean1')

%% 演化以及时间平均，时间涨落

O_mean2=zeros(validp,1);
O_fluct2=zeros(validp,1);

for k=1:validp
    psi0=product_states(proj_perms(k,:),N);
    O_exp=evolve(E,V,psi0,O2,T);
    O_mean2(k)=mean(O_exp);
    O_fluct2(k)=std(O_exp);
end

%%

save('O_fluct_szsz.mat','O_fluct2')
save('O_mean_szsz.mat','O_mean2')


%%

load('O_fluct_s2z.mat','O_fluct1')
load('O_fluct_szsz.mat','O_fluct2')
load('O_mean_s2z.mat','O_mean1')
load('O_mean_szsz.mat','O_mean2')

%%
k1=100; % 1 2 100  3


% Tt=linspace(0,8000,80000);
% Tt=Tt';
Tt=linspace(0,800,8000);
Tt=Tt';

%%
psi_sample=product_states(proj_perms(k1,:),N);

% initial_O=psi_sample*psi_sample';
c=V'*psi_sample;
psi_sample_rn=sparse(V*c/(abs(c'*c)^(1/2)));
initial_O=psi_sample_rn*psi_sample_rn';

O_exp=evolve(E,V,psi_sample,O1,Tt);
ovlp=evolve(E,V,psi_sample,initial_O,Tt);


%%

c = V' * psi_sample;
c = c / sqrt(c'*c);
c2 = abs(c) .^2;
Epsi=c2'*E;
Npsi=c2'*N_exp;

bet_alp=Bet_Alp(Epsi,Npsi,E,N_exp);
beta=bet_alp(1);
alpha=bet_alp(2);
betac=BetaC(Epsi,E);

O_diag=real(diag(V'*O1*V));

upC=0;
upG=0;
downC=0;
downG=0;
for k=1:validp
    exponC=betac*E(k);
    exponG=beta*E(k)-alpha*N_exp(k);
    upC=upC+exp(-exponC)*O_diag(k);
    upG=upG+exp(-exponG)*O_diag(k);
    downC=downC+exp(-exponC);
    downG=downG+exp(-exponG);
end
O1_c=upC/downC;
O1_g=upG/downG;


%%
fill_x = [0, 8000, 8000, 0];  % 矩形的x坐标
fill_y = [O_mean1(k1)-O_fluct1(k1), O_mean1(k1)-O_fluct1(k1), ...
          O_mean1(k1)+O_fluct1(k1), O_mean1(k1)+O_fluct1(k1)];  % 矩形的y坐标
% 在小图中也绘制填充区域
fill_x_inset = [0, 140, 140, 0];
fill_y_inset = [O_mean1(k1)-O_fluct1(k1), O_mean1(k1)-O_fluct1(k1), ...
                O_mean1(k1)+O_fluct1(k1), O_mean1(k1)+O_fluct1(k1)];


% 使用fill函数填充区域
plot(Tt,O_exp, 'Color',[0 0.6 0.8], 'LineWidth', 0.5)
hold on
fill(fill_x, fill_y, [0.8, 1, 0],'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on
plot([0,8000],[O_mean1(k1),O_mean1(k1)],'Color','g','LineWidth',2)
hold on
plot([0,8000],[O_mean1(k1)-O_fluct1(k1),O_mean1(k1)-O_fluct1(k1)],'Color', [0.8, 1, 0])
hold on
plot([0,8000],[O_mean1(k1)+O_fluct1(k1),O_mean1(k1)+O_fluct1(k1)],'Color', [0.8, 1, 0])
hold on
plot([0,8000],[O1_c,O1_c],'LineStyle','-.','Color','r','LineWidth',1.5)
hold on
plot([0,8000],[O1_g,O1_g],'LineStyle','-.','Color','b','LineWidth',1.5)
box on
xlim([0,800])
xticks(0:4000:8000)
% yticks(0.5:0.2:0.9)
% yticks(0:0.5:1)
yticks(0.2:0.4:1)
% xticklabels({'0','4','8'})
set(gca,'linewidth',1.1)
set(gca,'FontName','Times New Roman','FontSize',18)
set(gcf, 'Position', [500, 300, 550, 400]);
hold off

% 创建内嵌小图
% 设置小图的位置和大小 [左边界, 下边界, 宽度, 高度] (归一化坐标，0-1之间)
axes('Position', [0.2, 0.18, 0.45, 0.32])  % 调整这些值来改变小图位置和大小

% 绘制小图内容 (x轴范围0-40的细节)
plot(Tt, O_exp, 'Color',[0 0.6 0.8],'LineWidth', 1)
hold on
plot(Tt, ovlp, 'Color',[0.5 0.5 0.5],'LineWidth', 1,'LineStyle','-')
hold on
fill(fill_x_inset, fill_y_inset, [0.8, 1, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% 绘制小图中的水平线
plot([0, 140], [O_mean1(k1), O_mean1(k1)], 'Color', 'g', 'LineWidth', 1.5)
plot([0, 140], [O_mean1(k1)-O_fluct1(k1), O_mean1(k1)-O_fluct1(k1)], 'Color', [0.8, 1, 0], 'LineWidth', 0.5)
plot([0, 140], [O_mean1(k1)+O_fluct1(k1), O_mean1(k1)+O_fluct1(k1)], 'Color', [0.8, 1, 0], 'LineWidth', 0.5)
plot([0, 140], [O1_c, O1_c], 'LineStyle', '-.', 'Color', 'r', 'LineWidth', 1.2)
plot([0, 140], [O1_g, O1_g], 'LineStyle', '-.', 'Color', 'b', 'LineWidth', 1.2)
% 设置小图的显示范围
xlim([0, 140])
xticks(0:60:120)
ylim([0.7,0.735])
ax = gca;
% 调整x轴刻度标签与轴线的距离，正值增加距离，负值减少距离
ax.XRuler.TickLabelGapOffset = -5;
box on
set(gca, 'YTickLabel', [])
set(gca,'linewidth',1)
set(gca,'FontName','Times New Roman','FontSize',15)
set(gcf, 'Position', [500, 300, 550, 400]);

%%

Nbar=(N_exp+N_exp')/2;
Ebar=(E+E')/2;

%%
omega=zeros(validp,validp);
for k=1:validp-1
    for i=k+1:validp
        omega(k,i)=density_sts(Ebar(k,i),Nbar(k,i));
        omega(i,k)=omega(k,i);
    end
end

omega_reverse=omega .^ (-1);
omega_reverse(1:validp+1:end)=0;

save('Omega&Omega_reverse','omega','omega_reverse')

% load('Omega&Omega_reverse.mat')

%%

est_fluct=zeros(validp,1);
for k=1:validp
    psi1=product_states(proj_perms(k,:),N);
    ci=V'*psi1;
    ci=ci / sqrt(ci'*ci);
    c2=abs(ci) .^ 2;
    est_fluct(k)=c2'*(omega_reverse)*c2;
end

est_fluct=est_fluct .^ (1/2);
%% 
scatter(est_fluct,O_fluct1,60,'filled', 'MarkerEdgeColor', [0 0 1],'LineWidth',0.5,'MarkerFaceColor', [0 0.93 1] ,'MarkerFaceAlpha',0.3)
% hold on
% scatter(est_fluct(1),O_fluct1(1),60,'Marker','o','LineWidth',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [0 0.93 1])
% hold on
% scatter(est_fluct(2),O_fluct1(2),60,'Marker','o','LineWidth',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [0 0.93 1])
% hold on
% scatter(est_fluct(3),O_fluct1(3),60,'Marker','o','LineWidth',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [0 0.93 1])
% hold on
% scatter(est_fluct(100),O_fluct1(100),60,'Marker','o','LineWidth',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [0 0.93 1])
% xlim([2,20])
% ylim([0,0.13])
% xticks(2:8:18)
% yticks(0:0.06:0.12)
% box on
% set(gca,'linewidth',1)
% set(gca,'FontName','Times New Roman','FontSize',15)
% set(gcf, 'Position', [500, 300, 550, 400]);


%%

load('HSnorm2.mat','HSnorm')
HSnorm2_nd=HSnorm .^ 2;
HSnorm2_nd(1:validp+1:end)=0;

est_fluct_rho=zeros(validp,1);
for k=1:validp
    psi1=product_states(proj_perms(k,:),N);
    ci=V'*psi1;
    ci=ci / sqrt(ci'*ci);
    c2=abs(ci) .^ 2;
    est_fluct_rho(k)=c2'*(HSnorm2_nd)*c2;
end

est_fluct_rho=est_fluct_rho .^ (1/2);

%% 
scatter(est_fluct_rho,O_fluct1,50,'filled', 'MarkerEdgeColor', [0 0 1],'LineWidth',0.5,'MarkerFaceColor', [0 0.93 1],'MarkerFaceAlpha',0.3)
hold on
scatter(est_fluct_rho(1),O_fluct1(1),60,'Marker','o','LineWidth',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [0 0.93 1],'MarkerFaceAlpha',0.3)
hold on
scatter(est_fluct_rho(2),O_fluct1(2),60,'Marker','o','LineWidth',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [0 0.93 1],'MarkerFaceAlpha',0.3)
hold on
scatter(est_fluct_rho(3),O_fluct1(3),60,'Marker','o','LineWidth',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [0 0.93 1],'MarkerFaceAlpha',0.3)
hold on
scatter(est_fluct_rho(100),O_fluct1(100),60,'Marker','o','LineWidth',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor', [0 0.93 1],'MarkerFaceAlpha',0.3)
xlim([0,0.2])
ylim([0,0.13])
xticks(0:0.1:0.2)
yticks(0:0.06:0.12)
box on
set(gca,'linewidth',1)
set(gca,'FontName','Times New Roman','FontSize',15)
set(gcf, 'Position', [500, 300, 550, 400]);
