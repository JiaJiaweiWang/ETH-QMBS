load('E.mat','E')
load('V.mat','V')
load('sn.mat','sn')
load('tn.mat','tn')
load('N_exp.mat','N_exp')

validp=length(E);

%%

d=3;
N=10; % 格点数

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

O1=kron(kron(sz,sz),rdd2);

O1_diag=real(diag(V'*O1*V));

scatter3(E,N_exp,O1_diag)

%%

aaaaa=V'*V;
aaaaa=sum(aaaaa(:));

%%

O1_c=zeros(validp,1);
O1_g=zeros(validp,1);

for i=1:validp
    betaC=BetaC(E(i),E);
    x=Bet_Alp(E(i),N_exp(i),E,N_exp);
    beta=x(1);
    alpha=x(2);
    upC=0;
    vsubG=0;
    downC=0;
    downG=0;
    for k=1:validp
        exponC=betaC*E(k);
        exponG=beta*E(k)-alpha*N_exp(k);
        upC=upC+exp(-exponC)*O1_diag(k);
        vsubG=vsubG+exp(-exponG)*O1_diag(k);
        downC=downC+exp(-exponC);
        downG=downG+exp(-exponG);
    end
    O1_c(i)=upC/downC;
    O1_g(i)=vsubG/downG;
end

O1_c(isnan(O1_c))=O1_diag(isnan(O1_c));

%%


[Ene,sortnum]=sort(E);
O1_c_line=O1_c(sortnum);

Ocdeviat=abs(O1_c-O1_diag);

% scatter(E,O1_diag)
% hold on
% plot(Ene, O1_c_line)

scatter(E,O1_diag,40,Ocdeviat,'filled','MarkerEdgeColor', 'k' ,'LineWidth',0.2)
c = colorbar('Location', 'eastoutside'); % 显示颜色条
clim([min(Ocdeviat) max(Ocdeviat)]);
c.Ticks = linspace(0, 0.3, 2); % 设置颜色条刻度
colormap("cool")
hold on
plot(Ene,O1_c_line,'LineWidth',2,'color','r','LineStyle','-.')
hold on
scatter(E(sn),O1_diag(sn),60,'MarkerEdgeColor', 'r','LineWidth',2)

xlim([min(E)-0.2,max(E)+0.2])
ylim([min(O1_diag)-0.05,max(O1_diag)+0.05])
xticks(-5:5:5)
yticks(-0.1:0.3:0.5)
box on
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 300, 550, 400]);


%%
ft=fit([E,N_exp],O1_diag,'poly23');

Odeviat=abs(O1_diag-ft(E,N_exp));

[Eg,Ng]=meshgrid(linspace(min(E), max(E), 20), linspace(min(N_exp), max(N_exp), 20));
Og=ft(Eg,Ng);

scatter3(E,N_exp,O1_diag,40,Odeviat,'filled', 'MarkerEdgeColor', 'k')
xlim([min(E), max(E)])
ylim([min(N_exp), max(N_exp)])
zlim([-0.17,0.53])
c = colorbar('Location', 'eastoutside'); % 显示颜色条
clim([0 max(Odeviat)]);
c.Ticks = linspace(0 , 0.06,3); % 设置颜色条刻度
colormap("cool")
% colorbar off
hold on
surf(Eg, Ng, Og, 'FaceAlpha', 0.1, 'EdgeColor', 'b','FaceColor','black','EdgeAlpha', 0.5)
hold on
scatter3(E(sn),N_exp(sn),O1_diag(sn),60,'MarkerEdgeColor', 'r','LineWidth',1)
%  'MeshStyle', 'Row'
xticks(-5:5:5)
yticks(0:1:2)
zticks(-0.1:0.3:0.5)
box on
set(gca,'linewidth',1.1)
set(gca,'FontName','Times New Roman','FontSize',18)
set(gcf, 'Position', [500, 300, 550, 400]);



%%

scatter(E(tn),O1_g(tn)-O1_diag(tn),35,'LineWidth', 1, 'MarkerEdgeColor', [0 0 1],'MarkerFaceColor', [0 0.93 1])
hold on
scatter(E(sn),O1_g(sn)-O1_diag(sn),50,'LineWidth', 2, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [0 0.93 1])
hold on
scatter(E(sn),O1_c(sn)-O1_diag(sn),50,'marker','x', 'LineWidth', 1.2, 'MarkerEdgeColor', [1 0 0])
hold on
scatter(E(tn),O1_c(tn)-O1_diag(tn),30,'marker','x','LineWidth', 0.7, 'MarkerEdgeColor', [1 0.5 0],'MarkerEdgeAlpha', 0.8)
ylim([min(O1_c-O1_diag)-0.05,max(O1_c-O1_diag)+0.05])
xlim([min(E)-0.1,max(E)+0.1])
xticks(-8:4:8)
% yticks(0:0.5:1)
% xlabel('E','FontName','Arial')
% ylabel('log_{10}|\langle E_i|\psi(0)\rangle|^2','FontName','Arial')
box on
% colorbar off
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 500, 550, 400]);



%% 子系统的划分以及约化密度矩阵

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

%%
p=1;

schattenc=zeros(validp,1);
schatteng=zeros(validp,1);

for i=1:validp
    if E(i)==min(E)
        betaC=BetaC(E(i)+0.00001,E);
    elseif E(i)==max(E)
        betaC=BetaC(E(i)-0.00001,E);
    else
        betaC=BetaC(E(i),E);
    end
    x=Bet_Alp(E(i),N_exp(i),E,N_exp);
    beta=x(1);
    alpha=x(2);
    vsubG=zeros(d^A,d^A);
    vsubC=zeros(d^A,d^A);
    downC=0;
    downG=0;
    for k=1:validp
        exponC=betaC*E(k);
        exponG=beta*E(k)-alpha*N_exp(k);
        vsubC=vsubC+exp(-exponC)*vsub(:,:,k);
        vsubG=vsubG+exp(-exponG)*vsub(:,:,k);
        downC=downC+exp(-exponC);
        downG=downG+exp(-exponG);
    end
    vsubC=vsubC/downC;
    vsubG=vsubG/downG;
    schatteng(i)=sdistance(vsubG,vsub(:,:,i),p);
    schattenc(i)=sdistance(vsubC,vsub(:,:,i),p);
end

%%


scatter(E(tn),schatteng(tn),35,'LineWidth', 1, 'MarkerEdgeColor', [0 0 1],'MarkerFaceColor', [0 0.93 1])
hold on
scatter(E(sn),schatteng(sn),50,'LineWidth', 2, 'MarkerEdgeColor', [1 0 0], 'MarkerFaceColor', [0 0.93 1])
hold on
scatter(E(sn),schattenc(sn),50,'marker','x', 'LineWidth', 1.2, 'MarkerEdgeColor', [1 0 0])
hold on
scatter(E(tn),schattenc(tn),30,'marker','x','LineWidth', 0.7, 'MarkerEdgeColor', [1 0.5 0],'MarkerEdgeAlpha', 0.8)
ylim([0,0.35])
xlim([min(E)-0.1,max(E)+0.1])
xticks(-8:4:8)
yticks(0:0.1:3)
% xlabel('E','FontName','Arial')
% ylabel('log_{10}|\langle E_i|\psi(0)\rangle|^2','FontName','Arial')
box on
% colorbar off
set(gca,'linewidth',1.2)
set(gca,'FontName','Times New Roman','FontSize',20)
set(gcf, 'Position', [500, 500, 550, 400]);