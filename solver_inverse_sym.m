d=3;
N = 10; % 格点数
all_combinations = dec2base(0:3^N-1, 3, N) - '1'; % 转换为 -1, 0, +1

unique_perms_map = containers.Map('KeyType', 'char', 'ValueType', 'logical');

for i = 1:size(all_combinations, 1)
    current_perm = all_combinations(i,:);
    min_perm = current_perm;
    
    % 找到最小旋转表示
    for shift = 1:N-1
        rotated_perm = circshift(current_perm, [0, shift]);
        if is_less_than(rotated_perm, min_perm)
            min_perm = rotated_perm;
        end
    end
    
    % 如果最小表示未出现过，则加入
    key = mat2str(min_perm);
    if ~isKey(unique_perms_map, key)
        unique_perms_map(key) = true;
    end
end

% 提取唯一排列
unique_perms = cellfun(@(x) str2num(x), keys(unique_perms_map)', 'UniformOutput', false);
unique_perms = cell2mat(unique_perms);


proj_perms=zeros(size(unique_perms,1),N);
validp=0;
for i=1:size(unique_perms,1)
    symm_exit=1;
    for k=1:N
        rk=mod(k,N);
        if unique_perms(i,k)==1 && unique_perms(i,rk+1)==-1
            symm_exit=0;
        end
    end
    if symm_exit==1
        validp=validp+1;
        proj_perms(validp,:)=unique_perms(i,:);
    end
end
proj_perms=proj_perms(1:validp,:);

proj_perms_is = zeros(validp, N);
proj_perms_is(1,:) = proj_perms(1,:);
inverse_index = 1;

for i = 2:validp
    is_unique = true;
    
    % 计算当前排列的反演版本
    current_inverse = -flip(proj_perms(i,:));
    
    % 找到反演版本的最小表示
    min_inverse = current_inverse;
    for shift = 1:N-1
        rotated_inverse = circshift(current_inverse, [0, shift]);
        if is_less_than(rotated_inverse, min_inverse)
            min_inverse = rotated_inverse;
        end
    end
    
    % 检查最小反演版本是否已存在
    for k = 1:inverse_index
        if isequal(proj_perms_is(k,:), min_inverse)
            is_unique = false;
            break;
        end
    end
    
    if is_unique
        inverse_index = inverse_index + 1;
        proj_perms_is(inverse_index,:) = proj_perms(i,:);
    end
end

% 截断结果数组
proj_perms_is = proj_perms_is(1:inverse_index, :);

%%

nump=0;
numn=0;
eigenbasisp=zeros(d^N,inverse_index);
eigenbasisn=zeros(d^N,inverse_index);

for i=1:inverse_index
    current_inverse=-flip(proj_perms_is(i,:));
    min_inverse=current_inverse;
    for shift = 1:N-1
        rotated_inverse_perm = circshift(current_inverse, [0, shift]);
        if is_less_than(rotated_inverse_perm, min_inverse)
            min_inverse = rotated_inverse_perm;
        end
    end

    if isequal(min_inverse,proj_perms_is(i,:))
        nump=nump+1;
        permutat=proj_perms_is(i,:);
        for k=1:N
            eigenbasisp(:,nump)=eigenbasisp(:,nump)+product_states(permutat,N);
            permutat=circshift(permutat,1);
        end
        eigenbasisp(:,nump)=eigenbasisp(:,nump)/(abs(eigenbasisp(:,nump)'*eigenbasisp(:,nump))^(1/2));
    else
        nump=nump+1;
        numn=numn+1;

        permutat=proj_perms_is(i,:);
        for k=1:N
            eigenbasisp(:,nump)=eigenbasisp(:,nump)+product_states(permutat,N);
            permutat=circshift(permutat,1);
        end
        eigenbasisn(:,numn)=eigenbasisp(:,nump);
        permutat=-flip(permutat);
        for k=1:N
            eigenbasisp(:,nump)=eigenbasisp(:,nump)+product_states(permutat,N);
            eigenbasisn(:,numn)=eigenbasisn(:,numn)-product_states(permutat,N);
            permutat=circshift(permutat,1);
        end
        eigenbasisp(:,nump)=eigenbasisp(:,nump)/(abs(eigenbasisp(:,nump)'*eigenbasisp(:,nump))^(1/2));
        eigenbasisn(:,numn)=eigenbasisn(:,numn)/(abs(eigenbasisn(:,numn)'*eigenbasisn(:,numn))^(1/2));
    end
end

eigenbasisp=eigenbasisp(:,1:nump);
eigenbasisn=eigenbasisn(:,1:numn);


%% Hamiltonian

s1=sparse([0,1/(2^(1/2)),0;1/(2^(1/2)),0,0;0,0,0]);
s2=sparse([0,0,0;0,0,1/(2^(1/2));0,1/(2^(1/2)),0]);
l=sparse([0,0,0;0,0,0;0,0,1]);
r=sparse([1,0,0;0,0,0;0,0,0]);
sx=s1+s2;
id3=speye(d);
idd=speye(d^N);
H=sparse(d^N,d^N);
for i=1:N
    if(i==1)
        p1=sx;
    else
        p1=id3;
    end
    for j=2:N
        if(j==i)
            p1=kron(p1,sx);
        else
            p1=kron(p1,id3);
        end
    end
    H=H+p1;
end
for i=1:N
    if(i==1)
        p2=s1;
    else
        if(i==N)
            p2=l;
        else
            p2=id3;
        end
    end
    for j=2:N
        if(j==i)
            p2=kron(p2,s1);
        else
            if(j==i+1)
                p2=kron(p2,l);
            else
                p2=kron(p2,id3);
            end
        end
    end
    H=H-p2;
end
for i=1:N
    if(i==1)
        p3=r;
    else
        if(i==N)
            p3=s2;
        else
            p3=id3;
        end
    end
    for j=2:N
        if(j==i)
            p3=kron(p3,r);
        else
            if(j==i+1)
                p3=kron(p3,s2);
            else
                p3=kron(p3,id3);
            end
        end
    end
    H=H-p3;
end


%%


quasi_particle=sparse(d^N,d^N);
for i=1:N
    if(i==1)
        p4=spin(1);
        p5=spin(0);
    else
        if(i==N)
            p4=spin(0);
            p5=spin(-1);
        else
            p4=id3;
            p5=id3;
        end
    end
    for j=2:N
        if(j==i)
            p4=kron(p4,spin(1));
            p5=kron(p5,spin(0));
        else
            if(j==i+1)
                p4=kron(p4,spin(0));
                p5=kron(p5,spin(-1));
            else
                p4=kron(p4,id3);
                p5=kron(p5,id3);
            end
        end
    end
    psiy=(p4+p5)/(2^(1/2));
    quasi_particle=quasi_particle+psiy*psiy';
end

%%

H_p=eigenbasisp'*H*eigenbasisp;
H_n=eigenbasisn'*H*eigenbasisn;

H_p = (H_p + H_p')/2;  % 强制对称化
H_n = (H_n + H_n')/2;

% 使用eigs而不是eig以获得更高精度
options.tol = 1e-14;
options.maxit = 3000;
options.disp = 0;

[V_p, E_p] = eigs(H_p, size(H_p, 1), 'sa', options);  % 'sa' 表示最小特征值
[V_n, E_n] = eigs(H_n, size(H_n, 1), 'sa', options);

Vp=eigenbasisp * V_p;
Vn=eigenbasisn * V_n;

% 对特征向量进行排序
[Ep, idx_p] = sort(diag(E_p), 'ascend');
Vp = Vp(:, idx_p);
[En, idx_n] = sort(diag(E_n), 'ascend');
Vn = Vn(:, idx_n);



% [V_p,E_p]=eig(H_p);
% [V_n,E_n]=eig(H_n);
% 
% Vp=eigenbasisp * V_p;
% Vn=eigenbasisn * V_n;
% 
% Ep=real(diag(E_p));
% En=real(diag(E_n));

V=[Vp,Vn];
E=[Ep;En];

N_exp=real(diag(V'*quasi_particle*V));

scatter(E(1:nump),N_exp(1:nump))
hold on
scatter(E(nump+1:validp),N_exp(nump+1:validp))

%%

sn=convhull(E,N_exp);
snn=size(sn,1);
sn=sn(1:snn-1);

allIndices = 1:validp;
% 获取非凸包点的索引
tn=setdiff(allIndices, sn);

scatter(E(sn),N_exp(sn))
hold on
scatter(E(tn),N_exp(tn))

%%

save('E.mat','E')
save('V.mat','V')
save('N_exp','N_exp')
save('positiveSym','nump')
save('sn','sn')
save('tn','tn')
save('validp','validp')
