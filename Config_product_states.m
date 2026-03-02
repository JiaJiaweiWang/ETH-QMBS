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
    is_valid=1;
    for k=1:N
        rk=mod(k,N);
        if unique_perms(i,k)==1 && unique_perms(i,rk+1)==-1
            is_valid=0;
        end
    end
    if is_valid==1
        validp=validp+1;
        proj_perms(validp,:)=unique_perms(i,:);
    end
end
proj_perms=proj_perms(1:validp,:);

save('proj_perms','proj_perms')

% function less = is_less_than(a, b)
%     % 比较两个排列的字典序
%     less = false;
%     for k = 1:length(a)
%         if a(k) < b(k)
%             less = true;
%             return;
%         elseif a(k) > b(k)
%             return;
%         end
%     end
% end