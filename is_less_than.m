function less = is_less_than(a, b)
    % 比较两个排列的字典序
    less = false;
    for k = 1:length(a)
        if a(k) < b(k)
            less = true;
            return;
        elseif a(k) > b(k)
            return;
        end
    end
end