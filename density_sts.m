function [den] = density_sts(x, y)
    % 添加持久变量以避免重复加载
    persistent validp Ee Nn h1 h2
    
    if isempty(validp)
        % 加载您的数据文件
        load('N_exp.mat', 'N_exp');
        load('E.mat', 'E');  % 假设您的E数据保存在E.mat中
        load('H1mean.mat','H1mean')
        load('H2mean.mat','H2mean')

        Ee=E;
        Nn=N_exp;
        validp = length(E);
        h1=H1mean;
        h2=H2mean;
        fprintf('数据加载完成，validp = %d\n', validp);
    end
    
    den = 0;
    for i = 1:validp
        distance = (((x - Ee(i)) .^ 2) / (h1^2)) + (((y - Nn(i)) .^ 2) / (h2^2));
        den = den + exp(-distance / 2);
    end
    den = den / (2 * pi * validp * h1 * h2);
end