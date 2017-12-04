function [ result, z_output, k_output ] = flc1( Z, K )
%   flc1 fuzzy logic controlller 1 in offline stage
%   This reliability index is the output of a Fuzzy Logic Controller (FLC).
%   The input parameters are Zn and Kn, which are partitioned into three fuzzy sets 
%   (with triangular membership functions) corresponding to different degrees of similarity between each
%   parameter and its median value.
%   
%  Attention: the output should be Low, Medium or High. 
%  use 1 low, 2 medium, 3 high

%% 先使用最大最小化方法对Z和K分别做标准化
z_normalized = zeros(1, length(Z));
for i = 1:length(Z)
    z_normalized(1,i) = (Z(1,i) - min(Z)) / (max(Z)-min(Z));     
end

k_normalized = zeros(1, length(K));
for i = 1:length(K)
    k_normalized(1,i) = (K(1,i) - min(K)) / (max(K)-min(K));     
end

%% compute similarity between each parameter and its median value.
z_median = median(z_normalized);
k_median = median(k_normalized);

z_output = 1 + ones(1, length(Z));   % 初始化所有的ap的z值属于"medium"
for i = 1:length(Z)
    if z_normalized(1,i) <= z_median     % 如果第i个ap的z小于所有的z的median,并且该z与median的相似度大于百分之五十，则该ap的z值属于"low"
        tmp = (z_normalized(1,i) - z_median) / (z_median - min(z_normalized));
        if tmp < -0.5
            z_output(1,i) = 1;
        end
    end
    if z_normalized(1,i) > z_median      % 如果第i个ap的z大于所有的z的median,并且该z与median的相似度大于百分之五十，则该ap的z值属于"high"
        tmp = (z_normalized(1,i) - z_median) / (max(z_normalized) - z_median);
        if tmp > 0.5
            z_output(1,i) = 3;
        end
    end     
end

k_output = 1 + ones(1, length(K));   % 初始化所有的ap的k值属于"medium"
for i = 1:length(K)
   if k_normalized(1,i) <= k_median
        tmp = (k_normalized(1,i) - k_median) / (k_median - min(k_normalized));
        if tmp < -0.5
            k_output(1,i) = 1;
        end      
   end
   
   if k_normalized(1,i) > k_median
        tmp = (k_normalized(1,i) - k_median) / (max(k_normalized) - k_median);
        if tmp > 0.5
            k_output(1,i) = 3;
        end       
   end
end

%% compute reliability index for every access point
% rules_reliability_index from table II 
rules_reliability_index = [0.025 0.05 0.1; 0.3 0.5 0.6; 0.7 0.9 1];

result = zeros(1, length(K));
for i = 1:length(K)
   result(1, i) = rules_reliability_index(z_output(1,i), k_output(1,i));
end

end % function end

