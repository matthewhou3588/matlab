function [ sum_reliability] = pso_flc1_costfunction( input_args, z, k )
% function [ sum_reliability, node_reliability ] = pso_flc1_costfunction( input_args, z, k )
%pso_flc1_costfunction pso cost funtion for flc1
%   input_args is 1*12 matrix
%   input_args(1,1:6) cL bL aM cM aH bH for flc1's input z
%   input_args(1,7:12) cL bL aM cM aH bH for flc1's input k
%   z and k are the z and k for an anchor node

%cL bL aM cM aH bH
membership_args = ones(2,9);

% rules_reliability_index from table II 
rules_reliability_index = [0.025 0.05 0.1; 0.3 0.5 0.6; 0.7 0.9 1];

reliability = ones(1, length(z));

%% cL bL aM cM aH bH for z
membership_args(1,1) = 0;                % aL
membership_args(1,2) = input_args(1,2);  % bL
membership_args(1,3) = input_args(1,1);  % cL
    
membership_args(1,4) = input_args(1,3);  % aM
membership_args(1,5) = 50;               % bM
membership_args(1,6) = input_args(1,4);  % cM
    
membership_args(1,7) = input_args(1,5);  % aH
membership_args(1,8) = input_args(1,6);  % bH
membership_args(1,9) = 100;              % cH

%% %% cL bL aM cM aH bH for k
membership_args(2,1) = 0;                % aL
membership_args(2,2) = input_args(1,8);  % bL
membership_args(2,3) = input_args(1,7);  % cL
    
membership_args(2,4) = input_args(1,9);  % aM
membership_args(2,5) = 50;               % bM
membership_args(2,6) = input_args(1,10);  % cM
    
membership_args(2,7) = input_args(1,11);  % aH
membership_args(2,8) = input_args(1,12);  % bH
membership_args(2,9) = 100;              % cH

%% 
z_k = [z;k];   % just for 'for sentence'  convience

for i=1:length(z)
    
    fuzzy_linguistic = ones(2, 3);   
    
    for row = 1:2        % for z and k
       for index = 1:3   % for every low, medium, and high
          a = membership_args(row, index * 3 -2);
          b = membership_args(row, index * 3 -1);
          c = membership_args(row, index * 3); 
          
          fuzzy_linguistic(row, index) = membershipFunction(a, b, c, z_k(row, i));    
       end        
    end
    
    [~, z_index] = max(fuzzy_linguistic(1,:));
    [~, k_index] = max(fuzzy_linguistic(2,:));
    
    reliability(1,i) = rules_reliability_index(z_index, k_index);  % node i's reliability    
end

sum_reliability = sum(reliability);
node_reliability = reliability;

end

function fuzzy = membershipFunction(a, b, c, x)
    fuzzy = 0;    
    if x <= b && x > a
        fuzzy = (x - a) / (b - a);
    elseif x > b && x < c
        fuzzy = (c - x) / (c - b);
    end
end


