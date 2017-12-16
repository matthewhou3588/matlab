function [ reliability ] = flc2( reliability_index, proximity_index )
%% flc in the online stage
%  input argument : reliability_index , is the output of flc1 , and its
%  value is belongs to [0,1], so no need normalize
%  input argument : proximity_index, need normalize before call flc2

[~, arg1] = fuzzyfunction(reliability_index);
[~, arg2] = fuzzyfunction(proximity_index);

% INFERENCE RULES OF THE model total reliability index from table IV 
rules = [0.025 0.05 0.1; 0.3 0.5 0.6; 0.7 0.9 1];

reliability = rules(arg1, arg2);

end

function [ fuzzy_value, liguistic ] = fuzzyfunction( x )
%% return value : liguistic is 1 stands for low, 2 for meidum and 3 for high
%

%% low
low = 0;
if x > 0 && x < 0.5
    low = (0.5-x)/0.5;
end

%% medium
medium = 0;
if x>0 && x<0.5
    medium = x/0.5;
end 
if x>=0.5 && x<1
   medium = (1-x)/0.5; 
end

%% high
high = 0;
if x>0.5 && x<1
   high = (x-0.5)/0.5; 
end

[fuzzy_value, liguistic] = max([low, medium, high]);

end
