function [ output_args ] = pso_flc1_costfunction( input_args, z, k )
%pso_flc1_costfunction pso cost funtion for flc1
%   Detailed explanation goes here

aL = 0;
bM = 50;
cH = 100;
%cL bL aM cM aH bH
cL = input_args(1,1);
bL = input_args(1,2);
aM = input_args(1,3);
cM = input_args(1,4);
aH = input_args(1,5);
bH = input_args(1,6);

for i=1:length(z)
    % z low
    z_low = membershipFunction(aL, bL, cL, z(1,i));    
    % z medium
    z_medium = membershipFunction(aM, bM, cM, z(1,i));
    % z medium
    z_high = membershipFunction(aH, bH, cH, z(1,i));
    
    z_fuzzy = max([z_low z_medium z_high]);
    
    % k low
    k_low = membershipFunction(aL, bL, cL, z(1,i));    
    % z medium
    z_medium = membershipFunction(aM, bM, cM, z(1,i));
    % z medium
    z_high = membershipFunction(aH, bH, cH, z(1,i));
    
    z_fuzzy = max([z_low z_medium z_high]);
    
end
end

function fuzzy = membershipFunction(a, b, c, x)
    fuzzy = 0;    
    if x <= bL && x > aL
        fuzzy = (x - aL) / (bL - aL);
    elseif x > bL && x < cL
        fuzzy = (cL - x) / (cL - bL);
    end
end


