function best_variables = run_pso(zn, kn) 
tic
clc
%clear all
 close all
rng default

LB=[0 0 0 0 0 0 0 0 0 0 0 0];                           %lower bounds of variables
UB=[100 100 100 100 100 100 100 100 100 100 100 100];   %upper bounds of variables

% pso parameters values
m=12;        % number of variables
n=100;      % population size
inertia_weight = 0.73; 
wmax=0.9;   % inertia weight
wmin=0.4;   % inertia weight
c1=1.48;       % acceleration factor
c2=1.48;       % acceleration factor

% pso main program----------------------------------------------------start
maxite=20;    % set maximum number of iteration
maxrun=10;      % set maximum number of runs need to be
for run=1:maxrun
    %run
    % pso initialization----------------------------------------------start
    for i=1:n
        for j=1:m
            x0(i,j)=round(LB(j)+rand()*(UB(j)-LB(j)));  
        end
    end
    
    x=x0;       % initial population   % position of every particle in pso
    v=0.1*x0;   % initial velocity     % velocity of every particle in pso
    for i=1:n
        % f0(i,1)=ofun(x0(i,:));
        f0(i,1)=pso_flc1_costfunction(x0(i,:), zn, kn);
    end
    [fmax0,index0]=max(f0);
    pbest=x0;           % initial pbest
    gbest=x0(index0,:); % initial gbest
    % pso initialization------------------------------------------------end
    
    % pso algorithm---------------------------------------------------
    ite=1;
    tolerance=1;
    while ite<=maxite && tolerance>10^-12

        %w=wmax-(wmax-wmin)*ite/maxite; % update inertial weight
        w=inertia_weight;   % fix in paper

        % pso velocity updates
        for i=1:n
            for j=1:m
                v(i,j)=w*v(i,j)+c1*rand()*(pbest(i,j)-x(i,j))...
                +c2*rand()*(gbest(1,j)-x(i,j));      % eq. 7
            end
        end
        
        % pso position update
        for i=1:n
            for j=1:m
                x(i,j)=x(i,j)+v(i,j);
            end
        end
        
        % handling boundary violations
        for i=1:n
            for j=1:m
                if x(i,j)<LB(j)
                    x(i,j)=LB(j);
                elseif x(i,j)>UB(j)
                    x(i,j)=UB(j);
                end
            end
        end
         
        % evaluating fitness
        for i=1:n
            f(i,1)=pso_flc1_costfunction(x(i,:), zn, kn);
        end
        
        % updating pbest and fitness
        % in papaer: the aim is to maximize node reliability
        for i=1:n
            if f(i,1)>f0(i,1)
                pbest(i,:)=x(i,:);
                f0(i,1)=f(i,1);
            end
        end
        
        [fmax,index]=max(f0);   % finding out the best particle: maximize node reliability
        ffmax(ite,run)=fmax;    % storing best fitness
        ffite(run)=ite;         % storing iteration count

        % updating gbest and best fitness
        if fmax>fmax0
            gbest=pbest(index,:);
            fmax0=fmax;
        end
        
        % calculating tolerance
        if ite>100;
            tolerance=abs(ffmax(ite-100,run)-fmax0);
        end
        
        % displaying iterative results
%         if ite==1
%             disp(sprintf('Iteration Best particle Objective fun'));
%         end
%         disp(sprintf('%8g %8g %8.4f',ite,index,fmax0));
        ite=ite+1;
    end
    % pso algorithm-----------------------------------------------------end
    gbest;
    %fvalue=10*(gbest(1)-1)^2+20*(gbest(2)-2)^2+30*(gbest(3)-3)^2;
    fvalue=pso_flc1_costfunction(gbest, zn, kn);
    fff(run)=fvalue;
    rgbest(run,:)=gbest;
%     disp(sprintf('--------------------------------------'));
end
% pso main program------------------------------------------------------end
disp(sprintf('\n'));
disp(sprintf('*********************************************************'));
disp(sprintf('Final Results-----------------------------'));
[bestfun,bestrun]=max(fff)
best_variables=rgbest(bestrun,:)
disp(sprintf('*********************************************************'));
toc

%% PSO convergence characteristic
% plot(ffmax(1:ffite(bestrun),bestrun),'-k');
% xlabel('Iteration');
% ylabel('Fitness function value');
% title('PSO convergence characteristic')


%############################################################---------end

end
