clear all
initial_eqm

% Load the estimated productivity changes
load dT dTa dTn
dTa_est=dTa;
dTn_est=dTn;
inieqm=new;

%% Loop over all permutations (takes ~10 minutes)
experiments=perms([4 3 2 1]); % 1=productivity; 2=internal trade; 3=external trade; 4=migration costs
iterations=length(experiments);

loop=1;
warning off
for iteration=1:iterations %outer loop takes ~60 seconds
        % reset
        Cnijk=Cnijk2000;
        cijhat_target=ones(2*(N-1),2*(N-1));
        dTa=ones(N,1); dTn=ones(N,1);
        dni_ag=ones(N,N); dni_na=ones(N,N);
        
        % first change
        for change=1:4
            tic
            if experiments(iteration,change)==1
                dTa=dTa_est; dTn=dTn_est;
            elseif experiments(iteration,change)==2
                dni_ag(1:N-1,1:N-1)=dni_asym_ag(1:N-1,1:N-1); dni_na(1:N-1,1:N-1)=dni_asym_na(1:N-1,1:N-1);
            elseif experiments(iteration,change)==3
                dni_ag(N,:)=dni_asym_ag(N,:); dni_na(N,:)=dni_asym_na(N,:); dni_ag(:,N)=dni_asym_ag(:,N); dni_na(:,N)=dni_asym_na(:,N);
            elseif experiments(iteration,change)==4
                cijhat_target=Cnijk2005./Cnijk2000; 
            end
                cijhat=cijhat_target;
                for i=1:2*(N-1), cijhat(i,i)=0; end
                Cnijk=Cnijk2000.*cijhat;
                new=fsolve(@(X)main_simulate(X),ones(6*N,1),optimset('Display','off'));
            post_simulate
            
            gather_results(iteration,change)=realGDP_change;
           
            disp(char(sprintf('Completed %d out of %d Iterations',loop,iterations*4),...
                sprintf('Estimated Time Remaining (Based on Last Loop): %d Minutes',round((iterations*4-loop)*toc/60))))
            loop=loop+1;
        end
end
save decomposition_results

clear all
load decomposition_results
marginal_effects=[gather_results(:,1) gather_results(:,2:4)./gather_results(:,1:3)];
shares=log(marginal_effects)./log(prod(marginal_effects,2)*ones(1,4));
for change=1:4
    value(change,1)=mean(shares(experiments==change));
    growth(change,1)=mean(log(marginal_effects(experiments==change)));
    growth_SD(change,1)=std(log(marginal_effects(experiments==change)));
end

% Top Panel of Table 9
top_panel=[growth value growth_SD]

