clear all
initial_eqm

% Load the productivity
load dT dTa dTn % use dT_baseline for the homothetic preferences w/out IO links; use dT for the official data data; dT_homothetic_IOlinks for official data
dTa_est=dTa;
dTn_est=dTn;
inieqm=new;

%% Run the loops for only migration cost sub-components
experiments=perms([4:-1:1]); % 1=between-prov, within-nonag; 2=between-prov, within-ag; 3=between-prov, ag-nonag; 4=within-prov, ag-nonag
iterations=length(experiments);
basetypes=[kron([1;0],ones(4,1)) kron(ones(4,1),[0;1]) kron(ones(2,1),[0;0;1;1])];

% The various migration cost components
cijhat_full=Cnijk2005./Cnijk2000;
dummies_btwn_na=repmat([kron(ones(N-1,1),[0;0]) kron(ones(N-1,1),[0;1])],1,N-1).*(1-kron(eye(N-1,N-1),ones(2,2))); % between provinces within non-ag
dummies_btwn_ag=repmat([kron(ones(N-1,1),[1;0]) kron(ones(N-1,1),[0;0])],1,N-1).*(1-kron(eye(N-1,N-1),ones(2,2))); % between provinces within ag
dummies_btwn_btwn=repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1).*(1-kron(eye(N-1,N-1),ones(2,2))); % between provinces between sectors
dummies_wthn_btwn=repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1).*kron(eye(N-1,N-1),ones(2,2)); % within provinces between sectors

% Warning: This takes about 1 hour to run
loop=1;
warning off
for base=1:length(basetypes)

        % reset
        Cnijk=Cnijk2000;
        dTa=ones(N,1); dTn=ones(N,1);
        dni_ag=ones(N,N); dni_na=ones(N,N);
        
        % Different bases
        if basetypes(base,1)==1
            dTa=dTa_est; dTn=dTn_est;
        end
        if basetypes(base,2)==1
            dni_ag(1:N-1,1:N-1)=dni_asym_ag(1:N-1,1:N-1); dni_na(1:N-1,1:N-1)=dni_asym_na(1:N-1,1:N-1);
        end
        if basetypes(base,3)==1
            dni_ag(N,:)=dni_asym_ag(N,:); dni_na(N,:)=dni_asym_na(N,:); dni_ag(:,N)=dni_asym_ag(:,N); dni_na(:,N)=dni_asym_na(:,N);
        end
        new=fsolve(@(X)main_simulate(X),inieqm,optimset('Display','off'));
        baseguess=new;
        post_simulate
        realGDP_change_base=realGDP_change;
        
for iteration=1:iterations

    cijhat_target=ones(2*(N-1),2*(N-1));
    
        % Loop over migration cost changes, adding each component in
        % sequence. The marginal contributions are inferred later.
        for change=1:4
            tic
            
            if experiments(iteration,change)==1 % between-prov, within-nonag
                cijhat_target=cijhat_target.*(cijhat_full.*dummies_btwn_na+(1-dummies_btwn_na));
            elseif experiments(iteration,change)==2 % between-prov, within-ag
                cijhat_target=cijhat_target.*(cijhat_full.*dummies_btwn_ag+(1-dummies_btwn_ag));
            elseif experiments(iteration,change)==3 % between-prov, ag-nonag
                cijhat_target=cijhat_target.*(cijhat_full.*dummies_btwn_btwn+(1-dummies_btwn_btwn));
            elseif experiments(iteration,change)==4 % within-prov, ag-nonag
                cijhat_target=cijhat_target.*(cijhat_full.*dummies_wthn_btwn+(1-dummies_wthn_btwn));
            end
            
            cijhat=cijhat_target; % leave this on, always
            for i=1:2*(N-1), cijhat(i,i)=0; end % set diagonals to zero
            Cnijk=Cnijk2000.*cijhat;
            if change==1; guess=baseguess; else guess=new; end % this ensures a close guess is used, but resets each iteration
            new=fsolve(@(X)main_simulate(X),guess,optimset('Display','off'));
            post_simulate
            
%             for weight=[1:-0.25:0]; % replace with [1:-0.25:0] to gradually step to solution (helpful to ensure convergence; though takes longer)
%                 cijhat=(1-weight).*cijhat_target+weight.*ones(2*(N-1),2*(N-1)); % leave this on, always
%                 for i=1:2*(N-1), cijhat(i,i)=0; end
%                 Cnijk=Cnijk2000.*cijhat;
%                 if weight==1, guess=inieqm; else guess=new; end
%                 new=fsolve(@(X)simulate_effectivelabor(X),guess,optimset('Display','off'));
%             end
%             post_simulate

            gather_results(iteration,change,base)=realGDP_change./realGDP_change_base;
            
            disp(char(sprintf('Completed %d out of %d Iterations',loop,length(basetypes)*iterations*4),...
                sprintf('Estimated Time Remaining (Based on Last Loop): %d Minutes',round((length(basetypes)*iterations*4-loop)*toc/60))))
            loop=loop+1;
        end

end
end
save decomposition_results2

clear all
load decomposition_results2

for i=1:8
    marginal_effects(:,:,i)=[gather_results(:,1,i) gather_results(:,2:4,i)./gather_results(:,1:3,i)];
    shares(:,:,i)=log(marginal_effects(:,:,i))./log(prod(marginal_effects(:,:,i),2)*ones(1,4));
for change=1:4
    chunk=shares(:,:,i);
    chunk2=marginal_effects(:,:,i);
    MEs(:,i,change)=log(chunk2(experiments==change));
    value(change,1)=mean(chunk(experiments==change));
end
end

for change=1:4
    temp=MEs(:,:,change);
    growth(change,1)=mean(temp(:));
    growth_SD(change,1)=std(temp(:));
end

% Bottom panel of Table 9
bottom_panel=[growth growth./.571 growth_SD]


