clear all
initial_eqm

% Simulate all changes, and efficiency/productivity changes
dni_ag=dni_asym_ag;
dni_na=dni_asym_na;
cijhat=Cnijk2005./Cnijk2000;
for i=1:2*(N-1), cijhat(i,i)=0; end % set diagonals to one
Cnijk=Cnijk2000.*cijhat;
load dT dTa dTn
new=fsolve(@(X)main_simulate(X),new,optimset('Display','off'));
post_simulate

% Reset the new initial equilibrium to 2005, following all trade and migration cost changes and productivity changes
baseGDP=realGDP_change;
baseW=aggregate_welfare_change;
baseCnnjj=Cnnjj;
basemij=mij;
basedV=dV;
base_btwn_migshare=migrant_shareofpop_between_new;
guess=new;

%% Lower internal trade costs to Canada's level
dni_ag(1:N-1,1:N-1)=0.502.*dni_asym_ag(1:N-1,1:N-1); for i=1:N, dni_ag(i,i)=1; end
dni_na(1:N-1,1:N-1)=0.933.*dni_asym_na(1:N-1,1:N-1); for i=1:N, dni_na(i,i)=1; end
new=fsolve(@(X)main_simulate(X),guess,optimset('Display','off'));
guess2=new;
post_simulate
results_table(1,:)=[realGDP_change./baseGDP-1 aggregate_welfare_change./baseW-1];

%% Raise migration rates to U.S.'s level
dni_ag(1:N-1,1:N-1)=dni_asym_ag(1:N-1,1:N-1); for i=1:N, dni_ag(i,i)=1; end
dni_na(1:N-1,1:N-1)=dni_asym_na(1:N-1,1:N-1); for i=1:N, dni_na(i,i)=1; end
cijhat_target(1:2*(N-1),1:2*(N-1))=4.5; for i=1:2*(N-1), cijhat_target(i,i)=1; end 
cijhat=cijhat_target.*(1-kron(eye(N-1,N-1),ones(2,2)))+kron(eye(N-1,N-1),ones(2,2)); % This makes the changes only between provinces
Cnijk=Cnijk2005.*cijhat;
new=fsolve(@(X)main_simulate(X),guess,optimset('Display','off'));
post_simulate
results_table(2,:)=[realGDP_change./baseGDP-1 aggregate_welfare_change./baseW-1];

migrant_shareofpop_between_new % the new migrant stock share, should be one-third

%% Both changes together
dni_ag(1:N-1,1:N-1)=0.502.*dni_asym_ag(1:N-1,1:N-1); for i=1:N, dni_ag(i,i)=1; end
dni_na(1:N-1,1:N-1)=0.933.*dni_asym_na(1:N-1,1:N-1); for i=1:N, dni_na(i,i)=1; end
new=fsolve(@(X)main_simulate(X),guess2,optimset('Display','off'));
post_simulate
results_table(3,:)=[realGDP_change./baseGDP-1 aggregate_welfare_change./baseW-1]
