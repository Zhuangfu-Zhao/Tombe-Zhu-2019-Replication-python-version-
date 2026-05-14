close all
clear all
initial_eqm
inieqm=new;

%% Lower internal trade costs
dni_ag(1:N-1,1:N-1)=dni_asym_ag(1:N-1,1:N-1);
dni_na(1:N-1,1:N-1)=dni_asym_na(1:N-1,1:N-1);
new=fsolve(@(X)main_simulate(X),ones(6*N,1),optimset('Display','off'));
post_simulate
top_panel(1,:)=results;

% The No Migration Case
new=fsolve(@(X)simulate_fixedlabor(X),new(1:4*N),optimset('Display','off'));
post_simulate_fixedlabor
bottom_panel(1,:)=results_nomigration;

%% Lower external trade costs
dni_ag=ones(N,N);
dni_na=ones(N,N);
dni_ag(N,:)=dni_asym_ag(N,:);
dni_na(N,:)=dni_asym_na(N,:);
dni_ag(:,N)=dni_asym_ag(:,N);
dni_na(:,N)=dni_asym_na(:,N);
new=fsolve(@(X)main_simulate(X),ones(6*N,1),optimset('Display','off'));
post_simulate
top_panel(2,:)=results;

% Look with fixed labour
new=fsolve(@(X)simulate_fixedlabor(X),new(1:4*N),optimset('Display','off'));
post_simulate_fixedlabor
bottom_panel(2,:)=results_nomigration;

%% Lower all trade costs
dni_ag=dni_asym_ag;
dni_na=dni_asym_na;
new=fsolve(@(X)main_simulate(X),ones(6*N,1),optimset('Display','off'));
post_simulate
top_panel(3,:)=results;

% Look with fixed labour
new=fsolve(@(X)simulate_fixedlabor(X),new(1:4*N),optimset('Display','off'));
post_simulate_fixedlabor
bottom_panel(3,:)=results_nomigration;

%% Table 8 Output
% To simulate without intermediate inputs, include lines 70-74 in
% initial_eqm.m and then redo Table8.m from the start
top_panel
bottom_panel

