%% Lower all migration costs
close all
clear all
initial_eqm
dni_ag=ones(N,N);
dni_na=ones(N,N);
cijhat=Cnijk2005./Cnijk2000; % leave this on, always
for i=1:2*(N-1), cijhat(i,i)=0; end % this sets diagonals to one
Cnijk=Cnijk2000.*cijhat;
new=fsolve(@(X)main_simulate(X),new,optimset('Display','off'));
post_simulate
top_panel(1,:)=results

% To re-simulate with no land (eta=0), no housing (alpha=1), or high theta (we use 9999) 
% go to initial_eqm.m and change the necessary parameters. Then redo this section of code.

%% Agriculture to non-agriculture migration cost changes
close all
clear all
initial_eqm
cijhat=Cnijk2005./Cnijk2000; % leave this on, always
cijhat=cijhat.*repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1)+1-repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1);
for i=1:2*(N-1), cijhat(i,i)=0; end % set diagonals to one
Cnijk=Cnijk2000.*cijhat;
new=fsolve(@(X)main_simulate(X),new,optimset('Display','off'));
post_simulate
mid_panel(1,:)=results;

cijhat=Cnijk2005./Cnijk2000; % Reset the changes, to do just within-province
cijhat=cijhat.*repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1)+1-repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1);
cijhat=cijhat.*kron(eye(N-1,N-1),ones(2,2))+1-kron(eye(N-1,N-1),ones(2,2)); 
for i=1:2*(N-1), cijhat(i,i)=0; end % set diagonals to one
Cnijk=Cnijk2000.*cijhat;
new=fsolve(@(X)main_simulate(X),new,optimset('Display','off'));
post_simulate
mid_panel(2,:)=results;

cijhat=Cnijk2005./Cnijk2000; % Reset the changes, to do just between-province
cijhat=cijhat.*repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1)+1-repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1);
cijhat=cijhat.*(1-kron(eye(N-1,N-1),ones(2,2)))+kron(eye(N-1,N-1),ones(2,2));
for i=1:2*(N-1), cijhat(i,i)=0; end % set diagonals to one
Cnijk=Cnijk2000.*cijhat;
new=fsolve(@(X)main_simulate(X),new,optimset('Display','off'));
post_simulate
mid_panel(3,:)=results

%% Between provinces migration cost changes
close all
clear all
initial_eqm
cijhat=Cnijk2005./Cnijk2000; % leave this on, always
cijhat=cijhat.*(1-kron(eye(N-1,N-1),ones(2,2)))+kron(eye(N-1,N-1),ones(2,2));
for i=1:2*(N-1), cijhat(i,i)=0; end % set diagonals to one
Cnijk=Cnijk2000.*cijhat;
new=fsolve(@(X)main_simulate(X),new,optimset('Display','off'));
post_simulate
bottom_panel(1,:)=results;

cijhat=Cnijk2005./Cnijk2000; % Reset the changes, to do just within-agriculture
cijhat=cijhat.*(1-kron(eye(N-1,N-1),ones(2,2)))+kron(eye(N-1,N-1),ones(2,2));
cijhat=cijhat.*repmat([kron(ones(N-1,1),[1;0]) kron(ones(N-1,1),[0;0])],1,N-1)+1-repmat([kron(ones(N-1,1),[1;0]) kron(ones(N-1,1),[0;0])],1,N-1);
for i=1:2*(N-1), cijhat(i,i)=0; end % set diagonals to one
Cnijk=Cnijk2000.*cijhat;
new=fsolve(@(X)main_simulate(X),new,optimset('Display','off'));
post_simulate
bottom_panel(2,:)=results;

cijhat=Cnijk2005./Cnijk2000; % Reset the changes, to do just within-nonagriculture
cijhat=cijhat.*(1-kron(eye(N-1,N-1),ones(2,2)))+kron(eye(N-1,N-1),ones(2,2));
cijhat=cijhat.*repmat([kron(ones(N-1,1),[0;0]) kron(ones(N-1,1),[0;1])],1,N-1)+1-repmat([kron(ones(N-1,1),[0;0]) kron(ones(N-1,1),[0;1])],1,N-1);
for i=1:2*(N-1), cijhat(i,i)=0; end % set diagonals to one
Cnijk=Cnijk2000.*cijhat;
new=fsolve(@(X)main_simulate(X),new,optimset('Display','off'));
post_simulate
bottom_panel(3,:)=results
