clear all
initial_eqm

% Migration Weights in 2000
weight=(L0*ones(1,2*(N-1))).*mij2000.*(ones(2*(N-1),2*(N-1))-eye(2*(N-1),2*(N-1)));
% weight=weight.*repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1); % turn on for between ag and non-ag
% weight=weight.*(1-kron(eye(N-1,N-1),ones(2,2))); % turn on for just between-province migration
% weight=weight.*kron(eye(N-1,N-1),ones(2,2)); % turn on for just within-province migration
% weight=weight.*repmat([kron(ones(N-1,1),[1;0]) kron(ones(N-1,1),[0;0])],1,N-1); % turn on for within ag
% weight=weight.*repmat([kron(ones(N-1,1),[0;0]) kron(ones(N-1,1),[0;1])],1,N-1); % turn on for within non-ag

% Overall Migration Costs (Harmonic Means)
number=sum(sum(weight))./sum(L0);
weight=weight./sum(sum(weight));
avg_costs_2000=1/sum(sum(Cnijk2000.*weight)); 
avg_costs_2005=1/sum(sum(Cnijk2005.*weight));
change=avg_costs_2005./avg_costs_2000;

% Table 5, one row at a time (adjust the comments above)
[number avg_costs_2000 avg_costs_2005 change]
