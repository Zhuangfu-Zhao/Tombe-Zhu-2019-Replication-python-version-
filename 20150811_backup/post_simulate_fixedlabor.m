dwa=new(1:N);
dwn=new(N+1:2*N);
dPa=new(2*N+1:3*N);
dPn=new(3*N+1:4*N);
dLa=ones(N,1);
dLn=ones(N,1);
dL=[dLa dLn]'; dL=dL(:);
dw=[dwa dwn]'; dw=dw(:);

% New sectoral revenues
R_ag_new=dwa.*dLa.*R_ag;
R_na_new=dwn.*dLn.*R_na;

% New household final demand spending
Ir_new=(beta_ag+eta_ag).*R_ag_new./alpha;
dIr=Ir_new./Ir;
Iu_new=(beta_na+eta_na).*R_na_new./alpha;
dIu=Iu_new./Iu;
Da_new_r=alpha*epsilon.*Ir_new;
Da_new_u=alpha*epsilon.*Iu_new;
Dn_new_r=(alpha*(1-epsilon)).*Ir_new;
Dn_new_u=(alpha*(1-epsilon)).*Iu_new;
Ds_new_r=(1-alpha).*Ir_new;
Ds_new_u=(1-alpha).*Iu_new;

% Change in land prices
dra(:,1)=dwa.*dLa;
drn(:,1)=dwn.*dLn;

% New trade shares for agriculture and nonagriculture
inside_ag=dwa'.^beta_ag.*(dPa'.^gamma(1,1).*dPn'.^gamma(1,2)).^(1-beta_ag-eta_ag).*dra'.^eta_ag./(dTa'.^(1./theta));
denominator_ag=sum(pi_ag.*(dni_ag.*kron(inside_ag,ones(N,1))).^(-theta),2);
pi_ag_new=(pi_ag.*(dni_ag.*kron(inside_ag,ones(N,1))).^(-theta))./kron(denominator_ag,ones(1,N));
inside_na=dwn'.^beta_na.*(dPa'.^gamma(2,1).*dPn'.^gamma(2,2)).^(1-beta_na-eta_na).*drn'.^eta_na./(dTn'.^(1./theta));
denominator_na=sum(pi_na.*(dni_na.*kron(inside_na,ones(N,1))).^(-theta),2);
pi_na_new=(pi_na.*(dni_na.*kron(inside_na,ones(N,1))).^(-theta))./kron(denominator_na,ones(1,N));

% Sector expenditures
X_ag_new=inv(pi_ag_new')*R_ag_new;
X_na_new=inv(pi_na_new')*R_na_new;

% Change in real income per effective labour, by sector
dP=dPa.^epsilon.*dPn.^(1-epsilon);
dVa=dwa./(dP.^alpha.*dra.^(1-alpha));
dVn=dwn./(dP.^alpha.*drn.^(1-alpha));
dV=[dVa dVn]'; dV=dV(:);
dPa_vector=[dPa dPa]'; dPa_vector=dPa_vector(:);
dv_vector=[dIr./dLa dIu./dLn]'; dv_vector=dv_vector(:);
Vinew=dV(1:2*(N-1)).*Vi;

% Change in Real GDP
realGDP_change=sum(dV(1:2*(N-1)).*(mu.*L1base./sum(mu.*L1base)));

% Trade shares
Xi=X_ag+X_na;
Xi_new=X_ag_new+X_na_new;
pi=(pi_ag.*(X_ag*ones(1,N))+pi_na.*(X_na*ones(1,N)))./(X_ag*ones(1,N)+X_na*ones(1,N));
pi_new=(pi_ag_new.*(X_ag_new*ones(1,N))+pi_na_new.*(X_na_new*ones(1,N)))./(X_ag_new*ones(1,N)+X_na_new*ones(1,N));
international_import_share=sum((Xi(1:N-1)./sum(Xi(1:N-1))).*pi(1:N-1,N));
internal_import_share=sum((Xi(1:N-1)./sum(Xi(1:N-1))).*(sum(pi(1:N-1,1:N-1),2)-diag(pi(1:N-1,1:N-1))));
home_share=sum((Xi(1:N-1)./sum(Xi(1:N-1))).*diag(pi(1:N-1,1:N-1)));
international_import_share_new=sum((Xi_new(1:N-1)./sum(Xi_new(1:N-1))).*pi_new(1:N-1,N));
internal_import_share_new=sum((Xi_new(1:N-1)./sum(Xi_new(1:N-1))).*(sum(pi_new(1:N-1,1:N-1),2)-diag(pi_new(1:N-1,1:N-1))));
home_share_new=sum((Xi_new(1:N-1)./sum(Xi_new(1:N-1))).*diag(pi_new(1:N-1,1:N-1)));

% Welfare Change -- ignoring migration, equations 11 and 12 are the same thing
aggregate_welfare_change=realGDP_change;

% Main Results Table
results_nomigration=[100.*[internal_import_share_new-internal_import_share international_import_share_new-international_import_share] ...
    realGDP_change-1 aggregate_welfare_change-1];


