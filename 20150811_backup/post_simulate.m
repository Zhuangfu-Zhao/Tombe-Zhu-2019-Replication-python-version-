dwa=new(1:N);
dwn=new(N+1:2*N);
dPa=new(2*N+1:3*N);
dPn=new(3*N+1:4*N);
dLa=new(4*N+1:5*N);
dLn=new(5*N+1:6*N);
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

% Migration Shares
mnn1=diag(mij2000); mnn0=ones(2*(N-1),1);
while sum((mnn0-mnn1).^2)>1e-20
    mnn0=mnn1;
    M_mL_term=(L1base.*dL(1:2*(N-1),1))./(mnn0.*L0);
    temp=reshape(M_mL_term,2,N-1)';
    M_mL_term_ag=temp(:,1);
    M_mL_term_na=temp(:,2);
    Cnnjj_ag=1+(eta_ag+(1-alpha)*beta_ag).*M_mL_term_ag./(alpha.*beta_ag);
    Cnnjj_na=1+(eta_na+(1-alpha)*beta_na).*M_mL_term_na./(alpha.*beta_na);
    Cnnjj=[Cnnjj_ag Cnnjj_na]'; Cnnjj=Cnnjj(:);
    cij=Cnijk+eye(2*(N-1),2*(N-1)).*repmat(Cnnjj,1,2*(N-1));
    Vjmat=repmat((dV(1:2*(N-1)).*Vi)',2*(N-1),1);
    mij=(cij.*Vjmat).^kappa; 
    mij=mij./(sum(mij,2)*ones(1,2*(N-1)));
    mnn1=mean([diag(mij) mnn0],2);
end

% New labour allocation
L1=(L0'*mij)';

% Welfare Change
aggregate_welfare_change=sum(omega_welfare.*dV(1:2*(N-1)).*(Cnnjj./Cnnjj2000).*(diag(mij)./diag(mij_base)).^(-1/kappa)); % For derivation, see "MigrationModel - WithRedistribution.pdf"

% Various migration summary statistics
migrant_shareofpop_2000=sum(L0.*(1-diag(mij)))./sum(L0);
migrant_shareofpop_2005=sum(L0.*(1-diag(mij2005)))./sum(L0);
migrant_shareofpop_new=sum(L0.*(1-diag(mij)))./sum(L0);
migrants_model=migrant_shareofpop_new/migrant_shareofpop_2000;
migrants_data=migrant_shareofpop_2005/migrant_shareofpop_2000;

mij_between_new=mij.*(1-kron(eye(N-1,N-1),ones(2,2))); % turn on for just between-province migration
mij_within_new=mij.*kron(eye(N-1,N-1),ones(2,2)).*(1-eye(2*N-2,2*N-2)); % turn on for just within-province migration
mij_ruralurban_new=mij.*repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1); % turn on for between ag and non-ag
mij_urbanurban_new=mij_between_new.*repmat([kron(ones(N-1,1),[0;0]) kron(ones(N-1,1),[0;1])],1,N-1); % turn on for within non-ag
mij_ruralrural_new=mij_between_new.*repmat([kron(ones(N-1,1),[1;0]) kron(ones(N-1,1),[0;0])],1,N-1); % turn on for within ag

mij_urbanrural_new=mij_ruralurban_new; 
mij_urbanrural_new(:,[2:2:60])=0; % Just U to R
mij_ruralurban_new(:,[1:2:60])=0; % Just R to U

migrant_shareofpop_within_new=sum(L0.*sum(mij_within_new,2))./sum(L0);
migrant_shareofpop_between_new=sum(L0.*sum(mij_between_new,2))./sum(L0);
migrant_shareofpop_ruralurban_new=sum(L0.*sum(mij_ruralurban_new,2))./sum(L0);
migrant_shareofpop_urbanrural_new=sum(L0.*sum(mij_urbanrural_new,2))./sum(L0);
migrant_shareofpop_urbanurban_new=sum(L0.*sum(mij_urbanurban_new,2))./sum(L0);
migrant_shareofpop_ruralrural_new=sum(L0.*sum(mij_ruralrural_new,2))./sum(L0);

mij_between_2000=mij2000.*(1-kron(eye(N-1,N-1),ones(2,2))); % turn on for just between-province migration
mij_within_2000=mij2000.*kron(eye(N-1,N-1),ones(2,2)).*(1-eye(2*N-2,2*N-2)); % turn on for just within-province migration
mij_ruralurban_2000=mij2000.*repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1); % turn on for between ag and non-ag

mij_urbanrural_2000=mij_ruralurban_2000; 
mij_urbanrural_2000(:,[2:2:60])=0; % Just U to R
mij_ruralurban_2000(:,[1:2:60])=0; % Just R to U
mij_ruralrural_2000=mij_between_2000.*repmat([kron(ones(N-1,1),[1;0]) kron(ones(N-1,1),[0;0])],1,N-1); % turn on for within ag
mij_urbanurban_2000=mij_between_2000.*repmat([kron(ones(N-1,1),[0;0]) kron(ones(N-1,1),[0;1])],1,N-1); % turn on for within non-ag

migrant_shareofpop_within_2000=sum(L0.*sum(mij_within_2000,2))./sum(L0);
migrant_shareofpop_between_2000=sum(L0.*sum(mij_between_2000,2))./sum(L0);
migrant_shareofpop_ruralurban_2000=sum(L0.*sum(mij_ruralurban_2000,2))./sum(L0);
migrant_shareofpop_urbanrural_2000=sum(L0.*sum(mij_urbanrural_2000,2))./sum(L0);
migrant_shareofpop_urbanurban_2000=sum(L0.*sum(mij_urbanurban_2000,2))./sum(L0);
migrant_shareofpop_ruralrural_2000=sum(L0.*sum(mij_ruralrural_2000,2))./sum(L0);

migrants_model_within=migrant_shareofpop_within_new/migrant_shareofpop_within_2000;
migrants_model_between=migrant_shareofpop_between_new/migrant_shareofpop_between_2000;

migrants_model_ruralrural=migrant_shareofpop_ruralrural_new-migrant_shareofpop_ruralrural_2000;
migrants_model_ruralurban=migrant_shareofpop_ruralurban_new-migrant_shareofpop_ruralurban_2000;
migrants_model_urbanrural=migrant_shareofpop_urbanrural_new-migrant_shareofpop_urbanrural_2000;
migrants_model_urbanurban=migrant_shareofpop_urbanurban_new-migrant_shareofpop_urbanurban_2000;

mij_between_2005=mij2005.*(1-kron(eye(N-1,N-1),ones(2,2))); % turn on for just between-province migration
mij_within_2005=mij2005.*kron(eye(N-1,N-1),ones(2,2)).*(1-eye(2*N-2,2*N-2)); % turn on for just within-province migration
migrant_shareofpop_within_2005=sum(L0.*sum(mij_within_2005,2))./sum(L0);
migrant_shareofpop_between_2005=sum(L0.*sum(mij_between_2005,2))./sum(L0);

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

% Regional trade for the migration plots
intflow_new=((Xi_new(1:N-1)*ones(1,N-1)).*pi_new(1:N-1,1:N-1));
intexp_new=sum(intflow_new)'-diag(intflow_new);
intimp_new=sum(intflow_new,2)-diag(intflow_new);
intflow=((Xi(1:N-1)*ones(1,N-1)).*pi(1:N-1,1:N-1));
intexp=sum(intflow)'-diag(intflow);
intimp=sum(intflow,2)-diag(intflow);
exports_new=Xi_new(N).*pi_new(N,1:N-1)';
imports_new=Xi_new(1:N-1).*pi_new(1:N-1,N);
exports=Xi(N).*pi(N,1:N-1)';
imports=Xi(1:N-1).*pi(1:N-1,N);

% Change in Real GDP
realGDP_change=sum(dV(1:2*(N-1)).*dL(1:2*(N-1)).*(mu.*L1base./sum(mu.*L1base)));

% Share of labour by sector
Lnew=reshape(L1',2,N-1)';
Lold=reshape(L1base',2,N-1)';
La_new=Lnew(:,1);
Ln_new=Lnew(:,2);
ag_labour_share_newold=sum([La_data(1:N-1) La_new(1:N-1)])./sum(L1);

% Change in ag share of labour in the data
L1data=(L0'*mij2005)';
Lnew_data=reshape(L1data',2,N-1)';
La_new_data=Lnew_data(:,1);
Ln_new_data=Lnew_data(:,2);
ag_labour_share_newold_data=sum([La_data(1:N-1) La_new_data])./sum(L1);

% The main table of results
results=[100.*[internal_import_share_new-internal_import_share international_import_share_new-international_import_share] ...
    migrants_model_within-1 migrants_model_between-1 ...
    realGDP_change-1 aggregate_welfare_change-1];
