close all
clear all
initial_eqm
inieqm=new;

global va vn
va=Ir(1:N-1)./La;
vn=Iu(1:N-1)./Ln;

%% Original objects
chunk1r=(1-alpha).*Ir;
chunk2r=alpha.*(eta_ag./(eta_ag+beta_ag)).*Ir;
chunk1u=(1-alpha).*Iu;
chunk2u=alpha.*(eta_na./(eta_na+beta_na)).*Iu;
rS=[chunk1r+chunk2r chunk1u+chunk2u]'; rS=rS(:);
global v f w v_matrix
f=diag(rS(1:2*(N-1))./(mnn0.*L0));
w=[alpha.*beta_ag.*Ir(1:N-1)./(La.*(beta_ag+eta_ag)) alpha.*beta_na.*Iu(1:N-1)./(Ln.*(beta_na+eta_na))]'; w=w(:);
v=[Ir(1:N-1)./La Iu(1:N-1)./Ln]'; v=v(:);
vL=v.*L1;
v_matrix=repmat(w(1:2*(N-1))',2*(N-1),1)+f;

%% Let land income follow migrants
Cnijk=Cnijk2000;
dni_ag=ones(N,N);
dni_na=ones(N,N);
new=fsolve(@(X)simulate_landreform(X),ones(8*N,1),optimset('Display','off'));
dwa=new(1:N);
dwn=new(N+1:2*N);
dPa=new(2*N+1:3*N);
dPn=new(3*N+1:4*N);
dLa=new(4*N+1:5*N);
dLn=new(5*N+1:6*N);
dva=new(6*N+1:7*N); 
dvn=new(7*N+1:8*N);

dL=[dLa dLn]'; dL=dL(:);
dw=[dwa dwn]'; dw=dw(:);
dv=[dva dvn]'; dv=dv(:);

% New sectoral revenues
R_ag_new=dwa.*dLa.*R_ag;
R_na_new=dwn.*dLn.*R_na;

% Counterfactual land returns
chunk1r=(1-alpha).*Ir;
chunk2r=alpha.*(eta_ag./(eta_ag+beta_ag)).*Ir;
chunk1u=(1-alpha).*Iu;
chunk2u=alpha.*(eta_na./(eta_na+beta_na)).*Iu;
dra(:,1)=dva.*dLa.*chunk1r./(chunk1r+chunk2r)+dwa.*dLa.*chunk2r./(chunk1r+chunk2r);
drn(:,1)=dvn.*dLn.*chunk1u./(chunk1u+chunk2u)+dwn.*dLn.*chunk2u./(chunk1u+chunk2u);

% Land income from home, "f"
global v f w v_matrix
w_new=w.*dw(1:2*(N-1));
vL_new=[dva(1:N-1).*dLa(1:N-1).*Ir(1:N-1) dvn(1:N-1).*dLn(1:N-1).*Iu(1:N-1)]'; vL_new=vL_new(:);
wa_new=dwa(1:N-1).*(beta_ag.*R_ag(1:N-1)./La);
wn_new=dwn(1:N-1).*(beta_na.*R_na(1:N-1)./Ln);
wL_new=[wa_new.*dLa(1:N-1).*La wn_new.*dLn(1:N-1).*Ln]'; wL_new=wL_new(:);
firm_land_spending=[eta_ag.*R_ag_new eta_na.*R_na_new]'; firm_land_spending=firm_land_spending(:);
f_new=repmat(((1-alpha).*vL_new+firm_land_spending(1:2*(N-1)))./L0,1,2*(N-1)); % per-capita land earnings, source is row col is dest
v_new_matrix=repmat(w_new',2*(N-1),1)+f_new;
dP=dPa.^epsilon.*dPn.^(1-epsilon);
dcostliving=[(dP.^alpha.*dra.^(1-alpha)) (dP.^alpha.*drn.^(1-alpha))]'; dcostliving=dcostliving(:);
dV_new_matrix=(v_new_matrix./v_matrix)./repmat(dcostliving(1:2*(N-1))',2*(N-1),1);

% New migration matrix
mij_new=mij2000.*dV_new_matrix.^kappa; mij_new=mij_new./repmat(sum(mij_new,2),1,2*(N-1));

% New income
L1=(L0'*mij_new)';
reshape(vL_new',2,N-1)';
Ir_new=ans(:,1);
Iu_new=ans(:,2);
Ir_new(N,1)=Ir(N).*dwa(N).*dLa(N);
Iu_new(N,1)=Iu(N).*dwn(N).*dLn(N);

% New household final demand spending
Da_new_r=alpha*epsilon.*Ir_new;
Da_new_u=alpha*epsilon.*Iu_new;
Dn_new_r=(alpha*(1-epsilon)).*Ir_new;
Dn_new_u=(alpha*(1-epsilon)).*Iu_new;

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

% New labour allocation
L1=(L0'*mij_new)';

% Share of labour by sector ... does ag's share increase or decrease?
Lnew=reshape(L1',2,N-1)';
Lold=reshape(L1base',2,N-1)';
La_new=Lnew(:,1);
Ln_new=Lnew(:,2);
ag_labour_share_newold=sum([La_data(1:N-1) La_new(1:N-1)])./sum(L1);

% Welfare and GDP changes
dP=dPa.^epsilon.*dPn.^(1-epsilon);
dVa=dva./(dP.^alpha.*dra.^(1-alpha));
dVn=dvn./(dP.^alpha.*drn.^(1-alpha));
dV=[dVa dVn]'; dV=dV(:);
aggregate_welfare_change=sum(omega_welfare.*diag(dV_new_matrix).*(diag(mij_new)./diag(mij_base)).^(-1/kappa));
GDPshares=[(beta_ag+eta_ag).*R_ag (beta_na+eta_na).*R_na Ds_r Ds_u]./GDP;
GDPshares_r=[GDPshares(:,1)./(GDPshares(:,1)+GDPshares(:,3)) GDPshares(:,3)./(GDPshares(:,1)+GDPshares(:,3))];
dYa=sum(GDPshares_r.*[(R_ag_new./R_ag)./dPa (Ds_new_r./Ds_r)./dra],2);
GDPshares_u=[GDPshares(:,2)./(GDPshares(:,2)+GDPshares(:,4)) GDPshares(:,4)./(GDPshares(:,2)+GDPshares(:,4))];
dYn=sum(GDPshares_u.*[(R_na_new./R_na)./dPn (Ds_new_u./Ds_u)./drn],2);
dY=[dYa dYn]'; dY=dY(:);
realGDP_change=sum(dY(1:2*(N-1)).*Y./sum(Y)); % same as the first realGDP with the alphas

% Many different measures of migration patters, initially and in the counterfactual
migrant_shareofpop_2000=sum(L0.*(1-diag(mij)))./sum(L0);
mij=mij_new;
migrant_shareofpop_new=sum(L0.*(1-diag(mij)))./sum(L0);
mij_between_new=mij.*(1-kron(eye(N-1,N-1),ones(2,2))); % turn on for just between-province migration
mij_within_new=mij.*kron(eye(N-1,N-1),ones(2,2)).*(1-eye(2*N-2,2*N-2)); % turn on for just within-province migration
mij_ruralurban_new=mij.*repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1); % turn on for between ag and non-ag
mij_urbanrural_new=mij_ruralurban_new; 
mij_urbanrural_new(:,[2:2:60])=0; % Just U to R
mij_ruralurban_new(:,[1:2:60])=0; % Just R to U
migrant_shareofpop_within_new=sum(L0.*sum(mij_within_new,2))./sum(L0);
migrant_shareofpop_between_new=sum(L0.*sum(mij_between_new,2))./sum(L0);
migrant_shareofpop_ruralurban_new=sum(L0.*sum(mij_ruralurban_new,2))./sum(L0);
migrant_shareofpop_urbanrural_new=sum(L0.*sum(mij_urbanrural_new,2))./sum(L0);
mij_between_2000=mij2000.*(1-kron(eye(N-1,N-1),ones(2,2))); % turn on for just between-province migration
mij_within_2000=mij2000.*kron(eye(N-1,N-1),ones(2,2)).*(1-eye(2*N-2,2*N-2)); % turn on for just within-province migration
mij_ruralurban_2000=mij2000.*repmat([kron(ones(N-1,1),[0;1]) kron(ones(N-1,1),[1;0])],1,N-1); % turn on for between ag and non-ag
mij_urbanrural_2000=mij_ruralurban_2000; 
mij_urbanrural_2000(:,[2:2:60])=0; % Just U to R
mij_ruralurban_2000(:,[1:2:60])=0; % Just R to U
migrant_shareofpop_within_2000=sum(L0.*sum(mij_within_2000,2))./sum(L0);
migrant_shareofpop_between_2000=sum(L0.*sum(mij_between_2000,2))./sum(L0);
migrant_shareofpop_ruralurban_2000=sum(L0.*sum(mij_ruralurban_2000,2))./sum(L0);
migrant_shareofpop_urbanrural_2000=sum(L0.*sum(mij_urbanrural_2000,2))./sum(L0);
migrants_model_within=migrant_shareofpop_within_new/migrant_shareofpop_within_2000;
migrants_model_between=migrant_shareofpop_between_new/migrant_shareofpop_between_2000;

% Table 11
top_panel=[aggregate_welfare_change-1 realGDP_change-1 migrants_model_within-1 migrants_model_between-1]'
lower_panel=[ag_labour_share_newold;
    [migrant_shareofpop_urbanrural_2000 migrant_shareofpop_urbanrural_new];
    [migrant_shareofpop_ruralurban_2000 migrant_shareofpop_ruralurban_new]]
