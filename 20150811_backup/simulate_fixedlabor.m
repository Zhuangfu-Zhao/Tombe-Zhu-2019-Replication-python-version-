function F=simulate_fixedlabor(X)

global N theta beta_ag beta_na eta_ag eta_na dni_ag dni_na pi_ag pi_na alpha dTa dTn
global R_ag R_na Ir Iu epsilon gamma

dwa=X(1:N);
dwn=X(N+1:2*N);
dPa=X(2*N+1:3*N);
dPn=X(3*N+1:4*N);
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

% The system of equations
F=[
    
    dPa-denominator_ag.^(-1/theta);
    dPn-denominator_na.^(-1/theta);
    
    X_ag_new-(Da_new_r+Da_new_u+(1-beta_ag-eta_ag).*gamma(1,1).*R_ag_new+(1-beta_na-eta_na).*gamma(2,1).*R_na_new);
    X_na_new-(Dn_new_r+Dn_new_u+(1-beta_na-eta_na).*gamma(2,2).*R_na_new+(1-beta_ag-eta_ag).*gamma(1,2).*R_ag_new);

    ];
