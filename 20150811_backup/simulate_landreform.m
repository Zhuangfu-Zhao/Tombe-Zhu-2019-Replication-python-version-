function F=simulate_landreform(X)

global N theta beta_ag beta_na eta_ag eta_na dni_ag dni_na pi_ag pi_na L0 kappa alpha dTa dTn  mij2000 La Ln
global  R_ag R_na Ir Iu  epsilon La_data Ln_data L1base f gamma 

dwa=X(1:N);
dwn=X(N+1:2*N);
dPa=X(2*N+1:3*N);
dPn=X(3*N+1:4*N);
dLa=X(4*N+1:5*N);
dLn=X(5*N+1:6*N);
dva=X(6*N+1:7*N); 
dvn=X(7*N+1:8*N);

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

% The system of equations
F=[
    
    dPa-denominator_ag.^(-1/theta);
    dPn-denominator_na.^(-1/theta);

    L1(1:2*(N-1))./L1base(1:2*(N-1))-dL(1:2*(N-1));
    dLa(N).*La_data(N)+dLn(N).*Ln_data(N)-(La_data(N)+Ln_data(N));
    dwa(N)-dwn(N);
    
    X_ag_new-(Da_new_r+Da_new_u+(1-beta_ag-eta_ag).*gamma(1,1).*R_ag_new+(1-beta_na-eta_na).*gamma(2,1).*R_na_new);
    X_na_new-(Dn_new_r+Dn_new_u+(1-beta_na-eta_na).*gamma(2,2).*R_na_new+(1-beta_ag-eta_ag).*gamma(1,2).*R_ag_new);
    
    vL_new(1:2*(N-1))-(wL_new(1:2*(N-1))+sum(f_new.*repmat(L0,1,2*(N-1)).*mij_new)'); % equation 15
    dva(N)-dwa(N);
    dvn(N)-dwn(N);

    ];
