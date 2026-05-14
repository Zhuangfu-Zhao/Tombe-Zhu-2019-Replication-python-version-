close all
clear all
global N theta beta_ag beta_na alpha eta_ag eta_na gamma kappa dni_ag dni_na pi_ag pi_na mu L La_data Ln_data L1base drYL_data La Ln

% Merge in the mu and price data
global drYL_data
data=importdata('employment_realGDP_data.csv'); % use data.csv for the BTZ version; "data - official data.csv" for official data
[N vars]=size(data.textdata);
region=data.textdata(2:N,1);
rYLa_data=data.data(:,1);
rYLn_data=data.data(:,2);
La_data=data.data(:,3);
Ln_data=data.data(:,4);
l=(La_data+Ln_data);
L(1:N,1)=sum(La_data+Ln_data);
rYL_data=(rYLa_data.*La_data+rYLn_data.*Ln_data)./(La_data+Ln_data);
mu=[rYLa_data rYLn_data]'; 
mu=mu(:);
rYLa_data2005=data.data(:,5);
rYLn_data2005=data.data(:,6);
drYLa_data=data.data(:,7);
drYLn_data=data.data(:,8);
drYL_data=[drYLa_data drYLn_data]'; 
drYL_data=drYL_data(:);

% International values
L(N,1)=2102.978587;

% Load Initial Trade Shares
pi_ag=importdata('trade_ag.csv'); 
pi_ag=pi_ag./repmat(sum(pi_ag,2),1,N);
pi_na=importdata('trade_na.csv'); 
pi_na=pi_na./repmat(sum(pi_na,2),1,N);

% Load Changes in Trade Costs
data=importdata('tauhat.csv');
dni_ag_measured=reshape(data.data(:,1),N,N)';
dni_na_measured=reshape(data.data(:,2),N,N)';
dni_asym_ag=reshape(data.data(:,3),N,N)';
dni_asym_na=reshape(data.data(:,4),N,N)';

% Load Migration Shares
global cij L0 mij2000
data=importdata('mij2000.csv');
mij2000=reshape(data.data(:,2),60,60)'; % each row is a province-industry pair
data=importdata('mij2005.csv');
mij2005=reshape(data.data(:,2),60,60)'; % each row is a province-industry pair

% Initial Employment Vector
L1=[La_data Ln_data]';
L1=L1(:); 

% Hukou Registrations, by Province
global La Ln
L0=(L1'*inv(mij2000))';
reshape(L1',2,N-1)';
La=ans(:,1);
Ln=ans(:,2);

% Production Function Input Shares
beta_ag=0.287; 
beta_na=0.219; % These reflect Brandt, Restuccia, and Adamopolous's paper for ag's VA shares + the IO data
eta_ag=0.277; 
eta_na=0.025;

% Remove comments here for the "no land" case
% eta_ag=0; eta_na=0;
% beta_ag=0.287+0.277; beta_na=0.219+0.025;

% Intermediate Input Shares
gamma=[0.397 0.603;
    0.063 0.937];
% gamma=eye(2,2); % Comment out for IO links
% rescale=(beta_ag+eta_ag); % no IO links at all
% beta_ag=beta_ag./rescale; eta_ag=eta_ag./rescale; 
% rescale=(beta_na+eta_na); % no IO links at all
% beta_na=beta_na./rescale; eta_na=eta_na./rescale;

% Common Parameters
alpha=1-0.13; % set to 1 for no housing; otherwise 1-0.13
kappa=1.5; 
theta=4;
epsilon=0.09545; % agriculture's share of total end-use in the IO tables

% Initial total nominal expenditures and income (see write-up)
global X_ag X_na R_ag R_na I Ir Iu omega_ag omega_ag2 omega_na2 epsilon
I0=zeros(N,1); 
I=ones(N,1);
X0_ag=zeros(N,1); 
X_ag=ones(N,1);
X0_na=zeros(N,1); 
X_na=ones(N,1);
while sum((X0_ag-X_ag).^2)+sum((X0_na-X_na).^2)>1e-20
    X0_ag=X_ag;
    X0_na=X_na;
    R_ag=pi_ag'*X0_ag;
    R_na=pi_na'*X0_na;

    Ir=(beta_ag+eta_ag).*R_ag./alpha; % Rural Income
    Iu=(beta_na+eta_na).*R_na./alpha; % Urban Income
    I=Ir+Iu; % Total Provincial Income
    Ir=Ir./sum(I); 
    Iu=Iu./sum(I); % Normalization

    Da_r=alpha*epsilon.*Ir; 
    Da_u=alpha*epsilon.*Iu; 
    Dn_r=alpha.*(1-epsilon).*Ir; 
    Dn_u=alpha.*(1-epsilon).*Iu; 
    Ds_r=(1-alpha).*Ir;
    Ds_u=(1-alpha).*Iu;

    X_ag=(Da_r+Da_u+(1-beta_ag-eta_ag).*gamma(1,1).*R_ag+(1-beta_na-eta_na).*gamma(2,1).*R_na);
    X_na=(Dn_r+Dn_u+(1-beta_na-eta_na).*gamma(2,2).*R_na+(1-beta_ag-eta_ag).*gamma(1,2).*R_ag);

end

%%
% Solve for initial and 2005 real income per worker
global Vi dVi_data
Via(1:N-1,1)=rYLa_data; % real income ag, 2000
Vin(1:N-1,1)=rYLn_data; % real income na, 2000
Via2005(1:N-1,1)=rYLa_data.*drYLa_data; % real income ag, 2005
Vin2005(1:N-1,1)=rYLn_data.*drYLn_data; % real income na, 2005
Vi=[Via Vin]'; 
Vi=Vi(:);
Vi2005=[Via2005 Vin2005]'; 
Vi2005=Vi2005(:);
dVi_data=[Vi2005./Vi;1.055;1.055];
reshape(dVi_data',2,N)';
dVa_data=ans(:,1);
dVn_data=ans(:,2);

% Solve for initial migration costs (per effective labour)
Vimat=repmat(Vi,1,2*(N-1));
Vjmat=repmat(Vi',2*(N-1),1);
cij=(mij2000./repmat(diag(mij2000),1,2*(N-1))).^(1/kappa)./(Vjmat./Vimat);
cij2000=cij;
Vimat2005=repmat(Vi2005,1,2*(N-1));
Vjmat2005=repmat(Vi2005',2*(N-1),1);
cij2005=(mij2005./repmat(diag(mij2005),1,2*(N-1))).^(1/kappa)./(Vjmat2005./Vimat2005);
mij=(cij2000.*Vjmat).^kappa;
mij=mij./(sum(mij,2)*ones(1,2*(N-1)));
cijhat=ones(2*(N-1),2*(N-1));
L1base=(L0'*mij)';
L2005=(L0'*mij2005)';

% Land Rebate Adjustments to Migration Costs
global Cnnjj Cnnjj2005 Cnijk Cnijk2005
M_mL_term=L1base./(diag(mij).*L0);
temp=reshape(M_mL_term,2,N-1)';
M_mL_term_ag=temp(:,1);
M_mL_term_na=temp(:,2);
Cnnjj_ag=1+(eta_ag+(1-alpha)*beta_ag).*M_mL_term_ag./(alpha.*beta_ag);
Cnnjj_na=1+(eta_na+(1-alpha)*beta_na).*M_mL_term_na./(alpha.*beta_na);
Cnnjj=[Cnnjj_ag Cnnjj_na]'; 
Cnnjj=Cnnjj(:);
Cnijk=(1-eye(2*(N-1),2*(N-1))).*repmat(Cnnjj,1,2*(N-1)).*cij;
Cnijk2000=Cnijk;
Cnnjj2000=Cnnjj;

% As above, but for 2005
M_mL_term2005=L2005./(diag(mij2005).*L0);
temp=reshape(M_mL_term2005,2,N-1)';
M_mL_term_ag2005=temp(:,1);
M_mL_term_na2005=temp(:,2);
Cnnjj_ag2005=1+(eta_ag+(1-alpha)*beta_ag).*M_mL_term_ag2005./(alpha.*beta_ag);
Cnnjj_na2005=1+(eta_na+(1-alpha)*beta_na).*M_mL_term_na2005./(alpha.*beta_na);
Cnnjj2005=[Cnnjj_ag2005 Cnnjj_na2005]'; Cnnjj2005=Cnnjj2005(:);
Cnijk2005=(1-eye(2*(N-1),2*(N-1))).*repmat(Cnnjj2005,1,2*(N-1)).*cij2005;

% International initial labour distribution
Ln_data(N,1)=L(N)./(1+(beta_ag./beta_na).*(R_ag(N)./R_na(N)));
La_data(N,1)=L(N)-Ln_data(N);

% Average symmetric trade costs
tau_ag=((diag(pi_ag)*diag(pi_ag)')./(pi_ag.*pi_ag')).^(1./(2.*theta));
tau_na=((diag(pi_na)*diag(pi_na)')./(pi_na.*pi_na')).^(1./(2.*theta));

%% Get the Initial Equilibrium of the Model
global dTa dTn
dTa=ones(N,1);
dTn=ones(N,1);
dni_ag=ones(N,N);
dni_na=ones(N,N);
new=fsolve(@(X)main_simulate(X),ones(6*N,1),optimset('Display','off'));
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

% Initial Welfare Vector
global omega_welfare
omega_welfare=(L0./sum(L0)).*Vi.*Cnnjj2000.*diag(mij).^(-1/kappa)./sum((L0./sum(L0)).*Vi.*Cnnjj2000.*diag(mij).^(-1/kappa));

% Sets base values, for possible use later
global mij_base dL_base L1base
mij_base=mij;
dL_base=dL;
L1base=L1;
