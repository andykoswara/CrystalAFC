function [fn,Tw,Te,T,f,delta,Rf,C,Sigma,rest,L43,blockage,t,x,r,L] = ...
    encrustpfcV6(Nx,flow,xmax,ti,tf,Nt,Tc,Tw0,Te0,T0,Tfeed,C0,delta0,f0)
%""""""""
    % Simulate PFC-PBM-Encrustation ODEs with high-resolution finite 
        % volume (HRFV) 
    % Args:
      % Nx: grid resolution for reactor length [-]
      % flow: linear flow rate [m/min]
      % xmax: reactor length [m]
      % ti: initial time [s]
      % tf: final time [s]
      % Nt: grid resolution for time [-]
      % Tc: jacket temperature [^{o}C]
      % Tw0: initial wall temperature [^{o}C]
      % Te0: initial encrust temperature [^{o}C]
      % T0: initial tube temperature [^{o}C]
      % Tfeed: initial feed temperature [^{o}C]
      % C0: initial tube concentration []
      % delta0: initial encrust thickness [m]
      % f0: initial number of crystals [-]
    % Returns:
      % fn: ode output vector
      % Tw: wall temperature [^{o}C]
      % Te: encrust temperature [^{o}C]
      % T: tube temperature [^{o}C]O
      % f: number of crystals [-]
      % delta: encrust thickness [m]
      % Rf: flow radius [m]
      % C: tube concentration [kg crystal/kg solvent]
      % Sigma: level of supersaturation [-]
      % rest: residence time [-]
      % L43: crystal mean size [micron]
      % blockage: percentage of tube closing due to encrust [%]
      % t: time vector
      % x: reactor x-directional vector
      % r: reactor radial vector
      % L: reactor length [m]

    % NOTE:
    % V1: model starts with plug-flow crystallization with encrustation dynamics
      % but without mass and heat balance on the tube side
    % V2: includes encrustation with mass balance, but w/out heat balance
      % on the tube side
    % V3: includes encrustation with mass balance and heat flux across
      % tube wall and encrust. Heat conduction in axial direction is also added
      % for wall side and encrust
    % V4: PBM and mass balance correction, since x-section is not constant
      % Here mass balance is updated which takes into account the deposition
      % at the crystallizer wall
    % V5: includes encrustation AND decrustation with mass balance and
      % heat flux across tube wall and encrust, but without crystal formation.
    % V6: includes en/decrustation AND crystal growth and dissolution
      % with mass balance and heat flux across tube wall and encrust. Further
      % modified to take initial conditions as input.
  
%""""""""

%% Initialize parameters
Nr=20; % grid resolution for reactor radius
NL=20; % for crystal length
xmin=0; % [m]
delx=(xmax-xmin)/Nx;
x=delx/2:delx:xmax-delx/2;
Lmax=200e-6; % [m], maximum size of the crystals considered
Lmin=0; % minimum size
delL = (Lmax-Lmin)/NL; % [m], grid size
L=(delL/2:delL:Lmax-delL/2); % [m], grid range
L3 = L'.^3*ones(1,Nx); % grid volume range
e = 1e-10; % small constant to avoid division by zero in gradient
    % calculation
    
save('grid_params.mat', 'Nr', 'NL', 'delx', 'x', 'Lmax', 'Lmin',... 
    'delL', 'L', 'L3', 'e')

%% for robustness analysis
R_kr=1;
R_alpha=1;
R_kg=1;
R_kd=1;

%% tube-side properties
ID=1.50E-2; % [m], inner tube diameter
Ri=ID/2; % radius (7.5 mm)
R0=ID/2+1.2e-3; % Outer radius, 1.2 mm thick pyrex glass
Cp_w=753; % [J/(kg.K)]
K_w=1.005; % [W/(m.K)]
rho_w=2230; % [kg/m^3]
delR=(R0-Ri)/Nr;
R=(Ri+delR/2:delR:R0-delR/2)';
R=R*ones(1,Nx);
delr=1/Nr; % dimensionless r
r=(0+delr/2:delr:1-delr/2)';r=flipud(r); % flipping is used because
    % direction of dimensionless r is opposite
r=r*ones(1,Nx);

h=500*2; % [W/(m^2.K)], taken for tubular exchanger from internet
hT=500*2; % [W/(m^2.K)], tube side overall heat transfer coefficient
Cp=4185.5; % [J/(Kg.K)], slurry specific heat capacity
rho_l=1080; % [kg/m^3], density
K=0.58; % [W/(m.K)], conductivity
    % rho,Cp taken for water
eta=600e-6; % [Pa s], Viscosity [Pa s]
freq=2;
ampl=0.04; % [m]
w=2*pi*freq*ampl; % velocity used for mass transfer correlation
Reynolds=w*rho_l*ID/eta; % Oscillatory Reynolds Number

save('tube_params.mat', 'ID', 'Ri', 'R0', 'Cp_w', 'K_w', 'rho_w', 'delR',...
    'R', 'delr', 'r', 'h', 'hT', 'Cp', 'rho_l', 'K', 'eta', 'freq',...
    'ampl', 'w', 'Reynolds')

%% Crystal parameters (Potash-Alum system)
% Crystal geometric properties
phiV=0.62; % [-], volume shape factor
rho_c=1750; % [kg/m^3], crystal density
Rg=8.314; % [J/K.mol] ideal gas constant

% nucleation kinetics (primary and secondary)
ja=1.7e8; % primary
jb=5.64e6; % ------
kb=3.14e7; % secondary
j=1;% ------
b=1.32; % -------

% growth kinetics
Kg0=2.05e5; % growth rate constant
gamma=7.18e2; % variable related to size-dependent growth
beta=6.1e-5; % ---------------------------------------
delEg=5.77e4; % activation energy
g=1.42; % reaction order

% dissolution kinetics (from "Fines Removal..." paper)
Kd=10.7192e-9; % [m/s], dissolution rate constant
d=0.5122; % [], reaction order
q=0.3382; % [], another order

% solubility
a_sol=[0.0000458354514087806  0.000243183573604693 ...
    0.04633927075106140]; % potash-alum

save('crys_params.mat', 'phiV', 'rho_c', 'Rg', 'ja', 'jb', 'kb', 'j', ...
    'b', 'Kg0', 'gamma', 'beta', 'delEg', 'g', 'Kd', 'd', 'q', 'a_sol')

%% Encrustation parameters
Rf0=Ri-delta0; % flow radius
voidage=0.2;% [], from Matthias paper page 122
Cp_e=870; % J/(kg.K) sodium chloride
rho_e=(1-voidage)*rho_c+voidage*rho_l; % [Kg/m^3],
K_e=1.35; % [W/(m.K)], thermal conductivity (CaSO4)
dp=36e-6; % [m], mean crystal size of the fouling layer
Diff=1.57e-9; % [m2/s], diffusivity
Kt=1e-6; % Linear coefficient of expansion [1/K]
Kr0=2.36;% [kg/(m^2.s)] since our concentration is in kg/kg, not in kg/m^3,
    % This value is corresponding to 2.36e6 [m^4/(kg.s)]. Thus a higher value
    % was used
Grav=9.81; % [m/s^2], gravitational accl
alpha_e=1.27e-2;% % dissolution rate

save('encrust_params.mat', 'Rf0', 'voidage', 'Cp_e', 'rho_e', 'K_e',...
    'dp', 'Diff', 'Kt', 'Kr0', 'Grav', 'alpha_e')

%% simulation parameters
% volumetric flow rates
Vsol=flow/60*1e-6; % [m^3/sec], solvent flow rate
Af0=pi*Rf0.^2; % [m^2], tube side (effective) area
ux0=Vsol./Af0; % [m/s], tube flow rate

%% initial conditions
T0A=T0.*Af0;
C0A=C0.*Af0; %newly defined conserved variable (for mass balance)
f0A=f0.*(ones(NL,1)*Af0); % [#], number of crystals
Tw0=reshape(Tw0,Nr*Nx,1);
Te0=reshape(Te0,Nr*Nx,1);
f0A=reshape(f0A,NL*Nx,1);
fn0=[Tw0;Te0;T0A';f0A;delta0';C0A'];
tspan=ti:(tf-ti)/Nt:tf; % time span

save('HRFV_params.mat', 'Vsol', 'Af0', 'ux0', 'T0A', 'C0A', 'f0A', ...
     'Tw0', 'Te0', 'fn0', 'tspan')

options=odeset('OutputFcn',@odeprint);
[t,fn]=ode15s(@odefun,tspan,fn0,options); % solving odes

%% output processing
Tw=fn(end,1:Nr*Nx);Tw=reshape(Tw,Nr,Nx);
Te=fn(end,Nr*Nx+1:2*(Nr*Nx));Te=reshape(Te,Nr,Nx);
TA=fn(:,2*(Nr*Nx)+1:2*(Nr*Nx)+Nx);
fA=fn(end,2*(Nr*Nx)+Nx+1:2*(Nr*Nx)+Nx+(NL*Nx));fA=reshape(fA,NL,Nx);
delta=fn(:,2*(Nr*Nx)+(NL*Nx)+Nx+1:2*(Nr*Nx)+(NL*Nx)+2*Nx);
CA=fn(:,2*(Nr*Nx)+(NL*Nx)+2*Nx+1:2*(Nr*Nx)+(NL*Nx)+3*Nx);
Rf=Ri-delta;
Af=pi*Rf.^2;
T=TA./Af;
f=fA./(ones(NL,1)*Af(end,:));
C=CA./Af;
C_sat=polyval(a_sol,T);
Sigma=C./C_sat; %supersaturation
rest=sum(Af.*delx,2)/Vsol;
Nt=size(t,1);
L43=zeros(Nt,1);
for i=2:Nt
    fAtemp=fn(i,2*(Nr*Nx)+Nx+1:2*(Nr*Nx)+Nx+(NL*Nx));fAtemp=...
        reshape(fAtemp,NL,Nx);
    ftemp=fAtemp./(ones(NL,1)*Af(i,:));
    L43(i)=trapz(L,L'.^4.*ftemp(:,end))/trapz(L,L'.^3.*ftemp(:,end))*1e6;
        % [mu-m]
end
blockage=(pi*(ID/2)^2-pi*(ID./2-max(delta,[],2)).^2)./(pi*(ID/2)^2)*100;
    % [%]

%% encrustation ODE

    function df_dt=odefun(t,fn)

    %% assigning input variables
    Tw=fn(1:Nr*Nx);Tw=reshape(Tw,Nr,Nx);
    Te=fn(Nr*Nx+1:2*(Nr*Nx));Te=reshape(Te,Nr,Nx);
    TA=fn(2*(Nr*Nx)+1:2*(Nr*Nx)+Nx);TA=TA';
    fA=fn(2*(Nr*Nx)+Nx+1:2*(Nr*Nx)+Nx+(NL*Nx));fA=reshape(fA,NL,Nx);
    delta=fn(2*(Nr*Nx)+(NL*Nx)+Nx+1:2*(Nr*Nx)+(NL*Nx)+2*Nx);delta=delta';
    CA=fn(2*(Nr*Nx)+(NL*Nx)+2*Nx+1:2*(Nr*Nx)+(NL*Nx)+3*Nx);CA=CA';

    %% update values
    Rf=Ri-delta;
    Af=pi*Rf.^2; % Flow area
    C=CA./Af;
    T=TA./Af;
    f=fA./(ones(NL,1)*Af);
    ux=Vsol./Af; % fluid velocity

    if ~isreal(C)
        disp('Warning: C is complex !!')
        pause
    end

    %% wall temperature
    T_EW=(Te(Nr,:)+Tw(1,:).*delta*K_w/K_e*delr/delR)./...
        (1+delta*K_w/K_e*delr/delR);   %B.C. for Encrust-Wall interface

    fdiffy=diff(Tw);
    diffr=([2*(Tw(1,:)-T_EW);fdiffy;2*(Tc-Tw(Nr,:))]);
    diffr2=(diffr(1:end-1,:)+diffr(2:end,:))/2; % better-approx. of
        % first derivative
    diffr=diff(diffr); % second derivative

    diffx=[2*(Tw(:,1)-Tfeed*ones(Nr,1)) diff(Tw,1,2) zeros(Nr,1)];
    diffx2=diff(diffx,1,2);

    temp=K_w/(rho_w*Cp_w)*(diffr/delR^2+diffr2./R/delR + K_w*diffx2);
    df_dt=reshape(temp,Nr*Nx,1);

    %% encrust temp
    T_TE=(2*K_e./delta/delr.*Te(1,:)+hT*T)./(2*K_e./delta/delr+hT);
        % Tube side-Encrust interface temperature from B.C.
    fdiffy=diff(Te);
    diffr=([2*(Te(1,:)-T_TE);fdiffy;2*(T_EW-Te(Nr,:))]);
    diffr2=(diffr(1:end-1,:)+diffr(2:end,:))/2;
    diffr=diff(diffr);

    diffx=[2*(Te(:,1)-Tfeed*ones(Nr,1)) diff(Te,1,2) zeros(Nr,1)];
    diffx2=diff(diffx,1,2);
    Delta=ones(Nr,1)*delta;
    temp=K_e/(rho_e*Cp_e)*(diffr./(delr^2*Delta.^2)+diffr2./...
        (Ri-r.*Delta)./Delta/delr+diffx2); % Heat conduction in axial direction
    df_dt=[df_dt;reshape(temp,Nr*Nx,1)];

    %% tube temp
    fdiffx=diff(TA);
    Tflux=zeros(1,Nx+1);
    Tflux(1)=ux0(1)*Tfeed*Af0(1);
    Tflux(2:end)=ux.*TA;
    temp=-1/delx*diff(Tflux)+K/(rho_l*Cp*delx^2)*Af.*diff([(T(1)-Tfeed)*2 fdiffx 0]) ...
        +2*pi*Rf./(rho_l*Cp).*h.*(T_TE-T);
    df_dt=[df_dt;reshape(temp,Nx,1)];

    %% CSD kinetics via HRFV

    % CSD flux
    % x direction (along the PFR)
    fluxx=zeros(NL,Nx+1); % flux along PFC length
    Gxf = (ones(NL,1)*ux).*fA; % u*n
    fdiffx=diff(Gxf,1,2);
    rx=(fdiffx(:,1:end-1)+e)./(fdiffx(:,2:end)+e);
    absrx=abs(rx);
    phix = (absrx+rx)./(1+absrx); % phi(r_i)
    fluxx(:,1)=ux0(1).*f0(:,1).*Af0(1);
    fluxx(:,2)=(Gxf(:,1)+Gxf(:,2))/2;
    fluxx(:,3:end-1)=Gxf(:,2:end-1)+0.5*phix.*fdiffx(:,2:end);
    fluxx(:,Nx+1)=Gxf(:,Nx);

    % y direction (along crystal size)
    C_sat=polyval(a_sol,T);
    tempG1=1-exp(-gamma*(L+beta)); % growth kinetic expression 1
        % (size factor)
    tempG2=Kg0*exp(-delEg/Rg./(273.15+T)); % expression 2
    tempD=-Kd./(L.^q); % dissolution

    GD=zeros(NL,Nx); % crystal growth/dissolution rate (size-dependent)
    S=(C./C_sat); % [-], supersaturation
    sigma=S-1;% [-], relative supersaturation
    for k=1:Nx
        if C(k)>=C_sat(k) % in case of growth
            sigmaG=sigma(k)^g; % reaction order
            GD(:,k)= tempG1*tempG2(k)*sigmaG;
            GD(:,k)=R_kg*GD(:,k); % robustness analysis
        else % dissolution
            sigmaG=(-sigma(k))^d;
            GD(:,k)=tempD'*sigmaG;
            GD(:,k)=R_kd*GD(:,k); % robustness analysis
        end
    end

    %% nucleation (B.C for growth)
    Mt=rho_c*phiV*sum(f.*L3)*delL; % magma density, suspension density
    Jprim=ja*exp(-jb./(273.15+T).^3./(log(S)).^2); % nucleation kinetcis
        % (primary)
    Jsec=kb*Mt.^j.*(C-C_sat).^b; % secondary
    B0=Jprim+Jsec;

    Gyf = GD.*fA;
    fdiffy=diff(Gyf);
    fluxy=zeros(NL+1,Nx); % flux along crystal size
    for k=1:Nx
        if C(k)>=C_sat(k)
            fluxy(1,k)=B0(k)*Af(k);
            fluxy(2,k)=(Gyf(1,k)+Gyf(2,k))/2;
            ry = (fdiffy(1:end-1,k)+e)./(fdiffy(2:end,k)+e);
            absry = abs(ry);
            phiy = (absry+ry)./(1+absry); % phi(r_j)
            fluxy(3:end-1,k)=Gyf(2:end-1,k)+0.5*phiy.*fdiffy(2:end,k);
            fluxy(NL+1,k)=Gyf(NL,k);
        else
            fluxy(1,k)=Gyf(1,k);
            ry = (fdiffy(2:end,k)+e)./(fdiffy(1:end-1,k)+e);
            absry = abs(ry);
            phiy = (absry+ry)./(1+absry); % phi(r_j)
            fluxy(2:NL-1,k)=Gyf(2:end-1,k)-0.5*phiy.*fdiffy(1:end-1,k);
            fluxy(NL,k)=(Gyf(NL,k)+Gyf(NL-1,k))/2;
            fluxy(NL+1,k)=0;
        end
    end

    temp=-diff(fluxx,1,2)/delx-diff(fluxy)/delL;
    df_dt=[df_dt;reshape(temp,NL*Nx,1)];

    %% encrust thickness
    Tf=T+0.55*(T_TE-T); % film temperature
    Csat_f=polyval(a_sol,Tf);

    km=Diff./(2*Rf)*0.034*(Reynolds)^0.875*(eta/rho_l/Diff)^0.33;
        % mass transfer coeff
    kr=Kr0*exp(-37143./(8.314*(273.14+Tf))); % reaction rate, [brahim et al.]
    kr=R_kr*kr; % robustness analysis
    P_by_K6=83.2*w.^0.54; % [N]

    % for 2nd order reaction
    temp1=km./(rho_e*K_e).*(0.5.*(km./kr)+rho_l.*(C-Csat_f)-...
        (1/4*(km.^2./kr.^2)+km./kr.*rho_l.*(C-Csat_f) ).^0.5)...
        -1/P_by_K6*(1+Kt*(Te(1,:)-Te(end,:)))*dp*(rho_l^2*eta*Grav)^...
        (1/3).*w.^2.*delta./K_e; % multiplication by liquid
        % density rho_l  in the above expression is to convert the
        % concentration from kg/kg to kg/m^3

    temp2=alpha_e*(C-Csat_f); % decrustation rate
    temp2=R_alpha*temp2; % robustness analysis

    temp=zeros(1,Nx);
    for k=1:Nx
        if C(k)>=Csat_f(k) % encrustation
            temp(k)=K_e*temp1(k); % encrust thickness
        else % decrustation
            if delta(k)<=1e-6
%                 delta(k)=1e-6;
                temp(k)=0;% prevent negative delta
            else
                temp(k)=temp2(k);
            end
        end
    end

    df_dt=[df_dt;reshape(temp,Nx,1)];

    %% mass balance
    cflux=zeros(1,Nx+1);
    cflux(1)=ux0(1)*C0(1)*Af0(1);
    cflux(2:end)=ux.*CA;
    flowmass=-1/delx*diff(cflux); % due to in/out flow
    crysmass=-((3*phiV*rho_c*delL)*(L.^2*Gyf)*(258.21/474.39)/rho_l);
        % crystal growth/dissolution
    crustmass=-rho_e/rho_l*2*pi*(Ri-delta).*temp; % en/decrustation
    temp=flowmass+crysmass+crustmass;
    df_dt=[df_dt;reshape(temp,Nx,1)];

    end % ode function
end % Main function
