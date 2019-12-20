%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heteropolar Radial Magnetic Bearing
% Rotating Loss Calculation Program
%
% David Meeker, Ph.D.
% dmeeker@ieee.org
%
% 03Mar2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This program calculates rotoring losses in heterpolar radial
% magnetic bearings using a technique similar to that described
%
%   D. C. Meeker, A. V. Filatov, and E. H. Maslen, "Effect of 
%   Magnetic Hysteresis on Rotational Losses in Heteropolar 
%   Magnetic Bearings, IEEE Transactions on Magnetics, 
%   40(5):3302-3307, Sept 2004.
%   http://www.femm.info/dmeeker/pdf/01333140.pdf
%
% The paper used a magnetic scalar potential approach.
% This program translates the approach into one based on
% magnetic vector potential.  Since FEMM doesn't natively
% support the boundary conditions on harmonics described
% in the paper, these boundary conditions are applied
% iteratively.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Radial Magnetic Bearing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Useful Constants
  muo=pi*4.e-7;      % Permeability of free space, H/m
  RPM=pi/30.;		 % Conversion from RPM to radians/second
  deg=pi/180.;       % Conversion from degrees to radians
  inch=0.0254;       % Conversion from inches to meters

  % Stator Parameters
  Rri=1;             % Rotor Inner Radius, inches
  Rro=1.789;         % Rotor Outer Radius, inches
  nr=360;            % number of rotor harmonics considered
  g = 0.030;         % Rotor/stator gap, inches
  Rsi=6.150/2;       % Stator Inner Radius, inches
  Rso=7.725/2;       % Stator Outer Radius, inches
  wp=0.79;           % Pole width, inches
  depth=2*1.715;     % Bearing axial length, inches

  % Coil Properties
  hc=0.80;           % Coil height, inches
  wc=0.25;           % Coil thickness, inches
  turns=94;          % Number of turns per coil
  ic=[1 1 1 1];      % Vector of currents in each coil, Amps

  % Rotor Parameters
  mu=3460*muo;       % Magnetic permeability of rotor material
  sigma=7.46*10^6;   % Electrial Conductivity of rotor material
  phih=20*deg;       % Hysteresis lag angle
  fill=1;            % Lamination stacking factor
  d=0.014*inch;      % Thickness of an individual lamination
  omega=30000*RPM;   % Rotor speed

  % Description of air in which the rotor spins
  rhoAir = 1.204;   		% Density of air in kg/m^3
  nuAir  = 1.511*10^(-5);	% Kinematic Viscosity of Air in m^2/s
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Bearing Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('*** Heteropolar Radial Magnetic Bearing');
disp('*** Rotating Loss Calculation Program');
disp(' ');
disp('David Meeker, Ph.D.');
disp('dmeeker@ieee.org');
disp(' ');
disp('03Mar2006');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Definition

openfemm
newdocument(0)
mi_probdef(0,'inches','planar',10^(-8),depth,30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw Stator Geometry

np=8;                % Number of bearing poles

Rp=Rro+g;
ai=2*asin(wp/(2*Rp));
ao=2*asin(wp/(2*Rsi));

mi_drawarc(Rso,0,-Rso,0,180,5);
mi_drawarc(-Rso,0,Rso,0,180,5);

% list that will contain the locations of block labels
% that define the coil windings
labellist=zeros(2*np,1);

for k=0:(np-1),
  ab=2*k*pi/np - pi/np;
  theta1=ab+ai/2;
  theta0=ab-ai/2;

  % draw pole face %
  mi_drawarc(Rp*cos(theta0),Rp*sin(theta0),Rp*cos(theta1),Rp*sin(theta1),ai/deg,1);
  theta3=ab+ao/2;
  theta2=ab-ao/2;
  x0=Rp*cos(theta0); y0=Rp*sin(theta0);
  x2=Rsi*cos(theta2); y2=Rsi*sin(theta2);
  x1=Rp*cos(theta1); y1=Rp*sin(theta1);
  x3=Rsi*cos(theta3); y3=Rsi*sin(theta3);

  % Draw sides of pole %
  mi_drawline(x0,y0,x2,y2);
  mi_drawline(x1,y1,x3,y3);

  % Draw coil on one side of pole %
  n=j*((x3+j*y3)-(x1+j*y1)); n=n/abs(n); t=n/j;
  c=((x3+j*y3)+(x1+j*y1))/2;
  mi_drawline(real(c+t*hc/2),imag(c+t*hc/2),real(c+t*hc/2+n*wc),imag(c+t*hc/2+n*wc));
  mi_drawline(real(c-t*hc/2),imag(c-t*hc/2),real(c-t*hc/2+n*wc),imag(c-t*hc/2+n*wc));
  mi_drawline(real(c+t*hc/2+n*wc),imag(c+t*hc/2+n*wc),real(c-t*hc/2+n*wc),imag(c-t*hc/2+n*wc));
  labellist(2*k+1)=(c+n*wc/2);

  % Draw coil on other side of pole %
  n=j*(-(x2+j*y2)+(x0+j*y0)); n=n/abs(n);t=n/j;
  c=((x2+j*y2)+(x0+j*y0))/2;
  mi_drawline(real(c+t*hc/2),imag(c+t*hc/2),real(c+t*hc/2+n*wc),imag(c+t*hc/2+n*wc));
  mi_drawline(real(c-t*hc/2),imag(c-t*hc/2),real(c-t*hc/2+n*wc),imag(c-t*hc/2+n*wc));
  mi_drawline(real(c+t*hc/2+n*wc),imag(c+t*hc/2+n*wc),real(c-t*hc/2+n*wc),imag(c-t*hc/2+n*wc));
  labellist(2*k+2)=(c+n*wc/2);

  % Draw arc on stator inner diameter between poles %
  mi_drawarc(Rsi*cos(theta3),Rsi*sin(theta3),Rsi*cos(theta2+2*pi/np),Rsi*sin(theta2+2*pi/np),(2*pi/np-ao)/deg,5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Properties

% stator materials
mi_addmaterial('Stator Iron',2500,2500,0,0,0,0,0,1,0,0,0);
mi_addmaterial('Air',1,1,0,0,0,0,0,1,0,0,0);
mi_addmaterial('Coil',1,1,0,0,0,0,0,1,0,0,0);

% circuit properties

mi_addcircprop('i1',ic(1),1);
mi_addcircprop('i2',ic(2),1);
mi_addcircprop('i3',ic(3),1);
mi_addcircprop('i4',ic(4),1);

q=['i1';'i1';'i1';'i1';'i2';'i2';'i2';'i2';'i3';'i3';'i3';'i3';'i4';'i4';'i4';'i4'];

for k=1:length(labellist),
    mi_addblocklabel(real(labellist(k)),imag(labellist(k)));
    mi_selectlabel(real(labellist(k)),imag(labellist(k)));
    mi_setblockprop('Coil',0,(Rso-Rsi)/4,q(k,:),0,0,turns*sign(cos(pi*k/2-pi/4))); 
    mi_clearselected;
end

mi_addblocklabel(0,(Rsi+Rso)/2);  
mi_selectlabel(0,(Rsi+Rso)/2);
mi_setblockprop('Stator Iron',0,(Rso-Rsi)/4,'<None>',0,0,1); 
mi_clearselected;
mi_addblocklabel(0,(Rro+g/2)); 
mi_selectlabel(0,(Rro+g/2));
mi_setblockprop('Air',0,g/2,'<None>',0,1,1); 
mi_clearselected;
mi_addblocklabel(0,0);
mi_selectlabel(0,0);
mi_setblockprop('<No Mesh>',1,0,'<None>',0,1,1);
mi_clearselected;

mi_addboundprop('A=0',0,0,0,0,0,0,0,0,0);

mi_selectarcsegment(0,Rso);
mi_selectarcsegment(0,-Rso);
mi_setarcsegmentprop(5,'A=0',0,0);
mi_clearselected;

mi_zoomnatural;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Rotor

for k=1:(2*nr),
  mi_addboundprop(['rotor', mat2str(k)],0,0,0,0,0,0,0,0,2);
  mi_drawline(real(Rro*exp(k*j*pi/nr)),imag(Rro*exp(k*j*pi/nr)),real(Rro*exp((k-1)*j*pi/nr)),imag(Rro*exp((k-1)*j*pi/nr)));
  mi_selectsegment(real(Rro*exp((k-1/2)*j*pi/nr)),imag(Rro*exp((k-1/2)*j*pi/nr)));
  mi_setsegmentprop(['rotor',mat2str(k)],0,1,0,0);
  mi_clearselected;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Fourier Transform Matrices

M=zeros(nr-1,2*nr);
u=(0:(2*nr-1))*(pi/nr) + pi/(2*nr);
for k=1:(nr-1), M(k,:) =exp(-j*k*u)/nr; end
muh=mu*exp(-j*phih);

ro=Rro*inch;
ri=Rri*inch;
dz=depth*inch;
 
k=1:(nr-1); 
mun=(fill*muh)*tanh(sqrt(j*k*omega*sigma*muh)*d/2)./(sqrt(j*k*omega*sigma*muh)*d/2);
V=(-(k/ro).*(1./mun).*coth(k*log(ro/ri)));

Q=nr*M'*diag(V)*M;

% To get the c1 parameter applied to each facet, just 
% compute Q*A where A is the vector of rotor nodal potentials.

if (exist('OCTAVE_VERSION')==0)
	mi_saveas([cd,'\autobrg.fem']);
else
	mi_saveas('autobrg.fem');
end

dA=zeros(2*nr,1);
errnorm=1;
iter=0;
errtol=1e-6;
disp(sprintf('error goal = %g',errtol));

while (1>0)

  iter=iter+1;

  % Analyze geometry and update boundary condition on rotor
  mi_analyze;
  mi_loadsolution;
  A=mo_geta(Rro*real(exp(j*u)),Rro*imag(exp(j*u)));
  dAo=dA; 
  dA=real(Q*A);
  errnorm=sqrt((dA-dAo)'*(dA-dAo)/(dA'*dA));
  disp(sprintf('error = %g',errnorm));
  if (errnorm>errtol)
    for k=1:(2*nr),
      mi_modifyboundprop(['rotor',mat2str(k)],8,-dA(k));
    end
  else
    break;
  end
  
  % Compute zero speed flux distribution on first pass
  if (iter==1)
    bno=zeros(2*nr,1);
    for k=1:(2*nr),
      bno(k)=mo_getb(real(Rro*exp(j*u(k))),imag(Rro*exp(j*u(k))))*[real(exp(j*u(k)));imag(exp(j*u(k)))];
    end
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Rotating Losses

a=M*A;
pn=1:(nr-1);
for k=1:(nr-1),
  pn(k) = -(omega*pi*abs((k*a(k))/mun(k))^2*coth(k*log(ro/ri))*imag(mun(k))*dz);
end
ElectricalLoss=pn*ones(nr-1,1);
disp(sprintf('\nHysteresis + Eddy Current Loss = %g',ElectricalLoss));

disp('Dominant harmonics:')
for k=1:(nr-1),
  if (pn(k) > 0.001*ElectricalLoss)
    disp(sprintf('%i	%g',k,pn(k)));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimated Windage Loss
% from NASA TN D-4849
% http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19680027690.pdf

% Use a Newton iteration to solve the equation in D-4849 for Drag Coefficient
ReynoldsNumber = (Rro*inch)*(g*inch)*omega/nuAir;
Cdr = 1/(2*ReynoldsNumber);
while 1,  
  Z = 1/sqrt(Cdr) - (2.04 + 1.768*log(ReynoldsNumber*sqrt(Cdr)));
  if (Z < 1e-10) break; end;
  dZ = -1/(2*Cdr^(3/2)) - 0.884/Cdr;
  Cdr = (-Z + Cdr*dZ)/dZ;
end;

% evaluate the D-4849 expression for windage loss on a
% smooth cylinder rotating insida a concentric cylinder
WindageLoss=pi*Cdr*rhoAir*(Rro*inch)^4*omega^3*(depth*inch);
disp(sprintf('\nWindage Loss = %g',WindageLoss));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternative Loss Calc

% mo_groupselectblock(0);
% Torque=mo_blockintegral(22);
% RotorPower=Torque*omega;
% disp(sprintf('Weighted Stress Tensor Torque = %g',Torque));
% disp(sprintf('Power Derived from Torque = %g',RotorPower));
% mo_clearblock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot some results

bn=zeros(2*nr,1);
for k=1:(2*nr),
  bn(k)=mo_getb(real(Rro*exp(j*u(k))),imag(Rro*exp(j*u(k))))*[real(exp(j*u(k)));imag(exp(j*u(k)))];
end

plot(u/deg,bno,u/deg,bn);
xlabel('Position on Rotor Surface, Degrees');
ylabel('Normal Flux Density, Tesla');
title(sprintf('Zero Speed and %.0f RPM Normal Flux Density at Rotor Surface',omega/RPM ));

closefemm
