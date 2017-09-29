clc, clear, close all hidden 

%% symbolische VAriablen definieren
% physikalische und geometrische Parameter des Systems
% syms JB JF rB rF mB g
% Kräfte und Zustandsvariablen, Koeffizienten
syms F MF alpha alpha_d alpha_dd phiF_d phiF_dd
% Elemente des Wirkungsplans/Zustandsraumes
syms x x_d u u1 u2 A B Y
% ???
syms s wgb wgR
% ???
syms Rq km iq iq_d uq Lq

%% Testwerte
%Rücksubstitution
F=0;
MF=0;
modellNull = [10/180*pi;0;0;0];
% Beobachterstartpunkt
beobachterNull = [0;0;0;0];
Mmax = 1000000;
Mmin = -1*Mmax;
% physikalische und geometrische Parameter des Systems
g = 9.81;
rB = 0.06; % Rollkreisdurchmeser des Balles
rBreal = rB+0.005; %wahrer Durchmesser Ball bei 5mm "Einsinken"
rF = 0.3; % Felgendurchmesser
mB = 0.3; % Ballmasse
mF = 3.3; % Masse Ersatzscheibe(bis Trägheitsmoment bekannt ist)
% Dummyträgheitsmomente für vereinfachte Modelle
JB = 2/3*mB*rBreal^2; % Für dünnwandige Hohlkugel(Matroids)
% JF = 1/2*mF*rF^2; % Für Zylinderförmige Felge

% JF=0.1; 

  for m=1:10
   
 
   %%% Parameter stellen
   JF=0.05+m*0.01; 

   
   
   
Ta=1/50; % Abtastrate für alpha


%% Koeffizienten der DGL aufgrund mechanischen Modells (siehe

Nphi = JF*JB+mB*(JF*rB^2+JB*rF^2);
%c1 = (mB*rB^2+JB)/Nphi;
%c1 = ((2*rF)/(2*mB*rF^2 + 7*JF)); %b2 bei thorben
c1 = (mB*rB^2+JB)/(JF*mB*rB^2+JB*mB*rF^2+JB*JF);% b2 in VL
%c3 = (JB*rF*mB*g)/Nphi;
%c3 = ((2*g*mB*rF)/(2*mB*rF^2 + 7*JF)); %b1 bei thorben
c3 = mB*g*JB*rF/(JF*mB*rB^2+JB*mB*rF^2+JB*JF);% b1 in VL

Nalpha = (rB+rF)*Nphi;
%d1 = (JB*rF)/Nalpha;
%d1 = ((2*mB*rF^2 + 5*JF)/(mB*(rB + rF)*(2*mB*rF^2 + 7*JF))); %a2 bei thorben
d1 = JB*rF/((rF+rB)*(JF*mB*rB^2+JB*mB*rF^2+JB*JF)); %a2 in VL
%d3 = ((JF*rB^2 +JB*rF^2)*mB*g)/Nalpha;
%d3 = ((2*g*mB^2*rF^2 + 5*JF*g*mB)/(mB*(rB + rF)*(2*mB*rF^2 + 7*JF)));% a1 bei thorben
d3 = mB*g*(JF*rB^2+JB*rF^2)/((rF+rB)*(JF*mB*rB^2+JB*mB*rF^2+JB*JF)); %a1 bei VL

%phiF_dd = c1*MF + c2*F + c3*sin(alpha);
%alpha_dd = d1*MF + d2*F + d3*sin(alpha);

%% Linearisierung für kleine Auslenkung um alpha=0

A= [0 1 0 0; d3 0 0 0;0 0 0 1; c3 0 0 0];   %Übertragungsmatrix

B= [0 ; d1; 0 ; c1];   % Eingagsmatrix

C= [1 0 0 0; 0 0 1 0]; %Ausgangsmatrix für Werte alpha, phi(!beob!), phipunkt

D= 0;

%Steuerbarkeit, Beobachtbarkeit
strbr = rank(ctrb(A,B));
bbtbr = rank(obsv(A,C));

%% Beobachter optimal Rückführmatrix
Q    = eye(4);
R    = diag([10, 10]);
L    = lqr(A',C',Q, R)';

%% optimaler Regler + Vorfilter
Q_ctrl = diag([   1, ...
                 0.1, ...
                 1, ...
                   10]);
R_ctrl = 1;

% optimaler Regler
K = lqr(A,B,Q_ctrl,R_ctrl);

% Vorfilter (stationäres Führungsgrößenfilter):
M   = 1 / ( (C)/(B*K-A)*B );

%% Simulation
   sim modell
   
   
   MF_max(m)=max(abs(simout.get.Data));


  end

   for n=1:10
   
 
   %%% Parameter stellen
   JF_P(:,n)=0.05+n*0.01;

   end
   
   
  plot(JF_P,MF_max,'-o','color',[0,0.5,0])
  xlabel('Trägheitsmoment der Felge [m^4]')
  ylabel('Anfangsdrehmoment der Felge [Nm]')
  grid on
  %set(gca,'Box','On','GridAlpha',1)
  set(findall(gcf,'Type','axes'),'FontSize',20,'LineWidth',1,'XColor','black','YColor','black')


