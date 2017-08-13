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
alpha_0 = 10/180 *pi;
alpha_d= 0;      
phiF = 0;
phiF_d = 0;
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


wgR = 5;
wgB=  4*wgR;

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



%% Trafo in BNF 

SYS=ss(A,B,C,D);
[A,B,C,D]=ssdata(SYS);
[sysB,Tb]=canon(SYS,'companion');
[Ab,Bb,Cb,Db]=ssdata(sysB);
Tbinv=Tb^-1;

%% Fehlerhafte Trafo in RNF
% sysR=ss(Ab',Cb',Bb',Db');
% Ar=Ab';
% Br=Bb';
% Cr=Cb';



%% Trafo RNF
Su=[B A*B A^2*B A^3*B];

Suinv= inv(Su);

Sun=Suinv(end,:);

Tr= ([Sun; Sun*A; Sun*A^2; Sun*A^3]);

Trinv=inv(Tr);

Arr=Tr*A*Trinv;
Brr=Tr*B;
Crr=C*Trinv;


%%Koeffizienten auslesen


a(5)=1./Brr(4);

for i=1:4
   a(i)=-Arr(4,i).*a(5);
end


%Polynomvorgabe für Beobachter nach 4% Verfahren

% p0b = wgB^4;
% p1b = (3.02)*wgB^3;
% p2b = (4.36)*wgB^2;
% p3b = (3.02)*wgB;
% p4b = 1;

%Polynomvorgabe für Beobachter nach Butterworth

% p0b = wgB^4;
% p1b = sqrt(4+2*sqrt(2))*wgB^3;
% p2b = (2+sqrt(2))*wgB^2;
% p3b = sqrt(4+2*sqrt(2))*wgB;
% p4b = 1;

%Polynomvorgabe für Beobachter Binomial

p0b = wgB^4;
p1b = 4*wgB^3;
p2b = 6*wgB^2;
p3b = 4*wgB;
p4b = 1;

Pb = roots([p4b,p3b,p2b,p1b,p0b]);

% Beobachter Rückführmatrix

RbTr = place(Ab',Cb',Pb);
Rb= RbTr';
AM = Ab-Rb*Cb;
A_cl = AM;

%% Regler R

%Polynomvorgabe für Regler nach Butterworth 4% Verfahren 
% p0r = wgR^4;
% p1r = (3.02)*wgR^3;
% p2r = (4.36)*wgR^2;
% p3r = (3.02)*wgR;
% p4r = 1;

%Polynomvorgabe für Regler nach Butterworth

% p0r = wgR^4;
% p1r = sqrt(4+2*sqrt(2))*wgR^3;
% p2r = (2+sqrt(2))*wgR^2;
% p3r = sqrt(4+2*sqrt(2))*wgR;
% p4r = 1;
%Polynomvorgabe für Regler Binomial

p0r = wgR^4;
p1r = 4*wgR^3;
p2r = 6*wgR^2;
p3r = 4*wgR;
p4r = 1;

r(5)=1./Brr(end);

for k=1:4
   r(k)=-Arr(i,4).*r(5);
end


r1=p0r/p4r * r(5)-r(1);
r2=p1r/p4r * r(5)-r(2);
r3=p2r/p4r * r(5)-r(3);
r4=p3r/p4r * r(5)-r(4);

R = [r1 r2 r3 r4];

Pr = roots([p4r,p3r,p2r,p1r,p0r]);
R = place(A,B,Pr);
Ar_cl = A-B*R;


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


