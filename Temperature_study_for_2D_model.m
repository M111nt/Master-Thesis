
%%
% temperature study
%%

clear all;
clc;
%close all;

hbar = 1.054e-34; % in J.S
q=1.602e-19;
m0=9.3e-31;
e0=8.85e-12;
k=1.38e-23;
ms = 0.04*m0;

tw=13e-9;
tox=5e-9;
Lg=80e-9;
Lg1=[50]*1e-9;
% Lg1=80e-9;
Ef=0;
erox=18;
ertw=13;

Tr=0.75;
W = 10e-9;H=5e-9;

Cox = erox*e0/tox; % planar
%Cox = Lg*2*erox*e0*(W+H)/tox + 2.23*erox*e0; % rectangular wire and valid for tox << W,H, where W and H are width and height of the nanosheetg
%Cox= Cox/(W+H);

Cc = ertw*e0/(0.36*tw);% planar

%Cc = 6.94*erox*e0*(W^2+H^2)/(W*H); % rectangular GAA (COMSOL simulations) 

%Cc = Cc/(W+H);

Cg = (1/Cox+1/Cc)^-1;

lambda = 3/pi * (tw + tox*ertw/erox);% 2D MOSFET
%lambda = sqrt(ertw/((W/H)*erox) * tox*tw); % rectangualr GAA

etopvec = linspace(-0.5,1);

alpha=1;% non parabolicity 
%


    
% for ii=1:length(etopvec)
%     n0 = @(E) q*D2Dfun(E,etopvec(ii),ms,alpha).*fd(E,Ef,T);
%     n0_tot(ii) = integral(n0,etopvec(ii),0.7);
% end

% figure(111);plot(etopvec,n0_tot);

% pp =spline(etopvec,n0_tot);

% currenet-voltage characteristics

nbr_steps=11;

Vgvec=linspace(-0.1,1,nbr_steps);
%Vdvec = [ 0.05 0.5];
%Vdvec=0:0.05:1;
Vdvec=0.5;
T_loop=100:200:500;

VT=0.5;
for tt = 1:length(T_loop)
    
    T = T_loop(tt);% temperature study
    
    for ii=1:length(etopvec)
    n0 = @(E) q*D2Dfun(E,etopvec(ii),ms,alpha).*fd(E,Ef,T);
    n0_tot(ii) = integral(n0,etopvec(ii),0.7);
    end
    
    pp =spline(etopvec,n0_tot);

    
for kk=1:length(Lg1)
    Lg=Lg1(kk);
for jj=1:length(Vdvec)
    
    Vd=Vdvec(jj);
    
    for ii=1:length(Vgvec)
        
        Vg = Vgvec(ii);
        
        xtop = xtopapprox(0.2,Vd,lambda,Lg);
        
        etop = etopapprox(0.2,Vd,lambda,Lg);
        
        ag =(Ecfun(xtop,0.2,Vd,lambda,Lg) - Ecfun(xtop,0.2+10e-3,Vd,lambda,Lg) )/10e-3;
        
        %ag=1;
        Csigma = Cg/ag;
        
        ad = (Ecfun(xtop,0.2,Vd,lambda,Lg) - Ecfun(xtop,0.2,Vd+10e-3,lambda,Lg) )/10e-3;
 %       ad=0;
        Cdrain = ad*Csigma;        
        
        Efd = Ef-Vd;
        
        fzerofunc = @(x) -ag*(Vg-VT) - ad*Vd + q*((2-Tr)*fnval(pp,x)+Tr*fnval(pp,x+Vd))/Csigma - x;
                
        etop = fzero(fzerofunc,Vg/2);
               

        i0fun=@(E) sqrt(q)*2*q^2/hbar * Tr .* Mfun(E,etop,ms,alpha) .* (fd(E,Ef,T)-fd(E,Efd,T)) .* (E>etop);
        
        itot(jj,ii) = integral(i0fun,etop,etop+1);
        

        
    end

 
end
%%
% Transfer characteristics
figure(1);hold on;
title('Temperature study when Vds=0.5V, Vt=0.5V, Tr=0.75');
plot(Vgvec-VT,itot(1,:)*1e-3,'LineWidth',2);

legend('100K','300K','500K');
xlabel('V_{GS} (V)');
ylabel('I_{DS} (mA/\mu m)');


%gm
% kk = 1;
% for kk = 1:1:5
%     x=Vgvec;
%     qq = 1 + 4*kk;
%     y=itot(qq,:)*1e-3/(2*(W+H));
% %     y=itot(6,:)*1e-3/(2*(W+H));
%     x1=(length(Vgvec)-1)/2;
%     pp_DC = csaps(x,y, [1,ones(1,x1),repmat(5,1,x1)], [], ...
%             [ones(1,x1), repmat(5,1,x1), 0]);
%         
%     figure(2);hold on;
%     fnplt(fnder(pp_DC));
%     
%     
%     ylabel('g_{m}(mS/\mum)');xlabel('V_{GS} (V)');
%     legend('V_{DS} = 0.2', 'V_{DS} = 0.4', 'V_{DS} = 0.6','V_{DS} = 0.8','V_{DS} = 1.0');
%     title('The quasi 2D ballistic simulation for V_{GS} and g_{m} when Tr = 0.5');
% end 
        
        
% SS (mV/dec)   
% Inv_ss = [0,diff(log10(itot))./diff(Vgvec)];
% ss1 =(1./Inv_ss);
% ss(kk)=min(ss1(ss1>0))*1e3; 
% 


% 
% gm (mS/um) 
% 
% gm(kk)=peak([0,diff(itot)./diff(Vgvec)]);


%Output Characteristics

% figure(3);hold on;
% title('Output Characteristics');
% %legend('-0.1,0.01,0.12,0.23,0.34,0.45,0.56,0.67,0.78,0.89,1');
% for ii=1:length(itot(1,:))
% plot(Vdvec,itot(:,ii)*1e-3,'b','LineWidth',2);hold on; xlabel('V_{DS} (V)');ylabel('I_{DS} (mA/\mu m)');
% end


% Ron (ohm-um)
% Ron=[0;diff(itot(:,7))./diff(Vdvec')];
% figure(4);plot(Vdvec,Ron);
% gd (mS/um)
% plot(Vgvec,itot);hold on;
end
end
%%
% figure(100); scatter(Lg1,ss);% figure(200); scatter(W1,ss);
% figure(101); scatter(Lg1,gm);% figure(201); scatter(W1,gm);
% figure(102); scatter(Lg1,DIBL);% figure(202); scatter(W1,DIBL);
% figure(103); scatter(Lg1,gd);% figure(203); scatter(W1,gd);
% figure(104); scatter(Lg1,Ron);% figure(204); scatter(W1,Ron);


%%
% for kk=1:length(Vgvec)
% Etop(kk) = etopapprox(Vgvec(kk),0.5,lambda,Lg);
% 
% end
% 
% plot(Vgvec,Etop);
function x0 = xtopapprox(Vg,Vd,lambda,Lg)

VT=0.5;
ns=(VT-Vg).*(VT>Vg);
x0 = lambda/2*log(ns./(ns+Vd))+Lg/2;

end

function e0 = etopapprox(Vg,Vd,lambda,Lg)

VT=0.5;
ns=(VT-Vg).*(VT>Vg);
e0 = ns - 2*sqrt(ns.*(ns+Vd))*exp(-Lg/(2*lambda));
end

function M = Mfun(E,Ec,ms,alpha)

hbar = 1.054e-34;
M= sqrt(2*ms*(E-Ec).*(1+alpha*(E-Ec)))/(pi*hbar) .* (E>Ec);

end


function Ecn = Ecfun(xtop,Vg,Vd,lambda,Lg)

k1=1/lambda;
VT=0.5;
ns = VT-Vg;
Ecn= -(ns*sinh((Lg-xtop)*k1)+(ns+Vd)*sinh(xtop*k1))/sinh(Lg*k1)+ns;

end

function f = fd(E,Ef,T)

q=1.602e-19;
k=1.38e-23;
f = 1./(1+exp((E-Ef)/(k*T/q)));

end

function D=D2Dfun(E,En,ms,alpha)

hbar = 1.054e-34;
D = ms/(2*pi*hbar^2) * (1+2*alpha*(E-En)).* (E>En);

end

function yb = bound(x,b1,bu)

yb = min(max(x,b1),bu);

end

function [dy] = gradient_un (xx,yy)

pp =spline(xx,yy); dpp = fnder(pp); dy = fnval(dpp,xx);

end