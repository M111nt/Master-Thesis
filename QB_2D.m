clear all;
clc;
%close all;

hbar = 1.054e-34; % in J.S
q=1.602e-19;
m0=9.3e-31;
e0=8.85e-12;
k=1.38e-23;
T=300;

Lg1=[20 40 60 80 100]*1e-9;
% Lg1=30e-9;
Ef=0;
erox=18;
er_fun=@(x) 13.10E+00+1.100*x;% Permitivity of InGaAs changes with Indium composition
tox=5e-9;

Tr=100e-9./(100e-9+Lg1); % 100 nm is the mean free path



% need to change these parameters for study
c=1; % indium composition in the channel
ertw=er_fun(c);
tw=5e-9;
No_Of_Subs =3;% I think 3 are good

% capactiacnes

Cox = 2*erox*e0/(tox); % double gate

Cc = ertw*e0/(0.18*tw);% double gate

Cg = (1/Cox+1/Cc)^-1;

% screeening length

lambda = sqrt((1/8)*tw*(tw+4*(ertw*tox/erox)));% double gate

% Electronic paramerers of the InGaAs NS for different widths and
% composition at T=300K ... will try 10 K later
if (c==1)
    if(tw==3e-9)
        %   InAs H or tw=3e-9
        ms=[0.0522147135872787 0.149480903841096 0.250058268496825 0.328317010081031];
        E1= [1.06618321193555 1.85240716742656 2.71438169424118 3.66113243248657]-1.06618321193555;
        alpha=[1.4246 0.9707 0.7547 0.6604];
    elseif (tw==5e-9)
        %   InAs H or tw=5e-9
        ms=[0.0382479760730793 0.0893802672775955 0.150643027326401 0.212998450406476];
        E1 = [0.923531185354211 1.34413881179939 1.81770251077953 2.30595103757644]-0.923531185354211;
        alpha=[2.1007 1.3038 1.1343 1.082];
    elseif (tw==7e-9)
        %   InAs H or tw=7e-9
        ms=[0.0328055954275140 0.0656770507756168 0.107021969742892 0.151137277961047];
        E1=[0.867455199046153 1.14316992867206 1.46855834487474 1.80363928961151]-0.867455199046153;
        alpha=[2.25 1.6116 1.4721 1.0347];
    elseif (tw==10e-9)
        %   InAs H or tw=10e-9
        ms=[0.0289504471431968 0.0492522786928460 0.0760798861787328 0.105332642764020];
        E1=[0.828026518697966 1.00083104939028 1.21739531603394 1.44544617042573]-0.828026518697966;
        alpha=[2.7237 1.6710 1.5781 1.2947];
    end
    
elseif (c==0.8)
    if(tw==3e-9)
        %   InAs H or tw=3e-9
        ms=[0.0945137480181476 0.289986054218374 0.729533920995655 1.96057513135361];
        E1= [1.39290463567938 2.28816308287555 3.29129166473010 4.37507324988198]-1.39290463567938;
        alpha=[0.6738 0.3107 0.25 0.3792];
    elseif (tw==5e-9)
        %   InAs H or tw=5e-9
        ms=[0.060213618759712 0.139443941023532 0.270182899685620 0.484069629523636];
        E1= [1.13756791112703 1.62073087098082 2.16938307262117 2.75406563347775]-1.13756791112703;
        alpha=[1.0761 0.6767 0.3245 0.1622];
        
    elseif (tw==7e-9)
        %   InAs H or tw=7e-9
        ms=[0.048030507724735 0.094246596212338 0.163593633999310 0.261692555296921];
        E1= [1.03909787721269 1.35513490341652 1.72422664106341 2.11883971709856]-1.03909787721269;
        alpha=[1.3259 0.7802 0.5118 0.3988];
        
    elseif (tw==10e-9)
        %   InAs H or tw=10e-9
        ms=[0.040126642920527 0.066664986812252 0.104749192838261 0.154285661938180];
        E1= [0.973424777045629 1.16917475810491 1.40866767049677 1.66861562246663]-0.973424777045629;
        alpha=[1.5428 0.9481 0.671 0.5988];
        
    end
    
elseif(c==0.53)
    if(tw==3e-9)
        %   InAs H or tw=3e-9
        ms=[0.113955623746157 0.319238517534564 0.799727739994583 2.431623121318251];
        E1= [1.66168716372682 2.58063493475950 3.60234584230419 4.70022967513303]-1.66168716372682;
        alpha=[0.3895 0.1807 0.18 0];
        
    elseif (tw==5e-9)
        %   InAs H or tw=5e-9
        ms=[0.072066496998488 0.150036570588999 0.277469522660927 0.485705957220856];
        E1= [1.36753197367436 1.85646032949688 2.41434913673922 3.00948730580227]-1.36753197367436;
        alpha=[0.6523 0.3763 0.2225 0.1127];
        
    elseif (tw==7e-9)
        %   InAs H or tw=7e-9
        ms=[0.058344653182972 0.102333945349827 0.168179097859467 0.260751743574167];
        E1= [1.25985746112165 1.57210970886656 1.94363077776011 2.34375143647637]-1.25985746112165;
        alpha=[0.7701 0.5249 0.338 0.267];
        
    elseif (tw==10e-9)
        %   InAs H or tw=10e-9
        ms=[0.050081263201672 0.074285955569948 0.109539967708303 0.155494366178171];
        E1= [1.19225117124206 1.37815146778978 1.61386442407014 1.87415198124336]-1.19225117124206;
        alpha=[0.9738 0.6165 0.4278 0.3848];
        
    end
    
end


%scatter([1 2 3 4],E1);% just plotting the energy levels

% calcualting the carrier density

etopvec = linspace(-0.5,1);

for ii=1:length(etopvec)
    n0_dummy=0;
    for aa=1:No_Of_Subs
        n0 = @(E) q*D2Dfun(E,etopvec(ii)+E1(aa),ms(aa)*m0,alpha(aa)).*fd(E,Ef,T);
        n0_dummy= n0_dummy+integral(n0,etopvec(ii)+E1(aa),etopvec(ii)+E1(aa)+1);
    end
    n0_tot(ii)=n0_dummy;
end
% figure(111);plot(etopvec,n0_tot);

pp =spline(etopvec,n0_tot);

% currenet-voltage characteristics

nbr_steps=31;

Vgvec=linspace(-0.1,1.5,nbr_steps);

Vdvec = 0:0.05:1;

% Vdvec=0:0.01:1;

VT=0.5;

for kk=1:1%length(Lg1)
    
    Lg=Lg1(kk);
    
    for jj=1:length(Vdvec)
        
        Vd=Vdvec(jj);
        
        for ii=1:length(Vgvec)
            Qs_dummy=0;
            itot_dummy=0;
            itot_dummy1=0;
            for aa=1:No_Of_Subs % number of subbands
                
                Vg = Vgvec(ii);
                
                
                xtop = xtopapprox(0.2,Vd,lambda,Lg);
                
                etop = etopapprox(0.2,Vd,lambda,Lg);
                
                ag =(Ecfun(xtop,0.2,Vd,lambda,Lg) - Ecfun(xtop,0.2+10e-3,Vd,lambda,Lg) )/10e-3;
                
%                 ag=1;
                
                Csigma = Cg/ag;
                
                ad = (Ecfun(xtop,0.2,Vd,lambda,Lg) - Ecfun(xtop,0.2,Vd+10e-3,lambda,Lg) )/10e-3;
                
%                 ad=0;
                
                Cdrain = ad*Csigma;
                
                Efd = Ef-Vd;
                
                fzerofunc = @(x) -ag*(Vg-VT) - ad*Vd + q*((2-Tr(kk))*fnval(pp,x+E1(aa))+Tr(kk)*fnval(pp,x+Vd+E1(aa)))/Csigma - x;
                
                etop = fzero(fzerofunc,Vg/2);   
                
                % Current calcualtions
                
                i0fun=@(E) sqrt(q)*2*q^2/hbar * Tr(kk) .* Mfun(E,etop+E1(aa),ms(aa)*m0,alpha(aa)) .* (fd(E,Ef,T)-fd(E,Efd,T)) .* (E>(etop+E1(aa)));
                
                itot_dummy=itot_dummy+integral(i0fun,etop+E1(aa),etop+E1(aa)+1); % summing the current from all subbands at a particular Vg and Vd
                
                
                              
            end
            
            itot(jj,ii) = itot_dummy;

        end

    end
    
    
    %% Transfer characteristics
     
    figure(1);
    hold on;
    plot(Vgvec-VT,itot(11,:)*1e-3 );hold on;xlim([-0.5 1.2]);xlabel('V_{GS} (V)');ylabel('I_{DS} (mA/\mu m)'); 
    
    figure(2);hold on;
    plot(Vgvec-VT,[0,diff(itot(11,:)*1e-3)./diff(Vgvec)]);xlim([-0.1 1.2]);
    xlabel('V_{GS} (V)');ylabel('g_{m} (mS/\mu m)');
    gm1=[0,diff(itot(11,:)*1e-3)./diff(Vgvec)];
    gm(kk)=gm1(35); % consider this ar VGS-VT=0.3 V
    
    Inv_ss = [0;diff(log10(itot(11,:)))./diff(Vgvec)];
    ss1 =(1./Inv_ss);
    ss(kk)=min(ss1(ss1>0));
    
    %   Output Characteristics
        figure(3);hold on;
    
        for ii=1:length(itot(1,:))
            plot(Vdvec,itot(:,ii)*1e-3,'b','LineWidth',2);hold on; % plot(Vd_ext(:,end),itot(:,ii)*1e-3,'r','LineWidth',2);
            xlabel('V_{DS} (V)');ylabel('I_{DS} (mA/\mu m)');
        end
    
    %     figure(4);plot(Vdvec,[0;diff(itot(:,11)*1e-3)./diff(Vdvec')]);
    
    %     Ron (ohm-um)
    
    %     gd1(:,kk)=[0;diff(itot(:,7))./diff(Vdvec')];
    %
    %     % figure(1);hold on; plot(Vdvec,gd1(:,kk));
    %
    %     Ron(kk)=1/max(gd1(:,kk))*1e6; % in ohm-um
    %
    %     % gd (mS/um)
    
    %     gd(kk)=gd1(15,kk);
    %
end

%% figure(100); scatter(Lg1/1e-9,ss);ylabel('ss (mV/dec)');xlabel('L_g(nm)');
% figure(101); hold on;scatter(Lg1/1e-9,gm);ylabel('g_m (mS/\mum)');xlabel('L_g(nm)');
% % figure(102); scatter(Lg1,DIBL);
% figure(103); scatter(Lg1/1e-9,gd);ylabel('g_d (mS/\mum)');xlabel('L_g(nm)');
% figure(104); scatter(Lg1/1e-9,Ron);ylabel('On-resistance (\Omega-\mum)');xlabel('L_g(nm)');

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
D = ms/(2*pi*hbar^2) * (1+2*alpha*(E-En)).*(E>En);

end

function yb = bound(x,b1,bu)

yb = min(max(x,b1),bu);

end

function [dy] = gradient_un (xx,yy)

pp =spline(xx,yy); dpp = fnder(pp); dy = fnval(dpp,xx);

end