
%-------------------------------
% For 1D ballistic
%-------------------------------

clear;
close all;
global hbar q k m0 e0 T;

hbar = 6.626e-34; %planck's constant
q = 1.602e-19; %electron charge
k = physconst('Boltzmann'); %boltzmann's constant
m0 = 9.109e-31; %free electron mass
e0 = 8.854e-12; %vacuum permittivity
T = 300; %tempreture


tw=13e-9; %well thickness
tox=5e-9; %oxide thickness
Lg=80e-9; %gate length
Efs=0; %source fermi level



%bandgap, effective mass, and permitivitty.
Eg = 0.5; 
ms = 0.023*m0;
ertw=13; %well er
erox=18; %oxide er

%Nanosheet height and width
W=30e-9; 
H=5e-9;

% nonparabolicity factor
alpha = 1/Eg*(1-ms/m0)^2; 



%-----------------------------------------------------------
vg_loop = -0.3:0.01:2.0; 
vd_loop = 0.5:0.1:1.0;
%vs = 0;
vt = 0;

%------------------------------------------------------------
%calculte gate capacitance in ballistic model 
Cox=erox*e0/tox;
Cq=q^2*0.03*m0/(pi*hbar^2);
Cg = Cox * Cq / (Cox + Cq);

%Csigma = 6.0990e-4; %calculated by previous matlab code
Csigma = Cg;
Cd = 0;
Cs = 0;


% subbands values
No_subs = 100;

%define the size of each matrix, faster computing speed
gamma = zeros(10);
m_nm = zeros(10);
alpha_nm = zeros(10);
Enm =zeros(10);
store = zeros(100,4);

ctrl = 1;
%% calculate the energy level
for nn =1:sqrt(No_subs)
    
    for mm=1:sqrt(No_subs)
        
        %The rest part is too small compared to 1. 
        %The gamma cannot get a correct value.
        %We should calculate the Enm_p first.
        %gamma(nn,mm) = sqrt(1 + (2*alpha*hbar^2*pi^2/ms)*((nn/W)^2+(mm/H)^2) );
        
        %Enm_potential 
        Enm_p = (1/q)*((hbar/(2*pi))^2*pi^2/(2*ms))*((nn/W)^2+(mm/H)^2);
        
        %Calculate gamma by calculting Enm_p
        gamma(nn,mm) = sqrt(1+4*alpha*Enm_p);
        
        %effective sub-band mass
        m_nm(nn,mm) = ms*gamma(nn,mm)/m0;
        
        %an effective sub-band non-parabolicity factor
        alpha_nm(nn,mm) = alpha/gamma(nn,mm);
        
        %the energy of square well with infinte height
        Enm (nn,mm) = (gamma(nn,mm)-1)/(2*alpha);
        
        %store the parameters
        store(ctrl,:) = [Enm(nn,mm) alpha_nm(nn,mm) m_nm(nn,mm) gamma(nn,mm)];
        ctrl = ctrl + 1;
        
    end
end

x = linspace(1,10,10);
y = linspace(1,10,10);
[X,Y] = meshgrid(x,y);

% figure(1);hold on;
% contour(X,Y,Enm);
% xlabel('Index n of width in x direction');
% ylabel('Index m of width in y direction');
% title('The energy of square well with infinte height in 1D ballistic');

store_new = sortrows(store, 1);
Enm_new = store_new(:,1);
alpha_nm_new = store_new(:,2);
m_nm_new = store_new(:,3);
%gamma_new = store_new(:,4);

%Efd= Efs-vd;

%E0_new = zeros(1:131);
itot = zeros(1,6);
%I_total = zeros(1:131);

%% integral the total electron plot the Vgs VS Ids

vg_loop = -0.3:0.01:2.0; 
vd_loop = 0.5:0.1:1.0;

for jj = 1:length(vd_loop)
    
    vd = vd_loop(jj);
    Efd= Efs-vd;

    for ii = 1:length(vg_loop)
    
        vg = vg_loop(ii);
    
        for kk =1:6%number of subbands
        
            n0 = @(E,E0)  sqrt(q)*D1D(E, Enm_new(kk)+E0, m_nm_new(kk)*m0, alpha_nm_new(kk)).*(fd(E,Efs,T) + fd(E,Efd,T));
        
            n0_tot = @(E0) ( integral(@(E) n0(E,E0), E0+Enm_new(kk), 2) );
        
            %calculate Uscf
            E0_new(jj,ii) = fzero(@(E0) - E0 - (vg-vt) + q*n0_tot(E0)/(Csigma), vg/2);
        
            %h = 2*pi*hbar
            itot(kk) = 2*q*k*T/hbar * (log(1+exp((Efs-Enm_new(kk)-E0_new(jj,ii))/(k*T/q))) - log(1+exp((Efd-Enm_new(kk)-E0_new(jj,ii))/(k*T/q))) ) ;
                
        end
    
        I_total(jj,ii) = sum(itot);
    
    end 
    % plot
    
    %voltage characteristic
    figure(2); hold on; 
    plot(vg_loop, I_total(jj,:)*1e-3/(2*(W+H)),'LineWidth',2);
    xlabel('V_{GS} (V)');
    ylabel('I_{DS}(mA/\mum)');
    legend('V_{DS} = 0.5', 'V_{DS} = 0.6', 'V_{DS} = 0.7','V_{DS} = 0.8','V_{DS} = 0.9','V_{DS} = 1.0');
    title('The 1D ballistic simulation for V_{GS} and I_{DS}');

    %gm
    x=vg_loop;
    y=I_total(jj,:)*1e-3/(2*(W+H));
    x1=(length(vg_loop)-1)/2;
    pp_DC = csaps(x,y, [1,ones(1,x1),repmat(5,1,x1)], [], ...
            [ones(1,x1), repmat(5,1,x1), 0]);
        %fnplt(pp_DC,'b');        

    figure(3);hold on;

    fnplt(fnder(pp_DC));
    %fnplt(fnder(pp_DC),'b');
    legend('V_{DS} = 0.5', 'V_{DS} = 0.6', 'V_{DS} = 0.7','V_{DS} = 0.8','V_{DS} = 0.9','V_{DS} = 1.0');
    ylabel('g_{m}(mS/\mum)');xlabel('V_{GS} (V)');
    title('The 1D ballistic simulation for V_{GS} and g_{m}');
    
    
end


%% plot the Vgs VS lg(Id)
kk = 1;
I_new = log10(I_total);
figure(11); hold on; 
title('The 1D ballistic simulation for V_{GS} and log_{10}(I_{DS})');
xlabel('V_{GS}(V)');
ylabel('log_{10}(I_{DS})');
   
for jj = 1: length(vd_loop)
        
    plot(vg_loop,I_new(kk,:),'LineWidth',1);
        
    kk = kk + 1;
    
    legend('V_{DS} = 0.5', 'V_{DS} = 0.6', 'V_{DS} = 0.7','V_{DS} = 0.8','V_{DS} = 0.9','V_{DS} = 1.0');

end

%% Ron

% Ids=I_total;
% Ron=vd./I_total;
% figure(5); hold on;
% plot(Ids,Ron);
% xlabel('I_{DS}');
% ylabel('R_{on}');
% legend('V_{DS}=0.5V');











%% temperature study
vg_loop = -0.3:0.01:2.0;
vd = 0.5;
Efd= Efs-vd;
T_loop = 100:200:500;

for tt = 1:length(T_loop)
    
    T = T_loop(tt);
    
    for ii = 1:length(vg_loop)
        
        vg = vg_loop(ii);
    
        for kk =1:6%number of subbands
            
            %1D ballistic
            n0_t = @(E,E0)  sqrt(q)*D1D(E, Enm_new(kk)+E0, m_nm_new(kk)*m0, alpha_nm_new(kk)).*(fd(E,Efs,T) + fd(E,Efd,T));
        
            n0_tot_t = @(E0) ( integral(@(E) n0_t(E,E0), E0+Enm_new(kk), 2) );
        
            %calculate Uscf
            E0_new_t(1,ii) = fzero(@(E0) - E0 - (vg-vt) + q*n0_tot_t(E0)/(Csigma), vg/2);
        
            %h = 2*pi*hbar
            itot_t(kk) = 2*q*k*T/hbar * (log(1+exp((Efs-Enm_new(kk)-E0_new_t(1,ii))/(k*T/q))) - log(1+exp((Efd-Enm_new(kk)-E0_new(1,ii))/(k*T/q))) ) ;
                
        end
    
        I_total_t(1,ii) = sum(itot_t);
    end 
    
    %y(tt) = I_total_t(1,:)*1e-3/(2*(W+H));
    
    figure(4); hold on; 
    plot(vg_loop, I_total_t(1,:)*1e-3/(2*(W+H)));
    
    xlabel('V_{GS} (V)');
    ylabel('I_{DS}(mA/\mum)');
    title('The 1D ballistic simulation for voltage transfer characteristic in different temperature');
    %legend('100K', '200K', '300K', '400K', '500K');
    legend('100K', '300K', '500K');

end 


%% integral the total electron plot the Vds VS Ids

vg_loop = 0.5:0.1:1.0;
vd_loop = -0.3:0.01:2.0;

for jj = 1:length(vg_loop)
    
    vg = vg_loop(jj);
    

    for ii = 1:length(vd_loop)
    
        vd = vd_loop(ii);
        
        Efd= Efs-vd;
    
        for kk =1:6%number of subbands
        
            n0 = @(E,E0)  sqrt(q)*D1D(E, Enm_new(kk)+E0, m_nm_new(kk)*m0, alpha_nm_new(kk)).*(fd(E,Efs,T) + fd(E,Efd,T));
        
            n0_tot = @(E0) ( integral(@(E) n0(E,E0), E0+Enm_new(kk), 2) );
        
            %calculate Uscf
            E0_new(1,ii) = fzero(@(E0) - E0 - (vg-vt) + q*n0_tot(E0)/(Csigma), vg/2);
        
            %h = 2*pi*hbar
            itot(kk) = 2*q*k*T/hbar * (log(1+exp((Efs-Enm_new(kk)-E0_new(1,ii))/(k*T/q))) - log(1+exp((Efd-Enm_new(kk)-E0_new(1,ii))/(k*T/q))) ) ;
                
        end
    
        I_total(jj,ii) = sum(itot);
    
    end 
    % plot
    
    %voltage characteristic
    figure(12); hold on; 
    plot(vd_loop, I_total(jj,:)*1e-3/(2*(W+H)),'LineWidth',2);
    xlabel('V_{DS} (V)');
    ylabel('I_{DS}(mA/\mum)');
    legend('V_{GS} = 0.5', 'V_{GS} = 0.6', 'V_{GS} = 0.7','V_{GS} = 0.8','V_{GS} = 0.9','V_{GS} = 1.0');
    title('The 1D ballistic simulation for V_{DS} and I_{DS}');

  
end







%% ---------------------------------------------------------------------
function D = D1D(E,En,ms,alpha) % 1D density of states

hbar = 6.626e-34;

D = sqrt(2*ms) .* (1 + 2*alpha*(E-En) ./ (pi*hbar*sqrt((E-En).*(1+alpha*(E-En))))).*(E>En) ;

end

function f = fd(E,Ef,T) % fermi dirac function

k = physconst('Boltzmann');
q = 1.602e-19;
f = 1./(1+exp((E-Ef)/(k*T/q)));

end


