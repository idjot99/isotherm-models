%clear all;
clf(1)
clf(2)
%%
% 
% A full-range moisture sorption model for cellulose-based materials
% yielding consistent net isosteric heat of sorption
%
% Johan Tryding (a,b), Henrik Askfelt (b), Marcus Alexandersson (b), and
% Matti Ristinmaa (a)
% (a) Division of Solid Mechanics, Lund University, SE-221 00 Lund, Sweden ; 
% (b) Tetra Pak, Ruben Rausingsgata, SE-221 87 Lund, Sweden
%
%%
% Specific gas constant
R=461.5e-3; % kJ/kg/kappa
%
%%
%
% choose isotherm type 1 to 7
iso=1;
%
isotherm{1}="Hendersson I";
isotherm{2}="Hendersson II";
isotherm{3}="Oswin I";
isotherm{4}="Oswin II";
isotherm{5}="LW I";
isotherm{6}="LW II";
isotherm{7}="GAB I";
%
markerColors={'r','r','g','g','b','b','k'};

%%
% Calibrated bleached data from Table 2 
%
if isotherm{iso}=="Hendersson I"   
    omegaREF=0.818; c=0.661; kappaInf=0.0841; theta0=402; n=3.64; 
elseif isotherm{iso}=="Hendersson II" 
    omegaREF=0.818; c=0.586; kappaInf=0.0803; theta0=403; n=3.55; b=1.34;    
elseif isotherm{iso}=="Oswin I"  
    omegaREF=0.818; c=0.467; kappaInf=0.0639; theta0=416; n=3.02; 
elseif isotherm{iso}=="Oswin II" 
    omegaREF=0.818; c=0.577; kappaInf=0.0730; theta0=406; n=3.42; b=0.704;      
elseif isotherm{iso}=="LW I" 
    omegaREF=0.818; c=0.389; kappaInf=0.0548; theta0=429; n=2.60;  
elseif isotherm{iso}=="LW II" 
    omegaREF=0.818; c=0.572; kappaInf=0.0967; theta0=409; n=3.29; b=0.562;       
elseif isotherm{iso}=="GAB I" 
    omegaREF=0.818; C0=61.43; qC=-149.8;         K0=1.109;    qK=-23.209;            
end     
%    
% Alexandersson [2]
    theta1=310; thetaREF=296.15; 
%  
% Isotherm temperatures
    T=[25 80];       % C
    theta=273.15+T;  % K 
%%
% Water acivity function and Net isosteric heat of sorption 
%
% Temparture dependency of kappatheta and omegaFSP
%
kappa   =@(theta,kappaInf,theta0,n) (kappaInf.*exp(1/(n+1).*(theta0./theta).^(n+1))); 
omegaFSP=@(theta,theta1,thetaREF,omegaREF)(omegaREF*(1+ (thetaREF-theta)/theta1));
%
HisoMax =@(omega,theta,theta0,n) (1./c.*R.*theta0.*(theta0./theta).^(n)); 
eta     =@(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
          (1+omega/omegaREF.*(theta/theta0).^(n+2).*(theta0/theta1)./ ...
          (1+(thetaREF-theta)/theta1-omega/omegaREF)./(1+(thetaREF-theta)/theta1));
%      
if isotherm{iso}=="Hendersson I"   
    xi =@(c)(1./c.*(1-c).^(1-c).*exp(-(1-c)));
    XI =@(omega,kappa,c,omegaFSP)((omega./xi(c)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP)(1-exp(-XI(omega,kappa,c,omegaFSP)));  
    w =@(aw,kappa,c,omegaFSP)(1./(1./omegaFSP+1./(xi(c).*kappa.*(-log(1-aw)).^c)));
    %
    h  =@(aw)(-(1-aw)./aw.*log(1-aw));
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF))).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            ); 
elseif isotherm{iso}=="Hendersson II"           
    xi =@(c,b)(1./c.*(2*b.*(1-c)./(b+1)).^(1-b.*c).*exp(2*b.*(c-1)./(b+1)).*(1-exp(2*b.*(c-1)./(b+1))).^(b-1));
    XI =@(omega,kappa,c,omegaFSP,b)((omega./xi(c,b)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP,b)((1-exp(-(XI(omega,kappa,c,omegaFSP,b)).^(1./b))).^(b)); 
    w =@(aw,kappa,c,omegaFSP,b)(1./(1./omegaFSP+1./(xi(c,b).*kappa.*((-log(1-aw.^(1./b))).^b).^c)));    
    %
    h  =@(aw,b) ((-(1-aw.^(1./b))./aw.^(1./b).*log(1-aw.^(1./b)))); 
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF,b) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF),b),b).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            );     
elseif isotherm{iso}=="Oswin I"  
    xi =@(c)(c./4.*(1./c+1).^(1+c).*(1./c-1).^(1-c));
    XI =@(omega,kappa,c,omegaFSP)((omega./xi(c)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP)(XI(omega,kappa,c,omegaFSP)./(1+XI(omega,kappa,c,omegaFSP)));     
    w  =@(aw,kappa,c,omegaFSP)(1./(1./omegaFSP+1./(xi(c).*kappa.*(aw./(1-aw)).^c)));
    %
    h  =@(aw)(1-aw);  
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF))).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            );    
elseif isotherm{iso}=="Oswin II" 
    xi =@(c,b)(1./c.*(b.*c+1).^(b.*c+1).*(b-b.*c).^(b-b.*c)./(b+1).^(b+1));
    XI =@(omega,kappa,c,omegaFSP,b)((omega./xi(c,b)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP,b)(((XI(omega,kappa,c,omegaFSP,b).^(1./b))./(1+XI(omega,kappa,c,omegaFSP,b).^(1./b))).^b);   
    w =@(aw,kappa,c,omegaFSP,b)(1./(1./omegaFSP+1./(xi(c,b).*kappa.*((aw.^(1./b)./(1-aw.^(1./b))).^b).^c)));
    %
    h  =@(aw,b) ((1-aw.^(1./b))); 
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF,b) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF),b),b).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            );          
elseif isotherm{iso}=="LW I"  
    p  =@(c)((sqrt(4*c+5)-(2*c+1))./(1+c)./2);
    xi =@(c)(1./c.*(p(c)).^(1-c)./(p(c)+1).*exp(-(1+c).*p(c)));
    XI =@(omega,kappa,c,omegaFSP)((omega./xi(c)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP)(1-exp(-W(XI(omega,kappa,c,omegaFSP)))); 
    w  =@(aw,kappa,c,omegaFSP)(1./(1./omegaFSP+1./(xi(c).*kappa.*(-log(1-aw)./(1-aw)).^c)));
    %
    h  =@(aw)(-(1-aw)./aw.*log(1-aw)./(1-log(1-aw))); 
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF))).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            );   
elseif isotherm{iso}=="LW II" 
    p  =@(c,b)( (sqrt(4*b.*(b + c) + 1) -2*(2*c - 1)*b - 1)./ (sqrt(4*b.*(b + c) + 1) +2*(2*c + 1)*b + 1));
    xi =@(c,b)(1./(b.*c).*(p(c,b)).^(1-b.*c)./(p(c,b)+1).*exp(-(1+b.*c).*p(c,b)));
    XI =@(omega,kappa,c,omegaFSP,b)((omega./xi(c,b)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP,b)((1-exp(-W(XI(omega,kappa,c,omegaFSP,b).^(1./b)))).^b); 
    w  =@(aw,kappa,c,omegaFSP,b)(1./(1./omegaFSP+1./(xi(c,b).*kappa.*(((-log(1-aw.^(1./b))./(1-aw.^(1./b))).^b).^c))));  
    %
    h  =@(aw,b) (-(1-aw.^(1./b))./aw.^(1./b).*log(1-aw.^(1./b))./(1-log(1-aw.^(1./b)))); 
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF,b) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF),b),b).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...     
            );          
elseif isotherm{iso}=="GAB I" 
    K  =@(theta,K0,qK) (K0.*exp(qK/R./theta)); 
    C  =@(theta,C0,qC) (C0.*exp(qC/R./theta));     
    aw =@(omega,omegaFSP,C,K)((sqrt(((C-2).*K.*(omega/omegaFSP-1)+(C-1).*K.^2-1).^2+4*(C-1).*K.^2* ...
         (omega/omegaFSP).^2)+(C-1).*K.^2+(C-2).*K.*(omega/omegaFSP-1)-1)./(2*(C-1).*K.^2.*omega/omegaFSP));
    w  =@(aw,omegaFSP,C,K)(omegaFSP.*aw.*(1-K)./(1-K.*aw).*(1-K+C.*K)./(1-K.*aw+C.*K.*aw));
    %
    y =@(omega,omegaFSP) (omega./omegaFSP);
    A =@(C,K,y) ((C-2).*K.*(y-1)+(C-1).*K.^2-1);
    B =@(C,K,y) (2.*(C-1).*K.^2.*y);
    dy = @(omega,omegaFSP,theta1,omegaREF) (omega./omegaFSP./omegaFSP.*omegaREF./theta1);
    rdA = @(theta,C,K,y,qC,qK,dy) (-C.*K.*qC.*(K+y-1) - K.*qK.*((C-2).*(y-1) + 2.*(C-1).*K) + (C-2).*K.*dy.*R.*theta.^2);
    rdB = @(theta,C,K,y,qC,qK,dy) (-2.*C.*qC.*K^2.*y - 4.*(C-1).*K.^2.*qK.*y + 2.*(C-1).*K.^2.*dy.*R.*theta.^2);
    %
    Hiso=@(theta,A,B,y,rdA,rdB,dy) ((A.*rdA + rdB.*y + B.*dy.*R.*theta.^2 )./((sqrt(A.^2+2*B.*y)+A).*sqrt(A.^2+2*B.*y)) + (rdA)./((sqrt(A.^2+2*B.*y)+A)) - (rdB)./(B));      
end    
%
%%
%
% Water activity function
%
figure(1)
set(gca,'FontSize',14)
hold on;
%
markerLine={'-',':'};
%
for i=1:length(theta)
    omega{i}=0.0:0.001:omegaFSP(theta(i),theta1,thetaREF,omegaREF);
    if isotherm{iso}=="GAB I"  
        plot(omega{i},aw(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF),C(theta(i),C0,qC),K(theta(i),K0,qK)), ...          
             markerLine{i},'color',markerColors{iso},'LineWidth',1)         
    elseif strfind(isotherm{iso}, "II")>0 %(aw,kappa,c,omegaFSP,b)
        plot(omega{i},aw(omega{i},kappa(theta(i),kappaInf,theta0,n),c, ...
             omegaFSP(theta(1),theta1,thetaREF,omegaREF),b),...
             markerLine{i},'color',markerColors{iso},'LineWidth',2);         
    else      
        plot(omega{i},aw(omega{i},kappa(theta(i),kappaInf,theta0,n),c, ...
             omegaFSP(theta(1),theta1,thetaREF,omegaREF)),...
             markerLine{i},'color',markerColors{iso},'LineWidth',1); 
    end
end
%
    axis([0 0.3 0 1])
    xlabel('moisture ratio, \omega [-]')    
    ylabel('water activity, a_{\omega} [-]')
    legend('T= 25C','T= 80C','Location','southeast')   
%
%%
% Net isosteric heat of sorption
%
figure(2)
set(gca,'FontSize',14)
hold on;
%
for i=1:length(theta)
    if isotherm{iso}=="GAB I"
       plot(omega{i}, ...
        Hiso(theta(i), ...
         A(C(theta(i),C0,qC),K(theta(i),K0,qK),y(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF))),  ...
         B(C(theta(i),C0,qC),K(theta(i),K0,qK),y(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF))), ...
         y(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF)), ...
         rdA(theta(i),C(theta(i),C0,qC),K(theta(i),K0,qK),y(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF)), ...
             qC,qK,dy(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF),theta1,omegaREF)), ...
         rdB(theta(i),C(theta(i),C0,qC),K(theta(i),K0,qK),y(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF)), ...
             qC,qK,dy(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF),theta1,omegaREF)), ...
         dy(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF),theta1,omegaREF)), ...
            markerLine{i},'color',markerColors{iso},'LineWidth',1);
    elseif strfind(isotherm{iso}, "II")>0 
       plot(omega{i},Hiso(omega{i},theta(i),c,omegaREF,kappaInf,theta0,n,theta1,thetaREF,b), ...
            markerLine{i},'color',markerColors{iso},'LineWidth',2);    
    else
       plot(omega{i},Hiso(omega{i},theta(i),c,omegaREF,kappaInf,theta0,n,theta1,thetaREF), ...
            markerLine{i},'color',markerColors{iso},'LineWidth',1);
    end
end
    axis([0 0.3 0 2000])
    xlabel('moisture ratio, \omega [-]')
    ylabel('net heat of sorption, \Delta{H}_{iso} [kJ/kg]')
    legend('T= 25C','T= 80C','Location','northeast')   
%
%%
% SUBROUTINES
%
%%
% Haley's method to calculate the Lambert W-function W(x) for x>=0
% 
% Reference:
% [1] https://blogs.mathworks.com/cleve/2013/09/02/the-lambert-w-function/
%
% [2] Corless, R., Gonnet, G., Hare, D., Jeffrey, D., Knuth, Donald 
% "On the Lambert W function". Advances in Computational Mathematics 5 (1996) 329-359 
%
function y = W(x)
% Lambert W-function W(x) for x>=0.
%
    y    = ones(size(x));
    yInf = Inf*y;
    while any(abs(y - yInf)./abs(y) > 1.e-6)
       yInf = y;
       y = y - (y.*exp(y)-x)./((exp(y).*(y+1) - (y+2).*(y.*exp(y)-x)./(2*y+2)));
    end
end

 




