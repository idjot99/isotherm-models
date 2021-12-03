clear all;
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
% choose isotherm type 1 to 
iso=4;
%
isotherm{1}="Hendersson";
isotherm{2}="Hendersson II";
isotherm{3}="Oswin";
isotherm{4}="Oswin II";
isotherm{5}="LW";
isotherm{6}="LW II";
isotherm{7}="GAB";
%
markerColors={'r','r','g','g','b','b','k'};

%%
% Calibrated bleached data from Table 2 
%
if isotherm{iso}=="Hendersson"   
    omegaFSP=0.818; c=0.653; kappaInf=0.0771; theta0=412; n=3.36; 
elseif isotherm{iso}=="Hendersson II" 
    omegaFSP=0.818; c=0.429; kappaInf=0.0669; theta0=419; n=3.08; b=2.51;    
elseif isotherm{iso}=="Oswin"  
    omegaFSP=0.818; c=0.462; kappaInf=0.0585; theta0=427; n=2.80; 
elseif isotherm{iso}=="Oswin II" 
    omegaFSP=0.818; c=0.604; kappaInf=0.0698; theta0=413; n=3.28; b=0.632;      
elseif isotherm{iso}=="LW" 
    omegaFSP=0.818; c=0.385; kappaInf=0.0498; theta0=442; n=2.41;  
elseif isotherm{iso}=="LW II" 
    omegaFSP=0.818; c=0.630; kappaInf=0.1129; theta0=412; n=3.30; b=0.461;       
elseif isotherm{iso}=="GAB" 
    omegaFSP=0.818; C0=59.8897; qC=-145.7556; K0=1.1081;  qK=-23.0159;            
end     
%    
% Isotherm temperatures
    T=[25 80];       % C
    theta=273.15+T;  % K 
%%
% Water acivity function and Net isosteric heat of sorption 
%
% Temparture dependency of kappatheta 
kappa   =@(theta,kappaInf,theta0,n) (kappaInf.*exp(1/(n+1).*(theta0./theta).^(n+1))); 
%      
if isotherm{iso}=="Hendersson"   
    xi =@(c)(1./c.*(1-c).^(1-c).*exp(-(1-c)));
    XI =@(omega,kappa,c,omegaFSP)((omega./xi(c)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP)(1-exp(-XI(omega,kappa,c,omegaFSP)));     
    w =@(aw,kappa,c,omegaFSP)(1./(1./omegaFSP+1./(xi(c).*kappa.*(-log(1-aw)).^c)));
    HisoMax =@(omega,theta,theta0,n) (1./c.*R.*theta0.*(theta0./theta).^(n));     
    h  =@(aw)(-(1-aw)./aw.*log(1-aw));    
    Hiso =@(omega,theta,c,kappaInf,theta0,n,omegaFSP) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP)) ...
            ); 
elseif isotherm{iso}=="Hendersson II"           
    xi =@(c,b)(1./c.*(2*b.*(1-c)./(b+1)).^(1-b.*c).*exp(2*b.*(c-1)./(b+1)).*(1-exp(2*b.*(c-1)./(b+1))).^(b-1));
    XI =@(omega,kappa,c,omegaFSP,b)((omega./xi(c,b)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP,b)((1-exp(-(XI(omega,kappa,c,omegaFSP,b)).^(1./b))).^(b)); 
    w =@(aw,kappa,c,omegaFSP,b)(1./(1./omegaFSP+1./(xi(c,b).*kappa.*((-log(1-aw.^(1./b))).^b).^c)));
    %
    HisoMax =@(omega,theta,theta0,n) (1./c.*R.*theta0.*(theta0./theta).^(n));     
    h  =@(aw,b) ((-(1-aw.^(1./b))./aw.^(1./b).*log(1-aw.^(1./b))));     
    Hiso =@(omega,theta,c,kappaInf,theta0,n,omegaFSP,b) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP,b),b) ...
            );     
elseif isotherm{iso}=="Oswin"  
    xi =@(c)(c./4.*(1./c+1).^(1+c).*(1./c-1).^(1-c));
    XI =@(omega,kappa,c,omegaFSP)((omega./xi(c)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP)(XI(omega,kappa,c,omegaFSP)./(1+XI(omega,kappa,c,omegaFSP)));     
    w  =@(aw,kappa,c,omegaFSP)(1./(1./omegaFSP+1./(xi(c).*kappa.*(aw./(1-aw)).^c)));
    %
    HisoMax =@(omega,theta,theta0,n) (1./c.*R.*theta0.*(theta0./theta).^(n));     
    h  =@(aw)(1-aw);      
    Hiso =@(omega,theta,c,kappaInf,theta0,n,omegaFSP) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP)) ...
            );    
elseif isotherm{iso}=="Oswin II" 
    xi =@(c,b)(1./c.*(b.*c+1).^(b.*c+1).*(b-b.*c).^(b-b.*c)./(b+1).^(b+1));
    XI =@(omega,kappa,c,omegaFSP,b)((omega./xi(c,b)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP,b)(((XI(omega,kappa,c,omegaFSP,b).^(1./b))./(1+XI(omega,kappa,c,omegaFSP,b).^(1./b))).^b); 
    w =@(aw,kappa,c,omegaFSP,b)(1./(1./omegaFSP+1./(xi(c,b).*kappa.*((aw.^(1./b)./(1-aw.^(1./b))).^b).^c)));
    %
    HisoMax =@(omega,theta,theta0,n) (1./c.*R.*theta0.*(theta0./theta).^(n));     
    h  =@(aw,b) ((1-aw.^(1./b)));     
    Hiso =@(omega,theta,c,kappaInf,theta0,n,omegaFSP,b) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP,b),b) ...
            );          
elseif isotherm{iso}=="LW"  
    p  =@(c)((sqrt(4*c+5)-(2*c+1))./(1+c)./2);
    xi =@(c)(1./c.*(p(c)).^(1-c)./(p(c)+1).*exp(-(1+c).*p(c)));
    XI =@(omega,kappa,c,omegaFSP)((omega./xi(c)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP)(1-exp(-W(XI(omega,kappa,c,omegaFSP)))); 
    w  =@(aw,kappa,c,omegaFSP)(1./(1./omegaFSP+1./(xi(c).*kappa.*(-log(1-aw)./(1-aw)).^c)));
    %
    HisoMax =@(omega,theta,theta0,n) (1./c.*R.*theta0.*(theta0./theta).^(n));     
    h  =@(aw)(-(1-aw)./aw.*log(1-aw)./(1-log(1-aw))); 
    Hiso =@(omega,theta,c,kappaInf,theta0,n,omegaFSP) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP)) ...
            );   
elseif isotherm{iso}=="LW II" 
    p  =@(c,b)( (sqrt(4*b.*(b + c) + 1) -2*(2*c - 1)*b - 1)./ (sqrt(4*b.*(b + c) + 1) +2*(2*c + 1)*b + 1));
    xi =@(c,b)(1./(b.*c).*(p(c,b)).^(1-b.*c)./(p(c,b)+1).*exp(-(1+b.*c).*p(c,b)));
    XI =@(omega,kappa,c,omegaFSP,b)((omega./xi(c,b)./kappa./(1-omega./omegaFSP)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP,b)((1-exp(-W(XI(omega,kappa,c,omegaFSP,b).^(1./b)))).^b); 
    w  =@(aw,kappa,c,omegaFSP,b)(1./(1./omegaFSP+1./(xi(c,b).*kappa.*(((-log(1-aw.^(1./b))./(1-aw.^(1./b))).^b).^c))));
    %
    HisoMax =@(omega,theta,theta0,n) (1./c.*R.*theta0.*(theta0./theta).^(n));     
    h  =@(aw,b) (-(1-aw.^(1./b))./aw.^(1./b).*log(1-aw.^(1./b))./(1-log(1-aw.^(1./b)))); 
    Hiso =@(omega,theta,c,kappaInf,theta0,n,omegaFSP,b) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP,b),b) ...     
            );          
elseif isotherm{iso}=="GAB" 
    K  =@(theta,K0,qK) (K0.*exp(qK/R./theta)); 
    C  =@(theta,C0,qC) (C0.*exp(qC/R./theta));     
    aw =@(omega,omegaFSP,C,K)((sqrt(((C-2).*K.*(omega/omegaFSP-1)+(C-1).*K.^2-1).^2+4*(C-1).*K.^2* ...
         (omega/omegaFSP).^2)+(C-1).*K.^2+(C-2).*K.*(omega/omegaFSP-1)-1)./(2*(C-1).*K.^2.*omega/omegaFSP));
    w  =@(aw,omegaFSP,C,K)(omegaFSP.*aw.*(1-K)./(1-K.*aw).*(1-K+C.*K)./(1-K.*aw+C.*K.*aw));
    %
    y =@(omega,omegaFSP) (omega./omegaFSP);
    A =@(C,K,y) ((C-2).*K.*(y-1)+(C-1).*K.^2-1);
    B =@(C,K,y) (2.*(C-1).*K.^2.*y);
    rdA = @(theta,C,K,y,qC,qK) (-C.*K.*qC.*(K+y-1) - K.*qK.*((C-2).*(y-1) + 2.*(C-1).*K));
    rdB = @(theta,C,K,y,qC,qK) (-2.*C.*qC.*K^2.*y - 4.*(C-1).*K.^2.*qK.*y);
    %
    Hiso=@(theta,A,B,y,rdA,rdB) ((A.*rdA + rdB.*y)./((sqrt(A.^2+2*B.*y)+A).*sqrt(A.^2+2*B.*y)) + (rdA)./((sqrt(A.^2+2*B.*y)+A)) - (rdB)./(B));      
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
    omega{i}=0.0:0.0001:omegaFSP;
    if isotherm{iso}=="GAB"  
        plot(omega{i},aw(omega{i},omegaFSP,C(theta(i),C0,qC),K(theta(i),K0,qK)), ...          
             markerLine{i},'color',markerColors{iso},'LineWidth',2)         
    elseif strfind(isotherm{iso}, "II")>0 %(aw,kappa,c,omegaFSP,b)
        plot(omega{i},aw(omega{i},kappa(theta(i),kappaInf,theta0,n),c, ...
             omegaFSP,b),...
             markerLine{i},'color',markerColors{iso},'LineWidth',3);         
    else      
        plot(omega{i},aw(omega{i},kappa(theta(i),kappaInf,theta0,n),c, ...
             omegaFSP),...
             markerLine{i},'color',markerColors{iso},'LineWidth',2); 
    end
end
%
    axis([0 omegaFSP 0 1])
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
    if isotherm{iso}=="GAB"
       plot(omega{i}, ...
        Hiso(theta(i), ...
         A(C(theta(i),C0,qC),K(theta(i),K0,qK),y(omega{i},omegaFSP)),  ...
         B(C(theta(i),C0,qC),K(theta(i),K0,qK),y(omega{i},omegaFSP)), ...
         y(omega{i},omegaFSP), ...
         rdA(theta(i),C(theta(i),C0,qC),K(theta(i),K0,qK),y(omega{i},omegaFSP), ...
             qC,qK), ...
         rdB(theta(i),C(theta(i),C0,qC),K(theta(i),K0,qK),y(omega{i},omegaFSP), ...
             qC,qK)), ...
            markerLine{i},'color',markerColors{iso},'LineWidth',2);
    elseif strfind(isotherm{iso}, "II")>0 
       plot(omega{i},Hiso(omega{i},theta(i),c,kappaInf,theta0,n,omegaFSP,b), ...
            markerLine{i},'color',markerColors{iso},'LineWidth',3);    
    else
       plot(omega{i},Hiso(omega{i},theta(i),c,kappaInf,theta0,n,omegaFSP), ...
            markerLine{i},'color',markerColors{iso},'LineWidth',2);
    end
end
    axis([0 omegaFSP 0 1700])
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

 




