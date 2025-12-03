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
R=461.5e-3; % kJ/kg/K
%
%%
%
% choose isotherm type 1 to 7
iso=6;
%
isotherm{1}="Henderson I";
isotherm{2}="Henderson II";
isotherm{3}="Oswin I";
isotherm{4}="Oswin II";
isotherm{5}="LW I";
isotherm{6}="LW II";
isotherm{7}="GAB I";
isotherm{8}="GAB";
%
markerColors={'r','r','g','g','b','b','k','k'};
%
plotAlsoOmega="Yes"; % "Yes" or "No"
%
%%
% Models calibrated to bleached fiber data from Leuk et al. Drying
% Technology 34(5):563-573, 2016.
%
if isotherm{iso}=="Henderson I"
    p = [0.5786    0.0629  377.8934    4.7100];
    omegaREF=2.49; c=p(1); kappaInf=p(2); theta0=p(3); n=p(4);
elseif isotherm{iso}=="Henderson II"
    p = [0.4979    0.0519  375.3296    4.7748    1.4827];
    omegaREF=2.49; c=p(1); kappaInf=p(2); theta0=p(3); n=p(4); b=p(5);    
elseif isotherm{iso}=="Oswin I"  
    p = [0.4306    0.0541  364.0613    6.0171];
    omegaREF=2.49; c=p(1); kappaInf=p(2); theta0=p(3); n=p(4); 
elseif isotherm{iso}=="Oswin II" 
    p = [0.5837    0.0694  373.8213    5.0128    0.5539];
    omegaREF=2.49; c=p(1); kappaInf=p(2); theta0=p(3); n=p(4); b=p(5);       
elseif isotherm{iso}=="LW I" 
    p = [0.3644    0.0491  359.4227    6.7272];
    omegaREF=2.49; c=p(1); kappaInf=p(2); theta0=p(3); n=p(4); 
elseif isotherm{iso}=="LW II" 
    p = [0.5125    0.0617  370.8969    5.2068    0.6089];
    omegaREF=2.49; c=p(1); kappaInf=p(2); theta0=p(3); n=p(4); b=p(5);         
elseif isotherm{iso}=="GAB I" 
    p = [34.9718  -86.0959    1.0098   -4.1838];
    omegaREF=2.49;C0=p(1); qC=p(2);      K0=p(3);  qK=p(4);        
elseif isotherm{iso}=="GAB" 
    p = [0.4273  468.9266    0.1011  297.7358    0.0578];
                  C0=p(1); qC=p(2);      K0=p(3);  qK=p(4); mGAB=p(5);     
end     
%    
% Alexandersson [9]
    theta1=310; thetaREF=296.15; 
%  
% Isotherm temperatures
    T=[25 80];       % C
    theta=273.15+T;  % K 
% Isotherm dimensionless constant m
    m=3;
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
          (1+m.*(omega/omegaFSP(theta,theta1,thetaREF,omegaREF)).^m./ ...
          (1-(omega/omegaFSP(theta,theta1,thetaREF,omegaREF)).^m).* ...
          (omegaREF./omegaFSP(theta,theta1,thetaREF,omegaREF)).*(theta/theta0).^(n+2).*(theta0/theta1));      
%      
if isotherm{iso}=="Henderson I"
    xi =@(omega,kappa,c,omegaFSP)((omega./kappa./(1-(omega./omegaFSP).^m)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP)(1-exp(-xi(omega,kappa,c,omegaFSP)));
    %
    v =@(aw,kappa,c,omegaFSP) (omegaFSP./(3.*kappa.*(-log(1-aw)).^c));
    %
    h  =@(aw)(-(1-aw)./aw.*log(1-aw));
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF) ...
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF))).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            );
elseif isotherm{iso}=="Henderson II"   
    xi =@(omega,kappa,c,omegaFSP,b)((omega./kappa./(1-(omega./omegaFSP).^m)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP,b)((1-exp(-(xi(omega,kappa,c,omegaFSP,b)).^(1./b))).^(b)); 
    % 
    v =@(aw,kappa,c,omegaFSP,b) (omegaFSP./(3.*kappa.*((-log(1-aw.^(1./b))).^b).^c));       
    %
    h  =@(aw,b) ((-(1-aw.^(1./b))./aw.^(1./b).*log(1-aw.^(1./b)))); 
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF,b) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF),b),b).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            );     
elseif isotherm{iso}=="Oswin I"  
    xi =@(omega,kappa,c,omegaFSP)((omega./kappa./(1-(omega./omegaFSP).^m)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP)(xi(omega,kappa,c,omegaFSP)./(1+xi(omega,kappa,c,omegaFSP)));     
    %
    v =@(aw,kappa,c,omegaFSP) (omegaFSP./(3.*kappa.*(aw./(1-aw)).^c)); 
    %
    h  =@(aw)(1-aw);  
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF))).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            );    
elseif isotherm{iso}=="Oswin II" 
    xi =@(omega,kappa,c,omegaFSP,b)((omega./kappa./(1-(omega./omegaFSP).^m)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP,b)(((xi(omega,kappa,c,omegaFSP,b).^(1./b))./(1+xi(omega,kappa,c,omegaFSP,b).^(1./b))).^b);   
    %
    v =@(aw,kappa,c,omegaFSP,b) (omegaFSP./(3.*kappa.*((aw.^(1./b)./(1-aw.^(1./b))).^b).^c)); 
    %
    h  =@(aw,b) ((1-aw.^(1./b))); 
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF,b) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF),b),b).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            );          
elseif isotherm{iso}=="LW I"  
    xi =@(omega,kappa,c,omegaFSP)((omega./kappa./(1-(omega./omegaFSP).^m)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP)(1-exp(-W(xi(omega,kappa,c,omegaFSP)))); 
    %
    v =@(aw,kappa,c,omegaFSP) (omegaFSP./(3.*kappa.*(-log(1-aw)./(1-aw)).^c));     
    %
    h  =@(aw)(-(1-aw)./aw.*log(1-aw)./(1-log(1-aw))); 
    Hiso =@(omega,theta,c,omegaREF,kappaInf,theta0,n,theta1,thetaREF) ... 
           ( HisoMax(omega,theta,theta0,n).* ...
             h(aw(omega,kappa(theta,kappaInf,theta0,n),c,omegaFSP(theta,theta1,thetaREF,omegaREF))).* ...
             eta(omega,theta,theta0,n,theta1,thetaREF,omegaREF) ...
            );   
elseif isotherm{iso}=="LW II" 
    xi =@(omega,kappa,c,omegaFSP,b)((omega./kappa./(1-(omega./omegaFSP).^m)).^(1./c));
    aw =@(omega,kappa,c,omegaFSP,b)((1-exp(-W(xi(omega,kappa,c,omegaFSP,b).^(1./b)))).^b); 
    %
    v =@(aw,kappa,c,omegaFSP,b) (omegaFSP./(3.*kappa.*((-log(1-aw.^(1./b))./(1-aw.^(1./b))).^b).^c)); 
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
    %
    y =@(omega,omegaFSP) ((omega./omegaFSP));
    A =@(C,K,y) ((C-2).*K.*(y-1)+(C-1).*K.^2-1);
    B =@(C,K,y) (2.*(C-1).*K.^2.*y);
    dy = @(omega,omegaFSP,theta1,omegaREF) ((omega./omegaFSP)./omegaFSP.*omegaREF./theta1);
    rdA = @(theta,C,K,y,qC,qK,dy) (-C.*K.*qC.*(K+y-1) - K.*qK.*((C-2).*(y-1) + 2.*(C-1).*K) + (C-2).*K.*dy.*R.*theta.^2);
    rdB = @(theta,C,K,y,qC,qK,dy) (-2.*C.*qC.*K^2.*y - 4.*(C-1).*K.^2.*qK.*y + 2.*(C-1).*K.^2.*dy.*R.*theta.^2);
    %
    Hiso=@(theta,A,B,y,rdA,rdB,dy) ((A.*rdA + rdB.*y + B.*dy.*R.*theta.^2 )./((sqrt(A.^2+2*B.*y)+A).*sqrt(A.^2+2*B.*y)) + (rdA)./((sqrt(A.^2+2*B.*y)+A)) - (rdB)./(B)); 
elseif isotherm{iso}=="GAB"     
    % GAB-model 
    % Furmaniak, Terzyk and Gauden, Journal of Food Engineering 82 (2007)
    % 528-535
    %
    K  =@(theta,K0,qK) (K0.*exp(qK/R./theta)); 
    C  =@(theta,C0,qC) (C0.*exp(qC/R./theta));      
    aw =@(omega,mGAB,C,K)((sqrt(C.^2.*(mGAB-omega).^2+4*C.*mGAB.*omega)-C*(mGAB-omega)-2*omega)./(2*(C-1)*K*omega));
    %
    omegaFSPGAB=@(mGAB,C,K)(mGAB*C*K./(1-K)./(1+(C-1).*K));
    %
    Hiso=@(theta,C0,qC,K0,qK,omega,mGAB) (... 
       ((1-K(theta,K0,qK).*aw(omega,mGAB,C(theta,C0,qC),K(theta,K0,qK))).^2)./(1+(C(theta,C0,qC)-1).*K(theta,K0,qK).^2.*aw(omega,mGAB,C(theta,C0,qC),K(theta,K0,qK)).^2).*qC + qK);      
end    
%
%   Invers for m=3
    if strfind(isotherm{iso}, "II")>0
        w=@(aw,kappa,c,omegaFSP,b) (omegaFSP.*((0.5+sqrt(0.25+v(aw,kappa,c,omegaFSP,b).^3)).^(1/3)-v(aw,kappa,c,omegaFSP,b)./(0.5+sqrt(0.25+v(aw,kappa,c,omegaFSP,b).^3)).^(1/3)));            
    elseif isotherm{iso}=="GAB I" 
        w=@(aw,omegaFSP,C,K)(omegaFSP.*aw.*(1-K)./(1-K.*aw).*(1-K+C.*K)./(1-K.*aw+C.*K.*aw));
    elseif isotherm{iso}=="GAB" 
        w=@(aw,mGAB,C,K)(mGAB.*C.*K.*aw./(1-K.*aw)./(1-K.*aw+C.*K.*aw));      
    else
        w=@(aw,kappa,c,omegaFSP) (omegaFSP.*((0.5+sqrt(0.25+v(aw,kappa,c,omegaFSP).^3)).^(1/3)-v(aw,kappa,c,omegaFSP)./(0.5+sqrt(0.25+v(aw,kappa,c,omegaFSP).^3)).^(1/3)));        
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
    if isotherm{iso}=="GAB I"  
        omega{i}=0.0:0.001:omegaFSP(theta(i),theta1,thetaREF,omegaREF);
        plot(omega{i},aw(omega{i},omegaFSP(theta(i),theta1,thetaREF,omegaREF),C(theta(i),C0,qC),K(theta(i),K0,qK)), ...          
             markerLine{i},'color',markerColors{iso},'LineWidth',2)    
        if plotAlsoOmega=="Yes"
            awX=0:1.e-3:1; 
            plot(w(awX,omegaFSP(theta(i),theta1,thetaREF,omegaREF),C(theta(i),C0,qC),K(theta(i),K0,qK)),awX,'w--','LineWidth',1); 
        end            
    elseif isotherm{iso}=="GAB"  
        omegaGAB{i}=0.0:0.001:omegaFSPGAB(mGAB,C(theta(i),C0,qC),K(theta(i),K0,qK));
        plot(omegaGAB{i},aw(omegaGAB{i},mGAB,C(theta(i),C0,qC),K(theta(i),K0,qK)), ...          
             markerLine{i},'color',markerColors{iso},'LineWidth',2)     
        if plotAlsoOmega=="Yes"
            awX=0:1.e-3:1; 
            plot(w(awX,mGAB,C(theta(i),C0,qC),K(theta(i),K0,qK)),awX,'w--','LineWidth',1); 
        end         
    elseif strfind(isotherm{iso}, "II")>0 
        omega{i}=0.0:0.001:omegaFSP(theta(i),theta1,thetaREF,omegaREF);
        plot(omega{i},aw(omega{i},kappa(theta(i),kappaInf,theta0,n),c, ...
             omegaFSP(theta(i),theta1,thetaREF,omegaREF),b),...
             markerLine{i},'color',markerColors{iso},'LineWidth',2); 
        if plotAlsoOmega=="Yes"         
            awX=0:1.e-3:1; 
            plot(w(awX,kappa(theta(i),kappaInf,theta0,n),c,omegaFSP(theta(i),theta1,thetaREF,omegaREF),b),awX,'k','LineWidth',1);  
        end
    else     
        omega{i}=0.0:0.001:omegaFSP(theta(i),theta1,thetaREF,omegaREF);
        plot(omega{i},aw(omega{i},kappa(theta(i),kappaInf,theta0,n),c, ...
             omegaFSP(theta(i),theta1,thetaREF,omegaREF)),...
             markerLine{i},'color',markerColors{iso},'LineWidth',2); 
        if plotAlsoOmega=="Yes"
            awX=0:1.e-3:1; 
            plot(w(awX,kappa(theta(i),kappaInf,theta0,n),c,omegaFSP(theta(i),theta1,thetaREF,omegaREF)),awX,'k','LineWidth',1); 
        end
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
    elseif isotherm{iso}== "GAB" 
        plot(omegaGAB{i},Hiso(theta(i),C0,qC,K0,qK,omegaGAB{i},mGAB),markerLine{i},'color',markerColors{iso},'LineWidth',1);
    elseif strfind(isotherm{iso}, "II")>0 
       plot(omega{i},Hiso(omega{i},theta(i),c,omegaREF,kappaInf,theta0,n,theta1,thetaREF,b), ...
            markerLine{i},'color',markerColors{iso},'LineWidth',2);    
    else
       plot(omega{i},Hiso(omega{i},theta(i),c,omegaREF,kappaInf,theta0,n,theta1,thetaREF), ...
            markerLine{i},'color',markerColors{iso},'LineWidth',1);
    end
end
    axis([0 0.45 0 1700])
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
    while any(abs(y-yInf)./abs(y)>1.e-6)
       yInf = y;
       y = y-(y.*exp(y)-x)./((exp(y).*(y+1)-(y+2).*(y.*exp(y)-x)./(2*y+2)));
    end
end

 




