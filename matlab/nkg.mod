
//NK model with Government - Chapter 7 (UNDERSTANDING DSGE MODELS)

var Y IIPP I GG CCR CNR G
 KP KG L LR LNR
 R W CM P PI PIW Q RB LAMBDAR LAMBDANR
 T tau_c tau_k tau_l TRANS B
 A Sm SG SIG STRANS Stau_c Stau_l Stau_k U;

varexo e e_m e_G e_IG e_TRANS e_tau_c e_tau_l e_tau_k; 

parameters sigma phi alpha1 alpha2 alpha3 beta delta
 deltaG rhoa psi theta thetaW psiW
 rhoG rhoIG rhoTRANS rhotau_c rhotau_l rhotau_k rhom
 gammaG gammaIG gammaTRANS gammatau_c gammatau_l gammatau_k
 phiG phiIG phiTRANS phitau_c phitau_l phitau_k
 gammaR gammaPI gammaY tau_css tau_lss tau_kss phic omegaR Psi1 Psi2 chi; 

alpha1 = 0.3;
alpha2 = 0.65;
alpha3 = 0.05;
beta= 0.985;
delta= 0.025;
deltaG = 0.025;
theta = 0.75;
thetaW = 0.75;
sigma = 2;
phi = 1.5;
psi= 8;
psiW = 21;
tau_css = 0.16;
tau_lss = 0.17;
tau_kss = 0.08;
phic = 0.8;
omegaR = 0.5;
Psi2 = 1;
chi = 1;

//Fiscal Policy Parameters
gammaG = 0; 
gammaIG = 0.1;
gammaTRANS 0.1;
gammatau_c = 0;
gammatau_l = 0;
gammatau_k = 0;
phiG = O;
phiIG = -0.1;
phiTRANS = -0.1;
phitau_c = 0;
phitau_l = 0;
phitau_k = 0;

//Taylor's Rule Parameters
gammaR = 0.8;
gammaY = 0.5;
gammaPI = 1.5;

//Autoregressive Shock Parameters
rhoa = 0.9;
rhoG = 0.9; 
rhoIG = 0.9;
rhoTRANS = 0.9;
rhotau_c = 0.9;
rhotau_l = 0.9;
rhotau_k = 0.9;
rhom = 0.9; 
Psii = (i+tau_css)*((ilbeta)-(1-delta));

model(linear);

#phiB = 1;
#phi_TRANS = 0.01;
#phi_IG = 0.02;
#Uss = 1;
#Pss = 1;
#PIss= 1;
#RBss = 1/beta;

#Rss = Pss*((1+tau_css)/(1-tau_kss))*((1/beta)-(1-delta));

#CMss = ((psi-l)lpsi)*(l-beta*theta)*Pss;

#Wss = alpha2*((CMss*0.2-alpha3)-(1/alpha2))*((alpha1/Rss)-(alpha1/alpha2));

#Al= ((1-phic*beta)*((1-phic)^(-sigma))*(1-beta*thetaW)*((psiW-1)/psiW)*((1-tau_lss)/(1+tau_css))*(Wss/Pss)* (Wss/(alpha2*CMss))^phi)^(1/sigma);

#A2 = (((Rss*(Pss-tau_lss*(1-alpha1)*CMss)-tau_kss*(Rss-delta)*alpha1*CMss)/(Pss*Rss*(1+tau_css)))-(delta*alpha1*CMss/Rss)-(phiB/Pss) *((1/RBss)-1) + phi_TRANS);

#Yss = (A1/A2)-(sigma/(sigma+phi));

#Bss = phiB*Yss;

#Lss = alpha2*CMss*(Yss/Wss);

#LRss= Lss;

#LNRss = Lss;

#KPss = alphal*CMss*(Yss/Rss);

#IPss = delta*KPss;

#IGss = phi_IG*Yss;

#KGss = IGss/deltaG;

#Css = (1/(Yss-(phi/sigma)))*Al;
#CRss = Css;
#CNRss = Css;

#Gss = Yss - IPss - IGss - Css;

#TRANSss = phiTRANS*Yss;

#Tss = Pss*Gss + Pss*IGss + Pss*TRANSss - Bss*((l/RBss)-1);

#LAMBDARss = ((CRss^(-sigma))*(l-phic*beta)*(1-phic)^(-sigma))/((l+tau_css)*Pss);

#LAMBDANRss = ((CNRss^(-sigma))*(l-phic*beta)*(1-phic)^(-sigma))/((l+tau_css)*Pss);

#Qss = LAMBDARss*Pss*(l+tau_css);

//1-Ricardian Lagrangian household
LAMBDAR + P + (tau_cssl(l+tau_css))*tau_c = (sigma/((1-phic)*(1-phic*beta)))*(phic*beta*(CR(+l)-CR)-(CR-CR(-1)));

//2-Phillips equation for Ricardian household wages
PIW = beta*PIW(+1)+((1-thetaW)*(1-beta*thetaW)lthetaW)*(phi*LR-LAMBDAR+(tau_lss/(1-tau_lss))*tau_l);

//3-Gross wage inflations
PIW = W - W(-1);

//4-Ricardian household budget constraint
Pss*CRss*((P+CR)*(i+tau_css)+tau_css*tau_c) + Pss*IPss*((P+IP)*(l+tau_css)+tau_css*tau_c) + (BsslRBss)*(B-RB) = Wss*LRss*((W+LR)*(l-tau_lss)-tau_lss*tau_l)+Rss*KPss*((R+KP(-1))*(1-tau_kss)-tau_kss*tau_k) + Bss*B(-1) + omegaR*TRANSss*TRANS;

//5-Tobin's Q
(Qss/beta)*Q = (1-delta)*Qss*Q(+i) + LAMBDARss*Rss*Uss*(1-tau_kss)*(LAMBDAR(+i)+R(+i)+U(+i)-(tau_kss/(1-tau_kss)) *tau_k(+1)) - LAMBDARss*Pss*Uss*Psi1*U(+1);

//6-Demand for installed capacity 
(1-tau_kss)*(Rss/Pss)*(R-P-(tau_kss/(1-tau_kss))*tau_k)=Psi2*Uss*U;

//7-Demand for investments
(1+tau_css)*LAMBDARss*Pss*(LAMBDAR+P+(tau_css/(1+tau_css))*tau_c) - Qss*Q+chi*Qss*(IP-IP(-1)) = chi*beta*Qss*(IP(+1)-IP);

//8-Law of motion of private capital
KP = (1-delta)*KP(-1) + delta*IP;

//9-Euler's equation (Public bond)
LAMBDAR - RB = LAMBDAR(+1);

//10-Non-Ricardian household Lagrangian
LAMBDANR + P + (tau_css/(1+tau_css))*tau_c = (sigma/((1-phic)*(1-phic*beta)))*(phic*beta*(CNR(+1)-CNR)-(CNR-CNR(-1)));

//11-Phillips equation for non-Ricardian household wages
PIW = beta*PIW(+1)+((1-thetaW)*(1-beta*thetaW)/thetaW)*(phi*LNR-LAMBDANR+(tau_lss/(1-tau_lss))*tau_l);

//12-Aggregate consumption
Css*C = omegaR*CRss*CR + (1-omegaR)*CNRss*CNR;

//13-Aggregate labor
Lss*L = omegaR*LRss*LR + (1-omegaR)*LNRss*LNR;

//14-Production Function
Y = A + alphal*(U+KP(-1)) + alpha2*L + alpha3*KG(-1);

//15- Problem of the firm trade-off (MRS=Relative price)
L - U - KP(-1) = R - W;

//16-Marginal Cost
CM = alpha2*W + alpha1*R - A - alpha3*KG(-1);

//17-Phillips Equation
PI = beta*PI(+1) + ((1-theta)*(1-beta*theta)/theta)*(CM-P);

//18-Gross Inflation Rate
PI(+1) = P(+1) - P;

//19-Government budget constraint
(Bss/RBss)*(B-RB)-Bss*B(-1) + Tss*T = Pss*Gss*(P+G) + Pss*IGss*(P+IG) + Pss*TRANSss*(P+TRANS);

//20-Government tax revenue
Tss*T = tau_css*Pss*(Css*(C+P+tau_c)+IPss*(IP+P+tau_c)) + tau_lss*Wss*Lss*(W+L+tau_l) + tau_kss*KPss*(Rss*(R+KP(-1)+tau_k) - delta*(KP(-1)+tau_k));

//21-Rule for the movement of public capital
KG = (1-deltaG)*KG(-1) + deltaG*IG;

//22-Rule for the movement of public spending
G = gammaG*G(-1) + (1-gammaG)*phiG*(B(-1)-Y(-1)-P(-l))+SG;

//23-Rule for the movement of public investment
IG = gammaIG*IG(-1) + (1-gammaIG)*phiIG*(B(-1)-Y(-1)-P(-i))+SIG;

//24-Rule for the movement of transfer of income
TRANS= gammaTRANS*TRANS(-1) + (1-gammaTRANS)*phiTRANS*(B(-1)-Y(-1)-P(-i)) + STRANS;

//25-Rule for the movement of tax on consumption
tau_c = gammatau_c*tau_c(-1) = (1-gammatau_c)*phitau_c*(B(-1)-Y(-1)-P(-i))+Stau_c;

//26-Rule for the movement of tax on labor Income
tau_l = gammatau_l*tau_l(-1) + (1-gammatau_l)*phitau_l*(B(-1)-Y(-1)-P(-i))+Stau_l;

//27-Rule for the movement of tax on consumption
tau_k = gammatau_k*tau_k(-1) + (1-gammatau_k)*phitau_k*(B(-1)-Y(-1)-P(-i))+Stau_k;

//28-Taylor's rule
RB = gammaR*RB(-1)+(1-gammaR)*(gamrnaPI*PI + gammaY*Y)+Sm;

//29-Equilibrim condition
Yss*Y = Css*C + IPss*IP + IGss*IG + Gss*G;

//30-Productivity shock
A = rhoa*A(-1) e; 

//31 - Shock in Public Spending
SG = rhoG*SG(-1) + e_G; 

//32 - Shock in Public Investment
SIG = rhoIG*SIG(-1) + e_IG; 

//33 - Shock in Transfer of Income
TRANS = rhoTRANS*STRANS(-1) + e_TRANS; 

//34 - Shock in tax on Consumption
Stau_c = rhotau_c*Stau_c(-1) - e_tau_c;

//35 - Shock in tax on Labor Income
Stau_l = rhotau_l*Stau_l(-1) - e_tau_l;

//36 - Shock in tax on Capital Income
Stau_k = rhotau_k*Stau_k(-1) - e_tau_k;

//37- Monetary Shock
Sm = rhom*Sm(-1) - e_m;

end; 

steady;

check(qz_zero_threshold=le-20); 

shocks;
var e; stderr 0.01;
var e_G; stderr 0.01;
var e_IG; stderr 0.01;
var e_TRANS; stderr 0.01;
var e_tau_c; stderr 0.01;
var e_tau_l; stderr 0.01;
var e_tau_k; stderr 0.01;
var e_m; stderr 0.01;
end; 

stoch_simul(periods=1000,qz_zero_threshold=1e-20) Y IP IG CR CNR G KP KG LR LNR R W UU PI RB T BB tau_c tau_k tau_l TRANS A; 

// eof 




