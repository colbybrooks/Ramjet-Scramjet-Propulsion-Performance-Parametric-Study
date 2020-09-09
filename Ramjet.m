function [s12,s23,s13,s1e,s14,T1,T2,T3,Te,T4,no,T,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%% i) Module 1 - State Conditions

Ts = 288;               % K
Ps = 101.3;             % kPa
zstar = 8404;           % m

T1 = (z < 7958).*Ts.*(1-(y-1)/y.*z/zstar) + (z >= 7958).*210;
P1 = (z < 7958).*Ps.*(1-(y-1)/y.*z/zstar).^(y/(y-1))...
    +(z>=7958).*33.6.*exp(-(z-7958)./6605);

Tt1 = T1.*(1+(y-1)/2.*M1.^2);
Pt1 = P1.*(1+(y-1)/2.*M1.^2).^(y/(y-1));
a1 = sqrt(y*R*T1);
v1 = a1.*M1;

%% ii) Module 2 - Inlet/Diffuser

Tt2 = Tt1;
T2 = Tt2.*(1+(y-1)/2*M2.^2).^(-1);
Pt2 = P1.*(1+nd*(y-1)/2*M1.^2).^(y/(y-1));
P2 = Pt2.*(1+(y-1)/2*M2.^2).^(-y/(y-1));
Cp2 = a+b*T2;
s12 = Cp2.*log(Tt2./Tt1)-R*log(Pt2./Pt1);
a2 = sqrt(y*R*T2);
v2 = a2.*M2;

%% iii) Module 3 - Combustor

y = ymixture;

Tt3choked = Tt2.*(1/(2*(y+1))*1./M2.^2.*(1+y*M2.^2).^2.*(1+(y-1)/2*M2.^2).^(-1));

for i = 1:length(Tt3choked)
if Tt3choked(i) < Tt3max
    Tt3(i) = Tt3choked(i);
    M3(i) = 1;
else
    Tt3(i) = Tt3max;
    cc(i) = Tt3(i)./Tt2(i).*(1+(y-1)/2*M2.^2).*M2.^2./(1+y*M2.^2).^2;
    bb(i) = 2*cc(i)*y-1;
    aa(i) = cc(i)*y^2-(y-1)/2;
    M3n(i) = sqrt((-bb(i)-sqrt(bb(i).^2-4*aa(i).*cc(i)))./(2*aa(i)));
    M3p(i) = sqrt((-bb(i)+sqrt(bb(i).^2-4*aa(i).*cc(i)))./(2*aa(i)));
        if M3p(i) > M3n(i) && M3p(i) <= 1
            M3(i) = M3p(i);
        else
            M3(i) = M3n(i);
        end
end
end

q23 = a.*(Tt3-Tt2)+1/2*b.*(Tt3.^2-Tt2.^2);
T3 = Tt3.*(1+(y-1)/2*M3.^2).^(-1);
P3 = P2;
Pt3 = P3.*(1+(y-1)/2.*M3.^2).^(y/(y-1));
a3 = sqrt(y*R*T3);
v3 = a3.*M3;
Cp3 = a+b*T3;
s23 = Cp3.*log(Tt3./Tt2)-R*log(Pt3./Pt2);
s13 = s12+s23;

%% iv) Module 4 - Converging Nozzle

Mtest = (2/(y-1)*nn*(1-(P1./Pt3).^((y-1)/y))./...
    (1-nn*(1-(P1./Pt3).^((y-1)/y)))).^(1/2);

if Mtest < 1
    Me = Mtest;
    Pe = P1;
else
    Me = 1;
    Pe = Pt3.*(1-1./nn.*((y-1)./(y+1))).^(y/(y-1));   
end

Tte = Tt3;
Te = Tte.*(1+(y-1)/2*Me.^2).^(-1);
Pte = Pe.*(1+(y-1)/2.*Me.^2).^(y/(y-1));
ae = sqrt(y*R*Te);
ve = ae.*Me;
me = Pe./(R*Te).*ve.*Ae*1000;           % Kg/s
Cpe = a+b*Te;
s3e = Cpe.*log(Tte./Tt3)-R*log(Pte./Pt3);
s1e = s12+s23+s3e;

%% v) Module 5 - External Flow

if Mtest < 1
    n_ext = 1;
else
    n_ext = Mtest.^(-0.3);
end


Tt4 = Tte;
T4 = Tte.*(1-n_ext.*(1-(P1./Pte).^((y-1)/y)));
M4 = (2/(y-1)*(Tt4./T4-1)).^(1/2);
P4 = P1;
Pt4 = P4.*(1+(y-1)/2*M4.^2).^(y/(y-1));

a4 = sqrt(y*R*T4);
v4 = a4.*M4;
Cp4 = a+b*T4;
se4 = Cp4.*log(Tt4./Tte)-R*log(Pt4./Pte);
s14 = s12+s23+s3e+se4;

%% vi) Module 6 - Performance Parameter

mi = me./(1+q23./qf);
mf = me-mi;
f = mf./mi;

T = mi.*((1+f).*ve-v1)+(Pe-P1)*Ae*10^3;
TSFC = mf./T*3600;
g = 9.81;
Isp = T./(mf*g);
veq = ve+(Pe-P1).*Ae*10^3./me;
nth = ((me.*veq.^2/2)-(mi.*v1.^2/2))./(mi.*q23);
np = 2./(1+veq./v1);
no = nth.*np;
P = T.*v1;
end

