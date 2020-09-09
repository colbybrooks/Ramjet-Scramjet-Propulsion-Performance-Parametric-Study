%% Aircraft Propulsion Project

%% A) T-S Diagrams
close all
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Validation Case 7 First in Sections              M2 - 0.15
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z = 4300;               % m
M1 = 2.4;
nd = 0.92;
M2 = 0.15;
qf = 43.2*10^6;         % J/kg
Tt3max = 2400;
nn = 0.94;
Ae = 0.015;
y = 1.4;
ymixture = 1.3;
R = 286.9;              % J/(kg-K)
a = 986;
b = 0.179;

[s12,s23,s13,s1e,s14,T1,T2,T3,Te,T4] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);

figure
n = 1*10^(6);

s_12 = linspace(0,s12,n);
T12 = linspace(T1,T2,n);
plot(s_12,T12)
hold on

ds23 = s23/n;
T23 = T2;
for i = 1:n;
    dT23(i) = T23(i)/(a+b*T23(i))*ds23;
    T23(i+1) = T23(i)+dT23(i);
    if T23(i+1) >= T3;
        break
    end
end
s_23 = linspace(s12,s13,length(T23));
plot(s_23,T23)

s_3e = linspace(s13,s1e,n);
T3e = linspace(T3,Te,n);
plot(s_3e,T3e)

s_e4 = linspace(s1e,s14,n);
Te4 = linspace(Te,T4,n);
plot(s_e4,Te4)

ds41 = s14/n;
T41 = T4;
cp = 1004;
for i = 1:n;
    dT41(i) = T41(i)/cp*ds41;
    T41(i+1) = T41(i)-dT41(i);
    if T41(i+1) <= T1;
        break
    end
end
s_41 = linspace(s14,0,length(T41));
plot(s_41,T41)
hold off
xlim([-500 2500]) 

legend('1-2','2-3','3-e','e-4','4-1','location','northwest')
s = 12;
xlabel('delta s (J/kg-k)','fontsize',s)
ylabel('T (k)','fontsize',s)
title('M2 = 0.15')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Validation Case 8 Begins Post Section Cutter     M2 = 0.4
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

M2 = 0.4;
[s12,s23,s13,s1e,s14,T1,T2,T3,Te,T4] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);

figure
n = 1*10^(6);

s_12 = linspace(0,s12,n);
T12 = linspace(T1,T2,n);
plot(s_12,T12)
hold on

ds23 = s23/n;
T23 = T2;
for i = 1:n;
    dT23(i) = T23(i)/(a+b*T23(i))*ds23;
    T23(i+1) = T23(i)+dT23(i);
    if T23(i+1) >= T3;
        break
    end
end
s_23 = linspace(s12,s13,length(T23));
plot(s_23,T23)

s_3e = linspace(s13,s1e,n);
T3e = linspace(T3,Te,n);
plot(s_3e,T3e)

s_e4 = linspace(s1e,s14,n);
Te4 = linspace(Te,T4,n);
plot(s_e4,Te4)

ds41 = s14/n;
T41 = T4;
cp = 1004;
for i = 1:n;
    dT41(i) = T41(i)/cp*ds41;
    T41(i+1) = T41(i)-dT41(i);
    if T41(i+1) <= T1;
        break
    end
end
s_41 = linspace(s14,0,length(T41));
plot(s_41,T41)
hold off
xlim([-500 2500]) 
title('M2 = 0.4')
legend('1-2','2-3','3-e','e-4','4-1','location','northwest')
xlabel('delta s (J/kg-k)','fontsize',s)
ylabel('T (k)','fontsize',s)

%% B) Atmospheric Model vs International Standard Atmosphere

data = load('InternationalStandardAtmosphere.txt');
z = data(:,1);
T1_ISA = data(:,2);
P1_ISA = data(:,3);

Ts = 288;               % K
Ps = 101.3;             % kPa
zstar = 8404;           % m
T1 = (z < 7958).*Ts.*(1-(y-1)/y.*z/zstar) + (z >= 7958).*210;
P1 = (z < 7958).*Ps.*(1-(y-1)/y.*z/zstar).^(y/(y-1))...
    +(z>=7958).*33.6.*exp(-(z-7958)./6605);
figure
plot(z,P1,z,P1_ISA)
xlabel('Altitude z (m)','fontsize',s)
ylabel('Pressure P (kPA)','fontsize',s)
legend('Atm Model','Inter Std Atm')
figure
plot(z,T1,z,T1_ISA)
ylim([200 290])
xlabel('Altitude z (m)','fontsize',s)
ylabel('Temperature T (K)','fontsize',s)
legend('Atm Model','Inter Std Atm')

%% C) Flight Mach Number Trade Study

z = 4300;               % m
M1 = [0.8:0.01:5];
nd = 0.92;
M2 = 0.15;
qf = 43.2*10^6;         % J/kg
Tt3max = 2400;
nn = 0.94;
Ae = 0.015;
y = 1.4;
ymixture = 1.3;
R = 286.9;              % J/(kg-K)
a = 986;
b = 0.179;

[~,~,~,~,~,~,~,~,~,~,no,T,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
figure
plot(M1,no)
title('M2 = 0.15')
xlabel('M1','fontsize',s)
ylabel('Overall Efficiency n_o','fontsize',s)
figure
plot(M1,T)
title('M2 = 0.15')
xlabel('M1','fontsize',s)
ylabel('Thrust T (N)','fontsize',s)
figure
plot(M1,TSFC)
title('M2 = 0.15')
xlabel('M1','fontsize',s)
ylabel('Thrust Specific Fuel Consuption TSFC (kg/hr)/N','fontsize',s)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Validation Case 8 

M2 = 0.4;
[~,~,~,~,~,~,~,~,~,~,no,T,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
figure
plot(M1,no)
title('M2 = 0.4')
xlabel('M1','fontsize',s)
ylabel('Overall Efficiency n_o','fontsize',s)
figure
plot(M1,T)
title('M2 = 0.4')
xlabel('M1','fontsize',s)
ylabel('Thrust T (N)','fontsize',s)
figure
plot(M1,TSFC)
title('M2 = 0.4')
xlabel('M1','fontsize',s)
ylabel('Thrust Specific Fuel Consuption TSFC (kg/hr)/N','fontsize',s)

%% D) Altitude Trade Study
z = [2000:30000];
M1 = 2.4;
M2 = 0.15;

[~,~,~,~,~,~,~,~,~,~,no,T,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
figure
plot(z,no)
title('M2 = 0.15')
xlabel('Altitude z (m)','fontsize',s)
ylabel('Overall Efficiency n_o','fontsize',s)
figure
plot(z,T)
title('M2 = 0.15')
xlabel('Altitude z (m)','fontsize',s)
ylabel('Thrust T (N)','fontsize',s)
figure
plot(z,TSFC)
title('M2 = 0.15')
xlabel('Altitude z (m)','fontsize',s)
ylabel('Thrust Specific Fuel Consuption TSFC (kg/hr)/N','fontsize',s)


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Validation Case 8

M2 = 0.4;
[~,~,~,~,~,~,~,~,~,~,no,T,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);

figure
plot(z,no)
title('M2 = 0.4')
xlabel('Altitude z (m)','fontsize',s)
ylabel('Overall Efficiency n_o','fontsize',s)
figure
plot(z,T)
title('M2 = 0.4')
xlabel('Altitude z (m)','fontsize',s)
ylabel('Thrust T (N)','fontsize',s)
figure
plot(z,TSFC)
title('M2 = 0.4')
xlabel('Altitude z (m)','fontsize',s)
ylabel('Thrust Specific Fuel Consuption TSFC (kg/hr)/N','fontsize',s)

%% E) Optimal Flight Speed per varying Altitude

n = 1*10^(5);
z_vector = [2000:500:20000];
M2 = 0.15;
M1 = linspace(0.8,4,n);

for i = 1:length(z_vector);
    z = z_vector(i);
[~,~,~,~,~,~,~,~,~,~,no,~,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
no_matrix(i,:) = no;
TSFC_matrix(i,:) = TSFC;
end

[TSFC,TSFC_index] = min(TSFC_matrix');
[no,no_index] = max(no_matrix');

for i = 1:length(TSFC_index);
    M1_TSFC(i) = M1(TSFC_index(i));
    M1_no(i) = M1(no_index(i));
end
figure
plot(z_vector,M1_TSFC)
xlabel('Altitude z (m)')
ylabel('Optimal M1')
title('M2 = 0.15, Minimize TSFC')
figure
plot(z_vector,M1_no)
xlabel('Altitude z (m)')
ylabel('Optimal M1')
title('M2 = 0.15, Maximize Overall Efficiency n_o')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Validation Case 8

M2 = 0.4;

for i = 1:length(z_vector);
    z = z_vector(i);
[~,~,~,~,~,~,~,~,~,~,no,~,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
no_matrix(i,:) = no;
TSFC_matrix(i,:) = TSFC;
end

[TSFC,TSFC_index] = min(TSFC_matrix');
[no,no_index] = max(no_matrix');

for i = 1:length(TSFC_index);
    M1_TSFC(i) = M1(TSFC_index(i));
    M1_no(i) = M1(no_index(i));
end
figure
plot(z_vector,M1_TSFC)
xlabel('Altitude z (m)')
ylabel('Optimal M1')
title('M2 = 0.4, Minimize TSFC')
figure
plot(z_vector,M1_no)
xlabel('Altitude z (m)')
ylabel('Optimal M1')
title('M2 = 0.4, Maximize Overall Efficiency n_o')


%% F) Inlet/Diffuser Trade Study

z = 4300;
M1 = 2.4;
nd = [0.5:0.05:1];
M2 = 0.15;

[~,~,~,~,~,~,~,~,~,~,no,T,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
figure
plot(nd,no)
title('M2 = 0.15')
xlabel('Inlet/Diffuser Efficiency n_d ','fontsize',s)
ylabel('Overall Efficiency n_o','fontsize',s)
figure
plot(nd,T)
title('M2 = 0.15')
xlabel('Inlet/Diffuser Efficiency n_d ','fontsize',s)
ylabel('Thrust T (N)','fontsize',s)
figure
plot(nd,TSFC)
title('M2 = 0.15')
xlabel('Inlet/Diffuser Efficiency n_d ','fontsize',s)
ylabel('Thrust Specific Fuel Consuption TSFC (kg/hr)/N','fontsize',s)


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Validation Case 8

M2 = 0.4;
[~,~,~,~,~,~,~,~,~,~,no,T,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
figure
plot(nd,no)
title('M2 = 0.4')
xlabel('Inlet/Diffuser Efficiency n_d ','fontsize',s)
ylabel('Overall Efficiency n_o','fontsize',s)
figure
plot(nd,T)
title('M2 = 0.4')
xlabel('Inlet/Diffuser Efficiency n_d ','fontsize',s)
ylabel('Thrust T (N)','fontsize',s)
figure
plot(nd,TSFC)
title('M2 = 0.4')
xlabel('Inlet/Diffuser Efficiency n_d ','fontsize',s)
ylabel('Thrust Specific Fuel Consuption TSFC (kg/hr)/N','fontsize',s)

%% G) - Nozzle Trade Study

z = 4300;
M1 = 2.4;
M2 = 0.15;
nd = 0.92;
nn = [0.5:0.05:1];

[~,~,~,~,~,~,~,~,~,~,no,T,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
figure
plot(nn,no)
title('M2 = 0.15')
xlabel('Nozzle Efficiency n_n ','fontsize',s)
ylabel('Overall Efficiency n_o','fontsize',s)
figure
plot(nn,T)
title('M2 = 0.15')
xlabel('Nozzle Efficiency n_n ','fontsize',s)
ylabel('Thrust T (N)','fontsize',s)
figure
plot(nn,TSFC)
title('M2 = 0.15')
xlabel('Nozzle Efficiency n_n ','fontsize',s)
ylabel('Thrust Specific Fuel Consuption TSFC (kg/hr)/N','fontsize',s)


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Validation Case 8

M2 = 0.4;
[~,~,~,~,~,~,~,~,~,~,no,T,TSFC] =...
    Ramjet(z,M1,nd,M2,qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
figure
plot(nn,no)
title('M2 = 0.4')
xlabel('Nozzle Efficiency n_n ','fontsize',s)
ylabel('Overall Efficiency n_o','fontsize',s)
figure
plot(nn,T)
title('M2 = 0.4')
xlabel('Nozzle Efficiency n_n ','fontsize',s)
ylabel('Thrust T (N)','fontsize',s)
figure
plot(nn,TSFC)
title('M2 = 0.4')
xlabel('Nozzle Efficiency n_n ','fontsize',s)
ylabel('Thrust Specific Fuel Consuption TSFC (kg/hr)/N','fontsize',s)

%% H) Mach# at State 2 Trade Study

z = 4300;
M1 = 2.4;
nn = 0.94;
M2 = [0.1:0.1:2.5];

no = zeros(1,length(M2));
T = zeros(1,length(M2));
TSFC = zeros(1,length(M2));
for i = 1:length(M2);
[~,~,~,~,~,~,~,~,~,~,no(i),T(i),TSFC(i)] =...
    Ramjet(z,M1,nd,M2(i),qf,Tt3max,nn,Ae,y,ymixture,R,a,b);
end

figure
plot(M2,no)
xlabel('M2 ','fontsize',s)
ylabel('Overall Efficiency n_o','fontsize',s)
figure
plot(M2,T)
xlabel('M2 ','fontsize',s)
ylabel('Thrust T (N)','fontsize',s)
figure
plot(M2,TSFC)
xlabel('M2 ','fontsize',s)
ylabel('Thrust Specific Fuel Consuption TSFC (kg/hr)/N','fontsize',s)

%% I) Scramjet/Ramjet Optimal Thrust
close all

z = 90000 * 0.3048;            % m
M1 = 5;
M2 = [0.1:0.01:5];
Tt3 = [1300:2400];

Thrust = zeros(length(M2),length(Tt3));

for i = 1:length(M2);
    for j = 1: length(Tt3);
        [T]=ramjetOptT(z,M1,nd,M2(i),qf,Tt3,Tt3max,nn,Ae,y,ymixture,R,a,b);
        Thrust(i,:) = real(T);
    end
end

[index_M2,index_Tt3] = find(Thrust == max(max(Thrust)));
Optimal_M2 = M2(index_M2)
Optimal_Tt3 = Tt3(index_M2)

maxThrust = max(max(Thrust))

%% Code

type Ramjet
type ramjetOptT
