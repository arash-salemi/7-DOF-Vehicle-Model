clc,
clear,
close all;

%% Mass
Mf = 58;             % mass of both front wheels
Muf = Mf/2;          % mass of onr front wheel
Mr = 62;             % mass of both rear wheels
Mur = Mr/2;          % mass of one rear wheel
M_sprung = 1100;     % mass of sprung mass
M_curb = 1220;       % mass of vehicle

%% Inertia for sprung mass

Ix_s = 3.66*10^8; %kg.mm^2
Iy_s = 1.46*10^9; %kg.mm^2
Iz_s = 1.58*10^9; %kg.mm^2

%% Distances

X_f = 0.291;    % mm
X_r = 2648.784; % mm
Y_f = 766.247;  % mm
Y_r = 748.666;  % mm

tf = Y_f;   % track
tr = Y_r;   % track
Weel_base = abs(X_f - X_r);  % wheel_base
a = 0.4*Weel_base; % distance to front wheels
b = 0.6*Weel_base; % distance to rear wheels

%% parameters

kf = 21.58; % stiffness coefficient for frong suspension (N/mm)
kr = 26.49; % stiffness coefficient for rear suspension (N/mm)
kt = 220;   % stiffness coefficient for tire (N/mm)
cf = 1.560; % damping coefficient for frong suspension (N/mm.s)
cr = 0.896; % damping coefficient for rear suspension (N/mm.s)

%% Input: U CAN CHANGE THIS

v = 40;     % Km/h
jade = 500; % Meter

V = v * (5/18);      % convert to M/s
t_final = jade/ V;   % S
t = 1:0.00001:t_final; % time span

%% SINE INPUT

lambda = 800; % Meter
W = 2*pi*V/lambda;
x = [0:0.0001:jade];
Wx =(2*pi)/lambda;

% z_gf = 0.20 * sin(W*t); % Meter
z_gf = 0.20 * sin(Wx*x); % Meter
% z_gr = 20 * sin(W*(t-(Weel_base*10^-3/V))); % Meter
z_gr = 0.20 * sin(Wx*(x-(Weel_base))); % Meter
% plot(x,z_gf)
% hold on;
% plot(x,z_gr)
% hold on;
% title('Sine')
% xlabel('x (m)');
% ylabel('Z ground (m)');
% legend('Z front','Z rear')


%% SINE SWEEP INPUT

A = ((0.040-0.001)/jade)*x + 0.001; % Meter
lambda = ((10-0.02)/jade)*x + 0.02; % Meter
w = 2*pi*V./lambda;  % for time
wx = (2*pi)./lambda; % for x
% z_gf = A .* sin(w.*t); % Meter
% z_gr = A .* sin(w.*(t-(Weel_base*10^-3/V))); % Meter
z_gf = A .* sin(wx.*(x+a)); % Meter
z_gr = A .* sin(wx.*(x-b)); % Meter
% plot(x,z_gf)
% hold on;
% plot(x,z_gr)
% title('Sine Sweep')
% xlabel('x (m)');
% ylabel('Z ground (m)');
% legend('Z front','Z rear')

%% Spring & Damper for front wheels

V_f_input = 2; %

V_f = [0.05 0.1 0.3 0.6 1 1.2 1.5 4];
reb = [15.5 31 92 141 205 238 286 690];
comp = [12 20 51 71 103 120 143 345];
% plot (V_f, reb,'b--o')
% hold on;

p = polyfit(V_f,reb,1);
reb_output = polyval(p,V_f_input); % kgF
p = polyfit(V_f,comp,3);
comp_output = polyval(p,V_f_input); % kgF
% plot (V_f, reb_output)

%% Spring & Damper for rear wheels

V_r_input = 2;

V_r = [0.05 0.1 0.3 0.6 1 1.2 1.5 4];
expansion = [13.5 26.5 80 122 178 206 248 598];
contraction = [10 17.5 44 61 90 103 124 300];

p = polyfit(V_f,expansion,3);
expansion_output = polyval(p,V_r_input); % kgF
p = polyfit(V_f,contraction,3);
contraction_output = polyval(p,V_r_input); % kgF
% plot (V_f, expansion_output)

%% force by reggressions

A_f = (kt*z_gf - reb_output)/Muf;
A_r = (kt*z_gf - expansion_output)/Mur;

%% Equations 

syms Z1 Z2 Z3 Z4 Zs Phi Theta DZ1 DZ2 DZ3 DZ4 DZs DPhi DTheta
syms z_gr z_gf

M = [Muf;Muf;Mur;Mur;M_sprung;Ix_s;Iy_s];
C = [-cf 0 0 0 cf cf*tf -cf*a;
    0 -cf 0 0 cf -cf*tf -cf*a;
    0 0 -cr 0 cr cr*tr cr*b;
    0 0 0 -cr cr -cr*tr cr*b;
    cf cf cr cr -2*(cf+cr) 0 2*cf*a-2*cr*b;
    cf*tf -cf*tf cr*tr -cr*tr 0 -2*(cf*tf^2+ cr*tr^2) 0;
    -cf*a -cf*a cr*b cr*b 2*(cf*a - cr*b) 0 -2*(cf*a^2+ cr*b^2)];

K = [-kf-kt 0 0 0 kf kf*tf -kf*a;
    0 -kf-kt 0 0 kf -kf*tf -kf*a;
    0 0 -kr-kt 0 kr kr*tr kr*b;
    0 0 0 -kr-kt kr -kr*tr kr*b;
    kf kf kr kr -2*(kf+kr) 0 2*kf*a-2*kr*b;
    kf*tf -kf*tf kr*tr -kr*tr 0 -2*(kf*tf^2+ kr*tr^2) 0;
    -kf*a -kf*a kr*b kr*b 2*(kf*a - kr*b) 0 -2*(kf*a^2+ kr*b^2)];

Z = [Z1 Z2 Z3 Z4 Zs Phi Theta]';
DZ = [DZ1 DZ2 DZ3 DZ4 DZs DPhi DTheta]';
input = [z_gf z_gf z_gr z_gr 0 0 0]';

DDZ = (K*Z + C*DZ + kt*input)./M;

%% Solver ODE45

% z1 = x(1) , d-z1 = x(2) , z2 = x(3) , d-z2 = x(4) , z3 = x(5) , d-z3 = x(6) , z4 = x(7) , d-z4 = x(8)
% zs = x(9) , d-zs = x(10) , phi = x(11) , d-phi = x(12) , theta = x(13) , d-theta = x(14)

% tspan = [0 t_final];
tspan = linspace(0,t_final,1000000);
initial = zeros(14,1);
odefun = @(t,x)[x(2);x(4);x(6);x(8);x(10);x(12);x(14);
    (1/Muf) * (-kt * (x(1)-(20 * sin(W*t))) + kf * (x(9)+tf*x(11)-a*x(13)-x(1)) + cf * (x(10)+tf*x(12)-a*x(14)-x(2)));
    (1/Muf) * (-kt * (x(3)-(20 * sin(W*t))) + kf * (x(9)-tf*x(11)-a*x(13)-x(3)) + cf * (x(10)-tf*x(12)-a*x(14)-x(4)));
    (1/Mur) * (-kt * (x(5)-(20 * sin(W*t))) + kr * (x(9)+tr*x(11)+b*x(13)-x(5)) + cr * (x(10)+tr*x(12)+b*x(14)-x(6)));
    (1/Mur) * (-kt * (x(7)-(20 * sin(W*t))) + kr * (x(9)-tr*x(11)+b*x(13)-x(7)) + cr * (x(10)-tr*x(12)+b*x(14)-x(8)));
    (1/M_sprung) * (-kf * (2*x(9)-2*a*x(13)-x(1)-x(3)) -kr * (2*x(9)-2*b*x(13)-x(5)-x(7)) -cf * (2*x(10)-2*a*x(14)-x(2)-x(4)) -cr * (2*x(10)-2*b*x(14)-x(6)-x(8)));
    (1/Ix_s) * (kf*tf * (-2*tf*x(11)+x(1)-x(3)) - kr*tr * (-2*tr*x(11)+x(5)-x(7)) + cf*tf * (-2*tf*x(12)+x(2)-x(4)) - cr*tr * (-2*tr*x(12)+x(6)-x(8)));
    (1/Iy_s) * (kf*a * (2*x(9)-2*a*x(13)-x(1)-x(3)) -kr*b * (2*x(9)+2*b*x(13)-x(5)-x(7)) + cf*a * (2*x(10)-2*a*x(14)-x(2)-x(4)) - cr*b * (2*x(10)+2*b*x(14)-x(6)-x(8)))];
[t, Zsol] = ode45(odefun,tspan,initial);

%% accelerations

a1=(1/Muf) * (-kt * (Zsol(:,1)-(20 * sin(W*t))) + kf * (Zsol(:,9)+tf*Zsol(:,11)-a*Zsol(:,13)-Zsol(:,1)) + cf * (Zsol(:,10)+tf*Zsol(:,12)-a*Zsol(:,14)-Zsol(:,2)));
a2=(1/Muf) * (-kt * (Zsol(:,3)-(20 * sin(W*t))) + kf * (Zsol(:,9)-tf*Zsol(:,11)-a*Zsol(:,13)-Zsol(:,3)) + cf * (Zsol(:,10)-tf*Zsol(:,12)-a*Zsol(:,14)-Zsol(:,4)));
a3=(1/Muf) * (-kt * (Zsol(:,5)-(20 * sin(W*t))) + kr * (Zsol(:,9)+tf*Zsol(:,11)+b*Zsol(:,13)-Zsol(:,5)) + cr * (Zsol(:,10)+tr*Zsol(:,12)+b*Zsol(:,14)-Zsol(:,6)));
a4=(1/Muf) * (-kt * (Zsol(:,7)-(20 * sin(W*t))) + kr * (Zsol(:,9)-tf*Zsol(:,11)+b*Zsol(:,13)-Zsol(:,7)) + cr * (Zsol(:,10)-tr*Zsol(:,12)+b*Zsol(:,14)-Zsol(:,8)));
as=(1/M_sprung) * (-kf * (2*Zsol(:,9)-2*a*Zsol(:,13)-Zsol(:,1)-Zsol(:,3)) -kr * (2*Zsol(:,9)-2*b*Zsol(:,13)-Zsol(:,5)-Zsol(:,7)) -cf * (2*Zsol(:,10)-2*a*Zsol(:,14)-Zsol(:,4)-Zsol(:,2))-cr * (2*Zsol(:,10)-2*b*Zsol(:,14)-Zsol(:,6)-Zsol(:,8)));
aphi=(1/Ix_s) * (kf*tf * (-2*tf*Zsol(:,11)+Zsol(:,1)-Zsol(:,3)) - kr*tr * (-2*tr*Zsol(:,11)+Zsol(:,5)-Zsol(:,7)) + cf*tf * (-2*tf*Zsol(:,12)+Zsol(:,2)-Zsol(:,4)) - cr*tr * (-2*tr*Zsol(:,12)+Zsol(:,6)-Zsol(:,8)));
atheta=(1/Iy_s) * (kf*a * (2*Zsol(:,9)-2*a*Zsol(:,13)-Zsol(:,1)-Zsol(:,3)) -kr*b * (2*Zsol(:,9)+2*b*Zsol(:,13)-Zsol(:,5)-Zsol(:,7)) + cf*a * (2*Zsol(:,10)-2*a*Zsol(:,14)-Zsol(:,2)-Zsol(:,4))- cr*b * (2*Zsol(:,10)+2*b*Zsol(:,14)-Zsol(:,6)-Zsol(:,8)));


%% plot

% figure (1)
% subplot(2,2,1)
% plot(t, a1), title('Vertical Acceleration of m_1') % for ploting z1
% xlabel('t (s)')
% ylabel('ddz_1 (m/s^2)')
% subplot(2,2,2)
% plot(t, a2), title('Vertical Acceleration of m_2') % for ploting z2
% xlabel('t (s)')
% ylabel('ddz_2 (m/s^2)')
% subplot(2,2,3)
% plot(t, a3), title('Vertical Acceleration of m_3') % for ploting z3
% xlabel('t (s)')
% ylabel('ddz_3 (m/s^2)')
% subplot(2,2,4)
% plot(t, a4), title('Vertical Acceleration of m_4') % for ploting z4
% xlabel('t (s)')
% ylabel('ddz_4 (m/s^2)')
% 
% figure (2)
% subplot(3,1,1)
% plot(t, as), title('Vertical Acceleration of m_s') % for ploting zs
% xlabel('t (s)')
% ylabel('ddz_s (m/s^2)')
% subplot(3,1,2)
% plot(t, aphi), title('Angular Acceleration of m_s About X Axis')  % for ploting phi
% xlabel('t (s)')
% ylabel('ddphi (rad/s^2)')
% subplot(3,1,3)
% plot(t, atheta), title('Angular Acceleration of m_s About Y Axis') % for ploting theta
% xlabel('t (s)')
% ylabel('ddtheta (rad/s^2)')
% 
% figure(3)
% subplot(2,2,1)
% plot(t, (Muf*a1)), title('Force Acted bottom of Tire_1') % for ploting z1
% xlabel('t (s)')
% ylabel('F (N)')
% subplot(2,2,2)
% plot(t, (Muf*a2)), title('Force Acted bottom of Tire_2') % for ploting z2
% xlabel('t (s)')
% ylabel('F (N)')
% subplot(2,2,3)
% plot(t, (Mur*a3)), title('Force Acted bottom of Tire_3') % for ploting z3
% xlabel('t (s)')
% ylabel('F (N)')
% subplot(2,2,4)
% plot(t, (Mur*a4)), title('Force Acted bottom of Tire_4') % for ploting z4
% xlabel('t (s)')
% ylabel('F (N)')
% 
% figure(4)
% subplot(2,2,1)
% plot(t, abs(Zsol(:,1)-Zsol(:,9))), title('Relative Displacement of Tire_1 with respect to m_s') % for ploting z1
% xlabel('t (s)')
% ylabel('Z rel (mm)')
% subplot(2,2,2)
% plot(t, abs(Zsol(:,3)-Zsol(:,9))), title('Relative Displacement of Tire_2 with respect to m_s') % for ploting z2
% xlabel('t (s)')
% ylabel('Z rel (mm)')
% subplot(2,2,3)
% plot(t, abs(Zsol(:,5)-Zsol(:,9))), title('Relative Displacement of Tire_3 with respect to m_s') % for ploting z3
% xlabel('t (s)')
% ylabel('Z rel (mm)')
% subplot(2,2,4)
% plot(t, abs(Zsol(:,7)-Zsol(:,9))), title('Relative Displacement of Tire_4 with respect to m_s') % for ploting z4
% xlabel('t (s)')
% ylabel('Z rel (mm)')
% 
