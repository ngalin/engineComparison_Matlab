% Reciprocating Piston Single Cylinder Simulation of Power Stroke
%clear all;

% constants:
gas_constant = 8.3145; %[J/mol K]


% simulation parameters:
%piston_mass - assuming same mass as vanes of equivalent size:
r = 0.012904;
R = 0.038711;
d = R - r;
vane_angle = deg2rad(40);
volume_vanes = pi*r^2*d + d*(R^2-r^2)*vane_angle;
piston_mass = 7800 * volume_vanes * 2;

%now to calculate bore and stroke given compression ratio, and V_initial:
compression_ratio = 9;
V_initial = mLtom3(3); %[m^3]
V_final = V_initial + V_initial*compression_ratio; %[m3]
bore = sqrt(16*r^2/pi);
crank_radius = 2 * (V_final - V_initial) / (pi*bore^2);
stroke = 2 * crank_radius;
ratio_connecting_rod_length_2_crank_radius = 3; %heywood, p43 
connecting_rod_length = ratio_connecting_rod_length_2_crank_radius * crank_radius; %[m]
n = ratio_connecting_rod_length_2_crank_radius;
area = pi * (bore/2)^2;

% engine chemistry parameters:
ignition_pressure = 7739863.84; %[Pa]
ignition_temperature = 2579.954613; %[K]

% simulation:
% first initialise the system, at time t=0:
idx = 1;
polytropic_idx = 1.3;
pressure(idx) = ignition_pressure;
temperature(idx) = ignition_temperature;
V(idx) = V_initial;
Fp(idx) = pressure(idx) * area;
torque(idx) = 0;
theta(idx) = 0;
lever_arm(idx) = crank_radius * (sind(theta(idx)) + sind(2*theta(idx)) / (2*sqrt((connecting_rod_length/crank_radius)^2 - sind(theta(idx))^2)));
a(idx) = Fp(idx) / piston_mass;
x(idx) = 0;
v(idx) = 0;
work(idx) = 0;
time(idx) = 0;

idx = 2;
del_t = 1e-9; %[s]
for t = del_t:del_t:0.1 %some arbitrary end time, long enough for model to reach condition when theta>180deg.
    delta_speed = a(idx-1)*del_t;
    v(idx) = v(idx-1) + delta_speed;
    delta_x = v(idx-1)*del_t + a(idx-1)*del_t^2 / 2;
    x(idx) = x(idx-1) + delta_x;
    time(idx) = t;
    
    cos_theta = (x(idx)^2 + 2*crank_radius^2 - 2*x(idx)*crank_radius - 2*x(idx)*n*crank_radius + 2*n*crank_radius^2) / (2*crank_radius^2 + 2*n*crank_radius^2 - 2*x(idx)*crank_radius);
    theta(idx) = abs(acosd(cos_theta)); 
    
    %knowing the crank angle allows us to calculate cylinder volume:
    s = crank_radius*cosd(theta(idx)) + sqrt(connecting_rod_length^2 - crank_radius^2*sind(theta(idx))^2);
    V(idx) = V(1) + pi*bore^2/4*(connecting_rod_length + crank_radius - s);
    
    pressure(idx) = pressure(1) * (V(1)/V(idx))^polytropic_idx;
    temperature(idx) = temperature(1)*(V(idx)/V(1))^(1-polytropic_idx);
    Fp(idx) = pressure(idx-1) * area;
    a(idx) = Fp(idx) / piston_mass;

    lever_arm(idx) = crank_radius * (sind(theta(idx)) + sind(2*theta(idx)) / (2*sqrt((connecting_rod_length/crank_radius)^2 - sind(theta(idx))^2)));
    torque(idx) = Fp(idx) * lever_arm(idx);
    
    work(idx) = Fp(idx) * delta_x;
    
    if (theta(idx) > 180) 
        break;
    end
    
    idx = idx+1;
end
    
%remember that this is one stroke out of 4, the other 3 don't have any
%torque or work generation.
fprintf("reciprocating piston - power stroke stats: \n");
fprintf("work done by gases during expansion (ignition) [J]: %f\n", sum(work));
fprintf("stroke time [ms]: %f\n",time(end)*1000);
fprintf("average torque [Nm]: %f\n",mean(torque));

rp_torque = torque;
rp_time = time;
rp_lever_arm = lever_arm;
rp_pressure = pressure;
rp_temperature = temperature;
rp_theta = theta;