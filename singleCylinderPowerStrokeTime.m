% Reciprocating Piston Single Cylinder Simulation of Power Stroke
clear all;

% constants:
gas_constant = 8.3145; %[J/mol K]
steel_density = 7800; %[kg/m3]
polytropic_idx = 1.3;

%fudge factors:
piston_mass_fudge = 1.; %accounting for the fact that piston is not a solid lump of steel sized about Vd, but including effects of rod/crank weights
ignition_temp_fudge = 0.1;
ignition_pressure_fudge = 0.1; 

% INPUTS:
initial_pressure = 101.325e3; %[Pa, atmospheric pressure]
initial_temperature = 300; %[K] 
step_ignition_temperature = 2000; %[K]

bore = 35; %[mm]
stroke = 26.0; %[mm]
CR = 8; %[compression ratio, []]
displacement = 25; %[cc/cm^3]
crank_radius = 13; %[mm]
connecting_rod_length = 39; %[mm]
%piston_mass = 50; %[g]

% derived params:
Vd = (bore/2)^2*pi*stroke; %[mm^3]
Vc = Vd / (CR-1); %[mm^3]


%converting everything to SI units:
V_initial = mm3tom3(Vc);
V_final = mm3tom3(Vd);
bore = mmtom(bore);
crank_radius = mmtom(crank_radius);
stroke = mmtom(stroke);
connecting_rod_length = mmtom(connecting_rod_length);
n = connecting_rod_length / crank_radius; %[unitless]
area = pi * (bore/2)^2;
compression_ratio = CR;

if exist('piston_mass','var') == 0
    piston_mass = V_final * steel_density * piston_mass_fudge;
else
    piston_mass = gtokg(piston_mass);
end


% engine chemistry parameters:
% compression pressure and temperature:
compression_pressure = initial_pressure * (compression_ratio^polytropic_idx);
compression_temperature = initial_temperature * (compression_ratio ^ (polytropic_idx - 1));
ignition_temperature = compression_temperature + step_ignition_temperature;
ignition_pressure = compression_pressure * (ignition_temperature / compression_temperature);

ignition_temperature = ignition_temp_fudge * ignition_temperature;
ignition_pressure = ignition_pressure_fudge * ignition_pressure;

% simulation:
% first initialise the system, at time t=0:
idx = 1;
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

fprintf("RPM based on stroke time: %f\n",1/(time(end)*2)*60);
rp_torque = torque;
rp_time = time;
rp_lever_arm = lever_arm;
rp_pressure = pressure;
rp_temperature = temperature;
rp_theta = theta;