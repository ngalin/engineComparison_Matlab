% Galin Engine Single Cylinder Simulation of Power Stroke
clear all;

vane_angle = 40; %[deg]
r = 0.012904; %[m]
R = 3*r; %[m]
d = 2*r; %[m]
area = (R-r)*d; %[m^2]
lever_arm = (R+r)/2; %[m]
vane_moment_of_inertia = 0.000623; %[kgm^2]
compression_ratio = 9;
V_initial = mLtom3(3); %[m^3]
V_final = compression_ratio*V_initial; %[m^3]
ignition_pressure = 7739863.84; %[Pa]
ignition_temperature = 2579.954613; %[K]
polytropic_idx = 1.3;

min_angle = V_initial / (R^2 - r^2) / d * 2; %[rad]
max_angle = V_final / (R^2 - r^2) / d * 2; %[rad]
ave_distance_travelled = lever_arm * (max_angle - min_angle); %[m]

% simulation of just the expansion stroke:
% first initialise the system, at time t=0:
idx = 1;
pressure(idx) = ignition_pressure;
temperature(idx) = ignition_temperature;
V(idx) = V_initial;
Fp(idx) = pressure(idx) * area;
torque(idx) = Fp(idx) * lever_arm;

a(idx) = torque(idx) / vane_moment_of_inertia;
x(idx) = min_angle;
v(idx) = 0;
work(idx) = 0;
time(idx) = 0;

idx = 2;
del_t = 1e-9; %[s]
for t = del_t:del_t:0.1 %some arbitrary end time, long enough for model to reach condition when theta>180deg.
    delta_speed = a(idx-1)*del_t;
    v(idx) = v(idx-1) + delta_speed;
    delta_x = (v(idx-1)*del_t + a(idx-1)*del_t^2 / 2) * 2; %both vanes moving away from each other, hence x 2
    x(idx) = x(idx-1) + delta_x;
    time(idx) = t;
    
    V(idx) = d * (R^2 - r^2) * x(idx)/2;
    
    pressure(idx) = pressure(1) * (V(1)/V(idx))^polytropic_idx;
    temperature(idx) = temperature(1)*(V(idx)/V(1))^(1-polytropic_idx);
    Fp(idx) = pressure(idx-1) * area;
    torque(idx) = Fp(idx) * lever_arm;
    a(idx) = torque(idx) / vane_moment_of_inertia;
    
    work(idx) = torque(idx) * delta_x; %in calculating delta_x (above) we already multiplied by factor of 2 to take into account two vanes moving away from each other
    
    if (rad2deg(x(idx)) > 90) 
        break;
    end
    
    idx = idx+1;
end

fprintf("galin engine - power stroke stats: \n");
fprintf("work done by gases during expansion (ignition) [J]: %f\n", sum(work));
fprintf("stroke time [ms]: %f\n",time(end)*1000);
fprintf("average torque [Nm]: %f\n",mean(torque));

ge_torque = torque;
ge_time = time;
ge_lever_arm = lever_arm;
ge_pressure = pressure;
ge_temperature = temperature;
ge_theta = rad2deg(x);