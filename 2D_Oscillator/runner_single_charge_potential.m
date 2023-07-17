clear;
addpath(genpath([fileparts(pwd), '/../../../Yee/2D/particles/common/poisson/']));

rng(0);
% Physical Parameters
e = 1.602e-19;
AMU = 1.67377e-27;
c = 299792458;

second = 1;
meter = 1;
kg = 1;
coulomb = 1;
ampere = coulomb / second;
newton = kg * meter / second^2;
farad = (second^2 * coulomb^2) / (kg * meter^2);
volt = (kg * meter^2) / (second^3 * ampere);

scale_T = 1e-9 * second;
scale_L = c*scale_T * meter;
scale_M = .0005485799 * AMU * kg;
scale_C = e * coulomb;
scale_A = scale_C/scale_T * ampere;
scale_F = 8.854e-12 * scale_L * farad/meter * farad;
scale_N = 1.25e-6 * scale_A^2 * newton/ampere^2 * newton;
scale_V = scale_M^-1 * scale_L^-2 * scale_T^2 / scale_C^-1 * volt;

T = scale_T * second;
L = scale_L * meter;
M = scale_M * kg;
C = scale_C * coulomb;
A = scale_A * ampere;
F = scale_F * farad;
N = scale_N * newton;
V = scale_V * volt;

m_e = 1;
q_e = -1;

m_ion = 1837 * m_e;
q_ion = -q_e;

kappa = 1;
eps_0 = 1;
mu_0 = 1;

n_bar = 1;

N_ions = 128*128*100;
% N_electrons = N_ions/40;
N_electrons = 1;

w_ele = n_bar / N_electrons;
w_ion = n_bar / N_electrons;

w_pe = sqrt(-N_electrons*n_bar*q_e/(eps_0*m_e));

ax = -.5;
bx = ax + 1;
ay = -.5;
by = ay + 1;

Lx = bx - ax;
Ly = by - ay;

nodes = 64;
Nx = nodes + 1;
Ny = nodes + 1;
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
% dt = min(dx/6,2/w_pe);
dt = .1*dx/6;

T_fin = 100000 * dt;

x = ax:dx:bx;
y = ay:dy:by;

ion_x = ax + Lx/2 + randn(N_ions,1)*Lx/12;
ion_y = ay + Ly/2 + randn(N_ions,1)*Ly/12;

ele_x = ax + 3*Lx/8 + rand(N_electrons,1)*Lx/4;
ele_y = ay + 3*Ly/8 + rand(N_electrons,1)*Ly/4;

ele_vx = randn(N_electrons,1);
ele_vy = randn(N_electrons,1);

qs_ion = q_ion * ones(N_ions,1);
qs_e = q_e * ones(N_electrons,1);

ms_ion = m_ion * ones(N_ions,1);
ms_e = m_e * ones(N_electrons,1);

rho_ion = scatter_charge_2D(ion_x, ion_y, qs_ion, w_ion, x, y, Nx, Ny)/N_ions;
% rho_ele = scatter_charge_2D(ele_x, ele_y, qs_e  , w_ele, x, y, Nx, Ny);
% 
phi_ion = -poisson_solver(rho_ion,dx,dy);
phi = phi_ion;

phi_x = zeros(size(phi));
phi_y = zeros(size(phi));

phi_x(:,1) = (-3*phi(:,1) + 4*phi(:,2) - phi(:,3))/(2*dx);
phi_x(:,2:end-1) = (phi(:,3:end) - phi(:,1:end-2))/(2*dx);
phi_x(:,end) = (3*phi(:,end) - 4*phi(:,end-1) + phi(:,end-2))/(2*dx);


phi_y(1,:) = (-3*phi(1,:) + 4*phi(2,:) - phi(3,:))/(2*dx);
phi_y(2:end-1,:) = (phi(3:end,:) - phi(1:end-2,:))/(2*dx);
phi_y(end,:) = (3*phi(end,:) - 4*phi(end-1,:) + phi(end-2,:))/(2*dx);
% phi_ele = -poisson_solver(rho_ele,dx,dy);

% phi_full = phi_ion + phi_ele;
% phi_full = phi_ion;
% surf(x,y,phi_full);

% electrons = [ele_x'; ele_y'; ele_vx'; ele_vy'; qs_e'; ms_e']';
% e = [.2;.2;0;0;-1;1];

N_steps = floor(T_fin / dt);

H_hist = zeros(N_steps,1);

e_x = .2;
e_y = .2;
e_vx = 0;
e_vy = .1;

e = [[e_x,e_y,e_vx,e_vy,-1,1]];

for i = 1:N_steps
    % Velocity half step
%     e_vx = e_vx + dt/2*q_e/m_e*potential_x(e_x,e_y);
%     e_vy = e_vy + dt/2*q_e/m_e*potential_y(e_x,e_y);
%     e_vx = e_vx + dt/2*q_e/m_e*potential_x_numerical(e_x,e_y,phi_x,x,y);
%     e_vy = e_vy + dt/2*q_e/m_e*potential_y_numerical(e_x,e_y,phi_y,x,y);
    e = update_velocity(e, x, y, phi_ion, dx, dy, dt / 2);
    % Space full step
    e = update_location(e, dt);
%     e_x = e_x + dt*kinetic_vx(e_vx,e_vy);
%     e_y = e_y + dt*kinetic_vy(e_vx,e_vy);
    % Eliminate particles out of bounds
    assert(e_x > ax & e_x < bx & e_y > ay & e_y < by);
    % Velocity half step
    e = update_velocity(e, x, y, phi_ion, dx, dy, dt / 2);
%     e_vx = e_vx + dt/2*q_e/m_e*potential_x(e_x,e_y);
%     e_vy = e_vy + dt/2*q_e/m_e*potential_y(e_x,e_y);
%     e_vx = e_vx + dt/2*q_e/m_e*potential_x_numerical(e_x,e_y,phi_x,x,y);
%     e_vy = e_vy + dt/2*q_e/m_e*potential_y_numerical(e_x,e_y,phi_y,x,y);
    
%     H = hamiltonian(e_x,e_y,e_vx,e_vy);
    H = hamiltonian(e(1),e(2),e(3),e(4));
    H_hist(i) = H;
    
%     if mod(i,10) == 0
%         scatter(0,0,20,'MarkerFaceColor',[1,0,0]);
%         hold on;
%         scatter(e_x,e_y,10,'MarkerFaceColor',[0,1,1]);
%         hold off;
%         axis([ax bx ay by]);
%         title("Electron Locations, t = " + num2str(i*dt));
%         drawnow;
%         disp(num2str(i) + " " + num2str(H));
%     end
end
figure;
plot((1:N_steps)*dt,H_hist);
ylabel("H");
xlabel("t");
title("Hamiltonian Over Time");

function V = potential(x,y)
    V = -1/sqrt(x^2 + y^2);
end
function V_x = potential_x(x,y)
    V_x = x/(x^2+y^2)^(3/2);
end
function V_y = potential_y(x,y)
    V_y = y/(x^2+y^2)^(3/2);
end

function V = potential_numerical(x,y,phi,x_mesh,y_mesh)
    V = gather_field_2D(phi,[x,y],x_mesh,y_mesh);
end
function V_x = potential_x_numerical(x,y,phi_x,x_mesh,y_mesh)
    V_x = -gather_field_2D(phi_x,[x,y],x_mesh,y_mesh);
end
function V_y = potential_y_numerical(x,y,phi_y,x_mesh,y_mesh)
    V_y = -gather_field_2D(phi_y,[x,y],x_mesh,y_mesh);
end

function T = kinetic(vx,vy)
    T = .5*(vx^2 + vy^2);
end
function T_vx = kinetic_vx(vx,vy)
    T_vx = vx;
end
function T_vy = kinetic_vy(vx,vy)
    T_vy = vy;
end


function H = hamiltonian(x,y,vx,vy)
    H = kinetic(vx,vy) + potential(x,y);
end