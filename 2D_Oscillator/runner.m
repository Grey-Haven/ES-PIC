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

% n_bar = 1e6;

% n_bar_ions = n_bar;
% n_bar_eles = n_bar/8;

N_ions = 128*128*10;
N_electrons = N_ions/8;
% N_electrons = 1;

n_bar_ions = N_ions;
n_bar_eles = N_electrons;

w_ele = n_bar_eles / N_electrons;
w_ion = n_bar_ions / N_ions;

ax = 0;
bx = ax + 1;
ay = 0;
by = ay + 1;

Lx = bx - ax;
Ly = by - ay;

% del_x = Lx/4;
% del_y = Ly/4;
del_x = 1;
del_y = 1;
w_pe = sqrt(N_electrons*w_ele/(del_x*del_y)*q_e^2/(eps_0*m_e));

nodes = 64;
Nx = nodes + 1;
Ny = nodes + 1;
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
% dt = .1*min(dx/6,2/w_pe);
dt = 1/(10*w_pe);

T_fin = 1000000 * dt;

x = ax:dx:bx;
y = ay:dy:by;

ion_x = ax + Lx/2 + randn(N_ions,1)*Lx/12;
ion_y = ay + Ly/2 + randn(N_ions,1)*Ly/12;

ele_x = ax + 3*Lx/8 + rand(N_electrons,1)*Lx/4;
ele_y = ay + 3*Ly/8 + rand(N_electrons,1)*Ly/4;

ele_vx = randn(N_electrons,1)*.1;
ele_vy = randn(N_electrons,1)*.1;

qs_ion = q_ion * ones(N_ions,1);
qs_e = q_e * ones(N_electrons,1);

ms_ion = m_ion * ones(N_ions,1);
ms_e = m_e * ones(N_electrons,1);

rho_ion = scatter_charge_2D(ion_x, ion_y, qs_ion, w_ion, x, y, Nx, Ny);
rho_ele = scatter_charge_2D(ele_x, ele_y, qs_e  , w_ele, x, y, Nx, Ny);

A = construct_poisson_matrix(size(rho_ion,1));
A_inv = inv(A);

phi_ion = -poisson_solver(rho_ion,dx,dy,A_inv);
phi_ele = -poisson_solver(rho_ele,dx,dy,A_inv);

phi_full = phi_ion + phi_ele;
% phi_full = phi_ion;
% surf(x,y,phi_full);

electrons = [ele_x'; ele_y'; ele_vx'; ele_vy'; qs_e'; ms_e']';

max_w_pe = max_plasma_frequency(electrons,x,y,w_ele,eps_0);
dt = 1/(500*max_w_pe);

N_steps = floor(T_fin / dt);

H_hist = zeros(N_steps,1);
t_hist = zeros(N_steps,1);
dt_hist = zeros(N_steps,1);

t = 0;

for i = 1:N_steps
    % Velocity half step
    electrons = update_velocity(electrons, x, y, phi_full, dx, dy, dt / 2);
    % Space full step
    electrons = update_location(electrons, dt);
    % Eliminate particles out of bounds
    electrons = remove_particles(electrons,ax,bx,ay,by);
    % Scatter charge
    rho_ele = scatter_charge_2D(electrons(:,1), electrons(:,2), electrons(:,5), w_ele, x, y, Nx, Ny);
    % Compute potential
    phi_ele = -poisson_solver(rho_ele,dx,dy,A_inv);
    phi_full = phi_ion + phi_ele;
    % Velocity half step
    electrons = update_velocity(electrons, x, y, phi_full, dx, dy, dt / 2);
    
    H = compute_hamiltonian(electrons, x, y, phi_full);
    
    t = t+dt;
    
%     max_w_pe = max_plasma_frequency(electrons,x,y,w_ele,eps_0);
%     dt = 1/(50*max_w_pe);
    
    H_hist(i) = H;
    t_hist(i) = t;
    dt_hist(i) = dt;

    if mod(i,100) == 0
%         scatter(ion_x,ion_y,5,'MarkerFaceColor',[1,0,0]);
%         hold on;
%         scatter(electrons(:,1),electrons(:,2),2,'MarkerFaceColor',[0,1,1]);
%         hold off;
%         axis([ax bx ay by]);
%         title("Electron Locations, t = " + num2str(t));
%         drawnow;
%         disp(num2str(i) + " " + num2str(H));
        plot(t_hist(1:i),H_hist(1:i))
        xlabel("t");
        ylabel("H(t)");
        title("Hamiltonian over Time, const dt = " + num2str(dt));
        drawnow;
    end
end
figure;
plot((1:N_steps)*dt,H_hist);
ylabel("H");
xlabel("t");
title("Hamiltonian Over Time");