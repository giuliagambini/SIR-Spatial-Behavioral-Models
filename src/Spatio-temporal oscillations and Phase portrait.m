%% Figure 11 - Spatio-temporal oscillations (D_I=0,D_S=2,A=7.66)
close all; clc; clear;

% Parameters
beta = 35; mu = 1; nu = 1.8; q2 = 2; h = 2;
D_S = 2; D_I = 0.0; 
A = 7.66; 

sigma = 1 + nu/mu; 
coeff = [(beta + q2), -beta/sigma, 1]; 
roots_I = roots(coeff);
I_ee = max(roots_I); 
S_ee = 1 - sigma * I_ee; 

A_bar = A / S_ee; 

% Domain
Lx = 100; Ly = 100; 
Nx = 201; Ny = 201; 
dx = Lx / (Nx-1); dy = Ly / (Ny-1);
x_vec = linspace(0, Lx, Nx); y_vec = linspace(0, Ly, Ny);


rng(42); 

% Initial point 
S_start = 0.41;  
I_start = 0.201; 

% Random noise
rho_S = -0.01 + 0.02 * rand(Ny, Nx);
rho_I = -0.01 + 0.02 * rand(Ny, Nx);

S = S_start * (1 + rho_S);
I = I_start * (1 + rho_I);

% temporal noise
dt = 0.001; 
t_end = 200; 
n_steps = round(t_end / dt);

iN = [1, 1:Ny-1]; iS = [2:Ny, Ny]; jW = [1, 1:Nx-1]; jE = [2:Nx, Nx];

S_history = zeros(n_steps, 1);
I_history = zeros(n_steps, 1);

mid_Y = round(Ny/2); mid_X = round(Nx/2);

fprintf('Evaluating... ===\n');
for step = 1:n_steps
    lapS = (S(iN,:) + S(iS,:) + S(:,jW) + S(:,jE) - 4*S) / dx^2;
    lapI = (I(iN,:) + I(iS,:) + I(:,jW) + I(:,jE) - 4*I) / dx^2;

    S_E = (S + S(:, jE))/2; S_W = (S + S(:, jW))/2; S_S = (S + S(iS, :))/2; S_N = (S + S(iN, :))/2;
    
    dI_E = (I(:, jE) - I) / dx; dI_W = (I - I(:, jW)) / dx; dI_S = (I(iS, :) - I) / dy; dI_N = (I - I(iN, :)) / dy;
    
    cross_diff = A_bar * ( (S_E.*dI_E - S_W.*dI_W)/dx + (S_S.*dI_S - S_N.*dI_N)/dy );
    
    Psi = (beta * (I.^2) .* S) ./ (1 + q2 * (I.^h));
    
    S = S + dt * (D_S * lapS + cross_diff + mu*(1 - S) - Psi);
    I = I + dt * (D_I * lapI + Psi - (mu + nu) * I);
    
    S_history(step) = S(mid_Y, mid_X);
    I_history(step) = I(mid_Y, mid_X);
end

%% Plot
figure(11); clf; set(gcf, 'Position', [100, 200, 1100, 450]);

subplot(1, 2, 1); 
imagesc(x_vec, y_vec, S); 
set(gca, 'YDir', 'normal'); shading interp; colormap parula; colorbar; axis square;
title('S(x,y) at t = 200');
xlabel('x \rightarrow'); ylabel('y \rightarrow');

subplot(1, 2, 2); 
plot(S_history, I_history, 'LineWidth', 1.2, 'Color', [0 0.4470 0.7410]); 
grid on; axis square;
title('Phase Portrait at center (50,50)');
xlabel('S(50,50) \rightarrow'); ylabel('I(50,50) \rightarrow');

fprintf('=== Fatto! ===\n');
try sound(sin(1:2000)); catch; end