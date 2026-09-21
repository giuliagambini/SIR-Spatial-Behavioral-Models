%% Figure 10 - Spatio-parametric map(Ds vs A) 
close all; clc; clear;

% parameters
beta = 35; mu = 1; nu = 1.8; q2 = 2; h = 2; D_I = 1;

% Equilibrium evaluation
sigma = 1 + nu/mu; 
coeff = [(beta + q2), -beta/sigma, 1]; 
I_ee = max(roots(coeff));     
S_ee = 1 - sigma * I_ee; 

% Domain and gradient
Lx = 200; Ly = 200;  
Nx = 401; Ny = 401;  
dx = Lx / (Nx-1); dy = Ly / (Ny-1);

ds_vec = linspace(1, 20, Nx);
a_vec = linspace(0, 10, Ny);
[DS, A_mat] = meshgrid(ds_vec, a_vec); % Mappa bidimensionale dei parametri
A_bar_mat = A_mat ./ S_ee;             % Matrice A_bar locale

[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));

% fix the seed 
rng(42); 

% Uniform random noise in the interval (-0.01, 0.01)
rho_S = -0.01 + 0.02 * rand(Ny, Nx);
rho_I = -0.01 + 0.02 * rand(Ny, Nx);

% Initial conditions (=paper ones)
S = S_ee * (1 + rho_S);
I = I_ee * (1 + rho_I);

% temporal growth
dt = 0.002; 
t_end = 50; 
n_steps = round(t_end / dt);

iN = [1, 1:Ny-1];  iS = [2:Ny, Ny];
jW = [1, 1:Nx-1];  jE = [2:Nx, Nx];

fprintf('=== Calcolo Figura 10 in corso... ===\n');
for step = 1:n_steps
    lapS = (S(iN,:) + S(iS,:) + S(:,jW) + S(:,jE) - 4*S) / dx^2;
    lapI = (I(iN,:) + I(I_ee*0+iS,:) + I(:,jW) + I(:,jE) - 4*I) / dx^2;
    
    % Media di A_bar sulle interfacce dei volumi finiti
    A_bar_E = (A_bar_mat + A_bar_mat(:, jE))/2;
    A_bar_W = (A_bar_mat + A_bar_mat(:, jW))/2;
    A_bar_S = (A_bar_mat + A_bar_mat(iS, :))/2;
    A_bar_N = (A_bar_mat + A_bar_mat(iN, :))/2;
    
    S_E = (S + S(:, jE))/2;  S_W = (S + S(:, jW))/2;
    S_S = (S + S(iS, :))/2;  S_N = (S + S(iN, :))/2;
    dI_E = (I(:, jE) - I) / dx;  dI_W = (I - I(:, jW)) / dx;
    dI_S = (I(iS, :) - I) / dy;  dI_N = (I - I(iN, :)) / dy;
    
    cross_diff = ( (A_bar_E.*S_E.*dI_E - A_bar_W.*S_W.*dI_W)/dx + (A_bar_S.*S_S.*dI_S - A_bar_N.*S_N.*dI_N)/dy );
    
    Psi = (beta * (I.^2) .* S) ./ (1 + q2 * (I.^h));
    ReacS = mu*(1 - S) - Psi;
    ReacI = Psi - (mu + nu) * I;
    
    S = S + dt * (DS .* lapS + cross_diff + ReacS);
    I = I + dt * (D_I * lapI + ReacI);
end

%% Plot
figure(10); clf; set(gcf, 'Position', [100, 200, 1100, 450]);
subplot(1, 2, 1); imagesc(ds_vec, a_vec, S); set(gca, 'YDir', 'normal');
shading interp; colormap parula; colorbar; axis square;
title('Spatial patterns exhibited by S');
xlabel('$D_S \rightarrow$', 'Interpreter', 'latex'); ylabel('$A \rightarrow$', 'Interpreter', 'latex');
caxis([0.25, 0.95]);

subplot(1, 2, 2); imagesc(ds_vec, a_vec, I); set(gca, 'YDir', 'normal');
shading interp; colormap parula; colorbar; axis square;
title('Spatial patterns exhibited by I');
xlabel('$D_S \rightarrow$', 'Interpreter', 'latex'); ylabel('$A \rightarrow$', 'Interpreter', 'latex');
caxis([0.0, 0.6]);