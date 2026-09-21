%% Figure 3 - Visual summary of Spatial Turing patterns 
close all; clc; clear;

% Biological parameters
beta = 35; mu = 1; nu = 1.8;

% Evaluation of the EE
coeff = [beta*(mu+nu), -beta*mu, mu*(mu+nu)];
I_ee = max(roots(coeff));
S_ee = (mu+nu)/(beta*I_ee);

% Creazione della griglia in cui X = D_S e Y = D_I
Lx = 100; Ly = 100; 
Nx = 201; Ny=201;
dx = Lx / (Nx-1); 
dy = Ly / (Ny-1);

% Parameter vectors
ds_vec = linspace(1, 20, Nx);
di_vec = linspace(0.1, 5, Ny);

% Costruiamo le matrici 2D in cui la diffusione varia punto per punto
[DS, DI] = meshgrid(ds_vec, di_vec);
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));

% fix the seed 
rng(42); 

% Uniform random noise in the interval (-0.01, 0.01)
rho_S = -0.01 + 0.02 * rand(Ny, Nx);
rho_I = -0.01 + 0.02 * rand(Ny, Nx);

% Initial conditions (=paper ones)
S = S_ee * (1 + rho_S);
I = I_ee * (1 + rho_I);

% Temporal growth (paper dt=0.001, t_end=2000)
dt = 0.001; % small temporal step 
t_end = 2000; % big enough time to let the pattern grow
n_steps = round(t_end / dt);

% Indeces for Neumann BC
iN = [1, 1:Ny-1];  iS = [2:Ny, Ny];
jW = [1, 1:Nx-1];  jE = [2:Nx, Nx];

fprintf('=== Starting to evaluate ===\n');

for step = 1:n_steps
    % Laplacian evaluation for the diffusion part
    lapS = (S(iN,:) + S(iS,:) + S(:,jW) + S(:,jE) - 4*S) / dx^2;
    lapI = (I(iN,:) + I(iS,:) + I(:,jW) + I(:,jE) - 4*I) / dx^2;
    
    % Reaction part evaluation
    ReacS = mu*(1 - S) - beta * (I.^2) .* S;
    ReacI = beta * (I.^2) .* S - (mu + nu) * I;
    
    % Finally add Reaction + Diffusion
    S = S + dt * (DS .* lapS + ReacS);
    I = I + dt * (DI .* lapI + ReacI);
end

% Plot
figure(30); clf;
set(gcf, 'Position', [100, 200, 1100, 450]);

% Panel (a): Susceptiblel(S)
subplot(1, 2, 1);
imagesc(ds_vec, di_vec, S);
set(gca, 'YDir', 'normal'); % In order to set the 0
shading interp; colormap jet; colorbar; axis square;
title('Spatial patterns exhibited by S');
xlabel('$D_S \rightarrow$', 'Interpreter', 'latex');
ylabel('$D_I \rightarrow$', 'Interpreter', 'latex');
caxis([0.05, 0.55]);

% Panel (b): Infected (I)  
subplot(1, 2, 2);
imagesc(ds_vec, di_vec, I);
set(gca, 'YDir', 'normal');
shading interp; colormap jet; colorbar; axis square;
title('Spatial patterns exhibited by I');
xlabel('$D_S \rightarrow$', 'Interpreter', 'latex');
ylabel('$D_I \rightarrow$', 'Interpreter', 'latex');
caxis([0.1, 1.9]);

fprintf('\n=== Done!! :) ===\n');


