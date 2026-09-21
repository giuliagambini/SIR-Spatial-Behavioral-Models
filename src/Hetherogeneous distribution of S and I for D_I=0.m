%% Figure 8 - Extreme case without mobility of the infected (DI = 0)
close all; clc; clear;

% parameters
beta = 35; mu = 1; nu = 1.8; q2 = 2; h = 2;
D_S = 10; 
D_I = 0; % no infected movement

sigma = 1 + nu/mu; 
coeff = [(beta + q2), -beta/sigma, 1]; 
roots_I = roots(coeff);
I_ee = max(roots_I); 
S_ee = 1 - sigma * I_ee; 

% Domain and resolution 
Lx = 100; Ly = 100; 
Nx = 201; Ny = 201; 
dx = Lx / (Nx-1); dy = Ly / (Ny-1);
x_vec = linspace(0, Lx, Nx); y_vec = linspace(0, Ly, Ny);

% fix the seed 
rng(42); 

% Uniform random noise in the interval (-0.01, 0.01)
rho_S = -0.01 + 0.02 * rand(Ny, Nx);
rho_I = -0.01 + 0.02 * rand(Ny, Nx);

% Initial conditions (=paper ones)
S = S_ee * (1 + rho_S);
I = I_ee * (1 + rho_I);

% temporal growth
dt = 0.001; 
t_end = 50; 
n_steps = round(t_end / dt);
iN = [1, 1:Ny-1]; iS = [2:Ny, Ny]; jW = [1, 1:Nx-1]; jE = [2:Nx, Nx];

fprintf('Evaluating... \n');
for step = 1:n_steps
    
    % only S diffuses (Laplaciano solo per S)
    lapS = (S(iN,:) + S(iS,:) + S(:,jW) + S(:,jE) - 4*S) / dx^2;
    
    % PSi
    Psi = (beta * (I.^2) .* S) ./ (1 + q2 * (I.^h));
    
    % reaction part
    ReacS = mu*(1 - S) - Psi;
    ReacI = Psi - (mu + nu) * I;
    
    % sum (Nota: D_I * lapI è stato rimosso del tutto solo per efficienza)
    S = S + dt * (D_S * lapS + ReacS);
    I = I + dt * ReacI; 
end

% Plot 
figure(80); clf; 
set(gcf, 'Position', [100, 200, 1200, 450]); % bigger view

% Left panel: 2D map
subplot(1, 2, 1);
imagesc(x_vec, y_vec, I); 
set(gca, 'YDir', 'normal');
shading interp; colormap parula; colorbar; axis square;
xlabel('x $\rightarrow$', 'Interpreter', 'latex'); 
ylabel('y $\rightarrow$', 'Interpreter', 'latex');

% Right panel: mesh 3D
subplot(1, 2, 2);
mesh(x_vec, y_vec, I); 
colormap parula; axis square;
xlabel('x $\rightarrow$', 'Interpreter', 'latex'); 
ylabel('y $\rightarrow$', 'Interpreter', 'latex');
zlim([0 1.5]); 
view(-37.5, 30); % just for better visualizing

fprintf('=== Figura 8 Generata con Successo! ===\n');
