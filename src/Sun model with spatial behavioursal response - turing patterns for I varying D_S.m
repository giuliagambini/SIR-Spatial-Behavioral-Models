%% Figure 6 - Spatial patterns exhibited by I for A=0, h=2, q_2=2.(D_S=5,10,15)
close all; clc; clear;

% parameters and diffusion
beta = 35; mu = 1; nu = 1.8;
q2 = 2; h = 2; % social distancing parameters in the specific case
D_I = 1; 
DS_array = [5, 10, 15]; 

% Correction
sigma = 1 + nu/mu; % 2.8 
coeff = [(beta + q2), -beta/sigma, 1]; 
roots_I = roots(coeff);
I_ee = max(roots_I); % root of the equilibrium
S_ee = 1 - sigma * I_ee; % right value of S at equilibrium(~0.4183) 

% Domain and resolution
Lx = 100; Ly = 100;  
Nx = 201; Ny = 201;  
dx = Lx / (Nx-1); 
dy = Ly / (Ny-1);
x_vec = linspace(0, Lx, Nx);
y_vec = linspace(0, Ly, Ny);

% set seed
rng(42); 

% Uniform random noise in the interval (-0.01, 0.01)
rho_S = -0.01 + 0.02 * rand(Ny, Nx);
rho_I = -0.01 + 0.02 * rand(Ny, Nx);

% Initializing the matrices
S = zeros(Ny, Nx, 3);
I = zeros(Ny, Nx, 3);

% Apply to both S and I
for k = 1:3
    S(:,:,k) = S_ee * (1 + rho_S);
    I(:,:,k) = I_ee * (1 + rho_I);
end

% temporal growth
dt = 0.004;     
t_end = 400;  
n_steps = round(t_end / dt);

iN = [1, 1:Ny-1];  iS = [2:Ny, Ny];
jW = [1, 1:Nx-1];  jE = [2:Nx, Nx];

fprintf('Evaluating...\n');

for step = 1:n_steps
    for k = 1:3
        Sk = S(:,:,k);
        Ik = I(:,:,k);
        D_S = DS_array(k);
        
        % Laplacian
        lapS = (Sk(iN,:) + Sk(iS,:) + Sk(:,jW) + Sk(:,jE) - 4*Sk) / dx^2;
        lapI = (Ik(iN,:) + Ik(iS,:) + Ik(:,jW) + Ik(:,jE) - 4*Ik) / dx^2;
        
        % Local incidence term with social stop psi(I) 
        Psi = (beta * (Ik.^2) .* Sk) ./ (1 + q2 * (Ik.^h));
        
        % Reaction terms 
        ReacS = mu*(1 - Sk) - Psi; 
        ReacI = Psi - (mu + nu) * Ik; 
        
        % sum
        S(:,:,k) = Sk + dt * (D_S * lapS + ReacS);
        I(:,:,k) = Ik + dt * (D_I * lapI + ReacI);
    end
end

%% Plots
figure(); clf;
set(gcf, 'Position', [100, 200, 1400, 400]);

for k = 1:3
    subplot(1, 3, k);
    imagesc(x_vec, y_vec, I(:,:,k));
    set(gca, 'YDir', 'normal'); 
    shading interp; colormap parula; colorbar; axis square;
    title(sprintf('I(x,y) for D_S = %d ', DS_array(k)));
    xlabel('x $\rightarrow$', 'Interpreter', 'latex');
    ylabel('y $\rightarrow$', 'Interpreter', 'latex');
    % caxis([0.0, 0.5]); 
end

fprintf('\n Done! :) \n');