close all
clear

%% data preparation
%%% Real values
rho_out_length_real = 0.01; % unit:m
rho_inner_length_real = 0.005; % unit:m
theta_length_real = 2*pi; % unit:rad
% z_length_real = 0.01; % unit:m

M_real = 1e6; % unit:A/m, this is the "amplitude"
obs_height_rho_real = 0.002; % 0.002 m "above" the magnet

% obs_height2_real = -0.010; % 0.005 m below the magnet
% obs_x_length_real = 0.05; % unit:m
% obs_y_length_real = 0.05; % unit:m

%%% Normalization factors
M_star = 1e6; % unit: A/m, also for H
L_star = 0.01; % unit: m
phi_star = M_star * L_star; % unit: A

%%% Normalized values
rho_out_length = rho_out_length_real / L_star;
rho_inner_length = rho_inner_length_real / L_star;
theta_length = theta_length_real;   % do not normalize this
% y_length = y_length_real / L_star;
% z_length = z_length_real / L_star;
M = M_real ./ M_star;
obs_height_rho = obs_height_rho_real / L_star;
% obs_height2 = obs_height2_real / L_star;
% obs_x_length = obs_x_length_real / L_star;
% obs_y_length = obs_y_length_real / L_star;

%% information on the magnet 
magnet_discretization(M, rho_out_length, rho_inner_length, theta_length);

%% loop all observation surfaces
%%% observation_surface(height, grid_num)
grid_num = 2000;   % only in theta direction

% create the surface
[grid_rho_up, grid_rho_down, grid_theta, delta_h_rho, delta_h_theta] = ...
    observation_surface(rho_out_length, obs_height_rho, grid_num);

%% perform integral on each grid
phi_rho_up = integral(grid_rho_up);
phi_rho_down = integral(grid_rho_down);

phi_theta = integral(grid_theta);

%% obtain H: surface 1
H_rho = -(phi_rho_up - phi_rho_down)/delta_h_rho;
H_theta = -(1/grid_theta{1,1}(1))*diff(phi_theta,1,2)/delta_h_theta;

% Hz_1 = -(phi_1_up - phi_1_down)/delta_h;
% we need to convert back to real values
H_rho_real = H_rho * M_star;
H_theta_real = H_theta * M_star;
% Hz_1_real = Hz_1 * M_star;

%% visiualization: results on surface 1
% x-axis --> theta, range (0, 2*pi)
% double y-axes --> H_rho, H_theta
x = rad2deg(delta_h_theta/2: delta_h_theta: 2*pi);
figure
yyaxis left
plot(x, H_rho)
ylabel('rho-component of H1')
yyaxis right
plot(x, H_theta)
ylabel('theta-component of H1')
xlabel('theta')
xticks([0, 90, 180, 270, 360])
legend('rho-component', 'theta-component')

%% save results: H, and the grid. We need to packup all data
% Note: we do not pack two surfaces together.
H_real = zeros(size(H_rho_real,1),size(H_theta_real,2),2);
H_real(:,:,1) = H_rho_real;
H_real(:,:,2) = H_theta_real;
% H_real(:,:,3) = Hz_1_real;
save H_real H_real

% obs_coord = [x; y; z*ones(1,length(x))];
% save obs_coord obs_coord

%% for easier comparison...
% H_real_reshape = reshape(H_real,grid_num*grid_num,[]);
