function [grid_rho_up, grid_rho_down, grid_theta, delta_h_rho, delta_h_theta]...
    = observation_surface(rho_start, obs_height_rho, grid_num)
% input: the height of the observation surface w.r.t coordinate system; the x_start,
% x_length, y_start, y_length

% In this function, we should generate all phi grids for H computation,
% including phi_up, phi_down, phi_dx, and phi_dy. The discretization size
% should equal to the grid size of phi_dx or phi_dy.

delta_h_theta = 2*pi/grid_num;
delta_h_rho = 0.01;   % artificial

%% generate a surface for probing: phi_rho_up
    grid_rho_up = cell(1, grid_num);
    for j = 1:size(grid_rho_up, 2)
        for i = 1:size(grid_rho_up, 1)
            grid_rho_up{i,j} = [rho_start+obs_height_rho+0.5*delta_h_rho,...
                0.5*delta_h_theta+(j-1)*delta_h_theta];
            % e.g. test_grid{i,j} = [-20+ (i-1)*(50/200),-20+ (j-1)*(50/200),20];
        end
    end

%% generate a surface for probing: phi_rho_down
    grid_rho_down = cell(1, grid_num);
    for j = 1:size(grid_rho_down, 2)
        for i = 1:size(grid_rho_down, 1)
            grid_rho_down{i,j} = [rho_start+obs_height_rho-0.5*delta_h_rho,...
                0.5*delta_h_theta+(j-1)*delta_h_theta];
            % e.g. test_grid{i,j} = [-20+ (i-1)*(50/200),-20+ (j-1)*(50/200),20];
        end
    end
    
%% generate a surface for probing: phi_theta
    grid_theta = cell(1, grid_num+1);
    for j = 1:size(grid_theta, 2)
        for i = 1:size(grid_theta, 1)
            grid_theta{i,j} = [rho_start + obs_height_rho,...
                (j-1)*delta_h_theta ];  % + delta_h_theta
            % e.g. test_grid{i,j} = [-20+ (i-1)*(50/200),-20+ (j-1)*(50/200),20];
        end
    end
    
% % generate a surface for probing: phi_theta_neg
%     grid_theta_neg = cell(1, grid_num);
%     for j = 1:size(grid_theta_neg, 2)
%         for i = 1:size(grid_theta_neg, 1)
%             grid_theta_neg{i,j} = [rho_start+obs_height_rho,...
%                 (j-1)*delta_h_theta];
%             e.g. test_grid{i,j} = [-20+ (i-1)*(50/200),-20+ (j-1)*(50/200),20];
%         end
%     end

end