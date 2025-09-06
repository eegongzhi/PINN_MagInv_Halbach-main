function magnet_discretization(M,rho_out_length,rho_inner_length,theta_length)
% input: non-dimensional rho-theta size of the magnet, and the magnetization
% amplitude M

%% discretization step
%%% consider 2D magnet case: length=? in all directions
num_grid_rho = 10;
num_grid_theta = 320;
% num_grid_z = 20;
h_rho = (rho_out_length-rho_inner_length)/num_grid_rho;  % d_rho
h_theta = theta_length/num_grid_theta; % d_theta
% hz = z_length/num_grid_z; % z=0:h_z:10

% assume M_rho = 1, M_theta = 0 everywhere

%% generate grid for rho-direction derivatives
rho_derivative_grid = cell(num_grid_rho+1, num_grid_theta);
M_rho_grid = rho_derivative_grid;

for j = 1:size(rho_derivative_grid, 2)
    for i = 1:size(rho_derivative_grid, 1)
        % loop every point rho_derivative_grid{i,j} and record coordinate
        rho_derivative_grid{i,j} = [rho_inner_length+(i-1)*h_rho,...
            h_theta/2+(j-1)*h_theta];

        % calculate the magnitude
        mag = 1;

        if j < size(rho_derivative_grid, 2)/20 + 1
            M_rho_grid{i,j} = [mag,0];  % The M vector @ the grid for M_x, e.g. [0 0 1]
        elseif (j < 3*size(rho_derivative_grid, 2)/20 + 1) && (j > 2*size(rho_derivative_grid, 2)/20)
            M_rho_grid{i,j} = [mag,0];
        elseif (j < 5*size(rho_derivative_grid, 2)/20 + 1) && (j > 4*size(rho_derivative_grid, 2)/20)
            M_rho_grid{i,j} = [mag,0];
        elseif (j < 7*size(rho_derivative_grid, 2)/20 + 1) && (j > 6*size(rho_derivative_grid, 2)/20)
            M_rho_grid{i,j} = [mag,0];
        elseif (j < 9*size(rho_derivative_grid, 2)/20 + 1) && (j > 8*size(rho_derivative_grid, 2)/20)
            M_rho_grid{i,j} = [mag,0];
        elseif (j < 11*size(rho_derivative_grid, 2)/20 + 1) && (j > 10*size(rho_derivative_grid, 2)/20)
            M_rho_grid{i,j} = [mag,0];
        elseif (j < 13*size(rho_derivative_grid, 2)/20 + 1) && (j > 12*size(rho_derivative_grid, 2)/20)
            M_rho_grid{i,j} = [mag,0];
        elseif (j < 15*size(rho_derivative_grid, 2)/20 + 1) && (j > 14*size(rho_derivative_grid, 2)/20)
            M_rho_grid{i,j} = [mag,0];    
        elseif (j < 17*size(rho_derivative_grid, 2)/20 + 1) && (j > 16*size(rho_derivative_grid, 2)/20)
            M_rho_grid{i,j} = [mag,0]; 
        elseif (j < 19*size(rho_derivative_grid, 2)/20 + 1) && (j > 18*size(rho_derivative_grid, 2)/20)
            M_rho_grid{i,j} = [mag,0]; 
        else
            M_rho_grid{i,j} = [-mag,0];   % [-mag,0]
        end
    end
end


%% generate grid for theta-direction derivatives
theta_derivative_grid = cell(num_grid_rho, num_grid_theta);  % special treat here due to periodicity
M_theta_grid = theta_derivative_grid;

for j = 1:size(theta_derivative_grid, 2)
    for i = 1:size(theta_derivative_grid, 1)
        % loop every point y_derivative_grid{i,j,k}
        theta_derivative_grid{i,j} = [rho_inner_length+h_rho/2+(i-1)*h_rho,...
            (j-1)*h_theta];

        mag = 1;
        if j < size(theta_derivative_grid, 2)/20 + 1
            M_theta_grid{i,j} = [mag,0];  % The M vector @ the grid for M_x, e.g. [0 0 1]
        elseif (j < 3*size(theta_derivative_grid, 2)/20 + 1) && (j > 2*size(theta_derivative_grid, 2)/20)
            M_theta_grid{i,j} = [mag,0];
        elseif (j < 5*size(theta_derivative_grid, 2)/20 + 1) && (j > 4*size(theta_derivative_grid, 2)/20)
            M_theta_grid{i,j} = [mag,0];
        elseif (j < 7*size(theta_derivative_grid, 2)/20 + 1) && (j > 6*size(theta_derivative_grid, 2)/20)
            M_theta_grid{i,j} = [mag,0];
        elseif (j < 9*size(theta_derivative_grid, 2)/20 + 1) && (j > 8*size(theta_derivative_grid, 2)/20)
            M_theta_grid{i,j} = [mag,0];
        elseif (j < 11*size(theta_derivative_grid, 2)/20 + 1) && (j > 10*size(theta_derivative_grid, 2)/20)
            M_theta_grid{i,j} = [mag,0];
        elseif (j < 13*size(theta_derivative_grid, 2)/20 + 1) && (j > 12*size(theta_derivative_grid, 2)/20)
            M_theta_grid{i,j} = [mag,0];  
        elseif (j < 15*size(theta_derivative_grid, 2)/20 + 1) && (j > 14*size(theta_derivative_grid, 2)/20)
            M_theta_grid{i,j} = [mag,0];    
        elseif (j < 17*size(theta_derivative_grid, 2)/20 + 1) && (j > 16*size(theta_derivative_grid, 2)/20)
            M_theta_grid{i,j} = [mag,0];  
        elseif (j < 19*size(theta_derivative_grid, 2)/20 + 1) && (j > 18*size(theta_derivative_grid, 2)/20)
            M_theta_grid{i,j} = [mag,0];     
        else
            M_theta_grid{i,j} = [-mag,0];   % [-mag,0]
        end
    end
end


%% generate div(M) on the center grid
center_grid = cell(num_grid_rho,num_grid_theta);
M_center = center_grid;

for j = 1:size(center_grid, 2)
    for i = 1:size(center_grid, 1)
        center_grid{i,j} = [rho_inner_length+h_rho/2+(i-1)*h_rho,...
            h_theta/2+(j-1)*h_theta];
        

        dMrdr = (M_rho_grid{i+1,j}(1)+M_rho_grid{i,j}(1))/(2*center_grid{i,j}(1)) - ...
            (M_rho_grid{i+1,j}(1) - M_rho_grid{i,j}(1))/h_rho;
        
        % apply periodicity in theta direction
        if j+1 > num_grid_theta
            temp_M = M_theta_grid{i,1}(2);
        else
            temp_M = M_theta_grid{i,j+1}(2);
        end
        
        dMtdt = -(M_theta_grid{i,j}(2) - temp_M)/...
            (h_theta*center_grid{i,j}(1));
%         dMzdz = -(M_z_grid{i,j,k}(3) - M_z_grid{i,j,k+1}(3))/hz;

        divM = dMrdr + dMtdt;
        M_center{i,j} = divM;
    end
end


%% save grids
clear divM dMrdr dMtdt i j
save magnet_discretization


