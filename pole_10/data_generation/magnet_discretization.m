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


% %% generate grid for z-direction derivatives
% z_derivative_grid = cell(num_grid_x, num_grid_y, num_grid_z+1);
% M_z_grid = z_derivative_grid;
% for k = 1:size(z_derivative_grid, 3)
%     for j = 1:size(z_derivative_grid, 2)
%         for i = 1:size(z_derivative_grid, 1)
%             % loop every point z_derivative_grid{i,j,k}
%             z_derivative_grid{i,j,k} = [hx/2+(i-1)*hx,...
%                 hy/2+(j-1)*hy,...
%                 (k-1)*hz];
%             % calculate magnitude
%             valuez = M*exp(-((hx/2+(i-1)*hx-x_0)^2/(2*sigma_x^2)+...
%                 (hy/2+(j-1)*hy-y_0)^2/(2*sigma_y^2)));
%             valuey = exp(-((hx/2+(i-1)*hx-x_0)^2/(2*sigma_x^2)+...
%                 ((k-1)*hz-z_0)^2/(2*sigma_z^2)));
%             valuex = exp(-(((k-1)*hz-z_0)^2/(2*sigma_z^2)+...
%                 (hy/2+(j-1)*hy-y_0)^2/(2*sigma_y^2)));
% %             mag = valuez *valuey *valuex;
%             mag = 1;
%             dis_xy = sqrt(((i-1)*hx-x_0)^2 + (hy/2+(j-1)*hy-y_0)^2);
% %             mag = M * (1- dis_xy/(2*d_max));
% 
%             % calculate the theta
%             theta = theta_max * (dis_xy/d_max);
%             % recover x-y-z
%             x_comp = mag*sin(theta)*(hx/2+(i-1)*hx-x_0)/dis_xy;
%             y_comp = mag*sin(theta)*(hy/2+(j-1)*hy-y_0)/dis_xy;
%             z_comp = mag*cos(theta);
% 
%             M_z_grid{i,j,k} = [x_comp,y_comp,z_comp];   % The M vector @ the grid for M_z
%         end
%     end
% end

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

%% visualization: plot M_x_grid on x_derivative_grid
% % generate meshgrid from "x_derivative_grid"
% x_coord = zeros(1, size(x_derivative_grid, 1));
% y_coord = zeros(1, size(x_derivative_grid, 2));
% z_coord = zeros(1, size(x_derivative_grid, 3));
% for idx = 1:length(x_coord)
%     x_coord(idx) = x_derivative_grid{idx,1,1}(1);
% end
% for idy = 1:length(y_coord) 
%     y_coord(idy) = x_derivative_grid{1,idy,1}(2);
% end
% for idz = 1:length(z_coord) 
%     z_coord(idz) = x_derivative_grid{1,1,idz}(3);
% end
% [X,Y,Z] = meshgrid(x_coord, y_coord, z_coord);
% Mx_3D = zeros(size(M_x_grid));
% My_3D = zeros(size(M_x_grid));
% Mz_3D = zeros(size(M_x_grid));
% 
% for idx = 1:size(Mx_3D,1)
%     for idy = 1:size(Mx_3D,2)
%         for idz = 1:size(Mx_3D,3)
%             Mx_3D(idx,idy,idz) = M_x_grid{idx,idy,idz}(1);
%         end
%     end
% end
% for idx = 1:size(My_3D,1)
%     for idy = 1:size(My_3D,2)
%         for idz = 1:size(My_3D,3)
%             My_3D(idx,idy,idz) = M_x_grid{idx,idy,idz}(2);
%         end
%     end
% end
% for idx = 1:size(Mz_3D,1)
%     for idy = 1:size(Mz_3D,2)
%         for idz = 1:size(Mz_3D,3)
%             Mz_3D(idx,idy,idz) = M_x_grid{idx,idy,idz}(3);
%         end
%     end
% end
% figure
% quiver3(X,Y,Z,permute(Mx_3D,[2,1,3]),permute(My_3D,[2,1,3]),permute(Mz_3D,[2,1,3]))
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
