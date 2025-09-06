function phi = integral(grid)
    %% import dependencies
    % Note: the input "grid" is an observation grid.
    load magnet_discretization.mat
%     surface = hx * hy; % if we use different discretization on different surfaces, we should change this
%     volumn = hx * hy * hz;
    
    %% here,we should initiate test_value_volumn and so on
    test_value_volumn = zeros(size(grid,1), size(grid,2));
    test_value_front_back = test_value_volumn;

    %% integral in volumn (reduced to surface in 2D cases)
    for probe_j = 1:size(grid, 2)
        for probe_i = 1:size(grid, 1)
            % loop all volumn elements
            for vol = 1:numel(center_grid)

                d_theta = grid{probe_i, probe_j}(2) - center_grid{vol}(2);
                surface = h_rho * h_theta * center_grid{vol}(1);

                dis = sqrt(...
                    grid{probe_i, probe_j}(1)^2+center_grid{vol}(1)^2 - ...
                     2*grid{probe_i, probe_j}(1)*center_grid{vol}(1)*cos(d_theta));

                test_value_volumn(probe_i, probe_j) = test_value_volumn(probe_i, probe_j)+...
                    M_center{vol} * surface / dis;
            end
        end
    end
    test_value_volumn = test_value_volumn/(-4*pi);

    %% integral on surface: outer and inner surface, reduced to line in 2D cases
    dis_positive = [];    %%%%%%%%%
    dl_positive = [];    %%%%%%%%%
    for probe_j = 1:size(grid, 2)
        for probe_i = 1:size(grid, 1)
            % rho_derivative_grid, M_rho_grid{1,:}
            % rho_derivative_grid, M_rho_grid{end,:}
            for sur_front_i = 1:size(rho_derivative_grid,2)
                % positive rho
                d_theta = grid{probe_i, probe_j}(2) - rho_derivative_grid{end,sur_front_i}(2);
%                 d_l = h_theta * 0.5*(rho_derivative_grid{end,sur_front_i}(1)+rho_derivative_grid{end-1,sur_front_i}(1));
                d_l = h_theta * rho_derivative_grid{end,sur_front_i}(1);
                dis = sqrt(grid{probe_i, probe_j}(1)^2 + rho_derivative_grid{end,sur_front_i}(1)^2-...
                   2*grid{probe_i, probe_j}(1)*rho_derivative_grid{end,sur_front_i}(1)*cos(d_theta)...
                   );

                test_value_front_back(probe_i,probe_j) = test_value_front_back(probe_i,probe_j) +...
                    (M_rho_grid{end,sur_front_i}(1)) * d_l / dis;
                dis_positive = [dis_positive;dis];
                dl_positive = [dl_positive;d_l];

                % negative rho
                d_theta = grid{probe_i, probe_j}(2) - rho_derivative_grid{1,sur_front_i}(2);
%                 d_l = h_theta * 0.5*(rho_derivative_grid{1,sur_front_i}(1)+rho_derivative_grid{2,sur_front_i}(1));
                d_l = h_theta * rho_derivative_grid{1,sur_front_i}(1);
                dis = sqrt(grid{probe_i, probe_j}(1)^2 + rho_derivative_grid{1,sur_front_i}(1)^2-...
                   2*grid{probe_i, probe_j}(1)*rho_derivative_grid{1,sur_front_i}(1)*cos(d_theta)...
                   );

                test_value_front_back(probe_i,probe_j) = test_value_front_back(probe_i,probe_j) -...
                    (M_rho_grid{1,sur_front_i}(1)) * d_l / dis;
            end

        end
    end
    test_value_front_back =test_value_front_back/(4*pi);
       
    %% sum: Here we only solve phi. 
    test_value = test_value_volumn + test_value_front_back;
%     imagesc(test_value)
    phi = test_value;
end
