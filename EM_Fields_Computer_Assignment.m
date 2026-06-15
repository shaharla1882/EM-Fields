clear; 
clc;

% ------------------------ Functions -------------------------

% Receives a point in grid
% if the point is on C2 return 1v, if the point on C1 return -1, otherwise,
% return 0.Return C2 points.
function [C2_points, potential_on_boundary] = boundary_check(y, x, h_out, w_out, w_in,...
    start_gap_y, start_gap_x, delta, y_axis_C2, x_axis_C2, C2_points)
    potential_on_boundary = 0;
    % All boundary conditions on C1 - V = -1v
    boundary_Conds_C1 = y == 0 || x == h_out || y == w_out || x == 0;
    % All boundary conditions on C2 - V = 1v
    boundary_Conds_C2 = (y == (start_gap_y + w_in - delta) && start_gap_x <= x && ...
        x <= (start_gap_x + delta)) || ...    % Horizontal line at y=4.9, from x=1.6 to x=1.7
    (x == (start_gap_x + w_in) && (start_gap_y + w_in - delta) <= y ...
        && y <= (start_gap_y + w_in)) || ...    % Vertical line at x=2.4, from y=4.9 to y=5
    (y == (start_gap_y + (w_in)/2) && (start_gap_x + delta) <= x ...
        && x <= (start_gap_x + w_in - delta)) || ...    % Horizontal line at y=4.6, from x=1.7 to x=2.3
    (x == (start_gap_x + w_in - delta) && (start_gap_y + (w_in)/2) <= y ...
        && y <= (start_gap_y + w_in - delta)) || ...    % Vertical line at x=2.3, from y=4.6 to y=4.9
    (y == (start_gap_y + w_in - delta) && (start_gap_x + w_in - delta) <= x ...
        && x <= (start_gap_x + w_in)) || ...    % Horizontal line at y=4.9, from x=2.3 to x=2.4
    (x == (start_gap_x + delta) && (start_gap_y + (w_in)/2) <= y && y ...
        <=  (start_gap_y + w_in - delta)) || ...    % Vertical line at x=1.7, from y=4.6 to y=4.9
    (x == start_gap_x && (start_gap_y + w_in - delta) <= y && y <= ...
        (start_gap_y + w_in)) || ...    % Vertical line at x=1.6, from y=4.9 to y=5
    ((abs(x - y - start_gap_x + start_gap_y + w_in) < 1e-10) && start_gap_x <= x ...
        && x <= (start_gap_x + w_in/2) && (start_gap_y + w_in) <= y && y ...
        <= (start_gap_y + (3/2)*(w_in))) || ...  % line: x = y - 3.4, x from 1.6 to 2, y from 5 to 5.4
    ((abs(x + y - start_gap_y - start_gap_x - 2*w_in) < 1e-10) &&  (start_gap_x ...
        + w_in/2) <= x && x <= (start_gap_x + w_in) ...
        && (start_gap_y + w_in) <= y && y <= (start_gap_y + (1.5)*(w_in)));  % line: x = -y + 7.4, x from 6 to 8, y from 6 to 8
    % Inside C2 = 1v
    is_inside_C2 = inpolygon(y, x, y_axis_C2, x_axis_C2);
    
    if boundary_Conds_C1
        potential_on_boundary = -1;
    elseif boundary_Conds_C2 
        m = (x/delta) + 1;
        n = (y/delta) + 1;
        potential_on_boundary = 1;
        C2_points{end+1} = [int32(m), int32(n)];
    elseif is_inside_C2
        potential_on_boundary = 1;
    end
   
end

% Receives M, N, h_out, w_out, w_in, start_gap_y, start_gap_x, delta.
% Show the potential at each point on grid.
% if the point is on C2 place in (m,n) 1v, if the point on C1 place -1v, otherwise,
% place 0. 
% Returns the potential_matrix, an array of S_area_points, C1_points, C2_points
function [potential_matrix, S_area_points, C1_points, C2_points] = ...
    reset_potential_matrix(M, N, h_out, w_out, w_in, start_gap_y, ...
    start_gap_x, delta, y_axis_C2, x_axis_C2)
    % define potenial matrix MxN of zeros
    potential_matrix = zeros(M, N);
    % Extract importent points as list of lists
    S_area_points = {};
    C1_points = {};
    C2_points = {};
    % m - in y axis (horizontal axis)
    % n - in x axis (vertical axis)
    % Loop over all dots in grid and reset the potential 𝜙^(0)
    for m = 1:M
        x = (m-1) * delta;
        for n= 1:N
            y = (n-1) * delta;
            % if the current point not on the boundaries C1 and C2 calculate
            % the potential with points around it
            [C2_points, currnet_potential] = boundary_check(y,x, h_out, w_out, w_in,...
                start_gap_y, start_gap_x, delta, y_axis_C2, x_axis_C2, C2_points);
            potential_matrix(m, n) = currnet_potential;
            if currnet_potential == 0 
                S_area_points{end+1} =[m, n];
            elseif currnet_potential == -1
                C1_points{end+1} = [m, n];
            end
        end
    end
end

% Receives potential_matrix,delta, w_out, h_out,t, tol, fig_potentials
% Show the potential as a 2-D Map for each tol
function Potential_vs_tol(potential_matrix,delta, w_out, h_out,t, tol, fig_potentials)
    % Build the 2-D Maps of potential
    nexttile(t);
    imagesc(0:delta:w_out, 0:delta:h_out, potential_matrix);
    set(gca, 'YDir', 'normal');
    xlabel('Y [mm]');
    ylabel('X [mm]');
    title(['Tol = ', num2str(tol)]);
    axis equal tight;
    colormap(turbo);
    colorbar;
    movegui("southwest", fig_potentials);
end


% Receives eta, delta, w_out, h_out,c, tol, fig_eta
% Show eta as a 2-D Map for each tol
function eta_vs_tol(eta, tol, delta)
    figure;
    [X, Y] = meshgrid(0:delta:(size(eta,2)-1)*delta, 0:delta:(size(eta,1)-1)*delta); 
    scatter(X(:), Y(:), 20, eta(:), 'filled'); 
    colorbar;
    colormap(turbo);  
    xlabel('Y [mm]'); 
    ylabel('X [mm]');
    title(['\eta [C/(m^2)] (Tol = ', num2str(tol), ')']);
    axis equal tight;
    grid on;
end 

% Receives potential_matrix, delta, w_out, h_out, e, 
% fig_field, tol, M, N, S_area_points, C1_points, C2_points, 
% Shows the Electric field on 2-D map with quiver
function [Q_C1,Q_C2,C_Of_Capacitor]=Electric_Field(potential_matrix, delta, w_out, h_out, e, ...
    fig_field, tol, M, N, S_area_points, C1_points, C2_points)
    % Calculation of E
    Ex = zeros(M, N);
    Ey = zeros(M, N);
    % Loop on all S Area in order to build the Electric field matrix in it
    for m = 2:M-1
        for n = 2:N-1
            if ismember([m, n], cell2mat(S_area_points'), 'rows')
                Ex(m, n) = (potential_matrix(m-1, n)- ...
                    potential_matrix(m+1, n))/((2*delta)/1000);
                Ey(m, n) = (potential_matrix(m, n-1) - ...
                    potential_matrix(m, n+1))/((2*delta)/1000);
            end
        end
    end
    % Amount of arrows to show
    step = 3;
    [yq, xq] = meshgrid(0:delta:w_out, 0:delta:h_out);
    xq_sparse = xq(1:step:end, 1:step:end);
    yq_sparse = yq(1:step:end, 1:step:end);
    Ex_sparse = Ex(1:step:end, 1:step:end);
    Ey_sparse = Ey(1:step:end, 1:step:end);
    
    % Show  E(x,y)
    nexttile(e);
    imagesc(0:delta:w_out, 0:delta:h_out, potential_matrix);
    set(gca, 'YDir', 'normal');
    axis equal tight;
    colormap(turbo);
    colorbar;
    hold on;
    quiver(yq_sparse, xq_sparse, Ey_sparse,Ex_sparse, 0.8, 'k', 'LineWidth', 1);
    title(['Tol = ', num2str(tol)]);
    xlabel('Y [mm]');
    ylabel('X [mm]');
    axis equal tight;
    movegui("southeast", fig_field);
    % Calculate the surface charge density on C1 and C2
    [Q_C1,Q_C2,C_Of_Capacitor]=surface_charge_density(Ex, Ey, C1_points, C2_points, S_area_points, ...
        delta, M, N, tol);
end

% Receives Ex, Ey, C1_points, C2_points, S_area_points, delta, M, N, 
% tol
% Calculate the surface charge density on C1 and C2. 
% inside C2 and outside C1: E = 0 [V/m]
function [Q_C1,Q_C2,C_Of_Capacitor] = surface_charge_density(Ex, Ey, C1_points, C2_points, S_area_points, ...
    delta, M, N, tol)
    % mm to meter
    delta_m = delta/1000;
    % Permittivity of free space [F/m]
    epsilon_0 = (1/(36*pi))*1e-9;   
    % Charge per unit length on C1
    Q_C1 = 0;
    % Charge per unit length on C2
    Q_C2 = 0;
    
    eta = zeros(M, N);

    already_calculated_C1 = {};
    already_calculated_C2 = {};
    
    % Calculate Q per unit on C1 and C2, and eta on C1 and C2
    [Q_C1, eta] = C1_surface_charge_density(Q_C1, eta, C1_points, ...
     already_calculated_C1, S_area_points, delta_m, epsilon_0, Ex, Ey);
    [Q_C2, eta] = C2_surface_charge_density(Q_C2, eta, C2_points, ...
    already_calculated_C2, S_area_points, delta_m, epsilon_0, Ex, Ey);
   
    % Show eta as a 2-D Map for each tol
    eta_vs_tol(eta, tol,delta)
 
    % Calculate Capacitance per unit length 
    % |1V - (-1V)| = 2 [V]
    Delta_V = 2; 
    % Capacitance is always positive [F]
    C_Of_Capacitor = abs(Q_C1) / Delta_V;
end

% Receives Q_C1, eta, C1_points, already_calculated_C1, S_area_points, 
% delta_m, epsilon_0, Ex, Ey
% Returns the Q per unit on C1 and eta updated with C1 points
function [Q_C1, eta] = C1_surface_charge_density(Q_C1, eta, C1_points, ...
    already_calculated_C1, S_area_points, delta_m, epsilon_0, Ex, Ey)
    M = size(eta, 1);
    N = size (eta, 2);
    for place = C1_points
        m=place{1}(1);
        n=place{1}(2);
        if isempty(already_calculated_C1) || ...
                ~any(cellfun(@(v) isequal(v, [m, n]), already_calculated_C1))
            % Left edge (y=1, x)
            % Normal vector points into dielectric: +Y (rightwards from conductor)
            if ismember([m,n+1], cell2mat(S_area_points'), 'rows')
                eta(m, n) =  epsilon_0 * Ey(m, n+1);
                % Charge per unit length on C1
                Q_C1 = Q_C1 + eta(m, n)*delta_m;
            % down edge (y=0, x=1)
            % Normal vector points into dielectric: +X (upwards from conductor)
            elseif ismember([m+1, n], cell2mat(S_area_points'), 'rows')
                eta(m, n) =  epsilon_0 * Ex(m+1, n);
                % Charge per unit length on C1
                Q_C1 = Q_C1 + eta(m, n)*delta_m;
            % right edge (y=w_out, x)
            % Normal vector points into dielectric: -Y (lefttwards from conductor)
            elseif ismember([m, n-1], cell2mat(S_area_points'), 'rows')
                eta(m, n) =  epsilon_0 * (-Ey(m, n-1));
                % Charge per unit length on C1
                Q_C1 = Q_C1 + eta(m, n)*delta_m;
            % up edge (y=w_out, x=h_out)
            % Normal vector points into dielectric: -X (downwards from conductor)   
            elseif ismember([m-1, n], cell2mat(S_area_points'), 'rows')
                eta(m, n) =  epsilon_0 * (-Ex(m-1, n));
                % Charge per unit length on C1
                Q_C1 = Q_C1 + eta(m, n)*delta_m;
            end
            already_calculated_C1{end+1} = [m, n];
        end
    end   
        % right down Corner (y=w_out, x=1)
        % Normal vector points into dielectric: (-Y,+X)
        eta(1,N) = (eta(1, N-1) + eta(2, N))/2;
        Q_C1 = Q_C1 + eta(1,N) * delta_m;
        % Left up Corner (y=1, x=h_out)
        % Normal vector points into dielectric: (Y,-X)
        eta(M,1) = (eta(M-1,1) + eta(M,2))/2;
        Q_C1 = Q_C1 + eta(M,1) * delta_m ;
        % Right up Corner (y=w_out, x=h_out)
        % Normal vector points into dielectric: (-Y,-X)
        eta(M,N) = (eta(M-1,N) + eta(M,N-1))/2;
        Q_C1 = Q_C1 + eta(M,N) * delta_m;
        % Left down Corner (y=1, x=1)
        % Normal vector points into dielectric: (+Y,+X)
        eta(1,1) = (eta(1, 2) + eta(2, 1))/2;
        Q_C1 = Q_C1 + eta(1,1) * delta_m;
end

% Receives Q_C2, eta, C2_points, already_calculated_C2, S_area_points, 
% delta_m, epsilon_0, Ex, Ey
% Returns the Q per unit on C2 and eta updated with C2 points
function [Q_C2, eta] = C2_surface_charge_density(Q_C2, eta, C2_points, ...
    already_calculated_C2, S_area_points, delta_m, epsilon_0, Ex, Ey)
    M = size(eta, 1);
    N = size(eta, 2);
    for place = C2_points
        m = place{1}(1);
        n = place{1}(2);
        % Check if the current point has already been calculated
        if isempty(already_calculated_C2) || ...
            ~any(cellfun(@(v) isequal(v, [m, n]), already_calculated_C2))
            % Store neighbor direction vectors [dx, dy]
            neighbors_dirs = {};
            % Store neighbor E-field components [Ex, Ey]
            E_neighbors = [];      
            % Check 4 optional neighbors
            % Up neighbor so current point is below it - +X
            if ismember([m+1, n], cell2mat(S_area_points'), 'rows')
                % Neighbor up (+X)
                neighbors_dirs{end+1} = [1, 0];          
                E_neighbors(end+1, :) = [Ex(m+1, n), Ey(m+1, n)];
            end
            % Down neighbor so current point is above it - -X
            if ismember([m-1, n], cell2mat(S_area_points'), 'rows')
                % Neighbor down (-X)
                neighbors_dirs{end+1} = [-1, 0];         
                E_neighbors(end+1, :) = [Ex(m-1, n), Ey(m-1, n)];
            end
            % Right neighbor so current point is left to it - +Y
            if ismember([m, n+1], cell2mat(S_area_points'), 'rows')
                % Neighbor right (+Y)
                neighbors_dirs{end+1} = [0, 1];          
                E_neighbors(end+1, :) = [Ex(m, n+1), Ey(m, n+1)];
            end
            % Left neighbor so current point is right to it - -Y
            if ismember([m, n-1], cell2mat(S_area_points'), 'rows')
                % Neighbor left (-Y)
                neighbors_dirs{end+1} = [0, -1];         
                E_neighbors(end+1, :) = [Ex(m, n-1), Ey(m, n-1)];
            end
             % If no neighbors found, set eta = 0 and continue
            if isempty(neighbors_dirs)
                eta(m, n) = 0;
                continue;
            end
            % Sum all neighbor direction vectors to estimate normal vector
            normal_vec = [0, 0];
            for direction = 1:length(neighbors_dirs)
                normal_vec = normal_vec + neighbors_dirs{direction};
            end
            
            % Normalize the normal vector
            norm_mag = norm(normal_vec);
            if norm_mag > 1e-12
                normal_vec = normal_vec/norm_mag;
            end

            % Calculate average E-field of neighbors
            E_avg = mean(E_neighbors, 1); % average [Ex, Ey]
              % Project E-field onto normal vector
            E_normal = dot(E_avg, normal_vec);
            % Compute surface charge density eta = epsilon_0 * |E_normal|
            eta(m, n) = epsilon_0 * E_normal;
            if n > 0.5*N 
                segment_length = delta_m * sqrt(2);
            else
                segment_length = delta_m;
            end
                Q_C2 = Q_C2 + eta(m, n) * segment_length;
                
        end
            % Mark point as calculated
            already_calculated_C2{end+1} = [m, n];
    end
        
    % Remaining points
    if length(already_calculated_C2) ~= length(C2_points) 
        [lia, locb] = ismember(cell2mat(C2_points'), cell2mat(already_calculated_C2'), 'rows');
        left_points =  C2_points(~lia);
        for place = left_points
            m = place{1}(1);
            n = place{1}(2);
            % Upper shape
            if m > 0.5*M
                eta(m, n) = (eta(m+1, n) + eta(m, n-1))/2;
            % Down shape
            else
                eta(m, n) = (eta(m-1, n) + eta(m, n-1))/2;
            end
            Q_C2 = Q_C2 + eta(m, n) * delta_m;
            already_calculated_C2{end+1} = [m, n];
        end
    end 
end

% Receives tol_list, iteration_list
% Create a graph of iteration vs tolarance with given parameters
function iteration_vs_tol(tol_list, iteration_list)
    fig_graph = figure;
    figure(fig_graph);
    % Graph of Iterations To Converge vs Tolerance
    semilogx(tol_list, iteration_list, 'rO--', 'LineWidth', 2);
    xlabel('Tolerance');
    ylabel('Iterations To Converge');
    title('Iterations vs. Tolerance');
    grid on;
    movegui(fig_graph, "northwest")
end

% Receives potential_matrix, tol, S_area_points
% Returns the amount of iterations until potential matrix has converged and
% the potential_matrix after converge
function [counter,potential_matrix] = converge_of_potential(potential_matrix, ...
    tol, S_area_points)
    new_matrix = potential_matrix;
    counter = 0;
    % If to_stop is true- continue calculate 𝜙^(i).
    % Otherwise, to_stop will be false and all the potentials have already been
    % calculated.
    to_stop = true;
    % Save errors
    errors = [];
    while to_stop
        counter = counter + 1;
        for place=1:length(S_area_points)
            % Take the current point (n - X axis, m - Y axis)
            m = S_area_points{place}(1);
            n = S_area_points{place}(2);
            % Take all 4 points potential around (n,m)
            left_v = potential_matrix(m, n-1);
            right_v = potential_matrix(m, n+1);
            up_v = potential_matrix(m+1, n);
            down_v = potential_matrix(m-1, n);
            current_potential = (down_v + up_v + left_v + right_v)/4; 
            new_matrix(m, n) = current_potential; 
        end
        % calculate the difference between the old V matrix and the new one
        diff = new_matrix - potential_matrix;         
        norm_diff = norm(diff, 'F');   
        norm_new = norm(new_matrix, 'F');
        rel_err = norm_diff/norm_new;
        errors(end+1) = rel_err;
        if rel_err < tol        
            to_stop = false;
        end
        potential_matrix = new_matrix;
    end
    error_iteration_vs_tol(errors, counter, tol)
end

function error_iteration_vs_tol(errors, counter, tol)
    figure;
    semilogy(1:counter, errors, '-o','LineWidth',2)
    xlabel('Iteration number')
    ylabel('Relative error (||diff||_F / ||potential||_F)')
    title('Error Convergence of Iteration vs Tolerance=', num2str(tol))
    grid on
end
% Receives M, N, h_out, w_out, w_in, start_gap_y, start_gap_x, delta,
% y_axis_C2, x_axis_C2
% All calculations are in here
function Main_algoritem(M, N, h_out, w_out, w_in, start_gap_y, start_gap_x, ...
    delta, y_axis_C2, x_axis_C2)
    Q_C1_arr = [0,0,0,0,0];
    Q_C2_arr = [0,0,0,0,0];
    C_Of_Capacitor_arr = [0,0,0,0,0];
    % The list of tolerances
    tol_list = [1e-1, 1e-2, 1e-3, 1e-4,1e-5];
    % Will save the amount of iteration to converge in each tol
    iteration_list = [];
    % The potential as a Map
    fig_potentials = figure;
   
    t = tiledlayout(fig_potentials, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, 'Potential [V] Maps for given Tolerances', 'FontSize', 14);
    % The electric feild arrows graph in a map
    fig_field = figure;
    movegui(fig_field, 'southeast');
    e = tiledlayout(fig_field, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(e, 'Electric Field [V/m] for Each Tolerance', 'FontSize', 12);
    % Loop over all tols
    for pos=1:length(tol_list)
        % Reset the potential martix (MxN) according to C1, S and C2 potentials
        % S_area_points represents the points in S area as list of lists
        [potential_matrix, S_area_points, C1_points, C2_points] = ...
        reset_potential_matrix(M, N, h_out, w_out, w_in, start_gap_y, ...
        start_gap_x, delta, y_axis_C2, x_axis_C2);
        % Calculate converge of potential_matrix for each tol
        [counter,potential_matrix] = converge_of_potential(potential_matrix, ...
            tol_list(pos), S_area_points);

        % Add the current amount of iterations to iteration_list 
        iteration_list = [iteration_list, counter];
        % Show the potential as a 2-D Map for each tol
        Potential_vs_tol(potential_matrix,delta, w_out, h_out,t, tol_list(pos), ...
            fig_potentials)
        % Shows the Electric field on 2-D map with quiver
        [Q_C1,Q_C2,C_Of_Capacitor]=Electric_Field(potential_matrix, delta, w_out, h_out, e, ...
        fig_field, tol_list(pos), M, N, S_area_points, C1_points, C2_points);
        Q_C1_arr(pos) = Q_C1;
        Q_C2_arr(pos) = Q_C2;
        C_Of_Capacitor_arr(pos) = C_Of_Capacitor;
    end
    hold off
    % Create a graph of iteration vs tolarance with given parameters
    iteration_vs_tol(tol_list, iteration_list)
 
    % Create a graph of Q_C2_arr and Q_C1_arr vs tolarance with given parameters
    Q_C1_and_Q_C2_vs_tol(tol_list, Q_C1_arr, Q_C2_arr)
    % Create a graph of C_Of_Capacitor_arr vs tolarance with given parameters
    C_Of_Capacitor_arr_vs_tol(tol_list,C_Of_Capacitor_arr)

    
end

% Receives tol_list, Q_C1_arr, Q_C2_arr
% Create a graph of Q_C2_arr and Q_C1_arr vs tolarance with given parameters
function Q_C1_and_Q_C2_vs_tol(tol_list, Q_C1_arr, Q_C2_arr)
    figure;
    semilogx(tol_list, Q_C1_arr, 'bo--', 'LineWidth', 2); 
    hold on;
    semilogx(tol_list, Q_C2_arr, 'ro--', 'LineWidth', 2); 
    grid on;
    xlabel('Tolerance'); 
    ylabel('Q [C/m]');
    title('Q_{C1} and Q_{C2} vs. Tolerance');
    legend('Q_{C1}', 'Q_{C2}', 'Location', 'best');
    set(gca, 'XDir', 'reverse'); 
end

% Receives tol_list,C_Of_Capacitor_arr
% Create a graph of C_Of_Capacitor_arr vs tolarance with given parameters
function C_Of_Capacitor_arr_vs_tol(tol_list, C_Of_Capacitor_arr)
    figure;
    loglog(tol_list, C_Of_Capacitor_arr, 'gs--', 'LineWidth', 2);
    grid on;
    xlabel('Tolerance'); 
    ylabel('Capacitance [F/m]');
    title('Capacitance vs. Tolerance');
    set(gca, 'XDir', 'reverse');
end

% ------------------------ Parameters -------------------------

% C1 frame information
w_out = 10;
h_out = 4;
delta = 0.1;

%C2 frame information
w_in = 0.8; 
w_C2 = 2*w_in;
h_C2 = w_in;
start_gap_y = 4.2;
start_gap_x = 1.6;
% Sample dots in grid where phi(x,y) will be calculated
% for y (m index)
M = h_out/delta + 1;
% for x (n index)
N = w_out/delta + 1;

% ------------------------ Start Main -------------------------

% C1 frame  
[Y, X] = meshgrid(0:delta:w_out, 0:delta:h_out);

% Show Capacitor
fig1 = figure;
figure(fig1);
plot(Y, X, 'k.','MarkerSize', 1);
xlabel('Y [mm]');
ylabel('X [mm]');
title('The Capacitor');
axis equal tight;
movegui(fig1, "northeast")
% Place grid's dots acording to given delta
xticks(0:1:w_out);
yticks(0:1:h_out);
grid on;
hold on

% Build C1 frame Line 
y_axis_C1 = [0,w_out,w_out,0];
x_axis_C1 = [0,0,h_out,h_out];
% Show C1 frame
plot([y_axis_C1 y_axis_C1(1)], [x_axis_C1 x_axis_C1(1)], 'k-', 'LineWidth', 2);

% Build C2 frame with Id's ended with 1 and 5
y_axis_C2 = [...
    start_gap_y + (w_in)/2, ...        % point 1
    start_gap_y + w_in - delta, ...    % point 2
    start_gap_y + w_in - delta, ...    % point 3
    start_gap_y + w_in , ...           % point 4
    start_gap_y + (3/2)*(w_in), ...    % point 5
    start_gap_y + w_in , ...           % point 6
    start_gap_y + w_in - delta, ...    % point 7
    start_gap_y + w_in - delta, ...    % point 8
    start_gap_y + (w_in)/2];           % point 9

x_axis_C2 = [...
    start_gap_x + w_in - delta, ...    % point 1
    start_gap_x + w_in - delta, ...    % point 2
    start_gap_x + w_in , ...           % point 3
    start_gap_x + w_in , ...           % point 4
    start_gap_x + w_in/2, ...          % point 5
    start_gap_x, ...                   % point 6
    start_gap_x, ...                   % point 7
    start_gap_x + delta, ...           % point 8
    start_gap_x + delta];              % point 9

% Show C2 frame
plot([y_axis_C2 y_axis_C2(1)], [x_axis_C2 x_axis_C2(1)], 'k-', 'LineWidth', 2);

% Build --- frame Line 
y_axis = [start_gap_y + w_C2,start_gap_y + w_C2,start_gap_y, start_gap_y];
x_axis = [start_gap_x + w_in, start_gap_x, start_gap_x, start_gap_x + w_in];

% Show --- frame
plot([y_axis y_axis(1)], [x_axis x_axis(1)], 'k--', 'LineWidth', 1);

% Seperate line between two shapes
y_axis = [start_gap_y+w_in,start_gap_y+w_in];
x_axis = [start_gap_x + w_in, start_gap_x];
plot([y_axis y_axis(1)], [x_axis x_axis(1)], 'k--', 'LineWidth', 1);

% Show 'ε,μ' on grid
text(1.1,3.5,'ε,μ', 'FontSize', 14);

% Show 'C1' on grid
text(10.1,0.5,'C1', 'FontSize', 14);

% Show 'C2' on grid
text(5.4,2.2,'C2', 'FontSize', 10);

% Show 'S' on grid
text(8.3,2.5,'S', 'FontSize', 14);
hold off
Main_algoritem(M, N, h_out, w_out, w_in, start_gap_y, start_gap_x, ...
    delta, y_axis_C2, x_axis_C2)



