% Clear all variables from the workspace
clear all
% Close all open figure windows
close all
% Clear the command window
clc

% Set a seed for reproducibility (optional)
seed = 42;
rng(seed);

% Set ntwork filename
network_filename = 'network_data_300_24-02-27_19-32-50.xlsx';



% Start node in the network
startnode = 8;
% Destination node in the network
destination = 2;

% Number of particles in the swarm
nbr_parti = 50;
% Maximum number of iterations for the PSO algorithm
maxiter = 100;

% Choose fitness type between [dist, power, delay, all]
fitness_type = "dist";
% Flag to control whether to print information for all iterations
print_all = true;

% Execute the PSO algorithm
pso(network_filename, startnode, destination, nbr_parti, maxiter, fitness_type, print_all);


% PSO Main Function
function [execution_time, dist, power, delay, fitness] = pso(network_filename, startnode, destination, nbr_parti, maxiter, fitness_type, print_all)
    % Load data from Excel file
    [dimTable,~,nodesdata] = xlsread(network_filename);

    % Extract information about neighbors, coordinates, delays, and powers
    % Get the number of nodes from the data
    nodes = size(dimTable, 1);
    % Get neighbor nodes
    neighbors = nodesdata(:, 2);
    % Get node coordinates
    coordinates = nodesdata(:, 3);
    % Get delays associated with nodes
    delays = cell2mat(nodesdata(:, 4));
    % Get power values associated with nodes
    powers = cell2mat(nodesdata(:, 5));


    % Define a table (T) to save data
    varNames = {'MaxIter', 'Source', 'Destination', 'NumParticles', 'BestPath', 'ExecutionTime', 'Distance', 'Power', 'Delay', 'Fitness'};
    T = cell(0, numel(varNames));

    % Initialize iteration counter
    n = 1;

    % Initialize variables for global best solution
    gbest = [];
    fgbest = Inf;

    % Set maximum and minimum velocity bounds for particles
    Vmax = 1;
    Vmin = -Vmax;

    % Record the starting time using the `tic` function
    start_time = tic;

    % Initialize particles using a cell array for particle representations
    particles = cell(1, nbr_parti);
    
    for k=1:nbr_parti
        % Initialize particle position
        particles{k}.pos = zeros(1, nodes);
        particles{k}.pos(1) = startnode;
        
        % Initialize particle velocity
        particles{k}.vel = rand(1, nodes) * (Vmax - Vmin) + Vmin;

        % Initialize the particle's personal best position
        particles{k}.pbest = [];
        
        % Initialize the value of the objective function for the personal best
        particles{k}.fpbest = 0;
        
        % Initialize the value of the objective function for the particle
        particles{k}.valeurF = 0;
    end

    
    % PSO Algorithm
    % Set constants for the PSO algorithm
    c1 = 2.05;
    c2 = 2.05;
    w = 0.72;

    % Iterate through the PSO algorithm until the maximum number of iterations is reached
    while n <= maxiter
        % Iterate through particles in the swarm
        for i = 1:nbr_parti
            % Extract the current particle
            particle = particles{i};

            % Iterate through nodes in the particle's path
            for j = 2:nodes
                % Handle the case where the index is out of bounds
                if (j > length(particle.pbest)) || (j > length(gbest))
                    % Update velocity randomly
                    particle.vel(j) = rand() * (Vmax - Vmin) + Vmin;
                else
                    % Update particle velocity using PSO update rule
                    particle.vel(j) = w * particle.vel(j) + c1 * rand() * (particle.pbest(j) - particle.pos(j)) + c2 * rand() * (gbest(j) - particle.pos(j));
                end
                
                % Apply velocity bounds
                if particle.vel(j) > Vmax
                    particle.vel(j) = Vmax;
                elseif particle.vel(j) < Vmin
                    particle.vel(j) = Vmin;
                end
                
                % Update particle position based on the velocity
                s = 1 / (1 + exp(-particle.vel(j) / 3));

                % Extract the list of neighbors for the current node
                listNgbs = neighbors{particle.pos(j - 1)};

                % Check if neighbors are represented as numeric values
                if isnumeric(listNgbs)
                    newnode = listNgbs;
                else
                    % Split and convert neighbors represented as strings
                    newlist = str2double(strsplit(listNgbs, '*'));

                    % Update particle position based on the sigmoid function
                    if ismember(destination, newlist)
                        newnode = newlist(find(newlist == destination, 1));
                    elseif (rand() < s)
                        newnode = newlist(randi(length(newlist)));
                    else
                        newnode = 1;
                    end
                end
                
                % Update the particle's position
                particle.pos(j) = newnode;

                % Check if the destination node is reached
                if newnode == destination
                    particle.pos = [particle.pos, zeros(1, nodes - j)];
                    break;
                end
            end

            % Evaluate fitness for the updated particle position
            [~, ~, ~, particle.valeurF] = evaluate(filter_path(particle.pos), powers, coordinates, delays, fitness_type);
            
            % Update personal best if the fitness improves
            particle.fpbest = particle.valeurF;
            if (particle.valeurF < particles{i}.fpbest) && (particle.pos(1) == startnode) && (sum(particle.pos == destination) == 1)
                particles{i} = particle;
            end

            % Update global best based on fitness and path conditions
            if (particle.fpbest < fgbest) && (particle.pos(1) == startnode) && (sum(particle.pos == destination) == 1)
                gbest = particle.pos;
                fgbest = particle.fpbest;
            end
        end

        % Display iteration counter
        fprintf('Iteration %d\n', n);

        if print_all
            fprintf('Best Solution: %s\n', num2str(filter_path(gbest)));
            fprintf('Fitness: %.2f\n\n', fgbest);
        end

        % Increment iteration counter
        n = n + 1;
    end

    % Record the total execution time
    execution_time = toc(start_time);

    % Display execution time
    fprintf('Execution Time: %.2f\n', execution_time);
    
    % Filter out the path
    gbest = filter_path(gbest);
    
    % Evaluate fitness for the best path
    [dist, power, delay, fitness] = evaluate(gbest, powers, coordinates, delays, fitness_type);

    % Append results to a table
    T = [T; {maxiter, startnode, destination, nbr_parti, strjoin(cellstr(num2str(gbest)), '->'), execution_time, dist, power, delay, fitness}];

    % Saving results to an Excel file
    % Convert cell array to a table with specified variable names
    T = cell2table(T, 'VariableNames', varNames);

    % Generate a filename based on startnode and destination
    results_filename = ['pso_results_from_', num2str(startnode), '_to_', num2str(destination), '.xlsx'];

    % Write the table to an Excel file
    writetable(T, results_filename, 'WriteVariableNames', true);

    % Display a message indicating that results have been saved
    disp('Results Saved.');
end


% Fitness Evaluation Function
function [ds, pw, dl, f] = evaluate(posx, powers, coordinates, delays, fitness_type)
    % Initialize fitness components
    ds = 0; % Distance
    pw = 0; % Power
    dl = 0; % Delay
    
    try
        % Accumulate distance
        for i = 1:length(posx)-1
            % Extract coordinates of the current and next node
            pointXY1 = str2double(strsplit(coordinates{posx(i)}, '*'));
            pointXY2 = str2double(strsplit(coordinates{posx(i + 1)}, '*'));

            % Accumulate Euclidean distance
            ds = ds + norm(pointXY1 - pointXY2);
        end

        % Accumulate power and delay for each node
        for i = 1:length(posx)
            % Accumulate power and delay
            pw = pw + powers(posx(i));
            dl = dl + delays(posx(i));
        end

        % Choose fitness type based on the provided fitness_type parameter
        if fitness_type == "dist"
            f = ds;
        elseif fitness_type == "power"
            f = pw;
        elseif fitness_type == "delay"
            f = dl;
        else
            f = ds + pw + dl;
        end
    catch exception
        f = Inf;
    end
end


% Filter out the path
function [new_path] = filter_path(path)
    new_path = path(path > 0);
end