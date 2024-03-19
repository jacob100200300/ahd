% Initialize the min and max number of nodes, and the step
min_nodes = 50;
max_nodes = 300;
step = 50;

% Get current date and time
current_datetime = datestr(now, 'yy.mm.dd_HH.MM.ss');

% Construct the foldername for the Excel files
foldername = ['nets_', current_datetime];

% Create the folder to store the Excel files
mkdir(foldername);

% Set a seed for reproducibility
seed = 42;
rng(seed);

% Loop to generate multiple networks
for n = min_nodes:step:max_nodes
    generate_network(n, foldername);
end
disp('Done.');


% Generate a network
function [] = generate_network(n, foldername)
    % This function generates a network of n nodes with random positions,
    % adjacent nodes based on distance threshold, delays, and power consumption. 
    % It then saves the generated information to an Excel file.

    % Initialize the distance threshold, and separator
    distance_threshold = 20;
    separator = '*';

    % Create the filename with current datetime
    filename = [foldername, '/', 'net_', num2str(n), '.xlsx'];
    disp(['Generating ', filename, ' ...']);

    % Initialize minimum and maximum (x, y) positions
    min_position = 0;
    max_position = 100;

    % Initialize minimum and maximum delays
    min_delay = 30;
    max_delay = 100;

    % Initialize minimum and maximum power consumption
    min_power = 1;
    max_power = 100;

    % Initialize minimum and maximum capacity
    min_capacity = 10;
    max_capacity = 150;

    % Initialize adjacency matrix, positions matrix, delays, and power consumption
    adjacency_matrix = zeros(n, n);
    positions_matrix = zeros(n, 2);
    delays = randi([min_delay, max_delay], 1, n);
    power_consumption = randi([min_power, max_power], 1, n);
    capacity = randi([min_capacity, max_capacity], 1, n);

    % Generate random (x, y) positions for each node
    positions_matrix(:, 1) = randi([min_position, max_position], 1, n);
    positions_matrix(:, 2) = randi([min_position, max_position], 1, n);

    % Generate adjacency matrix based on distance threshold
    for i = 1:n
        for j = i+1:n
            % Calculate the Euclidean distance between nodes i and j
            distance = norm(positions_matrix(i, :) - positions_matrix(j, :));

            % Check if the distance is less than the specified threshold
            if (distance < distance_threshold) || (i == 1)
                % If the distance is less than the threshold, nodes i and j are adjacent
                adjacency_matrix(i, j) = j;
                adjacency_matrix(j, i) = i;
            end
        end
    end

    % Create a cell array of adjacent nodes
    adjacent_nodes = cell(n, 1);
    for i = 1:n
        adjacent_nodes{i} = strjoin(cellstr(num2str(adjacency_matrix(i, adjacency_matrix(i, :) > 0)')), separator);
        adjacent_nodes{i} = strrep(adjacent_nodes{i}, ' ', '');
    end

    % Combine x and y coordinates in one column
    node_position = cell(n, 1);
    for i = 1:n
        node_position{i} = [num2str(positions_matrix(i, 1)), separator, num2str(positions_matrix(i, 2))];
        node_position{i} = strrep(node_position{i}, ' ', '');
    end

    % Create a table with network information
    node_data = table((1:n)', adjacent_nodes, node_position, num2cell(delays'), num2cell(power_consumption'), num2cell(capacity'), 'VariableNames', {'Node', 'AdjacentNodes', 'Position', 'Delay', 'PowerConsumption', 'Capacity'});

    % Save the table to an Excel file
    writetable(node_data, filename, 'WriteVariableNames', false);
end