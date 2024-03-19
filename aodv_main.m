% Clear all variables from the workspace
clear all
% Close all open figure windows
close all
% Clear the command window
clc

% Declare global variables
global nodes  powers coordinates  delays

% Nodes
nodes = 50;
% Start node in the network
startnode = 8;
% Destination node in the network
destination = 2;

% Set network filename
network_filename = 'network_data_300_24-03-10_21-09-40.xlsx';
% Load data from Excel file
[~,~,nodesdata] = xlsread(network_filename);

% Choose fitness type between [dist, power, delay, all]
fitness_type = "dis";
% Flag to control whether to print information for all iterations
print_all = true;

% Set the number of iterations to repeat the AODV algorithm 
iter_num = 100;

% Extract information about neighbors, coordinates, delays, and powers
% Get neighbor nodes
neighbors=nodesdata(:,2);
% Get node coordinates
coordinates = nodesdata(:,3);
% Get delays associated with nodes
delays = cell2mat(nodesdata(:,4));
% Get power values associated with nodes
powers = cell2mat(nodesdata(:,5));

capcity =cell2mat(nodesdata(:,6));
% Get capcity values associated with nodes
% Use the given startnode variable as the source node
src_node1 = startnode; 
% Use the given destination variable as the destination node
dst_node = destination; 
% Initialize an empty connectivity matrix of size nodes x nodes
inrange = zeros(nodes); 

% Iterate through the neighbors list
for i = 1:length(neighbors)
        % Check if the neighbors list contains multiple neighbors (denoted by '*')
    if ~isempty(strfind(neighbors{i}, '*'))
        % Extract and convert the list of neighbors to numeric values
        neighbors_list = str2double(strsplit(neighbors{i}, '*'));
        
        % Update the in-range matrix for each neighbor in the list
        for n=1:length(neighbors_list)
            inrange(i, neighbors_list(n)) = 1;
        end
    else
        % Case when there is only one neighbor in the list
        neighbors_list = neighbors{i};
        
        % Update the in-range matrix for the single neighbor
        inrange(i, neighbors_list) = 1;
    end
    
    inrange(i, neighbors_list) = 1;
end

% Initialize an empty matrix for node coordinates
nodeloc = zeros(nodes, 2); % Initialize an empty matrix for node coordinates

% Iterate through the list of coordinates
for i = 1:length(coordinates)
    % Split the coordinate string into numeric values
    coord = str2double(strsplit(coordinates{i}, '*'));
    
    % Update the node location matrix with the coordinates
    nodeloc(i, :) = coord;
end

% Define the size, variable types, and variable names for the table (T)
sz = [iter_num 8];
varTypes = {'double','double','string', 'double','double', 'double','double', 'double'};
varNames = {'Source', 'Destination', 'BestPath', 'ExecutionTime', 'Distance', 'Power', 'Delay', 'Fitness'};
T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);


% Initialize variable for accumulating average execution time
avg_time = 0;

% Repeat the AODV algorithm for 10 iterations
for iter = 1:iter_num
    % Record the start time using tic()
    tic();
    
    % Call the AODV function to find the path and total distance
    [total_dist, path] = aodv(src_node1, dst_node, inrange, nodeloc);
    
    % Record the execution time using toc()
    exec_time = toc();
    
    % Calculate power, delay, and fitness
    [dist, power, delay, fitness_val] = check(cat(2,path,zeros(1,nodes-length(path))), fitness_type);
    
    % Display Information if 'print_all' is true
    if print_all
        % Display iteration counter
        fprintf('Iteration %d\n', iter);
        % Display the best solution (path)
        fprintf('Best Solution: %s\n', num2str(path));
        % Display the fitness value
        fprintf('Fitness: %.2f\n\n', fitness_val);
    end
    
    % Update the result table for the current iteration
    T(iter,:) = {src_node1, dst_node, {sprintf('[%s]', sprintf('%d ', path))}, exec_time, total_dist, power, delay, fitness_val};
    
    % Accumulate the total execution time for calculating the average
    avg_time = avg_time + exec_time;
end

% Calculate average execution time
avg_time = avg_time / iter_num;

% Display the average execution time
fprintf('Average Execution Time: %.2f\n', avg_time);

% Generate a filename based on startnode and destination
results_filename = ['aodv_results_from_', num2str(startnode), '_to_', num2str(destination), '.xlsx'];

% Write the result table 'T' to an Excel file
writetable(T, results_filename, 'FileType', 'spreadsheet')

% Display a message indicating that results have been saved
disp('Results Saved.');


% AODV Protocol Function
% Function to get the shortest path using AODV protocol
function [total_dist, path] = aodv(src_node1, dst_node, inrange, nodeloc)
    % Reset the src_node to the original source node after every iteration
    src_node = src_node1;
    
    % Initialize routing table
    rtngtble = src_node;
    % Initialize temporary table 1
    tble1 = src_node;
    % Initialize temporary table
    tble = src_node;
    
    % Initialize counters
    cnt = 1;
    cnt1 = 1;
    counter = 1;
    dimnsn(cnt) = numel(rtngtble);
    
    % Repeat until the destination node is reached
    while rtngtble~=dst_node
        % Iterate through the temporary table 1
        for ii=1:numel(tble1)
            % Set the source node from the temporary table 1
            src_node = tble1(ii);
            
            % Find neighbors in range for the current source node
            temp = find(inrange(src_node,:));
            
            % Exclude nodes already present in the temporary table
            temp = temp(find(ismember(temp,tble)==0));
            
            % Store source node and its reachable neighbors in a structure
            str{cnt1} = [src_node,temp];
            
            % Update the temporary table with new neighbors
            tble = [tble, temp];
            
            % Increment counter
            cnt1 = cnt1 + 1;
        end
        
        % Seprate nodes which are not present in routing table
        tble1 = tble(find(ismember(tble,rtngtble)==0));
        
        % Update the routing table with nodes from the temporary table
        rtngtble = [rtngtble, tble];
        
        % Remove the repeated node in table
        [any, index] = unique(rtngtble, 'first');
        rtngtble = rtngtble(sort(index));
        
        % Check if destination node is present in the routing table
        if ismember(dst_node, rtngtble)
            % Find the structure cell that has the destination node
            dst_cell = find(cellfun(@equal, str,repmat({dst_node},1,length(str)))); % find out whihch structre cell has destination node
            dst = dst_cell;
            nodtble = dst_node;
            frst_node = dst;
            
            % Reconstruct the AODV path from destination to source
            while frst_node ~= src_node1
                frst_node = str{dst(1)}(1);
                dst = find(cellfun(@equal, str,repmat({frst_node},1,length(str)))); 
                nodtble = [nodtble, frst_node];
            end
            % msgbox('path found')
            
            % Final routing table
            nodtble = fliplr(nodtble);
            
            % Save all AODV paths for each change in vehicle position into a structure
            route{counter} = nodtble; 
            
            % Take out the distance of nodes in routing table from each other
            for ii=1:numel(nodtble)-1
                distnc(ii)=sqrt((nodeloc(nodtble(ii+1),1)-nodeloc(nodtble(ii),1))^2+(nodeloc(nodtble(ii+1),2)-nodeloc(nodtble(ii),2))^2);
            end
            
            % Total distnace from source to destination
            total_dist = sum(distnc);
            % Total Distance between hops in AODV path
            distance(counter) = total_dist;
            % Increment counter
            counter = counter+1;
        end
        
        % Increment counter for the while loop
        cnt=cnt+1;
        dimnsn(cnt)=numel(rtngtble);
        
        % Check if there is only one node in the routing table
        if numel(rtngtble)==1           
            msgbox('1-No Node in range, Execute again')
            return
        end
        
        % Check if the counter exceeds the maximum threshold
        if cnt>=5
            % h8=msgbox('No path found');
            break
        end
    end
    
    % Check if the routing table exists
    if ~exist('nodtble','var')
        errordlg('Transmission Range is less,Kindly enhance it')
        return
    end
    
    % Return the final AODV path
    path = nodtble;
end


% Equal Function
function log = equal(inpu, dst)
    % dst=35;
    if ismember(dst, inpu)
        log = 1;
    else
        log = 0;
    end
end