% Fitness Evaluation Function

function  [ds, pw, dl, f] = check(posx, fitness_type)
    % Declare global variables    
    global nodes  powers coordinates  delays

    % Initialize fitness components
    ds = 0; % Distance
    pw = 0; % Power
    dl = 0; % Delay
    
    try
        % Accumulate distance
        for i=1:nodes-1
            if(posx(i)~=0 & posx(i+1)~=0)
                % Extract coordinates of the current and next node
                listDists= cell2mat(coordinates(posx(i)));
                newlistDist = split(listDists,'*');
                pointXY1=[str2num(strrep(newlistDist{1},',','.')) str2num(strrep(newlistDist{2},',','.'))];

                listDists= cell2mat(coordinates(posx(i+1)));
                newlistDist = split(listDists,'*');
                pointXY2=[str2num(strrep(newlistDist{1},',','.')) str2num(strrep(newlistDist{2},',','.'))];

                % Accumulate Euclidean distance
                ds = ds + norm(pointXY1-pointXY2);
            end
        end

        % Accumulate power and delay for each node
        for i=1:nodes
            if(posx(i)~=0)
                % Accumulate power and delay
                pw = pw + powers(posx(i));
                dl = dl + delays(posx(i));
            end
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