%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear all; clc; close all
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Load in values for E, V
E = dlmread('E.txt'); Edat = E;
V = dlmread('V.txt');

% Initiate the first loop
loops(1,:) = E(1,:);
E(1,:) = [];
loop_index = [1,2];
while norm(size(E)) > 2.5   % When deleting values norm([1,2]) = 2.236
    k = 2; % Iterator for loop matrix
    i = 1; % Iterator for iterating through E
    inloop = true; % Boolean conditional 
    while inloop
        if i == max(size(E)) + 1
            i = 1; % Restart the loop if ran through all the points
        end
        node2 = E(i,2); % The next node in question
        if loops(k-1, loop_index(1)) == node2 % If they are connected
            loops(k, loop_index) = E(i,:);  % Set values
            E(i,:) = [];    % Delete entries
            if loops(1, loop_index(2)) == loops(k, loop_index(1))
                inloop = false; % If loop has been completed, break
            elseif size(E)*[1;0] == 1   % Checks for the length of E
                loops(k+1, loop_index) = E(1,:); % Handles the last index
                inloop = false; % Removes itself from the 
            end
            k = k + 1;
        else
            i = i + 1;
        end        
    end
    if norm(size(E)) > 2.5  % When deleting values norm([1,2]) = 2.236
        loop_index = loop_index + 2; % Increase the loop index by 2
        loops(1:end, loop_index) = zeros(max(size(loops)), 2); % Pre-allocate
        loops(1,loop_index) = E(1,:); % Iniate
        E(1,:) = [];    % Delete entries
    end   
end

num_loops = size(loops)*[0;1]/2; % Determine number of loops
fprintf(['The number of loops is ', num2str(num_loops)])
loop_index = [7,8];
for i = 1:num_loops
    dat = loops(1:end, loop_index(1)); % Grab values
    dat(dat == 0) = []; % Delete zero entry
    
    idx = find(dat == min(dat), 1);
    
    % Print tabulated values 
    fprintf(['\n\nLoop Number: ',num2str(i), ', with ', num2str(max(size(dat))), ' unique nodes\n'])
    val_loop = true;
    k = 1;
    looptot = 1;
    while val_loop
        if idx == max(size(dat))
            idx = 1;
        end
        fprintf('%10d ', dat(idx))
        k = k + 1;
        idx = idx + 1;
        looptot = looptot + 1;
        if k == 11
           fprintf('\n') 
           k = 1;
        end
        if looptot == max(size(dat))
           val_loop = false ;
        end
    end
    loop_index = loop_index - 2;
end
%write_to_latex(loops)

loop_index = [7,8];
ax = axes;
ax.ColorOrder = [1 0 0; 0 0 1; 0 1 0; 1 0 1];
ax.LineStyleOrder = {'-','--', '-.', '.'};
hold on
for i = 1:num_loops
    indices1 = loops(1:end, loop_index(1)); indices1(indices1 == 0) = [];
    indices2 = loops(1:end, loop_index(2)); indices2(indices2 == 0) = [];
    
    node1 = V(indices1,:); 
    sz = max(size(node1));
    node1((sz+1):(sz+2),:) = [V(indices2(1),:); V(indices1(1),:)];
    
    plot(node1(1:end, 1), node1(1:end, 2), 'linewidth', 2)
    loop_index = loop_index - 2;
end
xlabel('X-Axis', 'fontsize', 16)
ylabel('Y-Axis', 'fontsize', 16)
title('Edge Connectivity ', 'fontsize', 16)
set(gcf, 'Color', 'w', 'Position', [200 200 800 400]);
%export_fig('big_loops.eps')
xlim([-0.1, 0.7])
ylim([-0.1, 0.1])
%export_fig('small_loops.eps')

function write_to_latex(loops)
    fid = fopen('loop_results', 'w');
    loop_index = [7,8];   
    num_loops = size(loops)*[0;1]/2;
    
    for i = 1:num_loops
        dat = loops(1:end, loop_index(1)); % Grab values
        dat(dat == 0) = []; % Delete zero entry
        idx = find(dat == min(dat), 1);
        
        string = append("Loop Number: ", num2str(i), ', with ', num2str(max(size(dat))), " unique nodes & & & & & & & & & \\ \hline");
        fprintf(fid, '\n %s \n', string);
        val_loop = true;
        k = 1;
        looptot = 1;
        while val_loop
            if idx == max(size(dat))
                idx = 1;
            end           
            if k == 10
                lbreak = "\\";
                fprintf(fid, '%10d %s', dat(idx), lbreak);
                fprintf(fid, '\n');
                k = 1;
            else
                fprintf(fid, '%10d &', dat(idx));    
                k = k + 1;
                idx = idx + 1;
            end
            looptot = looptot + 1;
            if looptot == max(size(dat))
               val_loop = false ;
            end
        end
        loop_index = loop_index - 2;
    end
end