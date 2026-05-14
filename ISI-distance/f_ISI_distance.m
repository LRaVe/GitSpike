%% ISI-distance computation with auxiliary boundary spikes and plotting 
% Author: Laure WOLFF
% Date: May 2026

function [dist_matrix,I, I_mean] = f_ISI_distance(spikes_trains, tmin, tmax)
% F_ISI_DISTANCE Calculates the ISI-distance between multiple spike trains
%
% Inputs:
%   spikes_trains : cell array containing spike timing vectors
%   tmin, tmax    : time boundaries (minimum and maximum)
%
% Outputs:
%   dist_matrix   : matrix of pairwise ISI-distances
%   I_mean        : average population ISI-distance
%   Generates several plots of the pairwise ISI-distances

    num_trains = length(spikes_trains);
    I = [];
    dist_matrix = zeros(num_trains, num_trains);
    all_t_events = [tmin, tmax];
    pair_data = {}; 
    
    % Edge correction
    spikes = cell(1,num_trains);
    for i = 1:num_trains
        s = unique(spikes_trains{i}); 
        spikes{i} = s(s > tmin & s < tmax); 
    end

    if num_trains >= 2
   compteur = 1;
   num_pairs = (num_trains * (num_trains - 1)) / 2;
   num_cols = 2; % On fixe 2 colonnes pour que ce soit lisible
   num_rows = ceil(num_pairs / num_cols);
   figure('Name', 'Pairwise ISI Distances');
   for i = 1:num_trains
       for j = i+1:num_trains
           compteur = compteur + 1;
           t_all = unique([tmin, spikes{i}, spikes{j}, tmax]);
           all_t_events = [all_t_events, t_all];
           Iij = 0;
           It_list = [];
           for k = 1 : length(t_all)-1
               t_mid = (t_all(k) + t_all(k+1)) / 2;

               %% Correction edge %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               if isempty(spikes{i}) || t_mid < spikes{i}(1)
                   val_x = spikes{i}(1) - tmin; % first interval
               elseif t_mid > spikes{i}(end)
                   val_x = tmax - spikes{i}(end); % last interval
               else
                   idx = find(spikes{i} <= t_mid, 1, 'last');
                   val_x = spikes{i}(idx+1) - spikes{i}(idx); % other interval
               end

               % Pour le train j
               if isempty(spikes{j}) || t_mid < spikes{j}(1)
                   val_y = spikes{j}(1) - tmin;
               elseif t_mid > spikes{j}(end)
                   val_y = tmax - spikes{j}(end);
               else
                   idy = find(spikes{j} <= t_mid, 1, 'last');
                   val_y = spikes{j}(idy+1) - spikes{j}(idy);
               end


               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               I_t = abs(val_x - val_y) / max(val_x, val_y);
               Iij = Iij+I_t * (t_all(k+1) - t_all(k));
               It_list = [It_list, I_t];
           end
           dist_matrix(i,j) = Iij / (tmax - tmin);
           dist_matrix(j,i) = Iij / (tmax - tmin);
           I = [I,(1/(tmax-tmin))*Iij];
           % display (It_list);
           % display (t_all);
           I_plot = [It_list, It_list(end)];

           pair_data{end+1}.t = t_all;
           pair_data{end}.It = It_list;

           subplot(num_rows, num_cols, compteur - 1); 
           stairs(t_all, I_plot, 'LineWidth', 1.5); 
           xlabel('Time');
           ylabel('I_t');
           xlim([0 tmax]);   
           ylim([0 1]);
           title(['Pair ', num2str(i), ' & ', num2str(j)]);  
           subtitle(['Distance: ', num2str(Iij/(tmax-tmin), '%.3f')]);
           grid on;
       end
   end   
       I_mean = mean(I);
       t_global = unique(all_t_events); 
       I_matrix = zeros(length(pair_data), length(t_global)-1);
       
       for p = 1:length(pair_data)
           t_p = pair_data{p}.t;
           It_p = pair_data{p}.It;
           for k = 1:length(t_global)-1
               t_mid = (t_global(k) + t_global(k+1)) / 2;
               idx = find(t_p(1:end-1) <= t_mid, 1, 'last');
               I_matrix(p, k) = It_p(idx);
           end
       end
       I_pop_mean = mean(I_matrix, 1);
       
       %display(I)
       
       figure;
       imagesc(dist_matrix); 
       colorbar;
       title('Matrix of the ISI-distance');
       
       figure;
       stairs(t_global, [I_pop_mean, I_pop_mean(end)]);
       xlabel('Time'); 
       ylabel('Average I_t');
       xlim([0 tmax]);   
       ylim([0 1]);
       title('Evolution of Population Average ISI distance');
       grid on;
    else
       I_mean = 0; % Cas avec moins de 2 trains
    end
end