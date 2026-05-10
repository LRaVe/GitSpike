clear all;

tmin=0;
tmax=10;

num_trains=3;
spikes_trains=cell(1,num_trains);
spikes_trains{1} = [0 1 2 4 7]; 
spikes_trains{2} = [3 4 6 10];
spikes_trains{3} = [2 5];

% Initialization of these constantes
I=[];
isi_values = cell(1,num_trains);
dist_matrix = zeros(num_trains, num_trains);
all_t_events = [tmin, tmax];
pair_data = {}; 

% Edge correction (soft correction without adding false spikes)
spikes = cell(1,num_trains);
for i = 1:num_trains
    s = unique(spikes_trains{i}); 
    spikes{i} = s(s > tmin & s < tmax); 
end

% To calculate ISI-distance
if num_trains == 1
    I = isi_values{1};
end

if num_trains >= 2
   compteur = 1;
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
                   idx = find(spikes{j} <= t_mid, 1, 'last');
                   val_y = spikes{j}(idx+1) - spikes{j}(idx);
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

           figure(compteur);
           stairs(t_all, I_plot); 
           xlabel('Spikes');
           ylabel('I_t');
           xlim([0 10]);   
           ylim([0 1]);
           title(['Evolution of ISI distance - Pair ', num2str(i), ...
                ' and ', num2str(j)]);  
           subtitle(['ISI-distance: ', num2str(Iij/(tmax-tmin))]);
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
   %display(I_matrix)
   I_pop_mean = mean(I_matrix, 1);
   
   figure(compteur + 1);
   stairs(t_global, [I_pop_mean, I_pop_mean(end)]);
   xlabel('spike'); 
   ylabel('Average It');
   xlim([tmin tmax]); 
   ylim([0 1]);
   title('Evolution of Population Average ISI distance');
   subtitle(['Global ISI-distance: ', num2str(I_mean)]);
end 

display(I)
display(I_mean)

figure(1);
imagesc(dist_matrix);
colorbar;
xticks(1:num_trains); 
yticks(1:num_trains);
xlabel('Spike_trains');
ylabel('Spike_trains');
title('Matrix of the ISI-distance');