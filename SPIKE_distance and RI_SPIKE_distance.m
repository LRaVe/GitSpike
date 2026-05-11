%% SPIKE-distance computation with auxiliary boundary spikes
% Author: Maxime BELTOISE
% Date: May 2026

% clear all;
% close all;


measures=6;               % +1:ISI,+2:SPIKE,+4:RI-SPIKE,+8:SPIKE-Synchro,+16:SPIKE-order,+32:Spike Train Order
showing=0;                % +1:Spike Trains,+2:Distance,+4:Profile,+8:Matrix
plotting=8;               % +1:Spike Trains,+2:Distance,+4:Profile,+8:Matrix


%% =========================
% INPUT SPIKE TRAINS
% =========================
spikes{1} = [12 16 76 80];
spikes{2} = [8 20 72 84];
spikes{3} = [10 14 84 92];
spikes{4} = [12 44 48 80];
spikes{5} = [8 52 56 84];
spikes{6} = [10 92];

% global time window
t_min = 0;
t_max = 100;





%% =====================================================
% DISPLAY
% =====================================================

[spikes, aux_begin, aux_end] = add_auxiliary_spikes(spikes, t_min, t_max);


if mod(measures,4)>1                                                        % SPIKE-distance
    [SPIKE_distance, SPIKE_distance_profile, SPIKE_distance_matrix] = SPIKE_dist_N(spikes, t_min, t_max, aux_begin, aux_end);
    if mod(showing,4)>1
        disp(SPIKE_distance);
    end
    if mod(showing,8)>3
        disp(SPIKE_distance_profile);
    end
    if mod(showing,16)>7
        disp(SPIKE_distance_matrix);
    end
    if mod(plotting,8)>3
        figure;
        area(SPIKE_distance_profile(:,1), SPIKE_distance_profile(:,2));
        xlabel('Time');
        ylabel('SPIKE distance');
        title(['SPIKE-distance = ', num2str(SPIKE_distance)]);
        xlim([t_min t_max]);
        ylim([0 1]);
        colororder([0.5 0.5 1]);
        grid on;
    end
    if mod(plotting,16)>7
        figure;
        imagesc(SPIKE_distance_matrix);
        colorbar;
        colormap jet;
        set(gca,'XTick',1:length(spikes),'YTick',1:length(spikes));
        xlabel('Spike trains');
        ylabel('Spike trains');
        title(['SPIKE-distance = ', num2str(SPIKE_distance)]);
    end
end


if mod(measures,8)>3                                                        % RI-SPIKE-distance
    [RI_SPIKE_distance, RI_SPIKE_distance_profile, RI_SPIKE_distance_matrix] = RI_SPIKE_dist_N(spikes, t_min, t_max, aux_begin, aux_end);
    if mod(showing,4)>1
        disp(RI_SPIKE_distance);
    end
    if mod(showing,8)>3
        disp(RI_SPIKE_distance_profile);
    end
    if mod(showing,16)>7
        disp(RI_SPIKE_distance_matrix);
    end
    if mod(plotting,8)>3
        figure;
        area(RI_SPIKE_distance_profile(:,1), RI_SPIKE_distance_profile(:,2));
        xlabel('Time');
        ylabel('RI-SPIKE distance');
        title(['RI-SPIKE-distance = ', num2str(RI_SPIKE_distance)]);
        xlim([t_min t_max]);
        ylim([0 1]);
        colororder([0.5 0.5 1]);
        grid on;
    end
    if mod(plotting,16)>7
        figure;
        imagesc(RI_SPIKE_distance_matrix);
        colorbar;
        colormap jet;
        set(gca,'XTick',1:length(spikes),'YTick',1:length(spikes));
        xlabel('Spike trains');
        ylabel('Spike trains');
        title(['RI-SPIKE-distance = ', num2str(RI_SPIKE_distance)]);
    end
end



%% =========================================================
% FUNCTION
% =========================================================


%% =========================================================
% SPIKE-distance between TWO spike trains (2x2 case)
% =========================================================
function [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end)


    %% =====================================================
    % INITIALIZE PROFILE
    %% =====================================================

    SPIKE_distance_profile = {};

    %% =====================================================
    % MAIN LOOP
    %% =====================================================

    for idx_1 = 1 : length(spikes1)

        if spikes2(1)>spikes1(idx_1)  %In case first spikes are not synchronised
            idx_2=1; 
        elseif spikes2(end)<=spikes1(idx_1) %In case last spikes are not synchronised
            idx_2=length(spikes2)-1;
        else
            idx_2 = find(spikes2(1:end) <= spikes1(idx_1), 1, 'last');
        end

        %% ---------------------------------------------
        % Train 2 contribution
        %% ---------------------------------------------
        ISI_dist_2 = spikes2(idx_2+1) - spikes2(idx_2);

        % nearest-neighbor distances
        delta_tp_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_begin);

        delta_tf_2 = auxiliary_delta(spikes2(idx_2+1), spikes2, spikes1, idx_2+1, aux2_end);

        xp_2 = spikes1(idx_1) - spikes2(idx_2);
        xf_2 = spikes2(idx_2+1) - spikes1(idx_1);

        S_2 = ((delta_tp_2*xf_2) + (delta_tf_2*xp_2)) / ISI_dist_2;

        if idx_1>1
            ISI_dist_1 = spikes1(idx_1) - spikes1(idx_1-1);

            S_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_end);

            S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));

            SPIKE_distance_profile{end+1} = [spikes1(idx_1) S];
        end

        if idx_1<length(spikes1)
            ISI_dist_1 = spikes1(idx_1+1) - spikes1(idx_1);

            S_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_begin);

            S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));

            SPIKE_distance_profile{end+1} = [spikes1(idx_1) S];
        end
    end

    for idx_2 = 1 : length(spikes2)
    
        if spikes1(1)>spikes2(idx_2)  %In case first spikes are not synchronised
            idx_1=1;
        elseif spikes1(end)<=spikes2(idx_2) %In case last spikes are not synchronised
            idx_1=length(spikes1)-1;
        else
            idx_1 = find(spikes1(1:end) <= spikes2(idx_2), 1, 'last');
        end
        

        %% ---------------------------------------------
        % Train 1 contribution
        %% ---------------------------------------------
        ISI_dist_1 = spikes1(idx_1+1) - spikes1(idx_1);

        % nearest-neighbor distances
        delta_tp_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_begin);

        delta_tf_1 = auxiliary_delta(spikes1(idx_1+1), spikes1, spikes2, idx_1+1, aux1_end);

        xp_1 = spikes2(idx_2) - spikes1(idx_1);
        xf_1 = spikes1(idx_1+1) - spikes2(idx_2);

        S_1 = ((delta_tp_1*xf_1) + (delta_tf_1*xp_1)) / ISI_dist_1;

        if idx_2>1
            ISI_dist_2 = spikes2(idx_2) - spikes2(idx_2-1);

            S_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_end);

            S = ((S_2*ISI_dist_1) + (S_1*ISI_dist_2)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));

            SPIKE_distance_profile{end+1} = [spikes2(idx_2) S];
        end

        if idx_2<length(spikes2)
            ISI_dist_2 = spikes2(idx_2+1) - spikes2(idx_2);

            S_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_begin);

            S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));

            SPIKE_distance_profile{end+1} = [spikes2(idx_2) S];
        end
    end

    
    %% =====================================================
    % CONVERT PROFILE
    %% =====================================================

    profile_mat = cell2mat(SPIKE_distance_profile');
    profile_mat = sortrows(profile_mat, 1);

    % ==========================================
    % keep only points inside interval
    % ==========================================
    for i = 1:size(profile_mat,1)
    
        % ==========================================
        % point before t_min
        % ==========================================
        if profile_mat(i,1) < t_min
    
            % first index whose abscissa >= t_min
            idx = find(profile_mat(:,1) >= t_min, 1, 'first');

            % linear interpolation
            profile_mat(i,2) = profile_mat(i,2) + ((profile_mat(idx,2) - profile_mat(i,2)) / (profile_mat(idx,1) - profile_mat(i,1))) * (t_min - profile_mat(i,1));
    
            % projection onto boundary
            profile_mat(i,1) = t_min;
    
        % ==========================================
        % point after t_max
        % ==========================================
        elseif profile_mat(i,1) > t_max
    
            % last index whose abscissa <= t_max
            idx = find(profile_mat(:,1) <= t_max, 1, 'last');
    
            % linear interpolation
            profile_mat(i,2) = profile_mat(idx,2) + ((profile_mat(i,2) - profile_mat(idx,2)) / (profile_mat(i,1) - profile_mat(idx,1))) * (t_max - profile_mat(idx,1));
    
            % projection onto boundary
            profile_mat(i,1) = t_max;
    
        end
    end

    % remove duplicates
    [~, idx] = unique(profile_mat,'rows','stable');
    profile_mat = profile_mat(idx,:);

    % sort
    profile_mat = sortrows(profile_mat,1);

    %% =====================================================
    % FINAL DISTANCE
    %% =====================================================

    t = profile_mat(:,1);
    S = profile_mat(:,2);

    SPIKE_distance_2x2 = trapz(t,S)/(t_max-t_min);

end



%% =========================================================
% ADD AUXILIARY SPIKES (for edge correction)
%% =========================================================
function [spikes, aux_begin, aux_end] = add_auxiliary_spikes(spikes, t_min, t_max)
    % Add auxiliary spike times at the beginning and end of the observation window

    % ===============================================
    % ====== Handle cell array of spike trains ======
    % ===============================================

    if iscell(spikes)
        aux_begin = zeros(1, length(spikes));
        aux_end = zeros(1, length(spikes));
        for i = 1:length(spikes)
            [spikes{i}, aux_begin(i), aux_end(i)] = process_single_train(spikes{i}, t_min, t_max);
        end
        return;
    end

    % =======================================
    % ====== Handle single spike train ======
    % =======================================

    [spikes, aux_begin, aux_end] = process_single_train(spikes, t_min, t_max);
end

function [train, aux_begin, aux_end] = process_single_train(train, t_min, t_max)
    aux_begin = 0;
    aux_end = 0;

    train = sort(unique(train));

    % Add auxiliary spike at the beginning if needed
    if train(1) > t_min

        if length(train) >= 2
            aux = train(1) - max(train(1)-t_min, train(2)-train(1));
        else
            aux = t_min;
        end

        train = [aux train];
        aux_begin = 1;
    end

    % Add auxiliary spike at the end if needed
    if train(end) < t_max

        if length(train) >= 2
            aux = train(end) + max(t_max-train(end), train(end)-train(end-1));
        else
            aux = t_max;
        end

        train = [train aux];
        aux_end = 1;
    end
end



%% =========================================================
% AUXILIARY DELTA MANAGEMENT
%
% For auxiliary spikes:
% - beginning auxiliary spike inherits delta from right neighbor
% - ending auxiliary spike inherits delta from left neighbor
%% =========================================================
function delta = auxiliary_delta(spike, own_train, other_train, idx, aux_idx)

    % standard nearest-neighbor distance
    delta_std = min(abs(spike - other_train(:)));
    delta = delta_std;

    % auxiliary at beginning
    if (idx == 1) && aux_idx
        delta = min(abs(own_train(2) - other_train(:)));
    end

    % auxiliary at end
    if (idx == length(own_train)) && aux_idx
        delta = min(abs(own_train(end-1) - other_train(:)));
    end
end



%% =========================================================
% SPIKE-distance between N spike trains
% =========================================================
function [D_global, profile_global, D_matrix] = SPIKE_dist_N(spikes, t_min, t_max, aux_begin, aux_end)

    N = length(spikes);

    D_matrix = zeros(N);

    profiles = {};
    idx_prof = 0;

    %% =====================================================
    % PAIRWISE DISTANCES
    %% =====================================================

    for i = 1:N
        for j = i+1:N

            [d, prof] = SPIKE_dist_2x2(spikes{i}, spikes{j}, t_min, t_max, aux_begin(i), aux_end(i), aux_begin(j), aux_end(j));

            D_matrix(i,j) = d;
            D_matrix(j,i) = d;

            idx_prof = idx_prof + 1;

            profiles{idx_prof} = prof;
        end
    end

    %% =====================================================
    % GLOBAL DISTANCE
    %% =====================================================

    D_global = mean(D_matrix(triu(true(N),1)));

    %% =====================================================
    % ALL TIME COORDINATES
    %% =====================================================

    t_all = [];

    for p = 1:length(profiles)
        t_all = [t_all ; profiles{p}(:,1)];
    end

    t_all = unique(sort(t_all));

    %% =====================================================
    % GLOBAL PROFILE
    %% =====================================================

    profile_global = [];

    for k = 1:length(t_all)

        t = t_all(k);

        vals_left  = [];
        vals_right = [];

        has_discontinuity = false;

        %% -------------------------------------------------
        % scan all pairwise profiles
        %% -------------------------------------------------

        for p = 1:length(profiles)

            P = profiles{p};

            idx = find(P(:,1) == t);

            %% =============================================
            % CASE 1 : discontinuity in this profile
            %% =============================================

            if length(idx) == 2

                has_discontinuity = true;

                vals_left(end+1)  = P(idx(1),2);
                vals_right(end+1) = P(idx(2),2);

            %% =============================================
            % CASE 2 : single point
            %% =============================================

            elseif length(idx) == 1

                vals_left(end+1)  = P(idx,2);
                vals_right(end+1) = P(idx,2);

            %% =============================================
            % CASE 3 : interpolation
            %% =============================================

            else

                idx_before = find(P(:,1) < t, 1, 'last');
                idx_after  = find(P(:,1) > t, 1, 'first');

                if ~isempty(idx_before) && ~isempty(idx_after)

                    t1 = P(idx_before,1);
                    t2 = P(idx_after,1);

                    S1 = P(idx_before,2);
                    S2 = P(idx_after,2);

                    %% avoid division by zero
                    if t2 ~= t1
                        S_interp = S1 + (S2-S1)*(t-t1)/(t2-t1);
                    else
                        S_interp = S1;
                    end

                    vals_left(end+1)  = S_interp;
                    vals_right(end+1) = S_interp;

                end
            end
        end

        %% -------------------------------------------------
        % averaging
        %% -------------------------------------------------

        if has_discontinuity
            profile_global = [ profile_global ; t mean(vals_left) ; t mean(vals_right)];
        else
            profile_global = [profile_global ; t mean(vals_left)];
        end
    end

    % sort
    profile_global = sortrows(profile_global,1);

end



%% =========================================================
% RI-SPIKE-distance between TWO spike trains (2x2 case)
% =========================================================
function [RI_SPIKE_distance_2x2, profile_mat] = RI_SPIKE_dist_2x2(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end)


    %% =====================================================
    % INITIALIZE PROFILE
    %% =====================================================

    SPIKE_distance_profile = {};

    %% =====================================================
    % MAIN LOOP
    %% =====================================================

    for idx_1 = 1 : length(spikes1)

        if spikes2(1)>spikes1(idx_1)  %In case first spikes are not synchronised
            idx_2=1; 
        elseif spikes2(end)<=spikes1(idx_1) %In case last spikes are not synchronised
            idx_2=length(spikes2)-1;
        else
            idx_2 = find(spikes2(1:end) <= spikes1(idx_1), 1, 'last');
        end

        %% ---------------------------------------------
        % Train 2 contribution
        %% ---------------------------------------------
        ISI_dist_2 = spikes2(idx_2+1) - spikes2(idx_2);

        % nearest-neighbor distances
        delta_tp_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_begin);

        delta_tf_2 = auxiliary_delta(spikes2(idx_2+1), spikes2, spikes1, idx_2+1, aux2_end);

        xp_2 = spikes1(idx_1) - spikes2(idx_2);
        xf_2 = spikes2(idx_2+1) - spikes1(idx_1);

        S_2 = ((delta_tp_2*xf_2) + (delta_tf_2*xp_2)) / ISI_dist_2;

        if idx_1>1
            ISI_dist_1 = spikes1(idx_1) - spikes1(idx_1-1);

            S_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_end);

            S = (S_1 + S_2) / (2*mean([ISI_dist_1 ISI_dist_2]));

            SPIKE_distance_profile{end+1} = [spikes1(idx_1) S];
        end

        if idx_1<length(spikes1)
            ISI_dist_1 = spikes1(idx_1+1) - spikes1(idx_1);

            S_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_begin);

            S = (S_1 + S_2) / (2*mean([ISI_dist_1 ISI_dist_2]));
            
            SPIKE_distance_profile{end+1} = [spikes1(idx_1) S];
        end
    end

    for idx_2 = 1 : length(spikes2)
    
        if spikes1(1)>spikes2(idx_2)  %In case first spikes are not synchronised
            idx_1=1;
        elseif spikes1(end)<=spikes2(idx_2) %In case last spikes are not synchronised
            idx_1=length(spikes1)-1;
        else
            idx_1 = find(spikes1(1:end) <= spikes2(idx_2), 1, 'last');
        end
        

        %% ---------------------------------------------
        % Train 1 contribution
        %% ---------------------------------------------
        ISI_dist_1 = spikes1(idx_1+1) - spikes1(idx_1);

        % nearest-neighbor distances
        delta_tp_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_begin);

        delta_tf_1 = auxiliary_delta(spikes1(idx_1+1), spikes1, spikes2, idx_1+1, aux1_end);

        xp_1 = spikes2(idx_2) - spikes1(idx_1);
        xf_1 = spikes1(idx_1+1) - spikes2(idx_2);

        S_1 = ((delta_tp_1*xf_1) + (delta_tf_1*xp_1)) / ISI_dist_1;

        if idx_2>1
            ISI_dist_2 = spikes2(idx_2) - spikes2(idx_2-1);

            S_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_end);

            S = (S_1 + S_2) / (2*mean([ISI_dist_1 ISI_dist_2]));

            SPIKE_distance_profile{end+1} = [spikes2(idx_2) S];
        end

        if idx_2<length(spikes2)
            ISI_dist_2 = spikes2(idx_2+1) - spikes2(idx_2);

            S_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_begin);

            S = (S_1 + S_2) / (2*mean([ISI_dist_1 ISI_dist_2]));

            SPIKE_distance_profile{end+1} = [spikes2(idx_2) S];
        end
    end

    
    %% =====================================================
    % CONVERT PROFILE
    %% =====================================================

    profile_mat = cell2mat(SPIKE_distance_profile');
    profile_mat = sortrows(profile_mat, 1);

    % ==========================================
    % keep only points inside interval
    % ==========================================
    for i = 1:size(profile_mat,1)
    
        % ==========================================
        % point before t_min
        % ==========================================
        if profile_mat(i,1) < t_min
    
            % first index whose abscissa >= t_min
            idx = find(profile_mat(:,1) >= t_min, 1, 'first');

            % linear interpolation
            profile_mat(i,2) = profile_mat(i,2) + ((profile_mat(idx,2) - profile_mat(i,2)) / (profile_mat(idx,1) - profile_mat(i,1))) * (t_min - profile_mat(i,1));
    
            % projection onto boundary
            profile_mat(i,1) = t_min;
    
        % ==========================================
        % point after t_max
        % ==========================================
        elseif profile_mat(i,1) > t_max
    
            % last index whose abscissa <= t_max
            idx = find(profile_mat(:,1) <= t_max, 1, 'last');
    
            % linear interpolation
            profile_mat(i,2) = profile_mat(idx,2) + ((profile_mat(i,2) - profile_mat(idx,2)) / (profile_mat(i,1) - profile_mat(idx,1))) * (t_max - profile_mat(idx,1));
    
            % projection onto boundary
            profile_mat(i,1) = t_max;
    
        end
    end

    % remove duplicates
    [~, idx] = unique(profile_mat,'rows','stable');
    profile_mat = profile_mat(idx,:);

    % sort
    profile_mat = sortrows(profile_mat,1);

    %% =====================================================
    % FINAL DISTANCE
    %% =====================================================

    t = profile_mat(:,1);
    S = profile_mat(:,2);

    RI_SPIKE_distance_2x2 = trapz(t,S)/(t_max-t_min);

end



%% =========================================================
% RI-SPIKE-distance between N spike trains
% =========================================================
function [D_global, profile_global, D_matrix] = RI_SPIKE_dist_N(spikes, t_min, t_max, aux_begin, aux_end)

    N = length(spikes);

    D_matrix = zeros(N);

    profiles = {};
    idx_prof = 0;

    %% =====================================================
    % PAIRWISE DISTANCES
    %% =====================================================

    for i = 1:N
        for j = i+1:N

            [d, prof] = RI_SPIKE_dist_2x2(spikes{i}, spikes{j}, t_min, t_max, aux_begin(i), aux_end(i), aux_begin(j), aux_end(j));

            D_matrix(i,j) = d;
            D_matrix(j,i) = d;

            idx_prof = idx_prof + 1;

            profiles{idx_prof} = prof;
        end
    end

    %% =====================================================
    % GLOBAL DISTANCE
    %% =====================================================

    D_global = mean(D_matrix(triu(true(N),1)));

    %% =====================================================
    % ALL TIME COORDINATES
    %% =====================================================

    t_all = [];

    for p = 1:length(profiles)
        t_all = [t_all ; profiles{p}(:,1)];
    end

    t_all = unique(sort(t_all));

    %% =====================================================
    % GLOBAL PROFILE
    %% =====================================================

    profile_global = [];

    for k = 1:length(t_all)

        t = t_all(k);

        vals_left  = [];
        vals_right = [];

        has_discontinuity = false;

        %% -------------------------------------------------
        % scan all pairwise profiles
        %% -------------------------------------------------

        for p = 1:length(profiles)

            P = profiles{p};

            idx = find(P(:,1) == t);

            %% =============================================
            % CASE 1 : discontinuity in this profile
            %% =============================================

            if length(idx) == 2

                has_discontinuity = true;

                vals_left(end+1)  = P(idx(1),2);
                vals_right(end+1) = P(idx(2),2);

            %% =============================================
            % CASE 2 : single point
            %% =============================================

            elseif length(idx) == 1

                vals_left(end+1)  = P(idx,2);
                vals_right(end+1) = P(idx,2);

            %% =============================================
            % CASE 3 : interpolation
            %% =============================================

            else

                idx_before = find(P(:,1) < t, 1, 'last');
                idx_after  = find(P(:,1) > t, 1, 'first');

                if ~isempty(idx_before) && ~isempty(idx_after)

                    t1 = P(idx_before,1);
                    t2 = P(idx_after,1);

                    S1 = P(idx_before,2);
                    S2 = P(idx_after,2);

                    %% avoid division by zero
                    if t2 ~= t1
                        S_interp = S1 + (S2-S1)*(t-t1)/(t2-t1);
                    else
                        S_interp = S1;
                    end

                    vals_left(end+1)  = S_interp;
                    vals_right(end+1) = S_interp;

                end
            end
        end

        %% -------------------------------------------------
        % averaging
        %% -------------------------------------------------

        if has_discontinuity
            profile_global = [ profile_global ; t mean(vals_left) ; t mean(vals_right)];
        else
            profile_global = [profile_global ; t mean(vals_left)];
        end
    end

    % sort
    profile_global = sortrows(profile_global,1);

end

