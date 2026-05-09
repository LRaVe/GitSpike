% Test multi-train spike synchronization with corrected algorithm

clear all;


dataset = 2;

if dataset == 1
    t_min=0;
    t_max=10;
    train1 = [0, 1.9, 3.9, 7, 10];
    train2 = [0, 2, 7.1, 9, 10];
    train3 = [0, 2.1, 4.1, 6.9, 10];
    spikes = {train1, train2, train3};  

    fprintf('Train 1: '); disp(train1);
    fprintf('Train 2: '); disp(train2);
    fprintf('Train 3: '); disp(train3);
    fprintf('\n');
elseif dataset == 2
    t_min=0;
    t_max=100;
    spikes=cell(1,2);
    spikes{1} = [12 16 28 32 44 48 60 64 76 80];
    spikes{2} = [8 20 24 36 40 52 56 68 72 84];
    fprintf('Train 1: '); disp(spikes{1});
    fprintf('Train 2: '); disp(spikes{2});

else   
    error('Invalid dataset selection. Please choose 1 or 2.');
end


% Call multi-train function
[C_matrix, C_global] = f_spike_synchro_multi(spikes, t_min, t_max);

fprintf('=== Pairwise Coincidence Matrix ===\n');
disp(C_matrix);

fprintf('\n=== Global SPIKE-Synchronization Index ===\n');
fprintf('C_global: %.4f\n', C_global);
