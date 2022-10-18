% Kylie Hoyt
close all
%% Load EEG Data
addpath('C:\Program Files\MATLAB\R2022a\toolbox\')
load('ErrP_data_HW1.mat', 'trainingEpochs')
load('ErrP_channels.mat')
Cz = find(strcmp({params.chanlocs.labels},'Cz'));
fs = params.fsamp;
win = [floor(-0.2*fs), ceil(0.8*fs)] + fs;
t = -0.2:1/fs:0.8;

%% Grand Average for Each Magnitude at Cz
figure();
hold on
gavg = zeros(win(2)-win(1),1);
plots = []*5;
plt = 1;
for mag = [0, 3, 6, 9, 12]
    nummags = 0;
    for tr = 1:size(trainingEpochs.rotation_data, 3)
        if trainingEpochs.magnitude(tr) == mag
            nummags = nummags+1;
            gavg = gavg + trainingEpochs.rotation_data(win(1):win(2)-1, Cz, tr);
        end
    end
    gavg = gavg./nummags;
    plots(plt) = plot(t, gavg,"DisplayName",int2str(mag));
    plt = plt + 1;
end
xlabel("Time (s)")
ylabel("Amplitude")
ERN_avg = mean([0.286, 0.312, 0.268, 0.282]);
Pe_avg = mean([0.427, 0.493, 0.408, 0.392]);
xline(ERN_avg, '--r', "Mean ERN");
xline(Pe_avg, '--b', "Mean Pe");
xline(0, '--', "Cue");
title("Grand Average ErrP at Cz for Different Magnitude Perturbations")
legend(plots)


%% Topological Plots w/o Spatial Filter
ERN = floor(ERN_avg*fs+fs);
GAv_ERN = mean(trainingEpochs.rotation_data(ERN, :, find(trainingEpochs.label == 1)), 3);
figure();
sgtitle("Grand Average Amplitude w/o Spatial Filter")
subplot(1, 2, 1);
subtitle("ERN")
topoplot(GAv_ERN, params.chanlocs,'maplimits', [-3 3]);
cbar('vert',0,[-3 3]);
Pe = floor(Pe_avg*fs+fs);
GAv_Pe = mean(trainingEpochs.rotation_data(Pe, :, find(trainingEpochs.label == 1)), 3);
subplot(1, 2, 2);
subtitle("Pe")
topoplot(GAv_Pe, params.chanlocs,'maplimits', [-3 3]);
cbar('vert',0,[-3 3]);


%% CAR Grand Average at Cz
CAR_sig = zeros(size(trainingEpochs.rotation_data));
for tr = 1:size(trainingEpochs.rotation_data, 3)
    CAR_sig(:, :, tr) = CAR(squeeze(trainingEpochs.rotation_data(:, :, tr)));
end

figure();
hold on
gavg_CAR = zeros(win(2)-win(1),1);
plots = []*5;
plt = 1;
for mag = [0, 3, 6, 9, 12]
    nummags = 0;
    for tr = 1:size(trainingEpochs.rotation_data, 3)
        if trainingEpochs.magnitude(tr) == mag
            nummags = nummags+1;
            gavg_CAR = gavg_CAR + CAR_sig(win(1):win(2)-1, Cz, tr);
        end
    end
    gavg_CAR = gavg_CAR./nummags;
    plots(plt) = plot(t, gavg_CAR,"DisplayName",int2str(mag));
    plt = plt + 1;
end
xlabel("Time (s)")
ylabel("Amplitude")
ERN_avg_CAR = mean([0.294, 0.308, 0.265, 0.302]);
Pe_avg_CAR = mean([0.430, 0.413, 0.378, 0.392]);
xline(ERN_avg_CAR, '--r', "Mean ERN");
xline(Pe_avg_CAR, '--b', "Mean Pe");
xline(0, '--', "Cue");
title("Grand Average ErrP at Cz for Different Magnitude Perturbations w/ CAR")
legend(plots)

%% Topological Plots w/ CAR
ERN_CAR = floor(ERN_avg_CAR*fs+fs);
GAv_ERN_CAR = mean(CAR_sig(ERN_CAR, :, trainingEpochs.label == 1), 3);
figure();
sgtitle("Grand Average Amplitude w/ CAR")
subplot(1, 2, 1);
subtitle("ERN")
topoplot(GAv_ERN_CAR, params.chanlocs, 'conv', 'on', 'maplimits', [-3 3]);
cbar('vert',0,[-3 3]);
Pe_CAR = floor(Pe_avg_CAR*fs+fs);
GAv_Pe_CAR = mean(CAR_sig(Pe, :, trainingEpochs.label == 1), 3);
subplot(1, 2, 2);
subtitle("Pe")
topoplot(GAv_Pe_CAR, params.chanlocs,'conv', 'on', 'maplimits', [-3 3]);
cbar('vert',0,[-3 3]);


%% CCA
% Single Trials
Err_Y = []; Corr_Y = [];
for tr = 1:size(trainingEpochs.rotation_data, 3)
    if trainingEpochs.label(tr) == 1
        Err_Y = cat(2, Err_Y, squeeze(trainingEpochs.rotation_data(win(1):win(2), :, tr)).');
    else
        Corr_Y = cat(2, Corr_Y, squeeze(trainingEpochs.rotation_data(win(1):win(2), :, tr)).');
    end
end
Y = cat(2, Err_Y, Corr_Y);
% Grand Average
Err_X = repmat(mean(trainingEpochs.rotation_data(win(1):win(2), :, find(trainingEpochs.label == 1)), 3).', 1, sum(trainingEpochs.label));
Corr_X = repmat(mean(trainingEpochs.rotation_data(win(1):win(2), :, find(trainingEpochs.label == 0)), 3).', 1, size(trainingEpochs.label, 1) - sum(trainingEpochs.label));
X = cat(2, Err_X, Corr_X);
[Wy,  r] = canoncorr(Y.', X.');
Wy = Wy.';
figure();
sgtitle("CCA Top Five Spatial Filters")
subplot(1, 5, 1);
subtitle("1")
topoplot(Wy(1, :), params.chanlocs, 'conv', 'on');
cbar('vert',0,[-1 1]*max(abs(Wy(1, :))));
subplot(1, 5, 2);
subtitle("2")
topoplot(Wy(2, :), params.chanlocs,'conv', 'on');
cbar('vert',0,[-1 1]*max(abs(Wy(2, :))));
subplot(1, 5, 3);
subtitle("3")
topoplot(Wy(3, :), params.chanlocs,'conv', 'on');
cbar('vert',0,[-1 1]*max(abs(Wy(3, :))));
subplot(1, 5, 4);
subtitle("4")
topoplot(Wy(4, :), params.chanlocs,'conv', 'on');
cbar('vert',0,[-1 1]*max(abs(Wy(4, :))));
subplot(1, 5, 5);
subtitle("5")
topoplot(Wy(5, :), params.chanlocs,'conv', 'on');
cbar('vert',0,[-1 1]*max(abs(Wy(5, :))));

function CAR_sig = CAR(sig)
    [ch, dim] = min(size(sig));
    CAR_sig = sig - mean(sig, dim);
end
