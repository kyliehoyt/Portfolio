% Kylie Hoyt 
close all
% Load EEG Data from session
addpath('C:\Program Files\MATLAB\R2022a\toolbox\')
load('selectedChannels.mat')
[s1,h1] = sload('Subject_003_TESS_Online__feedback__s001_r001_2021_07_06_153655.gdf');
[s2,h2] = sload('Subject_003_TESS_Online__feedback__s001_r002_2021_07_06_154358.gdf');
[s3,h3] = sload('Subject_003_TESS_Online__feedback__s001_r003_2021_07_06_155155.gdf');
[s4,h4] = sload('Subject_003_TESS_Online__feedback__s001_r004_2021_07_06_160446.gdf');
fs = h1.SampleRate;
mu = [8, 12];
Cz = find(strcmp([h1.Label],'CZ              '));

% Selected channels and non-zero recordings
s1 = s1(1:h1.EVENT.POS(end-1), 1:32);
s2 = s2(1:h2.EVENT.POS(end-1), 1:32);
s3 = s3(1:h3.EVENT.POS(end-1), 1:32);
s4 = s4(1:h4.EVENT.POS(end-1), 1:32);


%% Temporal Filter
s1_mu = butter_filt(s1, fs, 2, mu);
s2_mu = butter_filt(s2, fs, 2, mu);
s3_mu = butter_filt(s3, fs, 2, mu);
s4_mu = butter_filt(s4, fs, 2, mu);

% % clip before 1st trial
% s1_mu_short = s1_mu(h1.EVENT.POS(1):end, :);
% s2_mu_short = s2_mu(h2.EVENT.POS(1):end, :);
% s3_mu_short = s3_mu(h3.EVENT.POS(1):end, :);
% s4_mu_short = s4_mu(h4.EVENT.POS(1):end, :);
% t1 = 0:1/fs:(length(s1_mu_short)-1)/fs;
% t2 = 0:1/fs:(length(s2_mu_short)-1)/fs;
% t3 = 0:1/fs:(length(s3_mu_short)-1)/fs;
% t4 = 0:1/fs:(length(s4_mu_short)-1)/fs;
% 
% % plot
% figure();
% sgtitle("mu (8-13Hz)")
% subplot(2, 2, 1);
% [s1_mu_sp, t1] = Space_Channels(s1_mu_short, t1);
% plot(t1, s1_mu_sp)
% ylabel("Amplitude")
% xlabel("Time (s)")
% subtitle("Run 1")
% subplot(2,2,2);
% [s2_mu_sp, t2] = Space_Channels(s2_mu_short, t2);
% plot(t2, s2_mu_sp)
% ylabel("Amplitude")
% xlabel("Time (s)")
% subtitle("Run 2")
% subplot(2,2,3);
% [s3_mu_sp, t3] = Space_Channels(s3_mu_short, t3);
% plot(t3, s3_mu_sp)
% ylabel("Amplitude")
% xlabel("Time (s)")
% subtitle("Run 3")
% subplot(2,2,4);
% [s4_mu_sp, t4] = Space_Channels(s4_mu_short, t4);
% plot(t4, s4_mu_sp)
% ylabel("Amplitude")
% xlabel("Time (s)")
% subtitle("Run 4")


%% Spatial Filter - Large Laplacian
% Second neighbor for each electrode; index corresponds to
% selectedChannels.mat
LargeNeighbors = [3, 5, nan, nan; 
    6, nan, nan, nan; 
    1, 7, nan, nan; 
    5, nan, nan, nan; 
    1, 3, 6, 15;
    2, 5, 7, 15;
    3, 6, 8, 17;
    7, nan, nan, nan;
    10, 20, nan, nan;
    9, 11, 21, nan;
    10, 12, 22, nan;
    11, 23, nan, nan;
    14, 20, nan, nan;
    4, 13, 15, 24;
    5, 14, 15, 25;
    2, 15, 17, 26;
    7, 16, 18, 27;
    8, 17, 19, 28;
    18, 23, nan, nan;
    9, 21, nan, nan;
    10, 20, 22, nan;
    11, 21, 23, nan;
    12, 22, nan nan;
    14, 25, 30, nan;
    15, 24, 26, 30;
    16, 25, 27, 31;
    17, 26, 28, 32;
    18, 27, 32, nan;
    25, 27, 30, 32;
    24, 25, 32, nan;
    26, nan, nan, nan;
    27, 28, 30, nan];

% % Filter application
% s1_muLL = Large_Laplacian(LargeNeighbors, selectedChannels, s1_mu);
% s2_muLL = Large_Laplacian(LargeNeighbors, selectedChannels, s2_mu);
% s3_muLL = Large_Laplacian(LargeNeighbors, selectedChannels, s3_mu);
% s4_muLL = Large_Laplacian(LargeNeighbors, selectedChannels, s4_mu);

% % clip before 1st trial
% s1_muLL_short = s1_muLL(h1.EVENT.POS(1):end, :);
% s2_muLL_short = s2_muLL(h2.EVENT.POS(1):end, :);
% s3_muLL_short = s3_muLL(h3.EVENT.POS(1):end, :);
% s4_muLL_short = s4_muLL(h4.EVENT.POS(1):end, :);
% 
% %plot
% figure();
% sgtitle("mu (8-13Hz) + Large Laplacian")
% subplot(2, 2, 1);
% [s1_muLL_sp, ttemp] = Space_Channels(s1_muLL_short, t1);
% plot(t1, s1_muLL_sp)
% subtitle("Run 1")
% ylabel("Amplitude")
% xlabel("Time (s)")
% subplot(2, 2, 2);
% [s2_muLL_sp, ttemp] = Space_Channels(s2_muLL_short, t2);
% plot(t2, s2_muLL_sp)
% subtitle("Run 2")
% ylabel("Amplitude")
% xlabel("Time (s)")
% subplot(2, 2, 3);
% [s3_muLL_sp, ttemp] = Space_Channels(s3_muLL_short, t3);
% plot(t3, s3_muLL_sp)
% subtitle("Run 3")
% ylabel("Amplitude")
% xlabel("Time (s)")
% subplot(2, 2, 4);
% [s4_muLL_sp, ttemp] = Space_Channels(s4_muLL_short, t4);
% plot(t4, s4_muLL_sp)
% subtitle("Run 4")
% ylabel("Amplitude")
% xlabel("Time (s)")

%% Task Trials
% [sL_muLL, sR_muLL, L_cue, R_cue] = runs2trials_Cue({s1_muLL, s2_muLL, s3_muLL, s4_muLL}, {h1, h2, h3, h4});
% t_Ltasks = 0:1/fs:(size(sL_muLL,2)-1)/fs;
% t_Rtasks = 0:1/fs:(size(sR_muLL,2)-1)/fs;

%plot
% figure();
% sgtitle("First Ten Task Trials - Cz")
% subplot(1, 2, 1);
% [sL_muLL_sp, t_Ltasks] = Space_Channels(squeeze(sL_muLL(Cz, :, 1:10)), t_Ltasks);
% plot(t_Ltasks, sL_muLL_sp)
% subtitle("Left Hand")
% xline(L_cue/fs);
% ylabel("Amplitude")
% xlabel("Time (s)")
% subplot(1, 2, 2);
% [sR_muLL_sp, t_Rtasks] = Space_Channels(squeeze(sR_muLL(Cz, :, 1:10)), t_Rtasks);
% plot(t_Rtasks, sR_muLL_sp)
% xline(R_cue/fs);
% subtitle("Right Hand")
% ylabel("Amplitude")
% xlabel("Time (s)")

%% mu Power with Large Laplacian
% sL_muLL_pow = sL_muLL.^2;
% sR_muLL_pow = sR_muLL.^2;
% 
% %plot
% figure();
% sgtitle("First Ten Task Trials mu Power - Cz")
% subplot(1, 2, 1);
% [sL_muLL_sp, t_Ltasks] = Space_Channels(squeeze(sL_muLL_pow(Cz, :, 1:10)), t_Ltasks);
% plot(t_Ltasks, sL_muLL_sp)
% subtitle("Left Hand")
% xline(L_cue/fs);
% ylabel("Amplitude")
% xlabel("Time (s)")
% subplot(1, 2, 2);
% [sR_muLL_sp, t_Rtasks] = Space_Channels(squeeze(sR_muLL_pow(Cz, :, 1:10)), t_Rtasks);
% plot(t_Rtasks, sR_muLL_sp)
% xline(R_cue/fs);
% subtitle("Right Hand")
% ylabel("Amplitude")
% xlabel("Time (s)")


%% Grand Average mu Power Topoplots - last 0.5s
% Without Large Laplacian
[sL_mu, sR_mu] = runs2trials_End({s1_mu, s2_mu, s3_mu, s4_mu}, {h1, h2, h3, h4}, 0.5);
sL_mu_pow = sL_mu.^2;
sR_mu_pow = sR_mu.^2;
sL_mu_GApow = mean(sL_mu_pow, [1 3]);
sR_mu_GApow = mean(sR_mu_pow, [1 3]);
%topoplots
figure();
subplot(1, 2, 1);
topoplot(sL_mu_GApow, selectedChannels, 'maplimits', [6.5 10.5]);
subtitle("Left Hand")
cbar('vert',0,[6.5 10.5]);
subplot(1, 2, 2);
topoplot(sR_mu_GApow, selectedChannels, 'maplimits', [6.5 10.5]);
subtitle("Right Hand")
cbar('vert',0,[6.5 10.5]);
sgtitle("Grand Average mu Power w/o Spatial Filter")

% With Large Laplacian
sL_muLL = zeros(32, 256, 40);
sR_muLL = zeros(32, 256, 40);
% for tr = 1:size(sL_mu, 3) % reapply on last 0.5s windows
%     sL_muLL(:, :, tr) = (Large_Laplacian(LargeNeighbors, selectedChannels, squeeze(sL_mu(:, :, tr)))).';
%     sR_muLL(:, :, tr) = (Large_Laplacian(LargeNeighbors, selectedChannels, squeeze(sR_mu(:, :, tr)))).';
% end

sL_mu_CAR = sL_mu - mean(sL_mu, 2);
sL_mu_CAR_Pow = sL_mu_CAR .^ 2;
sL_mu_CAR_GAPow = mean(sL_mu_CAR_Pow, [1 3]);
sR_mu_CAR = sR_mu - mean(sR_mu, 2);
sR_mu_CAR_Pow = sR_mu_CAR .^ 2;
sR_mu_CAR_GAPow = mean(sR_mu_CAR_Pow, [1 3]);
 

%topoplots
figure();
subplot(1, 2, 1);
topoplot(sL_mu_CAR_GAPow, selectedChannels, 'maplimits' , [0 4], 'electrodes', 'ptslabels');
subtitle("Left Hand")
cbar('vert',0, [0 4]);
subplot(1, 2, 2);
topoplot(sR_mu_CAR_GAPow, selectedChannels, 'maplimits', [0 4]);
subtitle("Right Hand")
cbar('vert',0, [0 4]);
sgtitle("Grand Average mu Power w/ CAR")

%% Grand Average mu Power Plots - last 0.5s
sL_mu_CAR_GAPow_end = mean(sL_mu_CAR_Pow, 3);
sR_mu_CAR_GAPow_end = mean(sR_mu_CAR_Pow, 3);
tL_end = 1/fs:1/fs:0.5;
tR_end = 1/fs:1/fs:0.5;
figure();
sgtitle("Grand Average Mu Power - Last 0.5s")
subplot(1,2,1);
[sL_mu_CAR_GApow_sp, tL_end] = Space_Channels(sL_mu_CAR_GAPow_end, tL_end);
plot(tL_end, sL_mu_CAR_GApow_sp)
subtitle("Left Hand")
ylabel("Power")
xlabel("Time (s)")
subplot(1, 2, 2);
[sR_muCAR_GApow_sp, tR_end] = Space_Channels(sR_mu_CAR_GAPow_end, tR_end);
plot(tR_end, sR_muCAR_GApow_sp)
subtitle("Right Hand")
ylabel("Power")
xlabel("Time (s)")
legend(flip(h1.Label(1:32)));

function [sig_Ltr, sig_Rtr] = run2trials(sig_run, h)
    trigs = h.EVENT.TYP;
    pos = h.EVENT.POS;
    starts = pos(trigs==1000);
    cues =  trigs(find(trigs==1000)+2);
    ends = pos(find(trigs==1000)+4);
    sig_Ltr{1} = [];
    sig_Rtr{1} = [];
    for i = 1:length(starts)
        if cues(i) == 769
            sig_Ltr{end+1} = sig_run(starts(i):ends(i), :);
        else
            sig_Rtr{end+1} = sig_run(starts(i):ends(i), :);
        end
    end
    max_L = min(cellfun('length', sig_Ltr));
    for i = 2:length(sig_Ltr)
        sig_Ltr{1, i} = mypadarray(sig_Ltr{1, i}, max_L-length(sig_Ltr{1, i})).';
    end
    max_R = max(cellfun('length', sig_Rtr));
    for i = 2:length(sig_Rtr)
        sig_Rtr{1, i} = mypadarray(sig_Rtr{1, i}, max_R-length(sig_Rtr{1, i})).';
    end

    sig_Ltr = reshape(cell2mat(sig_Ltr), 32, max_L, []);
    sig_Rtr = reshape(cell2mat(sig_Rtr), 32, max_R, []);
end

function [sig_Ltr, sig_Rtr] = runs2trials(sig_runs, hs)
    sig_Ltr{1} = [];
    sig_Rtr{1} = [];
    for r = 1:length(sig_runs)
        s = sig_runs{1, r};
        h = hs{1, r};
        trigs = h.EVENT.TYP;
        pos = h.EVENT.POS;
        starts = pos(trigs==1000);
        cues =  trigs(find(trigs==1000)+2);
        ends = pos(find(trigs==1000)+4);
        for i = 1:length(starts)
            if cues(i) == 769
                sig_Ltr{end+1} = s(starts(i):ends(i), :);
            else
                sig_Rtr{end+1} = s(starts(i):ends(i), :);
            end
        end
    end
    max_L = max(cellfun('length', sig_Ltr));
    for i = 2:length(sig_Ltr)
        sig_Ltr{1, i} = mypadarray(sig_Ltr{1, i}, max_L-length(sig_Ltr{1, i})).';
    end
    max_R = max(cellfun('length', sig_Rtr));
    for i = 2:length(sig_Rtr)
        sig_Rtr{1, i} = mypadarray(sig_Rtr{1, i}, max_R-length(sig_Rtr{1, i})).';
    end

    sig_Ltr = reshape(cell2mat(sig_Ltr), 32, max_L, []);
    sig_Rtr = reshape(cell2mat(sig_Rtr), 32, max_R, []);
end

function [sig_Ltr, sig_Rtr] = runs2trials_End(sig_runs, hs, dur)
    len = dur*512;
    sig_Ltr = zeros(len, 32, 40);
    sig_Rtr = zeros(len, 32, 40);
    lt = 1; rt = 1;
    for r = 1:length(sig_runs)
        s = sig_runs{1, r};
        h = hs{1, r};
        trigs = h.EVENT.TYP;
        pos = h.EVENT.POS;
        for i = 1:length(trigs)
            if trigs(i) == 7692 || trigs(i) == 7693
                sig_Ltr(:,:,lt) = s(pos(i)-len+1:pos(i), :);
                lt = lt + 1;
            elseif trigs(i) == 7702 || trigs(i) == 7703
                sig_Rtr(:,:,rt) = s(pos(i)-len+1:pos(i), :);
                rt = rt + 1;
            end
            
        end
    end
end

function [sig_Ltr, sig_Rtr, min_rest_lenL, min_rest_lenR] = runs2trials_Cue(sig_runs, hs)
    sig_Ltr = [];
    sig_Rtr = [];
    ft_len_L = []; rest_len_L = [];
    ft_len_R = []; rest_len_R = [];
    for r = 1:length(sig_runs)
        s = sig_runs{1, r};
        h = hs{1, r};
        trigs = h.EVENT.TYP;
        pos = h.EVENT.POS;
        starts = pos(trigs==1000);
        class = trigs(find(trigs==1000)+2);
        cues = pos(find(trigs==1000)+2);
        ends = pos(find(trigs==1000)+4);
        for i = 1:length(starts)
            if class(i) == 769
                ft_len_L(end+1) = ends(i) - cues(i);
                rest_len_L(end+1) = cues(i) - starts(i);
                sig_Ltr{end+1} = s(starts(i):ends(i), :);
            else
                ft_len_R(end+1) = ends(i) - cues(i);
                rest_len_R(end+1) = cues(i) - starts(i);
                sig_Rtr{end+1} = s(starts(i):ends(i), :);
            end
        end
    end
    min_ft_lenL = min(ft_len_L);
    min_ft_lenR = min(ft_len_R);
    min_rest_lenL = min(rest_len_L);
    min_rest_lenR = min(rest_len_R);
    for i = 1:length(sig_Ltr)
        sig_Ltr{1, i} = (sig_Ltr{1, i}(rest_len_L(i)-min_rest_lenL+1:rest_len_L(i)+min_ft_lenL, :)).';
    end
    for i = 1:length(sig_Rtr)
        sig_Rtr{1, i} = (sig_Rtr{1, i}(rest_len_R(i)-min_rest_lenR+1:rest_len_R(i)+min_ft_lenR, :)).';
    end
    min_L = min(cellfun('length', sig_Ltr));
    min_R = min(cellfun('length', sig_Rtr));
    sig_Ltr = reshape(cell2mat(sig_Ltr), 32, min_L, []);
    sig_Rtr = reshape(cell2mat(sig_Rtr), 32, min_R, []);
end

function LL_sig = Large_Laplacian(Neighbors, chaninfo, sig)
    Filter = zeros(length(Neighbors));
    for r = 1:length(Filter)
        sum_dij = 0;
        for n = Neighbors(r,:)
            if ~isnan(n)
                sum_dij = sum_dij + distance(chaninfo, r, n);
                Filter(r, n) = -1*distance(chaninfo, r, n);
            end
        end
        if sum_dij ~= 0
            Filter(r, :) = Filter(r, :)./sum_dij;
        end
        Filter(r,r) = 1;
    end
    LL_sig = (Filter*(sig.')).';
end

function dij = distance(chaninfo, c1, c2)
    dij = ((chaninfo(c2).X-chaninfo(c1).X).^2+(chaninfo(c2).Y-chaninfo(c1).Y).^2+(chaninfo(c2).Z-chaninfo(c1).Z).^2).^0.5;
end

function [s_spaced, t] = Space_Channels(s, t)
    [m, dim] = min(size(s));
    if dim == 2
        s = s.';
    end
    space = max(range(abs(s), 2));
    offset = (0:space:((m-1)*space))';
    offset = repmat(offset, 1, size(s,2));
    s_spaced = s + offset;
end

function CAR_sig = CAR(sig)
    dim = min(size(sig));
    Filter = -1/(dim)*ones(dim) + eye(dim);
    CAR_sig = Filter*sig;
end


