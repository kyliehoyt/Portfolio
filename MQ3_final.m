clear all
close all
clc

%% PART 1: Read the music files you download
[filename, pathname] = uigetfile('*.mp3', 'Pick the music mp3 file you want to load')
file      = strcat(pathname,filename);

[y,fs] = audioread(file); % y is the sampled data and fs is sampling frequency
JB1 = y(:,1);
[filename, pathname] = uigetfile('*.wav','Pick the music file wav you want to load')
file      = strcat(pathname,filename);
[y2,fs] = audioread(file);
JB2 = y2(:,1);
% We can use the function 'sound(y,fs)' to pay the music in MATLAB. Please
% copy this function to you Command Window to enjoy this music. ##NOTE: enter
% 'clear sound' in your Command Window to stop playing.

% Now if we want Jingle Bells to be played twice as faster, how do we change the
% input of this 'sound' function?
sound(JB1,2*fs);
pause(7);
clear sound
sound(JB2,2*fs);
pause(7);
clear sound

% These two songs SOUND VERY different, but they are more related than you
% think. Use Fourier transforms to identify the difference (there is only
% one operation) between these two songs.  
%----------------------------------------------------------

%a) Plot the two Fourier transforms (real and imaginary parts for each song) with the proper spacing on the X-axis
%in Hz (not rad/s), the sampling frequency of the audio waveforms is the variable fs.

%b) What is the difference between these two songs

%Please upload the graphs from part a and your explanation for how these
%sonds are different.

%----------------------------------------------------------
%Here are some commands that might help you.
%fft (this an optimized version of the DFT and it ALWAYS used in signal
%processing)
%fftshift (this shifts the fft spectrum, so that the "0" frequency is in the center
%real (gets the real part of a vector)
%imag (gets the imaginary part of a vector)