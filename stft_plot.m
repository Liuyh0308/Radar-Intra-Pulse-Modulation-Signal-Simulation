% 加载 .mat 文件
data = load('Barker-snr09-no00439.mat');  % 请替换 'your_filename.mat' 为你的文件名

% 提取信号序列
% 假设信号序列的变量名为 'signal'
signal = data.wav;

window = hamming(256);  % 使用 512 点的汉明窗，你可以调整此值
overlap = length(window)/2;  % 50% 重叠
nfft = 256;  % FFT 点数


% 计算 STFT
%[s, f, t] = mystft(signal, window,  overlap,  nfft);
[s, f, t] = spectrogram(signal, window,  overlap,  nfft);


%s = 10*log10(abs(s));

wlen = length(window);
C = sum(window)/wlen;
s = abs(s)/wlen/C;
    
s = 20*log10(s + 1e-6);

    %figure(1)
surf(f, t, s')
shading interp;
axis tight;
view(0, 90);
%set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel('Time, s');
xlabel('Frequency, Hz');
title('Amplitude spectrogram of the signal');



