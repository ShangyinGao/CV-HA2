function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.
tic
% varargin
P = inputParser;
P.addOptional('do_plot', false, @islogical);
P.addOptional('window_length', 15, @isnumeric);
P.addOptional('min_corr', 0.99, @isnumeric);

P.parse(varargin{:});

do_plot = P.Results.do_plot;
window_length = P.Results.window_length;
min_corr = P.Results.min_corr;

% ergaezen n nullspalte und nullzeile, passen numberOfZero zuc
% windows_length. numberOfZero > windows_length-13
numberOfZero = window_length+10;
% Zahl der punkten in dem Bildsegment
N = window_length^2;
% Zahl der Merkalpunkten
numberP1 = size(Mpt1, 2);
numberP2 = size(Mpt2, 2);


% init Korrespondenzen
Korrespondenzen = zeros(4, 1);
% vergroessen images
I1Extend = zeros(numberOfZero*2+size(I1, 1), numberOfZero*2+size(I1, 2));
I1Extend(numberOfZero+1 : numberOfZero+size(I1, 1), numberOfZero+1 : numberOfZero+size(I1, 2)) = I1;
I2Extend = zeros(numberOfZero*2+size(I2, 1), numberOfZero*2+size(I2, 2));
I2Extend(numberOfZero+1 : numberOfZero+size(I2, 1), numberOfZero+1 : numberOfZero+size(I2, 2)) = I2;


%%
% Bias-Gain-Modell
windowsP1 = zeros(N, numberP1);
windowsP2 = zeros(N, numberP2);
% waehlen Bildersegmente um P1
for i = 1:numberP1
    windowsP1(:, i) = reshape(I1Extend(numberOfZero+Mpt1(2, i)-(window_length-1)/2 : numberOfZero+Mpt1(2, i)+(window_length-1)/2, ...
        numberOfZero+Mpt1(1, i)-(window_length-1)/2 : numberOfZero+Mpt1(1, i)+(window_length-1)/2), [], 1);
    % Bias norminieren
    windowsP1BiasNorm(:, i) = double(windowsP1(:, i)) - mean(windowsP1(:, i));
    % Gain norminieren
    windowsP1BiasGainNorm(:, i) = windowsP1BiasNorm(:, i)/std(windowsP1BiasNorm(:, i));
end

for j = 1:numberP2
    windowsP2(:, j) = reshape(I2Extend(numberOfZero+Mpt2(2, j)-(window_length-1)/2 : numberOfZero+Mpt2(2, j)+(window_length-1)/2, ...
        numberOfZero+Mpt2(1, j)-(window_length-1)/2 : numberOfZero+Mpt2(1, j)+(window_length-1)/2), [], 1);
    windowsP2BiasNorm(:, j) = double(windowsP2(:, j)) - mean(windowsP2(:, j));
    windowsP2BiasGainNorm(:, j) = windowsP2BiasNorm(:, j)/std(windowsP2BiasNorm(:, j));
end

% NCC rechnen
% vektorielles Skalarprodukt wie auf dem PDF
NCC = (windowsP1BiasGainNorm'*windowsP2BiasGainNorm)/(N-1);
[M, I] = max(NCC, [], 2);
for i = 1:numberP1
    if M(i) > min_corr
        Korrespondenzen = [Korrespondenzen, [Mpt1(:, i); Mpt2(:, I(i))]];
    end
end
% Entfernung der ersten 0 Spalte
Korrespondenzen = Korrespondenzen(:, 2:end);

%%
% ploten
if do_plot
    imshow([I1, ones(size(I1, 1), 10), I2]); hold on;
    plot(Korrespondenzen(1, :)', Korrespondenzen(2, :)','rx'); hold on;
    plot(Korrespondenzen(3, :)'+size(I1, 2)*ones(size(Korrespondenzen, 2), 1), Korrespondenzen(4, :)','rx'); hold on;
    plot([Korrespondenzen(1, :); Korrespondenzen(3, :)+size(I1, 2)*ones(1, size(Korrespondenzen, 2))+10], ... 
        [Korrespondenzen(2, :); Korrespondenzen(4, :)], 'y-');
    
end

%%
timeElapsed = toc
% print zahl der gefundenen Korrespondenzen
fprintf('Zahl der gefundenen Korrespondenzen ist %d\n', size(Korrespondenzen, 2));

end


