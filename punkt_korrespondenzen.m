function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.

% WICHTIG: windows_length (odd) darf nicht mehr als 11 sein
window_length = 5;
% Zahl der punkten in dem Bildsegment
N = window_length^2;
min_corr = 0.99;
% Zahl der Merkalpunkten
numberP1 = size(Mpt1, 2);
numberP2 = size(Mpt2, 2);

% init 
Korrespondenzen = zeros(4, 1);


%%
% Loesung I
windowsP1 = zeros(N, numberP1);
windowsP2 = zeros(N, numberP2);
% waehlen Bildersegmente um P1
for i = 1:numberP1
windowsP1(:, i) = reshape(I1(Mpt1(2, i)-(window_length-1)/2 : Mpt1(2, i)+(window_length-1)/2, ...
    Mpt1(1, i)-(window_length-1)/2 : Mpt1(1, i)+(window_length-1)/2), [], 1);
% Bias norminieren
windowsP1BiasNorm(:, i) = double(windowsP1(:, i)) - mean(windowsP1(:, 1));
% Gain norminieren
windowsP1BiasGainNorm(:, i) = windowsP1BiasNorm(:, i)/std(windowsP1BiasNorm(:, i));
% Ergaenzung
end

for j = 1:numberP2
windowsP2(:, j) = reshape(I2(Mpt2(2, j)-(window_length-1)/2 : Mpt2(2, j)+(window_length-1)/2, ...
    Mpt2(1, j)-(window_length-1)/2 : Mpt2(1, j)+(window_length-1)/2), [], 1);
windowsP2BiasNorm(:, j) = double(windowsP2(:, j)) - mean(windowsP2(:, j));
windowsP2BiasGainNorm(:, j) = windowsP2BiasNorm(:, j)/std(windowsP2BiasNorm(:, j));
end

% NCC rechnen
% vektorielles Skalarprodukt wie auf dem PDF
NCC = (windowsP1BiasGainNorm'*windowsP2BiasGainNorm)/(N-1);

[M, I] = max(NCC, [], 2);
for i = 1:max(numberP1, numberP2)
    if M(i) > min_corr
        Korrespondenzen = [Korrespondenzen, [Mpt1(:, i); Mpt2(:, I(i))]];
    end
end
% Entfernung der ersten 0 Spalte
Korrespondenzen = Korrespondenzen(:, 2:end);

%%
% ploten
imshow([I1, ones(size(I1, 1), 10), I2]); hold on;
plot(Korrespondenzen(2, :)', Korrespondenzen(1, :)','rx'); hold on;
plot(Korrespondenzen(4, :)'+size(I1, 2)*ones(size(Korrespondenzen, 2), 1), Korrespondenzen(3, :)','rx'); hold on;
plot([Korrespondenzen(2, :); Korrespondenzen(4, :)+size(I1, 2)*ones(1, size(Korrespondenzen, 2))+10], ... 
    [Korrespondenzen(1, :); Korrespondenzen(3, :)], 'y-');

end


