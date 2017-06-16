function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.

% windows length muss odd sein
window_length = 5;
% Zahl der punkten in dem Bildsegment
N = window_length^2;
min_corr = 0.95;
% Zahl der Merkalpunkten
numberP1 = size(Mpt1, 2);
numberP2 = size(Mpt2, 2);

% init 
Korrespondenzen = zeros(4, 1);
windowsP1 = zeros(N, numberP1);
windowsP2 = zeros(N, numberP2);

%%
% waehlen Bildersegmente um P1
for i = 1:numberP1
windowsP1(:, i) = reshape(I1(Mpt1(2, i)-(window_length-1)/2 : Mpt1(2, i)+(window_length-1)/2, ...
    Mpt1(1, i)-(window_length-1)/2 : Mpt1(1, i)+(window_length-1)/2), [], 1);
% Bias norminieren
windowsP1BiasNorm(:, i) = double(windowsP1(:, i)) - mean(windowsP1(:, 1));
% Gain norminieren
windowsP1BiasGainNorm(:, i) = windowsP1BiasNorm(:, i)/(sqrt((norm(windowsP1BiasNorm(:, i)))^2/(N-1)));
% Ergaenzung
end

%%
for j = 1:numberP2
windowsP2(:, j) = reshape(I2(Mpt2(2, j)-(window_length-1)/2 : Mpt2(2, j)+(window_length-1)/2, ...
    Mpt2(1, j)-(window_length-1)/2 : Mpt2(1, j)+(window_length-1)/2), [], 1);
windowsP2BiasNorm(:, j) = double(windowsP2(:, j)) - mean(windowsP2(:, j));
windowsP2BiasGainNorm(:, j) = windowsP2BiasNorm(:, j)/(sqrt((norm(windowsP2BiasNorm(:, j)))^2/(N-1)));
end

%%
% NCC rechnen
NCC = (windowsP1BiasGainNorm'*windowsP2BiasGainNorm)/(N-1);
[M, I] = max(NCC, [], 2);
for i = 1:max(numberP1, numberP2)
    if M(i) > min_corr
        Korrespondenzen = [Korrespondenzen, [Mpt1(:, i); Mpt2(:, I(i))]];
    end
end
Korrespondenzen = Korrespondenzen(:, 2:end);

%%
% ploten
% locationP1 = [Korrespondenzen(2, :); Korrespondenzen(1, :)]';
% locationP2 = [Korrespondenzen(4, :); Korrespondenzen(3, :)]';
locationP1 = [Mpt1(2, 1:20); Mpt1(1, 1:20)]';
locationP2 = [Mpt2(4, 1:20); Mpt2(3, 1:20)]';

matchedPointsP1 = cornerPoints(locationP1);
matchedPointsP2 = cornerPoints(locationP2);
% showMatchedFeatures(IGray1,IGray2,matchedPointsP1,matchedPointsP2);
figure(1);
imshow(I1); hold on;
plot(matchedPointsP1);
figure(2);
imshow(I2); hold on;
plot(matchedPointsP2);

end


