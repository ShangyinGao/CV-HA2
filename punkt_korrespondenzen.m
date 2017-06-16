function [Korrespondenzen] = punkt_korrespondenzen(I1,I2,Mpt1,Mpt2,varargin)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.

% WICHTIG: windows_length (odd) darf nicht mehr als 11 sein
window_length = 11;
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
% % Loesung II
% NCC_tmp = zeros(numberP2, 1);
% NCCMaxForEachP1 = zeros(numberP1, 2);
% einMatrix = ones(window_length, 1);
% for i = 1:numberP1
%     windowsP1_tmp = I1(Mpt1(2, i)-(window_length-1)/2 : Mpt1(2, i)+(window_length-1)/2, ...
%         Mpt1(1, i)-(window_length-1)/2 : Mpt1(1, i)+(window_length-1)/2);
%     windowsP1BiasNorm_tmp = double(windowsP1_tmp) - mean(mean(windowsP1_tmp));
% %     windowsP1BiasGainNorm_tmp = windowsP1BiasNorm_tmp./(sqrt((norm(windowsP1BiasNorm_tmp))^2/(N-1)));
%     windowsP1BiasGainNorm_tmp = windowsP1BiasNorm_tmp./std(reshape(windowsP1BiasNorm_tmp, [], 1));
%     
% 
%     for j = 1:numberP2
%         windowsP2_tmp = I2(Mpt2(2, j)-(window_length-1)/2 : Mpt2(2, j)+(window_length-1)/2, ...
%             Mpt2(1, j)-(window_length-1)/2 : Mpt2(1, j)+(window_length-1)/2);
%         windowsP2BiasNorm_tmp = double(windowsP2_tmp) - mean(mean(windowsP2_tmp));
% %         windowsP2BiasGainNorm_tmp = windowsP2BiasNorm_tmp./(sqrt((norm(windowsP2BiasNorm_tmp))^2/(N-1)));
%         windowsP2BiasGainNorm_tmp = windowsP2BiasNorm_tmp./std(reshape(windowsP2BiasNorm_tmp, [], 1));
%         NCC_tmp(j) = trace(windowsP1BiasGainNorm_tmp'*windowsP2BiasGainNorm_tmp)/(N-1);
% 
%     end
%     [M, I] = max(NCC_tmp);
%     NCCMaxForEachP1(i, :) = [M, I];
%     
% end

%%
% ploten
locationP1 = [Korrespondenzen(2, :); Korrespondenzen(1, :)]';
locationP2 = [Korrespondenzen(4, :); Korrespondenzen(3, :)]';
% % locationP1 = [Mpt1(2, :); Mpt1(1, :)]';
% % locationP2 = [Mpt2(2, :); Mpt2(1, :)]';
% 
% matchedPointsP1 = cornerPoints(locationP1);
% matchedPointsP2 = cornerPoints(locationP2);
% % showMatchedFeatures(IGray1,IGray2,matchedPointsP1,matchedPointsP2);
% figure(1);
% imshow(I1); hold on;
% plot(matchedPointsP1);
% figure(2);
% imshow(I2); hold on;
% plot(matchedPointsP2);

figure(1),imshow(I1),hold on,
plot(Korrespondenzen(2, :)', Korrespondenzen(1, :)','ys');
figure(2),imshow(I2),hold on,
plot(Korrespondenzen(4, :)', Korrespondenzen(3, :)','ys');


end


