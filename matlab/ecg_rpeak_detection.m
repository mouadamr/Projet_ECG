%% =====================================================================
%  Projet 1 : Filtrage d'un signal ECG et detection des pics R
%             pour l'estimation de la frequence cardiaque
%  Module : FPGA, VHDL et Traitement de Signal Numerique - Master IISE
%  Methode : Pan-Tompkins (version simplifiee, chaine causale -> FPGA)
%  Donnees : MIT-BIH Arrhythmia Database, enregistrement 100 (fs = 360 Hz)
% =====================================================================
clear; close all; clc;

%% 1. Chargement du signal ------------------------------------------------
% Dossiers du projet (calcules a partir de l'emplacement de ce script)
scriptDir = fileparts(mfilename('fullpath'));
projDir   = fileparts(scriptDir);
dataFile  = fullfile(projDir, 'data', 'ecg_100.mat');
figDir    = fullfile(projDir, 'figures');
simDir    = fullfile(projDir, 'simulink');
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(simDir, 'dir'), mkdir(simDir); end

% Fichier prepare a partir de MIT-BIH record 100 (derivation MLII) :
%   ecg  : signal en mV     fs : 360 Hz     annR : positions des battements
%                                                  annotees par les cardiologues
D = load(dataFile);
ecg  = D.ecg(:);
fs   = D.fs;
annR = D.annR(:);
haveAnn = true;
N = numel(ecg);
t = (0:N-1)' / fs;
fprintf('Signal : %d echantillons, fs = %d Hz, duree = %.1f s\n', N, fs, N/fs);

%% 2. Filtrage de visualisation (0.5 - 40 Hz) -------------------------------
% Supprime la derive de ligne de base (<0.5 Hz), le secteur (60 Hz ici)
% et le bruit musculaire haute frequence.
[bv, av] = butter(2, [0.5 40] / (fs/2), 'bandpass');
ecg_f = filtfilt(bv, av, ecg);                     % phase nulle (hors ligne)

%% 3. Chaine de detection Pan-Tompkins (causale, implementable sur FPGA) ---
% 3.1 Passe-bande FIR 5-15 Hz (energie du complexe QRS)
Nfir = 64;
b_bp = fir1(Nfir, [5 15] / (fs/2));
x_bp = filter(b_bp, 1, ecg);
dBP  = Nfir/2;                                     % retard du FIR (phase lineaire)

% 3.2 Derivee a 5 points
b_der = [1 2 0 -2 -1] * (fs/8);
x_der = filter(b_der, 1, x_bp);
dDer  = 2;

% 3.3 Elevation au carre
x_sq = x_der .^ 2;

% 3.4 Integration par fenetre glissante (150 ms)
W = round(0.150 * fs);
x_mwi = filter(ones(1, W) / W, 1, x_sq);
dMWI  = round((W-1)/2);

totalDelay = dBP + dDer + dMWI;

%% 4. Detection des pics R -------------------------------------------------
minDist = round(0.25 * fs);                        % periode refractaire 250 ms
[pks_all, ~] = findpeaks(x_mwi, 'MinPeakDistance', minDist);
thr = 0.3 * prctile(pks_all, 90);                  % seuil fixe (a rendre adaptatif)
[~, locs] = findpeaks(x_mwi, 'MinPeakDistance', minDist, 'MinPeakHeight', thr);

% Recherche du vrai maximum R dans l'ECG filtre autour de la position estimee
rLocs = zeros(size(locs));
for k = 1:numel(locs)
    c  = locs(k) - totalDelay;
    i1 = max(1, c - round(0.20*fs));
    i2 = min(N, c + round(0.10*fs));
    [~, im] = max(ecg_f(i1:i2));
    rLocs(k) = i1 + im - 1;
end
rLocs = unique(rLocs);
rLocs = rLocs([true; diff(rLocs) > minDist]);      % supprime les doublons

%% 5. Frequence cardiaque ---------------------------------------------------
RR = diff(rLocs) / fs;                             % intervalles RR (s)
HR = 60 ./ RR;                                     % bpm
fprintf('Battements detectes : %d\n', numel(rLocs));
fprintf('Frequence cardiaque moyenne : %.1f bpm (min %.1f, max %.1f)\n', ...
        mean(HR), min(HR), max(HR));

%% 6. Validation par rapport aux annotations ------------------------------
if haveAnn
    tol = round(0.05 * fs);                        % tolerance 50 ms
    TP = sum(arrayfun(@(r) any(abs(annR - r) <= tol), rLocs));
    FP = numel(rLocs) - TP;
    FN = numel(annR) - TP;
    Se  = 100 * TP / (TP + FN);
    PPV = 100 * TP / (TP + FP);
    fprintf('Validation : TP=%d  FP=%d  FN=%d  Se=%.2f%%  PPV=%.2f%%\n', ...
            TP, FP, FN, Se, PPV);
end

%% 7. Etude de la quantification (preparation FPGA) -----------------------
nbits = 16;                                        % coefficients en Q1.15
b_q = round(b_bp * 2^(nbits-1)) / 2^(nbits-1);
[H,  f] = freqz(b_bp, 1, 2048, fs);
[Hq, ~] = freqz(b_q,  1, 2048, fs);
x_bp_q = filter(b_q, 1, ecg);
SQNR = 10*log10(sum(x_bp.^2) / sum((x_bp - x_bp_q).^2));
fprintf('Quantification %d bits : SQNR = %.1f dB\n', nbits, SQNR);

% Export des coefficients entiers pour le code VHDL
coeffs_int = round(b_bp * 2^(nbits-1));
dlmwrite(fullfile(simDir, 'fir_coeffs_q15.txt'), coeffs_int(:));   % compatible anciennes versions

%% 8. Export pour Simulink ----------------------------------------------------
% Bloc "From Workspace" -> variable simin ; comparer la sortie a x_mwi
simin = [t ecg];
save(fullfile(simDir, 'ecg_simulink.mat'), 'simin', 'fs', 'b_bp', 'b_der', 'W', 'x_mwi');

%% 9. Figures ---------------------------------------------------------------
idx = t < 10;                                      % 10 premieres secondes

figure('Name', 'Signal brut vs filtre');
subplot(2,1,1); plot(t(idx), ecg(idx)); grid on;
title('ECG brut'); xlabel('Temps (s)'); ylabel('mV');
subplot(2,1,2); plot(t(idx), ecg_f(idx)); grid on;
title('ECG filtre (0.5-40 Hz)'); xlabel('Temps (s)'); ylabel('mV');

figure('Name', 'Spectres');
[Pr, fr] = pwelch(ecg,   4096, [], [], fs);
[Pf, ~ ] = pwelch(ecg_f, 4096, [], [], fs);
plot(fr, 10*log10(Pr), fr, 10*log10(Pf)); grid on; xlim([0 100]);
legend('Brut', 'Filtre'); xlabel('Frequence (Hz)'); ylabel('DSP (dB/Hz)');
title('Densite spectrale de puissance (Welch)');

figure('Name', 'Spectrogramme');
spectrogram(ecg(t < 60), hamming(512), 384, 1024, fs, 'yaxis');
ylim([0 60]); title('Spectrogramme ECG brut (60 s)');

figure('Name', 'Etapes Pan-Tompkins');
stages = {ecg, x_bp, x_der, x_sq, x_mwi};
names  = {'ECG brut', 'Passe-bande 5-15 Hz', 'Derivee', 'Carre', 'Integration 150 ms'};
for k = 1:5
    subplot(5,1,k); plot(t(idx), stages{k}(idx)); grid on; title(names{k});
end
xlabel('Temps (s)');

figure('Name', 'Pics R detectes');
plot(t(idx), ecg_f(idx)); hold on; grid on;
r10 = rLocs(rLocs <= find(idx, 1, 'last'));
plot(t(r10), ecg_f(r10), 'rv', 'MarkerFaceColor', 'r');
title('Detection des pics R'); xlabel('Temps (s)'); ylabel('mV');

figure('Name', 'Frequence cardiaque');
plot(t(rLocs(2:end)) / 60, HR); grid on;
title('Frequence cardiaque instantanee'); xlabel('Temps (min)'); ylabel('bpm');

figure('Name', 'Quantification FIR');
plot(f, 20*log10(abs(H)), f, 20*log10(abs(Hq)), '--'); grid on;
legend('Double precision', sprintf('Q1.%d', nbits-1));
xlabel('Frequence (Hz)'); ylabel('|H| (dB)');
title('Effet de la quantification des coefficients');

%% 10. Sauvegarde automatique des figures (PNG) dans le dossier "figures" --
figs = findobj('Type', 'figure');
for k = 1:numel(figs)
    nom = regexprep(get(figs(k), 'Name'), '[^\w]', '_');
    saveas(figs(k), fullfile(figDir, [nom '.png']));
end
fprintf('Figures enregistrees dans : %s\n', figDir);
