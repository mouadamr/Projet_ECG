PROJET ECG - Filtrage et détection des pics R (Projet n°1)
=========================================================

Folders
-------
data\      ecg_100.mat        -> the heart signal (MIT-BIH record 100, 30 min, 360 Hz)
matlab\    ecg_rpeak_detection.m -> the main program (open it and press Run)
simulink\  created by the program: ecg_simulink.mat + fir_coeffs_q15.txt
           (later: ecg_model.slx, the Simulink block model)
figures\   created by the program: all graphs saved as PNG for the report
rapport\   the report (later)

How to run
----------
1. Open MATLAB.
2. Open matlab\ecg_rpeak_detection.m  (File > Open, or double-click it).
3. Click the green "Run" button (or press F5).
4. Read the results in the Command Window. 7 graphs open and are saved in figures\.

Expected results (record 100)
-----------------------------
Beats detected      : 2272  (reference: 2273)
Mean heart rate     : about 75.8 bpm
Sensitivity Se      : about 99.9 %
Precision  PPV      : about 99.9 %
Quantization 16 bits: SQNR about 78 dB
