# C_save_roi_per_participant

Subjects: 22 (from manifest.mat)
Cells: 48
Output rows: 1056

Each row: participant, band, measure, decomposition, reference, roi_value (high - low, ROI average).
roi_value = mean over ROI pairs of (high_condition_value - low_condition_value).
Per-subject contrasts computed from EEG files (not battery_persubject.csv).

Output: roi_per_participant.csv
