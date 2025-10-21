% NSVG v1.0 — NifTI Sub Volume Generator (made by Aayush Goud)
% Manual 3D lesion/ROI tracer with live 3D preview + preview/confirm +
% native file dialogs for opening and saving. Works with MNI images. You can
% toggle between the different cross sections to find the cross section you
% want to use. NOTE: once you start editing on one cross section, editing
% on another cross section will likely cause the program to crash.
%
% Utilization: trace out regions of interest (ROIs) on MRI images; whether
% it be brain bleeds, glioblastomae, or lesions caused by FUS thalamotomy
% or traumatic brain injury.
%
% Output: .nii/.nii.gz mask of a traced ROI that can be overlayed in NIfTI
% scenes to visualize otherwise non-native ROIs in relation to other brain
% regions. Initially designed to create .nii files to interface with the
% 3D scenes created by the leadDBS application (Horn et al).
%
% To run:  NSVG;   (in the MATLAB Command Window)
%
% Requirements: Image Processing Toolbox (drawfreehand), MRI NIfTI file in
% MNI coordinates (native/ACPC coords implementation under development).
%
% Not currently cleared for clinical use.
