classdef NSVG < handle
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------

    properties (Access=private)
        % Data
        info                 % niftiinfo of reference image
        V0   double          % original, display-normalized 3D volume (scanner order)
        V    double          % 3D volume in display space
        BW   logical         % 3D mask in display space
        orientation char     % 'axial'|'coronal'|'sagittal'
        dim  uint8 = 3       % slicing dim in display space (always 3 here)
        nSlices double
        k double             % current slice idx

        % UI handles
        f                    % main figure
        ax                   % 2D slice axes
        imH                  % 2D image handle
        ttl                  % title handle
        ax3                  % 3D preview axes
        p3                   % 3D preview patch
        sld                  % slider
        btnAccept
        btnUndo
        tglPrev
        btnFinish
        btnOpen
        ddOrient
        note                 % status line
        tb                   % toolbar
        rois                 % array of drawfreehand handles

        % last-used paths
        lastOpenDir char = ''
        lastSaveMask char = ''
        lastSaveSurf char = ''

        % mode state
        mode char = 'idle'   % 'idle'|'edit'|'preview'|'finish'|'abort'
    end

    methods
        function self = NSVG()
            self.buildUI();
        end
    end

    %================== UI Construction ==================%
    methods (Access=private)
        function buildUI(self)
            % Big main window
            self.f = figure( ...
                'Name','NSVG — NifTI Sub Volume Generator  |  Open → Trace → Finish', ...
                'Color','w','NumberTitle','off', ...
                'KeyPressFcn',@self.onKey, ...
                'WindowScrollWheelFcn',@self.onScroll, ...
                'Units','normalized','Position',[0.02 0.05 0.96 0.88]);

            % Axes: 2D slice (left)
            self.ax  = axes('Parent',self.f,'Position',[0.05 0.12 0.55 0.83]); %#ok<LAXES>
            self.imH = imshow(zeros(10,10),[],'Parent',self.ax); hold(self.ax,'on');
            colormap(self.ax, gray(256));
            self.ttl = title(self.ax,'Open a NIfTI to start');

            % Axes: live 3D (right)
            self.ax3 = axes('Parent',self.f,'Position',[0.65 0.12 0.30 0.83]);
            hold(self.ax3,'on'); axis(self.ax3,'equal'); axis(self.ax3,'vis3d');
            title(self.ax3,'Live 3D preview (MNI/native mm)');
            xlabel(self.ax3,'x'); ylabel(self.ax3,'y'); zlabel(self.ax3,'z');
            box(self.ax3,'on'); grid(self.ax3,'on');
            self.p3 = [];

            % Controls row
            self.btnOpen = uicontrol('Style','pushbutton','Units','normalized','Position',[0.05 0.04 0.10 0.04], ...
                'String','Open NIfTI…','Callback',@(~,~) self.onOpen());

            self.ddOrient = uicontrol('Style','popupmenu','Units','normalized','Position',[0.16 0.04 0.10 0.04], ...
                'String',{'axial','coronal','sagittal'}, 'Value',1, ...
                'TooltipString','Display orientation', ...
                'Callback',@(~,~) self.onChangeOrientation());

            self.sld = uicontrol('Style','slider','Units','normalized','Position',[0.27 0.04 0.18 0.04], ...
                'Min',1,'Max',2,'Value',1,'Enable','off','Callback',@(h,~) self.setK(round(get(h,'Value'))));

            self.btnAccept = uicontrol('Style','pushbutton','Units','normalized','Position',[0.46 0.04 0.10 0.04], ...
                'String','Accept slice','Enable','off','Callback',@(~,~) self.acceptSlice());

            self.btnUndo = uicontrol('Style','pushbutton','Units','normalized','Position',[0.57 0.04 0.10 0.04], ...
                'String','Undo all','Enable','off','Callback',@(~,~) self.clearSlice());

            self.tglPrev = uicontrol('Style','togglebutton','Units','normalized','Position',[0.68 0.04 0.12 0.04], ...
                'String','Prev→Template','Enable','off', ...
                'TooltipString','Copy previous slice contours as template', ...
                'Callback',@(h,~) setappdata(self.f,'usePrevTemplate',get(h,'Value')));

            self.btnFinish = uicontrol('Style','pushbutton','Units','normalized','Position',[0.81 0.04 0.14 0.04], ...
                'String','Finish & Save','Enable','off','Callback',@(~,~) self.finishPressed());

            % Status line
            self.note = annotation(self.f,'textbox',[0.05 0.005 0.90 0.03],'String','', ...
                'EdgeColor','none','HorizontalAlignment','center');

            % Toolbar: freehand button
            self.tb = uitoolbar(self.f);
            icon = self.local_poly_icon();
            uipushtool(self.tb,'CData',icon,'TooltipString','Draw freehand ROI', ...
                'ClickedCallback',@(~,~) self.drawNewFreehand());

            self.rois = gobjects(0);
            self.updateNote('Open a NIfTI (patient T1 in MNI or native).');
        end
    end

    %================== Open / Init ==================%
    methods (Access=private)
        function onOpen(self)
            if ~isempty(self.lastOpenDir) && isfolder(self.lastOpenDir)
                startDir = self.lastOpenDir;
            else
                startDir = pwd;
            end
            [f,p] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, 'Choose input NIfTI', startDir);
            if isequal(f,0), return; end
            inPath = fullfile(p,f);
            self.lastOpenDir = p;

            % Load NIfTI
            self.info = niftiinfo(inPath);
            Vraw = double(niftiread(self.info));
            if ndims(Vraw) > 3
                Vraw = Vraw(:,:,:,1); % first volume if 4D
            end
            Vraw = Vraw - min(Vraw(:));
            mx = max(eps, max(Vraw(:)));
            Vraw = Vraw / mx;
            Vraw(~isfinite(Vraw)) = 0;  % guard

            % Keep original; derive display volume from this every time
            self.V0 = Vraw;

            % Build display-space for selected orientation
            ori = self.getOrient();
            self.applyOrientation(ori);

            self.mode = 'edit';
            self.updateNote('Use toolbar to draw ROIs. Double-click to close path. Press A to Accept.');

            fprintf('[NSVG] Loaded %s  size=%s  class=%s\n', ...
                self.info.Filename, mat2str(size(self.V0)), class(self.V0));
            fprintf('[NSVG] Display orientation: %s  nSlices=%d\n', self.orientation, self.nSlices);
        end

        function onChangeOrientation(self)
            if isempty(self.V0), return; end
            ori = self.getOrient();
            self.applyOrientation(ori);
            delete(self.rois(ishandle(self.rois))); self.rois = gobjects(0);
            self.updateNote('Orientation changed. Mask cleared to avoid confusion.');
        end

        function applyOrientation(self, ori)
            % Recompute V/BW/k/nSlices from the original V0 for a given orientation
            switch lower(ori)
                case 'axial'
                    self.V = self.V0;                  self.orientation='axial';
                case 'coronal'
                    self.V = permute(self.V0, [3 1 2]); self.orientation='coronal';
                case 'sagittal'
                    self.V = permute(self.V0, [3 2 1]); self.orientation='sagittal';
                otherwise
                    self.V = self.V0; self.orientation='axial';
            end
            self.BW      = false(size(self.V));
            self.nSlices = size(self.V,3);
            self.k       = max(1, round(self.nSlices/2));
            self.updateSlider();
            self.recreateSliceImage();
            self.refreshOverlay();
            self.update3D();
            self.updateTitle();
        end

        function updateSlider(self)
            if self.nSlices > 1
                step = [1/(self.nSlices-1) min(10/(self.nSlices-1), 1)];
            else
                step = [1 1];
            end
            set(self.sld,'Min',1,'Max',self.nSlices,'Value',self.k,'SliderStep',step,'Enable','on');
            set([self.btnAccept,self.btnUndo,self.tglPrev,self.btnFinish], 'Enable','on');
        end

        function recreateSliceImage(self)
            if isgraphics(self.imH), delete(self.imH); end
            cla(self.ax);
            self.imH = imshow(self.getImg(self.k), [], 'Parent', self.ax);
            colormap(self.ax, gray(256));
            hold(self.ax,'on');
        end

        function s = getOrient(self)
            strs = get(self.ddOrient,'String'); s = strs{get(self.ddOrient,'Value')};
        end
    end

    %================== Interaction ==================%
    methods (Access=private)
        function setK(self, kNew)
            if isempty(self.V), return; end
            self.k = max(1, min(self.nSlices, kNew));
            if isgraphics(self.sld), set(self.sld,'Value',self.k); end
            % Update displayed slice; recreate if size mismatch (paranoia)
            try
                set(self.imH, 'CData', self.getImg(self.k));
            catch
                self.recreateSliceImage();
            end
            delete(self.rois(ishandle(self.rois))); self.rois = gobjects(0);
            % "Prev→Template"
            if self.k>1 && isappdata(self.f,'usePrevTemplate') && getappdata(self.f,'usePrevTemplate')
                prev = self.getSliceBW(self.k-1);
                if any(prev(:))
                    B = bwboundaries(prev,'noholes');
                    for b = 1:numel(B)
                        pos = self.roibounds2pos(B{b});
                        self.rois(end+1) = drawfreehand(self.ax,'Position',pos,'Closed',true); %#ok<AGROW>
                    end
                end
            end
            self.refreshOverlay(); self.updateTitle(); self.updateNote();
        end

        function onScroll(self, ~, evt)
            if isempty(self.V), return; end
            self.setK(self.k + sign(-evt.VerticalScrollCount));
        end

        function onKey(self, ~, evt)
            if isempty(self.V), return; end
            switch evt.Key
                case {'rightarrow','uparrow'}, self.setK(self.k+1);
                case {'leftarrow','downarrow'}, self.setK(self.k-1);
                case 'a', self.acceptSlice();
                case 'u', self.clearSlice();
                case 'q', self.finishPressed();
                case 'escape', self.mode='abort'; uiresume(self.f);
            end
        end

        function drawNewFreehand(self)
            if isempty(self.V), return; end
            self.rois(end+1) = drawfreehand(self.ax,'Closed',true); %#ok<AGROW>
        end

        function acceptSlice(self)
            if isempty(self.V), return; end
            M2 = false(size(self.getImg(self.k)));
            for r = 1:numel(self.rois)
                if isgraphics(self.rois(r))
                    M2 = M2 | createMask(self.rois(r), self.imH);
                end
            end
            self.putSliceBW(self.k, M2);
            delete(self.rois(ishandle(self.rois))); self.rois = gobjects(0);
            self.refreshOverlay(); self.updateNote(); self.update3D();
        end

        function clearSlice(self)
            if isempty(self.V), return; end
            self.putSliceBW(self.k, false(size(self.getImg(self.k))));
            delete(self.rois(ishandle(self.rois))); self.rois = gobjects(0);
            self.refreshOverlay(); self.updateNote(); self.update3D();
        end

        function finishPressed(self)
            if isempty(self.V), return; end
            self.acceptSlice();              % commit current edits
            self.mode = 'preview';
            % Simple preview loop
            keepEditing = true;
            while ishandle(self.f) && keepEditing
                choice = self.previewAndConfirm();
                switch choice
                    case "back"
                        keepEditing = false; % return to UI (still open)
                    case "discard"
                        close(self.f); return;
                    case "save"
                        [maskChosen, surfChosen, wasCancelled] = self.promptSavePaths();
                        if wasCancelled
                            keepEditing = false; % back to edit
                        else
                            self.lastSaveMask = maskChosen;
                            self.lastSaveSurf = surfChosen;
                            % Close UI and write outputs
                            if ishandle(self.f), close(self.f); end
                            self.writeOutputs(self.lastSaveMask, self.lastSaveSurf);
                            keepEditing = false;
                            return;
                        end
                end
            end
        end
    end

    %================== Rendering helpers ==================%
    methods (Access=private)
        function img = getImg(self, k), img = squeeze(self.V(:,:,k)); end

        function M2 = getSliceBW(self, k)
            M2 = squeeze(self.BW(:,:,k));
        end

        function putSliceBW(self, k, M2)
            self.BW(:,:,k) = logical(M2);
        end

        function refreshOverlay(self)
            delete(findobj(self.ax,'Type','Image','Tag','overlay'));
            M2 = self.getSliceBW(self.k);
            if any(M2(:))
                O = cat(3, zeros(size(M2)), ones(size(M2)), zeros(size(M2))); % green
                him = imshow(O,'Parent',self.ax); set(him,'AlphaData',0.3*M2,'Tag','overlay');
            end
            drawnow;
        end

        function T = getWorldAffine(self)
            % Returns a 4x4 voxel->world (mm) transform
            if isfield(self.info,'Transform') && ~isempty(self.info.Transform) && isfield(self.info.Transform,'T')
                T = self.info.Transform.T;
            else
                sp = self.info.PixelDimensions;
                T = diag([sp(:).', 1]);
            end
        end

        function update3D(self)
            if ~isempty(self.p3) && isgraphics(self.p3), delete(self.p3); self.p3 = []; end
            if isempty(self.V) || nnz(self.BW)==0, drawnow; return; end
            try
                [faces, verts] = isosurface(self.BW, 0.5);
                T = self.getWorldAffine();
                verts_h = [verts, ones(size(verts,1),1)];
                verts_w = verts_h * T';
                verts_mm = verts_w(:,1:3);
                self.p3 = patch(self.ax3,'Faces',faces,'Vertices',verts_mm, ...
                       'FaceColor',[0.85 0.2 0.2],'EdgeColor','none','FaceAlpha',0.85);
                camlight(self.ax3,'headlight'); lighting(self.ax3,'gouraud'); view(self.ax3,3);
            catch ME
                warning('LivePreview:Failed','%s', ME.message);
            end
            drawnow;
        end

        function updateTitle(self)
            set(self.ttl,'String',sprintf('%s slice %d/%d', self.capitalize(self.orientation), self.k, self.nSlices));
        end

        function updateNote(self, msg)
            if nargin==2
                set(self.note,'String',msg); return;
            end
            if isempty(self.info), set(self.note,'String',''); return; end
            vv = prod(self.info.PixelDimensions);
            txt = sprintf('Accepted voxels (slice %d): %d   |   Total voxels: %d   (Volume %.1f mm^3 / %.3f mL)', ...
                self.k, nnz(self.getSliceBW(self.k)), nnz(self.BW), nnz(self.BW)*vv, (nnz(self.BW)*vv)/1000);
            set(self.note,'String',txt);
        end

        function pos = roibounds2pos(~, B) % [row col] -> [x y]
            pos = [B(:,2), B(:,1)];
        end
    end

    %================== Preview & Save ==================%
    methods (Access=private)
        function choice = previewAndConfirm(self)
            % Returns "save"|"back"|"discard"
            choice = "back";
            pf = figure('Name','NSVG Preview (Accept & Write / Back / Discard)', ...
                        'Color','w','NumberTitle','off', ...
                        'Units','normalized','Position',[0.20 0.10 0.60 0.80], ...
                        'WindowStyle','modal');

            ax = axes('Parent',pf); hold(ax,'on'); axis(ax,'equal'); axis(ax,'vis3d'); view(ax,3);
            title(ax,'Lesion preview (world mm)'); xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z'); box(ax,'on'); grid(ax,'on');

            if nnz(self.BW) > 0
                [faces, verts] = isosurface(self.BW, 0.5);
                T = self.getWorldAffine();
                verts_h = [verts, ones(size(verts,1),1)];
                verts_w = verts_h * T';
                verts_mm = verts_w(:,1:3);
                patch(ax,'Faces',faces,'Vertices',verts_mm, ...
                      'FaceColor',[0.85 0.2 0.2],'EdgeColor','none','FaceAlpha',0.9);
                camlight(ax,'headlight'); lighting(ax,'gouraud');
                voxVol = prod(self.info.PixelDimensions); vol_mm3 = nnz(self.BW)*voxVol;
                annotation(pf,'textbox',[0.02 0.92 0.96 0.06],'String', ...
                    sprintf('Voxels: %d   Volume: %.1f mm^3 (%.3f mL)   Faces: %d', nnz(self.BW), vol_mm3, vol_mm3/1000, size(faces,1)), ...
                    'EdgeColor','none','HorizontalAlignment','center');
            else
                text(ax,0.5,0.5,0.5,'Mask is empty','Units','normalized','HorizontalAlignment','center','FontSize',12);
            end

            uicontrol('Style','pushbutton','Units','normalized','Position',[0.18 0.02 0.18 0.06], ...
                'String','Accept & Write','FontWeight','bold','Callback',@onAccept);
            uicontrol('Style','pushbutton','Units','normalized','Position',[0.41 0.02 0.18 0.06], ...
                'String','Back to Edit','Callback',@onBack);
            uicontrol('Style','pushbutton','Units','normalized','Position',[0.64 0.02 0.18 0.06], ...
                'String','Discard & Exit','Callback',@onDiscard);

            uiwait(pf);
            if isappdata(pf,'ch'), choice = string(getappdata(pf,'ch')); end
            if ishandle(pf), close(pf); end

            function onAccept(~,~), setappdata(pf,'ch','save');    uiresume(pf); end
            function onBack(~,~),   setappdata(pf,'ch','back');    uiresume(pf); end
            function onDiscard(~,~),setappdata(pf,'ch','discard'); uiresume(pf); end
        end

        function [maskPath, surfPath, cancelled] = promptSavePaths(self)
            cancelled = true;
            % sensible defaults near the input (if we can infer)
            dfltMask = self.lastSaveMask; dfltSurf = self.lastSaveSurf;
            try
                [refDir, refBase, ~] = fileparts(self.info.Filename);
                if isempty(dfltMask), dfltMask = fullfile(refDir,[refBase '_lesionMask.nii.gz']); end
                if isempty(dfltSurf), dfltSurf = fullfile(refDir,[refBase '_lesionSurface.mat']); end
            catch
                if isempty(dfltMask), dfltMask = fullfile(pwd,'lesionMask.nii.gz'); end
                if isempty(dfltSurf), dfltSurf = fullfile(pwd,'lesionSurface.mat'); end
            end

            [mf, md] = uiputfile( ...
                {'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'; '*.nii','NIfTI (*.nii)'; '*.nii.gz','NIfTI GZip (*.nii.gz)'}, ...
                'Save lesion mask as…', dfltMask);
            if isequal(mf,0) || isequal(md,0), maskPath=''; surfPath=''; return; end
            maskPath = NSVG.normalizeMaskPath(fullfile(md,mf));

            [sf, sd] = uiputfile({'*.mat','MAT-file (*.mat)'}, 'Save lesion surface as…', dfltSurf);
            if isequal(sf,0) || isequal(sd,0), maskPath=''; surfPath=''; return; end
            surfPath = fullfile(sd,sf);

            cancelled = false;
        end

        function writeOutputs(self, outMaskNifti, outSurfMat)
            % Build header like ref, but uint8 label data
            maskInfo = self.info;
            maskInfo.ImageSize    = size(self.BW);
            maskInfo.Datatype     = 'uint8';
            maskInfo.BitsPerPixel = 8;
            if isfield(maskInfo,'raw')
                maskInfo.raw.dim(2:4) = size(self.BW);
                maskInfo.raw.datatype = 2;
                maskInfo.raw.bitpix   = 8;
            end

            % Ensure dirs
            [md,~,~] = fileparts(outMaskNifti); if ~isempty(md) && ~exist(md,'dir'), mkdir(md); end
            [sd,~,~] = fileparts(outSurfMat);   if ~isempty(sd) && ~exist(sd,'dir'), mkdir(sd); end

            % Choose compression by extension
            try
                isGz = endsWith(lower(outMaskNifti),'.nii.gz');
                niftiwrite(uint8(self.BW), outMaskNifti, maskInfo, 'Compressed', isGz);
            catch ME
                % Retry uncompressed .nii
                alt = regexprep(outMaskNifti,'\.nii\.gz$','.nii','ignorecase');
                niftiwrite(uint8(self.BW), alt, maskInfo, 'Compressed', false);
                outMaskNifti = alt;
                warning('NSVG:WriteCompressedFailed','%s\nWrote uncompressed NIfTI instead: %s', ME.message, alt);
            end

            % Surface (world mm)
            faces = []; verts_mm = [];
            if nnz(self.BW) > 0
                [faces, verts] = isosurface(self.BW, 0.5);
                T = self.getWorldAffine();
                verts_h = [verts, ones(size(verts,1),1)];
                verts_w = verts_h * T';
                verts_mm = verts_w(:,1:3);
            end
            save(outSurfMat, 'faces', 'verts_mm');

            % expose to base
            assignin('base','lesion_surface_fv',struct('faces',faces,'vertices',verts_mm));

            fprintf('[NSVG] Saved mask: %s\n', outMaskNifti);
            fprintf('[NSVG] Saved surf: %s\n', outSurfMat);
        end
    end

    %================== small utils ==================%
    methods (Access=private, Static)
        function s = capitalize(t), s = lower(t); if ~isempty(s), s(1)=upper(s(1)); end; end
        function icon = local_poly_icon()
            try
                d = load(fullfile(matlabroot,'toolbox','matlab','icons','tool_polygon.mat'));
                icon = d.cdata;
            catch
                icon = zeros(16,16,3); icon(2:15,2:15,:) = 0.2;
            end
        end
        function outPath = normalizeMaskPath(inPath)
            % Ensure exactly one of: .nii  or .nii.gz
            p = inPath;
            low = lower(p);

            % already correct?
            if endsWith(low,'.nii.gz') || endsWith(low,'.nii')
                outPath = p; return;
            end

            % if user typed just ".gz" (e.g., "test.gz"), convert to ".nii.gz"
            if endsWith(low,'.gz')
                p = regexprep(p,'\.gz$','', 'ignorecase');
                outPath = [p '.nii.gz'];
                return;
            end

            % if user typed some other extension, replace it with ".nii.gz"
            [folder, base, ~] = fileparts(p);
            outPath = fullfile(folder, [base '.nii.gz']);
        end
    end
end
