function report = cvi_test_import(workspaceFile,dicomRoot)
%CVI_TEST_IMPORT Command-line diagnostic for cvi42 contour import.
%
% Example:
%   ffc_setup;
%   report = cvi_test_import('patient.cvi42ws','D:\ClinicalDicoms');
%
% This does not need an FCI dataset.  It verifies workspace decoding,
% DICOM matching, patient-coordinate conversion, cardiac-phase grouping and
% surface construction before the contours are drawn by FFCProcessing.

    cvi = cvi_import_contours(workspaceFile,dicomRoot);
    report = struct();
    report.MatchedContours = numel(cvi.contours);
    report.UnmatchedContours = cvi.unmatchedContourCount;
    report.NumberOfSeries = numel(cvi.series);
    report.DefaultSeries = cvi.defaultSeriesIndex;
    report.Series = repmat(struct('Description','','ContourCount',0, ...
        'NumberOfFrames',0,'NumberOfPhases',0,'PhaseLabels',{{}}, ...
        'SurfaceCount',0,'ClinicalSliceIndices',[],'ClinicalStackNormalPatient',[NaN NaN NaN]),numel(cvi.series),1);

    fprintf('\ncvi42 import diagnostic\n');
    fprintf('  Matched contours:   %d\n',report.MatchedContours);
    fprintf('  Unmatched contours: %d\n',report.UnmatchedContours);
    fprintf('  DICOM source groups:%d\n',report.NumberOfSeries);
    for s = 1:numel(cvi.series)
        series = cvi.series(s);
        report.Series(s).Description = series.Description;
        report.Series(s).ContourCount = series.ContourCount;
        report.Series(s).NumberOfFrames = series.NumberOfFrames;
        report.Series(s).NumberOfPhases = series.NumberOfPhases;
        report.Series(s).PhaseLabels = series.PhaseLabels;
        report.Series(s).SurfaceCount = numel(series.Surfaces);
        report.Series(s).ClinicalSliceIndices = series.ClinicalSliceIndices;
        report.Series(s).ClinicalStackNormalPatient = series.ClinicalStackNormalPatient;
        marker = ' ';
        if s == cvi.defaultSeriesIndex, marker = '*'; end
        fprintf('%s Series %d: %s\n',marker,s,series.Description);
        fprintf('    contours=%d, frames=%d, phases=%d, surfaces=%d\n', ...
            series.ContourCount,series.NumberOfFrames,series.NumberOfPhases, ...
            numel(series.Surfaces));
        fprintf('    available contour phases: %s\n',strjoin(series.PhaseLabels,', '));
        fprintf('    clinical slice indices: %s\n',mat2str(series.ClinicalSliceIndices));
        fprintf('    clinical stack normal: %s\n',mat2str(series.ClinicalStackNormalPatient,4));
    end
    fprintf('  * = default series selected by the app\n\n');
end
