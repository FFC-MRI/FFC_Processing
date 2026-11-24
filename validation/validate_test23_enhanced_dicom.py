from pathlib import Path
r=Path(__file__).resolve().parents[1]
w=(r/'dicom_handling/ffc_export_enhanced_mr.m').read_text(errors='ignore')
a=(r/'FFCProcessingCVI.m').read_text(errors='ignore')
s=(r/'ffc_start_cvi.m').read_text(errors='ignore')
checks={
    'uses prescription basis': 'geometry.prescriptionBasisPatient' in w,
    'no RAH/LPS reflection': 'rahToDicomLps' not in w and "diag([-1 -1 1])" not in w,
    'IOP direct prescription row/column': "iop = [readDirection(:).' phaseDirection(:).']" in w,
    'phase reversed only': 'frame = flip(frame,2);' in w and 'flip(flip(frame,1),2)' not in w,
    'source slice order reversed': 'sourceSliceIndex = nSlices + 1 - sliceIndex' in w,
    '2D gap parameter used': "numericParameter(scan,'SPACING_BETWEEN_SLICE',NaN)" in w,
    '2D centre spacing thickness plus gap': 'rawSliceThicknessM + rawSliceGapM' in w,
    'plane positions use prescription normal': 'sliceCenter = centre + sliceOffset*sliceDirection' in w,
    'TEST23 implementation version': 'ABDNFCI_TEST23' in w,
    'TEST23 app': 'TEST 23' in a,
    'TEST23 launcher': 'TEST23' in s,
}
failed=[k for k,v in checks.items() if not v]
print(f'{len(checks)-len(failed)}/{len(checks)} checks passed')
if failed:
    print('FAILED:',*failed,sep='\n - ')
    raise SystemExit(1)
