from pathlib import Path
root=Path(__file__).resolve().parents[1]
w=(root/'dicom_handling/ffc_export_enhanced_mr.m').read_text()
a=(root/'FFCProcessingCVI.m').read_text()
s=(root/'ffc_start_cvi.m').read_text()
checks={
 'enhanced MR SOP': '1.2.840.10008.5.1.4.1.1.4.1' in w,
 'direct IOP': "iop = [readDirection(:).' phaseDirection(:).'];" in w,
 'direct slice indexing': 'sourceSliceIndex = sliceIndex;' in w,
 'both in-plane axes reversed': 'frame = flip(flip(frame,1),2);' in w,
 'no reverse slice source': 'nSlices + 1 - sliceIndex' not in w,
 'patient XYZ centre retained': 'centre = double(geometry.centerPatientMm(:));' in w,
 'per-frame patient position': "hex2dec('0032')" in w and 'imagePosition' in w,
 '2D gap plus thickness': 'rawSliceThicknessM + rawSliceGapM' in w,
 'TEST24 implementation version': 'ABDNFCI_TEST24' in w,
 'TEST24 app': 'TEST 24' in a,
 'TEST24 launcher': 'TEST24' in s,
}
failed=[k for k,v in checks.items() if not v]
print({'checks_passed':sum(checks.values()),'checks_total':len(checks),'failed':failed})
raise SystemExit(1 if failed else 0)
