from pathlib import Path
r=Path(__file__).resolve().parents[1]
w=(r/'dicom_handling/ffc_export_enhanced_mr.m').read_text(errors='ignore')
m=(r/'dicom_handling/ffc_extract_dicom_metadata.m').read_text(errors='ignore')
a=(r/'FFCProcessingCVI.m').read_text(errors='ignore')
checks={
'RAH to LPS matrix': "diag([-1 -1 1])" in w,
'centre converted': 'rahToDicomLps * centreRAH' in w,
'basis converted': 'rahToDicomLps * basisRAH' in w,
'pixel flip both axes': 'flip(flip(frame,1),2)' in w,
'IOP from corrected directions': "iop = [readDirection(:).' phaseDirection(:).']" in w,
'no fake Patient ID': "patientID = 'FCI'" not in w,
'no fake Patient Name': "patientName = 'FCI^RECON'" not in w,
'raw XML reread': 'readRawXmlIdentity(scan)' in m,
'patient block ID support': "inPatientBlock && strcmp(nk,'ID')" in m,
'DOB extraction': "field='PatientBirthDate'" in m,
'sex extraction': "field='PatientSex'" in m,
'StudyID raw': 'rawMeta.StudyID' in w,
'raw SeriesID default': 'rawMeta.SeriesID' in w,
'AcquisitionDate tag': "hex2dec('0022')" in w,
'AcquisitionTime tag': "hex2dec('0032')" in w,
'TEST22 app': 'TEST 22' in a,
}
failed=[k for k,v in checks.items() if not v]
print(f'{len(checks)-len(failed)}/{len(checks)} checks passed')
if failed:
 print('FAILED:',*failed,sep='\n - '); raise SystemExit(1)
