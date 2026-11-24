from pathlib import Path
root=Path(__file__).resolve().parents[1]
app=(root/'FFCProcessingCVI.m').read_text(errors='ignore')
exp=(root/'dicom_handling'/'ffc_export_enhanced_mr.m').read_text(errors='ignore')
meta=(root/'dicom_handling'/'ffc_extract_dicom_metadata.m').read_text(errors='ignore')
reader=(root/'file_handling'/'read_raw_data_Cameleon.m').read_text(errors='ignore')
core=(root/'ImageReconCore.m').read_text(errors='ignore')
launcher=(root/'ffc_start_cvi.m').read_text(errors='ignore')
checks={
 'TEST20 launcher':'TEST20' in launcher,
 'TEST20 title':'TEST 20' in app,
 'batch choice':'All participant sequences' in app,
 'raw-only tooltip':'raw dataset metadata only' in app,
 'shared study UID':'StudyInstanceUID' in app and 'studyUID' in app,
 'shared frame UID':'FrameOfReferenceUID' in app and 'frameUID' in app,
 'serie metadata parser':'getSerieMetadataFromSerieXML' in reader,
 'reader retains dataset path':'par.datasetPath = data_path' in reader,
 'core retains source metadata':'source_metadata = struct' in core,
 'metadata helper':'ffc_extract_dicom_metadata' in exp,
 'patient name written':"hex2dec('0010'),hex2dec('0010')" in exp,
 'patient id written':"hex2dec('0010'),hex2dec('0020')" in exp,
 'birth date written':"hex2dec('0010'),hex2dec('0030')" in exp,
 'sex written':"hex2dec('0010'),hex2dec('0040')" in exp,
 'participant fallback':'participantLabel' in meta,
 'DOB candidates':'DATE_OF_BIRTH' in meta,
 'sex candidates':'SUBJECT_GENDER' in meta,
 'stable study key':'StudyKey' in meta,
 'no runtime dicom dependency':'dicominfo' not in meta.lower() and 'dicomread' not in meta.lower(),
 'implementation version': 'ABDNFCI_TEST20' in exp,
}
failed=[k for k,v in checks.items() if not v]
print({'checks_passed':sum(checks.values()),'checks_total':len(checks),'failed':failed})
raise SystemExit(1 if failed else 0)
