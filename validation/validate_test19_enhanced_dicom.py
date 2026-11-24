from pathlib import Path

root = Path(__file__).resolve().parents[1]
writer = (root / 'dicom_handling' / 'ffc_export_enhanced_mr.m').read_text()
app = (root / 'FFCProcessingCVI.m').read_text()
setup = (root / 'ffc_setup.m').read_text()
launcher = (root / 'ffc_start_cvi.m').read_text()

checks = {
    'Enhanced MR SOP Class UID': "1.2.840.10008.5.1.4.1.1.4.1" in writer,
    'Explicit VR Little Endian': "1.2.840.10008.1.2.1" in writer,
    'Shared Functional Groups': "hex2dec('5200'),hex2dec('9229')" in writer,
    'Per-frame Functional Groups': "hex2dec('5200'),hex2dec('9230')" in writer,
    'Pixel Measures FG': "hex2dec('0028'),hex2dec('9110')" in writer,
    'Frame Content FG': "hex2dec('0020'),hex2dec('9111')" in writer,
    'Plane Position FG': "hex2dec('0020'),hex2dec('9113')" in writer,
    'Plane Orientation FG': "hex2dec('0020'),hex2dec('9116')" in writer,
    'Frame Anatomy FG': "hex2dec('0020'),hex2dec('9071')" in writer,
    'Anatomic Region sequence': "hex2dec('0008'),hex2dec('2218')" in writer,
    'MR Image Frame Type FG': "hex2dec('0018'),hex2dec('9226')" in writer,
    'Pixel Value Transformation FG': "hex2dec('0028'),hex2dec('9145')" in writer,
    'Dimension Organization': "hex2dec('0020'),hex2dec('9221')" in writer,
    'Dimension Index Sequence': "hex2dec('0020'),hex2dec('9222')" in writer,
    'Acquisition Context': "hex2dec('0040'),hex2dec('0555')" in writer,
    'Research content qualification': "'RESEARCH'" in writer,
    'Unspecified rescale units (US)': "hex2dec('1054'),'LO','US'" in writer,
    'No obsolete RELATIVE rescale type': "'RELATIVE'" not in writer,
    'Field/time/slice dimension indices': "[fieldIndex temporalIndex sliceIndex]" in writer,
    'Exact FCI conditions in Frame Comments': "FCI evolution field=%s mT; evolution time=%s ms" in writer,
    'Native geometry helper': 'ffc_scan_spatial_geometry(scan)' in writer,
    'Magnitude pixel source': 'scan.magimage' in writer,
    'GUI button property': 'ExportEnhancedDICOMButton' in app,
    'GUI callback invokes writer': 'ffc_export_enhanced_mr(scan,outputFile)' in app,
    'DICOM handling folder on path': "'dicom_handling'" in setup,
    'TEST19 launcher': 'TEST19' in launcher,
}

failed = [name for name, ok in checks.items() if not ok]
for name, ok in checks.items():
    print(('PASS' if ok else 'FAIL') + ': ' + name)
print(f'\n{sum(checks.values())}/{len(checks)} static checks passed')
if failed:
    raise SystemExit('Failed checks: ' + ', '.join(failed))
