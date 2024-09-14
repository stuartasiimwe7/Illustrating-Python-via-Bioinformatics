import re
# Parse spectral file and peptide sequence file
def parse_input_files(spectral_file, peptide_file):
    spectra = {}
    with open(spectral_file, 'r') as f:
        spectrum_lines = f.read().splitlines()
    i = 0
    while i < len(spectrum_lines):
        if spectrum_lines[i].startswith('BEGIN IONS'):
            scan_num = None
            precursor_mz = None
            charge = None
            peaks = []
            i += 1
            while i < len(spectrum_lines) and not spectrum_lines[i].startswith('END IONS'):
                line = spectrum_lines[i]
                if line.startswith('TITLE='):
                    scan_num = int(re.findall(r'Scan (\d+)', line)[0])
                elif line.startswith('PEPMASS='):
                    precursor_mz = float(line.split('=')[1])
                elif line.startswith('CHARGE='):
                    charge = int(line.split('=')[1].replace('+', ''))
                elif line.count('\t') == 1:
                    mz, intensity = line.split('\t')
                    peaks.append((float(mz), float(intensity)))
                i += 1
            spectra[scan_num] = {'precursor_mz': precursor_mz, 'charge': charge, 'peaks': peaks}
    
    peptides = []
    with open(peptide_file, 'r') as f:
        peptide_lines = f.read().splitlines()
    for line in peptide_lines:
        spec_file, scan_num, precursor_mz, charge, sequence = line.split('\t')
        scan_num = int(scan_num)
        precursor_mz = float(precursor_mz)
        charge = int(charge)
        peptides.append({'spec_file': spec_file, 'scan_num': scan_num, 'precursor_mz': precursor_mz, 'charge': charge, 'sequence': sequence})
    
    return spectra, peptides

# Generate theoretical peak list for a given peptide
def generate_theoretical_peak_list(peptide, ion_type, fragment_mass_tolerance):
    aa_masses = {
        'G': 57.02146,
        'A': 71.03711,
        'S': 87.03203,
        'P': 97.05276,
        'V': 99.06841,
        'T': 101.04768,
        'C': 103.00919,
        'L': 113.08406,
        'I': 113.08406,
        'X': 113.08406, # Unknown amino acid, treated as Leu
        'N': 114.04293,
        'D': 115.02694,
        'Q': 128.05858,
        'K': 128.09496,
        'E': 129.04259,
        'M': 131.04049,
        'H': 137.05891,
        'F': 147.06841,
        'R': 156.10111,
        'Y': 163.06333,
