'''
The code below, tries not only to save the theoretical Peaks lists and Annotated Peak lists 
but also tries to do that in the specified in the assignnment. 
So to get rid of the error - this should be the end of the task - FIGHTING!!
'''

valid_aa_codes = set('ARNDCQEGHILKMFPSTWYV')
fragment_tol = 0.5
ion_types = ('b', 'y')

#generate all possible m/z for fragments of types ion_types and of charges from 1 to maxcharge
def fragments(peptide, ion_types=('b', 'y'), maxcharge=1):
    for i in range(1, len(peptide)):
        for ion_type in ion_types:
            for charge in range(1, maxcharge+1):
                if ion_type[0] in 'abc':
                    yield mass.fast_mass(peptide[:i], ion_type=ion_type, charge=charge)
                else:
                    yield mass.fast_mass(peptide[i:], ion_type=ion_type, charge=charge)

#Our Output files
theoretical_peak_list_file = open('output/Theoretical_peak_list.txt', 'w')
psm_annotated_peak_list_file = open('output/PSM_annotated_peak_list.txt', 'w')

#Reading PSM list file and MGF spectrum file
with mgf.read('data/hw2_test.mgf') as spectra, open('data/hw2_psmlist_test_v2.txt', 'r') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    next(reader, None)
    for row in reader:
        peptide = row[4]
        if not all(aa in valid_aa_codes for aa in peptide):
            continue

        # Generating Theoretical Peak Lists
        fragment_mz_values = list(fragments(peptide, ion_types=ion_types, maxcharge=2))
        #print(fragment_mz_values)
        #theoretical peak list to file
        theoretical_peak_list_file.write('BEGIN\n')
        theoretical_peak_list_file.write(f'PEPTIDE={peptide}\n')
        theoretical_peak_list_file.write(f'TITLE={row[1]}\n')
        
        for i, mz in enumerate(fragment_mz_values):
            ion_label = ion_types[i % len(ion_types)] + str((i // 2) + 1) + "+" * (i % 2 + i // 2) #generate a label for the ion
            theoretical_peak_list_file.write(f'{ion_label}\t{mz}\n')
        theoretical_peak_list_file.write('END\n\n')

        #Find matching peaks in the spectrum > write to annotated peak list file
        psm_annotated_peak_list_file.write('BEGIN\n')
        psm_annotated_peak_list_file.write(f'PEPTIDE={peptide}\n')
        psm_annotated_peak_list_file.write(f'TITLE={row[1]}\n')
        
        for mz, intensity in zip(spectra[0]['m/z array'], spectra[0]['intensity array']):
            if any(abs(mz - frag_mz) <= fragment_tol for frag_mz in fragment_mz_values):
                #ion_label = f'{ion_types[fragment_mz_values.index(frag_mz)%2]}{fragment_mz_values.index(frag_mz)//2+1}+'
               if mz in fragment_mz_values:
                    ion_label = f'{ion_types[fragment_mz_values.index(mz) % 2]}{fragment_mz_values.index(mz) // 2 + 1}+'
                    psm_annotated_peak_list_file.write(f'{ion_label}\t{mz}\t{intensity}\n')
        psm_annotated_peak_list_file.write('END\n\n')

theoretical_peak_list_file.close()
psm_annotated_peak_list_file.close()


valid_aa_codes = set('ARNDCQEGHILKMFPSTWYV')
tolerance = 0.5

with mgf.read('data/hw2_test.mgf') as spectra, open('data/hw2_psmlist_test_v2.txt', 'r') as tsvfile, open('theoretical_peak_list.txt', 'w') as theoretical_file, open('annotated_peak_list.txt', 'w') as annotated_file:
    spectrum = next(spectra)
    reader = csv.reader(tsvfile, delimiter='\t')
    next(reader, None)
    for row in reader:
        peptide = row[4]
        if not all(aa in valid_aa_codes for aa in peptide):
            continue
        #Gen Theoretical Peak Lists
        fragment_mz_values = fragments(peptide)

        theoretical_peaks = []
        for mz in fragment_mz_values:
            theoretical_file.write(str(mz) + '\n')
            theoretical_peaks.append(mz)

        annotated_peaks = []
        for i in range(len(spectrum['m/z array'])):
            peak = spectrum['m/z array'][i]
            for theoretical_peak in theoretical_peaks:
                if abs(peak - theoretical_peak) <= tolerance:
                    annotated_file.write(str(peak) + '\t' + str(spectrum['intensity array'][i]) + '\n')
                    annotated_peaks.append(peak)
                    break

        #Remove annotated peaks from theoretical peaks
        for peak in annotated_peaks:
            if peak in theoretical_peaks: