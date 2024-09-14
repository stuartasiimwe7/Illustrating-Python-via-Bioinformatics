from pyteomics import mgf, parser
from Bio.SeqUtils import ProtParam


# ###  Read in the Spectral File


spectral_file = "data/spectral_file.mgf"

with mgf.read(spectral_file) as spectra:
    for spectrum in spectra:
        # extract relevant information from each spectrum
        scan_num = spectrum['params']['title']
        precursor_mz = spectrum['params']['pepmass'][0]
        precursor_charge = spectrum['params']['charge'][0]
        peaks = spectrum['m/z array']
        intensities = spectrum['intensity array']
        # process peaks and intensities...


# ###  Read in the peptide sequence file 


peptide_file = 'data/psmlist.txt'

# Ion types to consider
ion_types = ['b', 'y']

# Fragment charge states to consider
charge_states = [1, 2]

# Fragment mass tolerance (in Da)
fragment_tol = 0.5

peptide_dict = {}
#Mass of proton is 1.0072766 a.m.u. or 1.6726 x 10^-27 kg.
mass_proton = 1.0072766


with open(peptide_file,'r') as f:
    next(f) # skip header row
    for line in f:
        line = line.strip()
        if not line:
            continue
        # split tab-separated line into columns
        SpecFile, ScanNum, PrecursorMZ, Charge, Peptide=line.split('\t')
        peptide_dict[Peptide] = (SpecFile, ScanNum, PrecursorMZ, Charge)
        print(peptide_dict)


for Peptide, (SpecFile, ScanNum, PrecursorMZ, Charge) in peptide_dict.items():
   # Print some information about the current Peptide
    print(f"Processing Peptide: {Peptide}, SpecFile: {SpecFile}, ScanNum: {ScanNum}")


from pyteomics import mass
# Generate theoretical peak list for each ion type and charge state
for ion_type in ion_types:
    for c in charge_states:
        ion_series = ion_type + str(c) + '+'
        ion_masses = mass.fast_mass(
            ion_series,
            seq=Peptide,
            ion_type=ion_type,
            Charge=c,
            aa_mass=mass.std_aa_mass,
            terminus_mass=mass.std_terminus_mass,
            ion_masses={'y': mass.fast_mass('y', seq=Peptide, ion_type='y', Charge=c, aa_mass=mass.std_aa_mass, terminus_mass=mass.std_terminus_mass)}
        )

        # Filter theoretical peaks based on fragment mass tolerance
        for i, ion_mass in enumerate(ion_masses):
            for neutral_loss in [0.0, mass.water, mass.ammonia]:
                ion_mass_nl = ion_mass - neutral_loss
                if abs(ion_mass_nl - PrecursorMZ) <= fragment_tol:
                    # Output annotated peak for current ion type, charge state, and neutral loss
                    print(f"{ion_series}\t{ion_mass_nl:.4f}\t{ion_mass_nl - PrecursorMZ:.4f}")


# ### Generate the theoretical peak list for each peptide:


def get_fragment_masses(Peptide, ion_type):
    """
    Calculate the expected masses of fragment ions for a given Peptide sequence and ion type
    """
    prot_param = ProtParam.ProteinAnalysis(str(Peptide))
    aa_masses = prot_param.monoisotopic_counts

    if ion_type == 'y':
        ion_masses = [sum(aa_masses[i:]) + 19.0178 for i in range(len(aa_masses))]  # add mass of H2O
    elif ion_type == 'b':
        ion_masses = [sum(aa_masses[:i]) + 1.0078 for i in range(len(aa_masses))]  # add mass of H

    return ion_masses


def get_peak_list(Peptide, ion_type, fragment_tol):
    """
    Generate a theoretical peak list for a given Peptide sequence, ion type, and fragment mass tolerance
    """
    fragment_masses = get_fragment_masses(Peptide, ion_type)
    peak_list = []

    for i, ion_mass in enumerate(fragment_masses):
        if i == 0:
            continue
        diff = ion_mass - fragment_masses[i-1]
        if abs(diff - 1.0078) <= fragment_tol:
            peak_list.append({'ion_type': ion_type, 'ion_num': i, 'mass': ion_mass, 'intensity': 1.0})

    return peak_list


theoretical_peak_lists = []
for Peptide in Peptide:
    peak_list = get_peak_list(Peptide['Peptide_seq'], ion_type, fragment_tol)
    theoretical_peak_lists.append({'scan_num': Peptide['scan_num'], 'peak_list': peak_list})



# ### Annotate the peaks in the spectral file:


from pyteomics import mass

# Iterate over each peptide in the peptide_dict
for peptide, (spec_file, scan_num, precursor_mz, precursor_charge) in peptide_dict.items():

    # Convert precursor_mz and precursor_charge to floats
    precursor_mz = float(precursor_mz)
    precursor_charge = int(precursor_charge)

    # Calculate the precursor mass from m/z and Charge
    precursor_mass = precursor_mz * precursor_charge - precursor_charge * mass.proton

    # Iterate over each ion type in ion_types
    for ion_type in ion_types:

        # Iterate over each charge state in charge_states
        for c in charge_states:

            # Generate ion series and masses
            ion_series = ion_type + str(c) + '+'
            ion_masses = mass.fast_mass(
                ion_series,
                seq=peptide,
                ion_type=ion_type,
                Charge=c,
                aa_mass=mass.std_aa_mass,
                terminus_mass=mass.std_terminus_mass,
                ion_masses={'y': mass.fast_mass('y', seq=peptide, ion_type='y', Charge=c, aa_mass=mass.std_aa_mass, terminus_mass=mass.std_terminus_mass)}
            )

            # Iterate over each ion mass and calculate neutral loss masses
            for ion_mass in ion_masses:
                for neutral_loss in [0.0, mass.water, mass.ammonia]:
                    ion_mass_nl = ion_mass - neutral_loss

                    # Filter theoretical peaks based on fragment mass tolerance
                    if abs(ion_mass_nl - precursor_mass) <= fragment_tol:
                        # Output annotated peak
                        print(f"{ion_series}\t{ion_mass_nl:.4f}\t{ion_mass_nl - precursor_mass:.4f}")



