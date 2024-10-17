# Define the mass table for the amino acids
mass_table = {
    "A": 71,
    "C": 103,
    "D": 115,
    "E": 129,
    "F": 147,
    "G": 57,
    "H": 137,
    "I": 113,
    "K": 128,
    "L": 113,
    "M": 131,
    "N": 114,
    "P": 97,
    "Q": 128,
    "R": 156,
    "S": 87,
    "T": 101,
    "V": 99,
    "W": 186,
    "Y": 163,
}


# Function to generate the spectrum
def generate_spectrum(peptide):
    prefix_masses = [0]  # Start with mass 0
    for i in range(len(peptide)):
        prefix_masses.append(prefix_masses[-1] + mass_table[peptide[i]])

    # Create the spectrum by considering both linear and circular subpeptides
    peptide_mass = prefix_masses[-1]  # The total mass of the peptide
    spectrum = [0]

    # Linear subpeptides
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            spectrum.append(prefix_masses[j] - prefix_masses[i])

    # Circular subpeptides
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide)):
            spectrum.append(peptide_mass - (prefix_masses[j] - prefix_masses[i]))

    return sorted(spectrum)


peptide = "AMGY"
spectrum = generate_spectrum(peptide)
print(spectrum)
