from ete3 import NCBITaxa

ncbi = NCBITaxa()

# Read serotype names from file
input_file = "dir_serotypes.txt"
output_file = "serotype_taxids.txt"

with open(input_file, "r") as f:
    serotypes = [line.strip() for line in f if line.strip()]

# Always prepend "Salmonella enterica subsp. enterica "
species_names = [f"Salmonella enterica subsp. enterica serovar {serotype}" for serotype in serotypes]

# Get taxonomic IDs
taxid_mapping = ncbi.get_name_translator(species_names)

# Save results
with open(output_file, "w") as f:
    for species, taxids in taxid_mapping.items():
        serotype_name = species.replace("Salmonella enterica subsp. enterica serovar ", "")
        f.write(f"{taxids[0]}\t{serotype_name}\n")
        print(f"{species}: {taxids[0]}")

print(f"TaxIDs saved to {output_file}")
