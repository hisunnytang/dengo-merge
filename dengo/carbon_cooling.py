from reaction_classes import Species, ion_cooling_rate, species_registry

for i in range(7):
    ion_state = i + 1
    speciesName = "c_%s" % ion_state
    # Check if the species already exists
    # in the species registry, if it does
    # we don't want to create it again
    if (speciesName in species_registry) == False:
        s = Species(speciesName, 12, i)
    else:
        s = species_registry[speciesName]
    ion_cooling_rate(s)
