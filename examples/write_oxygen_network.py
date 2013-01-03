from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
import dengo.oxygen_rates, dengo.oxygen_cooling
from dengo.chemistry_constants import tiny, kboltz, mh

# If only a subset of species are wanted put them here
# and change the commented lines below
want = ('OI', 'OII', 'OIII', 'de', 'ge')

oxygen = ChemicalNetwork()
oxygen.add_energy_term()

for ca in cooling_registry.values():
    # The following line can be used to specify a subset of species
    #if all(sp.name in want for sp in ca.species):
    if ca.name.startswith("O"):
       oxygen.add_cooling(ca)

for s in reaction_registry.values():
    # The following line can be used to specify a subset of species
    # if all(sp.name in want for sp in s.species):
    if s.name.startswith("O"):
        oxygen.add_reaction(s)

# This defines the temperature range for the rate tables
oxygen.init_temperature((1e0, 1e8))

# Set to false if you don't want intermediate solution output
oxygen.write_intermediate_solutions = True

# Generate initial conditions (switch to False to disable this)
generate_initial_conditions = True
start_neutral = False
initial_state = 'OI' # if using start_neutral = True, this should be the _1 state
                      # if starting in a close to CIE state this doesn't really matter

if generate_initial_conditions:
    import numpy as np
    NCELLS = 32
    density = 1
    init_array = np.ones(NCELLS)
    init_values = dict()
    init_values['density'] = density * init_array
    init_values[initial_state] = init_array.copy() # use conservation to set this below
           
    # set up initial temperatures values used to define ge
    temperature = np.logspace(4, 6.7, NCELLS)
    #temperature[:] = 1e7; # need to remove this line for the above one to matter
    init_values['T'] = temperature

    if start_neutral:
        X = 1e-6
 
        # populate initial fractional values for the other species
        for s in sorted(oxygen.required_species):
            if s.name != 'ge' and s.name != initial_state:
                if s.name == 'de':
                    continue
                else:
                    init_values[s.name] = X * init_array
                init_values[initial_state] -= init_values[s.name]
        init_values['de'] = init_array * 0.0
        for s in sorted(oxygen.required_species):
            if s.name in ("ge", "de"): continue
            init_values['de'] += init_values[s.name] * s.free_electrons
            print "Adding %0.5e to electrons from %s" % (
                (init_values[s.name] * s.free_electrons)[0], s.name)
        print "Total de: %0.5e" % (init_values['de'][0])

    else:
        # Unless neutral is desired, we'll start in a CIE solution
        import chianti.core as ch
        import chianti.util as chu

        # Note: ChiantiPy can return a full array of CIE values
        # for a given species, while that would minimize the number
        # of times we would have to call upon ChiantiPy, which would
        # would make the IC generation faster, it might get tricky
        # for a combined network.  However, with careful thought
        # a solution could be found.
        # This would use code that looks something like:
        #     ion_frac = ch.ioneq(s.number, init_values['T'])
        #     ion_frac = ion_frac.Ioneq
        #     init_values[s.name] = ion_frac[s.free_electrons] * init_array

        # populate initial fractional values for the species
        for s in sorted(oxygen.required_species):
            if s.name != 'ge':
                if s.name == 'de':
                    continue
                else:
                    ion_name = chu.zion2name(s.number, s.free_electrons + 1)
                    ion = ch.ion(ion_name, temperature=init_values['T'])
                    ion.ioneqOne()
                    ion_frac = ion.IoneqOne
                    init_values[s.name] = ion_frac * init_array
                
                # # add some random noise
                # init_values[s.name] += 0.5 * init_values[s.name] * np.random.randn(NCELLS)
                
                # in case something is negative or super small:
                init_values[s.name][init_values[s.name] < tiny] = tiny
                print init_values[s.name]
                
                # # conservation...
                # init_values[initial_state] -= init_values[s.name]
        
        init_values['de'] = init_array * 0.0
        for s in sorted(oxygen.required_species):
            if s.name in ("ge", "de"): continue
            init_values['de'] += init_values[s.name] * s.free_electrons
            print "Adding %0.5e to electrons from %s" % (
                (init_values[s.name] * s.free_electrons)[0], s.name)
        print "Total de: %0.5e" % (init_values['de'][0])

    # convert to masses to multiplying by the density factor and the species weight
    for s in oxygen.required_species:
        if s.name == 'ge': continue
        if s.name == 'de':
            init_values[s.name] *= (density) # this is still just number density
        else:
            init_values[s.name] *= (density * s.weight)

    #compute new total density and number density
    density = init_array * 0.0
    number_density = init_array * 0.0
    for s in sorted(oxygen.required_species):
        if s.name == 'ge': continue
        number_density += init_values[s.name]/s.weight
        if s.name == 'de': continue
        density += init_values[s.name]
    init_values['density'] = density

    # calculate ge (very crudely)
    gamma = 5.e0/3.e0
    init_values['ge'] = ((temperature * number_density * kboltz)
                         / (density * mh * (gamma - 1)))

    # Write the initial conditions file
oxygen.write_solver("oxygen", output_dir = ".",
                    init_values=init_values)
