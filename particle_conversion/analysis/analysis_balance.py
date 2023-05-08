import numpy as np
import cantera as ct

#========#
#========#

def read_data(filename):

    output = []

    with open(filename) as f:
        for line in f:

            tmp = [n for n in line.split()]
            output.append(tmp)

    return output

#========#
#========#

#####################################################
## START
#####################################################


gas = ct.Solution('gri30.xml')
avail_species = gas.species_names
elements = gas.element_names
#elements.remove('Ar')
#elements.append('AR')

species_dict = {}
data_dict = {}
heat_dict = {}
threshold = 1e-6
pressure0 = 20e5


# Run 1 with data filling data_dict and species_dict
#========================================================#
data = read_data('output-refactored.out')

for i,val in enumerate(data):

    val[1] = val[1].upper()
    val[2] = float(val[2])

    if val[2] < threshold:
        val[2] = 0.0

    if val[1] == 'CH2<S>':
        val[1] = 'CH2(S)'

    if val[0] in species_dict.keys() and val[1] in avail_species:
        species_dict[val[0]].update({val[1]:val[2]})
    elif val[1] in avail_species:
        species_dict.update({val[0]:{}})
        species_dict[val[0]].update({val[1]:val[2]})
    #-------
    elif val[0] in data_dict.keys() and val[1] not in avail_species:
        data_dict[val[0]].update({val[1]:val[2]})
    elif val[1] not in avail_species:
        data_dict.update({val[0]:{}})
        data_dict[val[0]].update({val[1]:val[2]})
#========================================================#


# Run 2 reading gas mass flows
#========================================================#
data = read_data('output-mass-balance.out')

for i,val in enumerate(data):

    if len(val) > 0:
        if val[0] in data_dict.keys():
            data_dict[val[0]].update({'MFLOW':float(val[1])})
#========================================================#


# Run 2.1 reading
#========================================================#
data = read_data('output-heat-balance.out')

for i,val in enumerate(data):

    if len(val) == 2 and not val[0].startswith('"') and not val[0].startswith('---'):
        heat_dict.update({val[0]:float(val[1])})
#========================================================#


# Run 3 Getting DPM Data
#========================================================#
data = read_data('output-dpm.out')
flag_outlet = False
flag_coal = False
tmp = []
dpm_outlet_dict = {'C':0}   # MFLOW kg/s
coal_dict = {
    'C':0.837,          # daf
    'H':0.048,          # daf
    'O':0.091,          # daf
    'N':0.024,          # daf
    'FC':0.5428,        # daf
    'VOL':0.4572,       # daf
    'ASH':0.0974,       # dry
    'MOIST':0.045,      # ar
    'MFLOW':0}          # ar

coal_moist = {'TEMP':300,'H2O':1}
# Moist

# Change ulti and proxy to (ar) status
coal_dict['C']  = coal_dict['C']    * (1-coal_dict['ASH']) * (1-coal_dict['MOIST'])
coal_dict['H']  = coal_dict['H']    * (1-coal_dict['ASH']) * (1-coal_dict['MOIST'])
coal_dict['O']  = coal_dict['O']    * (1-coal_dict['ASH']) * (1-coal_dict['MOIST'])
coal_dict['N']  = coal_dict['N']    * (1-coal_dict['ASH']) * (1-coal_dict['MOIST'])
coal_dict['FC']  = coal_dict['FC']  * (1-coal_dict['ASH']) * (1-coal_dict['MOIST'])
coal_dict['VOL'] = coal_dict['VOL'] * (1-coal_dict['ASH']) * (1-coal_dict['MOIST'])
coal_dict['ASH'] = coal_dict['ASH'] * (1-coal_dict['MOIST'])

# Basically you may assume that all materials, accept the unreacted carbon is put into the system
for i,val in enumerate(data):

    if len(val) > 0:

        if 'Mass' in val and 'Flow' in val:
            flag_coal = True

        if len(val) <= 8 and flag_coal == True:

            if val[0].startswith('Escaped'):
                for j, stringdpm in enumerate(val):
                    tmp.append(stringdpm)

                coal_dict['MFLOW'] = float(tmp[-3])
                flag_coal = False


        if 'Species' in val and 'Content' in val:
            flag_outlet=True

        if len(val) <= 8 and flag_outlet==True:

            if val[0].startswith('Escaped'):
                for j,stringdpm in enumerate(val):
                    if stringdpm.startswith('fc'):
                        dpm_outlet_dict['C'] += float(val[j+2])

# Generating a H2O source from Moist
data_dict.update({'coal_moist':{'H2O':1,'TEMP':300,'MFLOW':coal_dict['MOIST']*coal_dict['MFLOW']}})
species_dict.update({'coal_moist':{'H2O':1}})

# Generating a CHON source from Vol
vol_C = coal_dict['C'] - coal_dict['FC']
vol_H = coal_dict['H']
vol_O = coal_dict['O']
vol_N = coal_dict['N']

data_dict.update({'coal_vol':{'C':vol_C,'H':vol_H,'O':vol_O,'N':vol_N,'TEMP':300,'MFLOW':coal_dict['VOL']*coal_dict['MFLOW']}})
species_dict.update({'coal_vol':{'C':vol_C,'H':vol_H,'O':vol_O,'N':vol_N}})

# Generating a CHON source from Char
char_C = coal_dict['FC']

data_dict.update({'coal_char':{'C':char_C,'TEMP':300,'MFLOW':coal_dict['FC']*coal_dict['MFLOW']}})
species_dict.update({'coal_char':{'C':char_C}})


# Run 4 Extrating Difference (coal input required)
#========================================================#
element_dict = {}
element_diff = {}
element_diff_rel = {}

for i,valA in enumerate(species_dict.keys()):

    # todo gas species in mass fraction please
    element_dict.update({valA:{}})
    gas.TPY = data_dict[valA]['TEMP'], pressure0, species_dict[valA]

    #-------------------------------------------------------------#
    for k,valB in enumerate(elements):
        # Inlets and Outlets
        element_dict[valA].update({valB:gas.elemental_mass_fraction(valB) * data_dict[valA]['MFLOW']})
    #-------------------------------------------------------------#

        if i == 0:
            element_diff.update({valB:0})
            element_diff_rel.update({valB:0})
            coal_dict.update({valB:0})



for k, valB in enumerate(elements):
    for i,valA in enumerate(species_dict.keys()):
        element_diff[valB] += element_dict[valA][valB]

    if valB == 'C':
        element_diff[valB] = element_diff[valB] - dpm_outlet_dict[valB]

    element_diff_rel[valB] = element_diff[valB]/element_dict['outlet'][valB]



#coal_dict.update({'Ar':0})

#gas.TPY = 300,1e5,coal_dict
#gas()

#gas.TPX = 1273.15,20e5,species_dict['inlet_gas']
#gas.TPX = 1273.15,20e5,species_dict['inlet_coal']
#gas.TPX = 1273.15,20e5,species_dict['outlet']

#print(element_dict['outlet'])
#print(element_diff)

sum = 0
for i,val in enumerate(element_diff):
    sum += element_diff[val]

for i,val in enumerate(element_diff_rel):
    print(val,element_diff_rel[val])

print('Heat-Error-rel',heat_dict['Net']/heat_dict['outlet'])
print('Mass-tot',sum)
print('Mass-rel',-sum/data_dict['outlet']['MFLOW'])
print('--finished--')

