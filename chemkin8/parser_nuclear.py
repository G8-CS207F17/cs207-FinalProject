import xml.etree.ElementTree as ET
from chemkin8.chemkin import *
import numpy as np

def parseXML(file):
    """
    Parses the XML file containing reactions information. Returns reaction coefficients and reaction rates.
    
    INPUTS
    ========
    file: string, required
          path of the XML file
          
    RETURNS
    ========
    species_lst: list of species to fix order
    reactions_dict: dictionary containing all relevant information for calculating reaction rates
    """
    
    tree = ET.parse(file)
    root = tree.getroot()

    # initialize variables
    reactions_dict = {'reactants': {}, 'products': {}, 'rates': []}
    coeffs_dict = []
    species_lst = []
    reversible_lst = []

    for species in root.find('phase').find('speciesArray').text.split():
        species_lst.append(species)
        reactions_dict['reactants'][species] = []
        reactions_dict['products'][species] = []

    for child in root.find('reactionData').findall('reaction'):

        # Record whether reaction is reversible
        reversible_lst.append(child.get('reversible'))

        # Parse reactants and products
        species = list(reactions_dict['reactants'].keys())
        for reactant in child.find('reactants').text.split():
            key, value = reactant.split(':')
            reactions_dict['reactants'][key].append(float(value))
            species.remove(key)
        for reactant in species:
            reactions_dict['reactants'][reactant].append(0.0)

        species = list(reactions_dict['products'].keys())
        for product in child.find('products').text.split():
            key, value = product.split(':')
            reactions_dict['products'][key].append(float(value))
            species.remove(key)
        for product in species:
            reactions_dict['products'][product].append(0.0)

        # Parse reaction coefficients
        coeffs = child.find('rateCoeff')
        if coeffs.find('Nuclear'):
            new_dict = {}
            new_dict['type'] = 'Nuclear'
            try:
                new_dict['halfLife'] = float(coeffs.find('Nuclear').find('halfLife').text)
            except:
                raise ValueError('Missing halfLife. Please check your XML file.')
            reactions_dict['rates'].append(new_dict)

    return species_lst, reactions_dict, reversible_lst


file = 'tests/rxns_nuclear.xml'
species_lst, reactions_dict, reversible = parseXML(file)
print(species_lst)
print(reactions_dict)
print(reversible)
