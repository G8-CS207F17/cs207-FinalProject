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
    reactions_dict = []

    for child in root.find('reactionData').findall('reaction'):
        new_dict = {}
        new_dict['id'] = child.get('id')
        new_dict['reactants'] = []
        new_dict['r_mass'] = []
        new_dict['products'] = []
        new_dict['p_mass'] = []
        for reactant in child.find('reactants').text.split():
            key, value = reactant.split(':')
            new_dict['reactants'].append(key)
            new_dict['r_mass'].append(float(value))

        for product in child.find('products').text.split():
            key, value = product.split(':')
            new_dict['products'].append(key)
            new_dict['p_mass'].append(float(value))
        
        # Parse reaction coefficients
        coeffs = child.find('rateCoeff')
        if coeffs.find('Nuclear'):
            try:
                new_dict['halfLife'] = float(coeffs.find('Nuclear').find('halfLife').text)
            except:
                raise ValueError('Missing halfLife. Please check your XML file.')

        reactions_dict.append(new_dict)
    return reactions_dict


file = 'tests/rxns_nuclear.xml'
reactions_dict = parseXML(file)
print(reactions_dict)
