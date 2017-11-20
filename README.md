# CS207 G8 Final Project: Chemical Kinetics

[![Build Status](https://travis-ci.org/G8-CS207F17/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/G8-CS207F17/cs207-FinalProject)
[![Coverage Status](https://coveralls.io/repos/github/G8-CS207F17/cs207-FinalProject/badge.svg?branch=master)](https://coveralls.io/github/G8-CS207F17/cs207-FinalProject?branch=master)


Introduction
------------
Our package computes the progress rates of of a system chemical reactions. The package extracts parameters for each chemical reaction from an `xml` file. We also provide functions to calculate three types of reaction rate coefficients (Constant, Arrhenius and Modified Arrhenius) as well as progress rate for a reaction set. Based on the reaction rate coefficients and progress rates, the reaction rate of a system of reactions can be obtained.

Our package can also handle an arbitrary number of species and reactions. It is compatible with elementary and irreversible reactions, and can be easily extended to reversible or non-elementary reactions.




Installation
------------
Install the package by running `pip install chemkin8` in terminal.

To run the test suite, either navigate to the folder where this package is installed or download this repo and navigate to the root folder, then run `pytest` in terminal.

For example, if the python distribution is installed with anaconda, then the package will be installed at `/Users/myid/anaconda/lib/python3.6/site-packages/chemkin8` where `myid` is the user name of the laptop. Then, open terminal, run `cd /Users/myid/anaconda/lib/python3.6/site-packages/chemkin8`, followed by `pytest`. If your default python version is not python3.6, for example, if your default python version is 3.5, then cd into your python3.5 folder instead.




Basic Usage and Examples
------------------------
Let us have a set of reactions:
```
H + O2 => OH + O
H2 + O => OH + H
H2 + OH => H2O + H
```

We can add information about these reactions in an xml file of the form:
```
<?xml version="1.0"?>

<ctml>

    <phase>
        <speciesArray> H O OH H2 H2O O2 </speciesArray>
    </phase>

    <reactionData id="test_mechanism">
        <!-- reaction 01  -->
        <reaction reversible="no" type="Elementary" id="reaction01">
            <equation>H + O2 =] OH + O</equation>
            <rateCoeff>
                <Arrhenius>
                    <A>3.52e+10</A>
                    <E>7.14e+04</E>
                </Arrhenius>
            </rateCoeff>
            <reactants>H:1 O2:1</reactants>
            <products>OH:1 O:1</products>
        </reaction>

        <!-- reaction 02 -->
        <reaction reversible="no" type="Elementary" id="reaction02">
            <equation>H2 + O =] OH + H</equation>
            <rateCoeff>
                <modifiedArrhenius>
                    <A>5.06e-2</A>
                    <b>2.7</b>
                    <E>2.63e+04</E>
                </modifiedArrhenius>
            </rateCoeff>
            <reactants>H2:1 O:1</reactants>
            <products>OH:1 H:1</products>
        </reaction>

        <!-- reaction 03 -->
        <reaction reversible="no" type="Elementary" id="reaction03">
            <equation>H2 + OH =] H2O + H</equation>
            <rateCoeff>
                <Constant>
                    <k>1.0e+03</k>
                </Constant>
            </rateCoeff>
            <reactants>H2:1 OH:1</reactants>
            <products>H2O:1 H:1</products>
        </reaction>
    </reactionData>

</ctml>
```

We import the chemical kinetics library as follows:
```
from chemkin8 import chemkin
```

We can initialize the object and parse the XML file by passing in the path to the XML file:
```
testcase1 = chemkin.chemkin('path-to-xml')
```

We input the concentrations of each species and the temperature at which reactions occur. The reaction rates of each species in the system can be given as:
```
x = [1,1,1,1,1,1]
T = 1500
testcase1.reaction_rates(x, T)
```

The result is given as:
```
[-227364086.53073898, 227364586.53073898, 231985198.37073097, -2311055.9199959813, 500.0, -229675142.45073497]
```

Reversible Reactions:
```
Use same functions and methods as basic reactions.
The tag <reaction reversible = ""> in your xml file determines which functions will run in the background.

If tag <reaction reversible="yes"...> : calculations will be for reversible reactions.
If tag <reaction reversible="no"...>: calculations will be for non-reversible reactions.
```

Future Features
---------------
**Motivation:**
The aim of including these features is to expand the types of reactions that our `chemkin` library can process and generate information about.

**Suggested Features:**
1. Support for nuclear reactions
2. Support for radioactive entities
3. Generate half-life graphs for species

**Description of Features and Steps to add them:**

1. Nuclear reactions - We want to detect the type of nuclear reaction and generate the decay sequence until a stable product is obtained.:
    - Maintain database of radioactive elements
    - Check against database for atomic weight
        - Database contains atomic weight and atomic number
        - Also contains nature of nuclei - stable or unstable
          If product is stable, decay is halted
          else, continue reducing to stable products (decay series - more in step 3)

    - Detect type of reaction out of
        - alpha decay
        - beta decay
        - positron emission
        - electron capture
        - gamma emission
        - spontaneous fission
          which produces stable nuclei. Generate complete reactions, graphs, etc.

    - Decay series: produce set of reaction series generated to reach stable radioactive nuclei
    - Generate reactions and visualisation of decay, half live graphs of radioactive entities in products/series etc

2. Radioactive entities - If simply a radioactive element is provided, we generate the decay sequence as described for nuclear reactions.
    - If simply a radioactive entity is provided, its decomposition reaction/series can be generated
    - Reactions and graphs at each step are illustrated.

3. Visualisation - We want to generate graphs showing the decay rate of each reacting specie. This function/class can be called for the existing reaction sets also.
    - Calculate reaction rates for each entity
    - Generate decay graphs for products and reactants for reversible and irreversible reactions
    - Generate half-life graph for radioactive particles
