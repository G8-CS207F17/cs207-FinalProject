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

For example, if the python distribution is installed with anaconda, then the package will be installed at `/Users/myid/anaconda/lib/python3.6/site-packages` where `myid` is the user name of the laptop. Then, open terminal, run `cd /Users/myid/anaconda/lib/python3.6/site-packages`, followed by `pytest`.




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
