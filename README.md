# cs207-FinalProject
Final Project Repo for G8

[![Build Status](https://travis-ci.org/G8-CS207F17/cs207-FinalProject.svg?branch=master)](https://travis-ci.org/G8-CS207F17/cs207-FinalProject)
[![Coverage Status](https://coveralls.io/repos/github/G8-CS207F17/cs207-FinalProject/badge.svg?branch=master)](https://coveralls.io/github/G8-CS207F17/cs207-FinalProject?branch=master)

Introduction
------------
Our package is to compute the reaction rate of chemical reactions. We take inputs of chemical reactions from an .xml file. We extract the parameters of each chemical reactions. We have functions to calculate three types of reaction rate coefficients: Constant reaction rate coefficients, Arrhenius reaction rate coefficients, Modified Arrhenius (a.k.a Kooij) reaction rate coefficients. Our package can handle elementary reactions as well as irreversible reactions. It can also handle an arbitary number of species and reactions.


Installation
------------
Download from https://github.com/G8-CS207F17/cs207-FinalProject. 

Packages need to be installed: numpy, xml, pytest (pip install numpy xml pytest)

run unit test suites: pytest --cov --cov-report term-missing

run doctest suites: pytest --doctest-modules --cov --cov-report term-missing



Basic Usage and Examples
------------------------
Provide a few examples on using your software in some common situations.  You may want to show how the code works with a small set of reactions.
