---
description: Derived from Kristian Knudsen's PyEIS library.
---

# README

## KYLER'S SPEED CIRCUIT

Kyler's Speed Circuit is an Electrochemical Impedance Spectroscopy Analyzer that is tailored to a specific circuit\(Kyler's Circuit\) and fitted with said circuit. Through a customizeable process, whether one wants everything to be automated, or if one wants to look closely at a specific dataset or graph, Kyler' Speed Circuit can cover all essential functionality. It is derived from Kristian Knudsen's PyEIS library, which has a plethora of different circuits, most of which we do not need. 

## Update to come

We are working to secede from the PyEIS Library in order to maintain our own dependencies. We are also working to also adjust functioning in the Libraries to customize our own needs.

## Update as of 3/29
Currently, we are using a different version of the fitter, we moved over to a C# version application on .NET and are implementing all the functions there. We will still be using this version to fit simple values but will be moving our codebase over to a different repository.

Python notebooks should still be functioning and guiding. 
New Scripts have been written for the C# Application to use, as we wrapped the python code in a C# document to call the terminal.

