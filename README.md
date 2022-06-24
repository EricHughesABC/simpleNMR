# simpleNMR

simpleNMR is a python/QT5/matplotlib program. It runs on windows and macos. The aim of the software is to present the information derived from processing NMR spectra using commercial software such as MestReNova (chemical shifts, integrals, coupling patterns, and correlations from standard 2D experiments) in a way that makes it easier for the user to understand and keep track of the information while verifying a proposed structure.

The program displays the data  in an interactive manner. The idea behind the program was to position it as an analysis tool between the raw data and the use of pencil and paper to analyse NMR data and commercial, complete black box solutions that can provide an answer with very little manual intervention. The software takes a high level interactive approach to display the information from the NMR experiments in order quickly check if a proposed structure is consistent with the NMR data.

The screenshot below shows the main user interface of the program. It consists of two panels. THe left-hand panel shows the proton and carbon 1-D spectra of the molecule. The right hand side panel shows an image of the proposed molecule provided by the user on top of which are placed the carbon atoms of the molecule. COSY data is shown as links between carbon atoms.  The carbon atoms/nodes can be moved over the background so that they align with the corresponding atoms in the proposed molecule. HMBC corelations are displayed on the molecule when the cursor is positioned over a carbon atom. The corresponding proton and carbon peaks are highlighted in the spectra panel.

Furthermore, when a peak is highlighted in the spectra panel, the peak is highlighted in red, if the peak has a corresponding carbon/proton peak it is  highlighted. HMBC correlations are shown by highlighting further peaks in different colours and showing the links in the molecule panel. Information on what the chemical shift might correspond to in terms of functional group is shown up in a pop-up.

![simpleNMR](Screenshot_simpleNMR.png)



## Installation

The program comes with a requirements file. Due to the use of the module RDKIT, I believe the best way to install the software is to use a CONDA installation. RDKIT can be installed using pip, but requires a python version of 3.7

## Running the program

The program can be run from the command line by typing 

```python simpleNMR.py```

Then an  example excel file can be opened  from the **File** dropdown menu by clicking the **open** item.

A problem directory can be opened directly for the command line by starting the program with the path of the example directory on the commandline

```python simpleNMR.py exampleProblems\ethyleneDiamine```

## Example problems

There are a number of example problems in the **exampleProblems** folder and they are listed below. Some are from real data and others are taken from extracting the data from examples in the book:

```Guide to NMR Spectral Interpretation A Problem-Based Approach to Determining the Structures of Small Organic Molecules. Antonio Randazzo, 2018, Loghia Publishing. ISBN: 978-88-95122-40-3```

 - ch9_013
 - ch9_016
 - ch9_021
 - ch9_025
 - ch9_041
 - ch9_048
 - ch9_053
 - ch9_069
 - ch9_092
 - ch9_102
 
 - ethyleneDiamine
 - glycidyl_methacrylate
 - Rotenone
 - nowick
 - Problem74
 - Problem79
 - Problem90
 - Problem90a
 


