# simpleNMR

simpleNMR is a python/QT5/matplotlib program. It runs on windows and macos. The aim of the software is to present the information derived from processing NMR spectra using commercial software such as MestReNova (chemical shifts, integrals, coupling patterns, and correlations from standard 2D experiments) in a way that makes it easier for the user to understand and keep track of the information while verifying a proposed structure.

The program displays the data  in an interactive manner. The idea behind the program was to position it as an analysis tool between the raw data and the use of pencil and paper to analyse NMR data and commercial, complete black box solutions that can provide an answer with very little manual intervention. The software takes a high level interactive approach to display the information from the NMR experiments in order quickly check if a proposed structure is consistent with the NMR data.

The screenshot below shows the main user interface of the program. It consists of two panels: THe right-hand panel shows the proton and carbon 1-D spectra of the molecule. The left-hand side panel shows an image of the proposed molecule provided by the user on top of which are placed the carbon atoms of the molecule. COSY data is shown as red links between carbon atoms.  The carbon atoms/nodes can be moved over the background molecule so that they align with the corresponding atoms in the proposed molecule. HMBC corelations are displayed on the molecule when the cursor is positioned over a carbon atom. The corresponding proton and carbon peaks are highlighted in the spectra panel.

Furthermore, when a peak is highlighted in the spectra panel, the peak is highlighted in red, if the peak has a corresponding carbon/proton peak it is  highlighted. HMBC correlations are shown by highlighting further peaks in different colours and showing the links in the molecule panel. Information on what the chemical shift might correspond to in terms of functional group is shown up in a pop-up.

![simpleNMR](Screenshot_simpleNMR.png)

The initial positions of the carbon atom nodes over the background molecule image can be set randomly, from a previously saved session, or by using the JAVA HOSE code that can be found on [nmrshiftdb](https://nmrshiftdb.nmr.uni-koeln.de/) website which has been incorporated into the program.

## Installation

The program comes with a requirements file. Due to the use of the module RDKIT, I believe the best way to install the software is to use a CONDA installation. RDKIT can be installed using pip, but requires a python version of 3.7, but I have not been able to get it to work so for now I would stick to a CONDA installation.

### JAVA

In the repo, there are JDK binaries for windows, macos and linux in the folder **jre**. The binaries were downloaded from amazon [amazon corretto](https://docs.aws.amazon.com/corretto/latest/corretto-8-ug/downloads-list.html). The program checks to see if java is installed by the user and if not, defaults to these binaries to run the java code.

Sometimes the java code has been compiled for a later version and will not run. An error message will appear in the background terminal.  If this happens run the appropriate bat files to recompile the java code and restart the program simpleNMR.py

 - ```CompileJavaMacLinux.bat```
 - ```CompileJavaWindows10.bat```
 
 ## PyInstaller
 
 An executable of the program can be be created using pyinstaller for windows or macos.  A .spec file for creating a single executable has been created. The command to create an executable is 
 
 ```pyinstaller simpleNMR_with_includes-F.spec```
 
 Before running the command it is best practise to copy the following files over to a new directory and run the aboove command from there.
 
 - ```about_dialog.py``` 
 - ```excelheaders.py```
 - ```moleculePlot.py```
 - ```nmrProblem.py```
 - ```nmrmol.py```
 - ```nx_pylab.py```
 - ```qt5_tabs_001.py```
 - ```simpleNMR.py```
 - ```spectraPlot.py```
 - ```xy3_dialog.py```
 - ```jre```
 - ```csTables```
 - ```cdk-2.7.1.jar```
 - ```predictorc```
 - ```NewTest.class```
 
 In the ```jre``` folder keep only the java runtime environment that matches the operating system the pyinstaller command is run on.
 
 Make sure that pyinstaller is ran from the correct python environment that simpleNMR.py works on.
  
## Running the program

The program can be run from the command line by typing 

```python simpleNMR.py```

Then an  example excel file can be opened  from the **File** dropdown menu by clicking the **open** item.

A problem directory can be opened directly from the command line by starting the program with the path of the example directory on the commandline

```python simpleNMR.py exampleProblems\ethyleneDiamine```

## Requirements

 - matplotlib==3.5.1
 - networkx==2.4
 - nmrglue==0.8
 - numpy==1.21.5
 - pandas==1.3.4
 - Pillow==9.2.0
 - PyQt5==5.15.6
 - PyYAML==6.0
 - rdkit==2020.09.10
 - scipy==1.7.3

## Example problems

There are a number of example problems in the **exampleProblems** folder and they are listed below. Some are from real data and others are taken from extracting the data from examples in the book:

```Guide to NMR Spectral Interpretation A Problem-Based Approach to Determining the Structures of Small Organic Molecules. Antonio Randazzo, 2018, Loghia Publishing. ISBN: 978-88-95122-40-3```

 
