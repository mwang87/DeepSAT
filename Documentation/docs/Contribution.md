
## Instructions for submitting NMR data to the Moliverse (the training set for SMART)
Please read all of the instructions below before you fill in the data collection files.

### Definitions
Folder = a folder on your computer for organizing files; File = an excel file in .xlsx or .csv format; 

Sheet = a sheet within an excel file; Training database = the database used to train the AI; 

Test dataset = the inquiry dataset awaiting the AI to generate structural hypotheses.

### Submitting structures to the Moliverse library
Thank you for contributing to the SMART training database (AKA the Moliverse Library). For this submission, we ask you to fill out some information in the metadata file (only one metadata file for all compounds being submitted), and a second file with the tabulated NMR data as described below (the NMR data for each compound in its own file). 

Connection between the two files will be by a 6 digit code number.  The starting letter/number combination for your submitted compounds will be assigned by the SMART operators, and then you should increment up for each ensuing compound by a value of 1 (e.g. AAA001 is the assigned code, the next compound in your series is AAA002, and so forth). Although compound replicates need to be avoided, please don’t worry as the SMART can cope with some replicates. You can continue adding more entries on the pre-existing metadata file with the already given code pattern. 

In the Metadata file and the NMR Tables, only the items with an orange colored background are required; the green fields are desirable and will help us with developing more accurate and informative future versions of SMART. If for any compound, a specific data entry is not available, please leave it as blank.

As noted above, in the Metadata file, only the SMILES entry is mandatory. Only when the isomeric SMILES is not available for a compound, canonical SMILES can be entered instead. The SMILES can be found on the PubChem website for majority of published compounds. If you cannot find the SMILES on PubChem, the 3D structure in ChemDraw can be converted to SMILES by the program.

### Convert Structure to SMILES
To convert a ChemDraw structure to SMILES, you must:
1) Select the structure using the selection tool; 
2) From the Edit menu, point to Copy As, and then choose SMILES;
3) Paste the string in the target cell of the excel metadata file sheet.

### Minimum information for NMR tables
For compiling the NMR tables, please follow the examples in "NMR Tables" file in this same folder as the format for your data. NOTE: proton shifts should be rounded to exactly 2 decimal digits, and 13C shifts should be rounded to exactly a single decimal digit.
For the NMR tables, the atom numbering is not needed, but just be certain that the proton shift matches the carbon shift in every row of the table.

In the NMR table files, wherever there are diastereotopic protons on a methylene carbon (i.e., CH2 with two distinct proton shifts), please add a separate entry for both the carbon and proton:

| Column A | Column B | Column C               |   
|----------|----------|------------------------|
| 13C      | 1H       | Carbon type (optional) |
| 43.2     | 3.12     | CH2                    |
| 43.2     | 3.40     | CH2                    |

In the NMR table files, please put 13C shifts in Column A of the sheet and proton shifts in Column B of the sheet, as shown above.
If you have any questions, please contact Chen Zhang (chz023 at-sign ucsd dot edu) or Raphael Reher (rreher at-sign ucsd dot edu)

Please download data submission package [here](https://tinyurl.com/vee67qk). After you complete, please email your package back to chz023 at-sign ucsd dot edu.

## References, Q&A

To reference the system please cite the papers listed below. To ask questions, please post on the user [forum](https://groups.google.com/forum/#!forum/smartnmr). We will respond you as soon as we can.

1. Zhang C\*, Idelbayev Y\*, Roberts N, Tao Y, Nannapaneni Y, Duggan BM, Min J, Lin EC, Gerwick EC, Cottrell GW, Gerwick WH. Small Molecule Accurate Recognition Technology (SMART) to Enhance Natural Products Research. *Scientific Reports.* 2017, **7(1)**, 14243. [DOI: 10.1038/s41598-017-13923-x](https://doi.org/10.1038/s41598-017-13923-x) *\*These authors contributed equally to this work.*

2. Reher R\*, Kim H\*, Zhang C\*, Mao HH, Wang M, Nothias LF, Caraballo-Rodriguez AM, Glukhov E, Teke B, Leao T, Alexander KL, Duggan BM, Van Everbroeck EL, Dorrestein PC, Cottrell GW, Gerwick WH. A Convolutional Neural Network-based approach for the Rapid Characterization of Molecularly Diverse Natural Products. *Journal of the American Chemical Society.* 2020, **142(9)**, 4114-4120.  [DOI: 10.1021/jacs.9b13786](https://doi.org/10.1021/jacs.9b13786) *\*These authors contributed equally to this work.*   

3. Coming soon...

End of Instructions.

<a href="https://smart.ucsd.edu/classic" target="_blank"><img src="https://user-images.githubusercontent.com/20175888/70386594-ecd8dc00-194e-11ea-8378-ba1929e90ae4.png" alt="SMART_LOGO" title="USE SMART" align="right" width="250" height="100" border="10" /></a>
