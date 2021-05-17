## User's Guide for SMART 3.0 Analysis

Welcome to use the SMART 3.0 to test your compound(s) @ [SMART 3.0](https://smart3.ucsd.edu/).

- The current version of SMART 3.0 as of 05/16/2021 consists of 2D NMR spectra from 144,254 natural products. 
- One SMART 3.0 analysis should take < 20 seconds.
- If your results are dissatisfying please try to process your data again manually (go to **How to process a raw HSQC spectrum to a NMR table** and then delete noise and duplicate annotations, add peaks missed by auto-peak picking etc.) 

## How to process a raw HSQC spectrum to a NMR table with MestreNova (version 12 and newer)

1. Open your raw HSQC spectrum in MestreNova (preferences: modern view)
    - Drag&Drop your HSQC file (for Bruker data you find your spectrum under: pdata/1/2rr)
    - Depending on purity and concentration of your sample and aquisition time your
      spectrum looks more or less clean and may need additional processing (see 2.)
    ![image](https://user-images.githubusercontent.com/57916837/70353294-f153a680-1821-11ea-96fb-85f3fea4f2b4.png)
2. For processing your HSQC spectrum click on 'Processing' tab
    - click on 'Auto Phase Correction' (optional: correct manually)
    - click on 'Auto Baseline Correction'
    - click on 'More Processing' --> click on 'Reduce t1 noise'
    You should see a clean spectrum now.
    ![image](https://user-images.githubusercontent.com/57916837/70353422-2b24ad00-1822-11ea-90a9-424dd2619d1a.png)
3.  Annotate HSQC spectrum with chemical shifts (1H,13C)
    - click on 'Analysis' tab
    - click on Auto Peak Picking (Important: Check by manually adding missed peaks and removing duplicated, nonsense and solvent peak         annotations).
    Now each peak should be annotated with two numbers separated by comma (1H, 13C chemical shifts)
    ![HSQC_annotation](https://user-images.githubusercontent.com/57916837/70353559-7939b080-1822-11ea-84ea-87cc07946574.png)
    
4. Generate NMR table from annotated HSQC spectrum
    - click on 'Analysis' tab
    - click on 'NMR Peaks Table'
    - right click on table, setup report, setup table
    - customize table by unchecking every value but for f2 (change visible name to 1H), f1 (change visible name to 13C), and Intensity (optional, but essential for Edited-HSQC)
    - copy all (ctrl+A)
    - click on 'copy peaks' and choose 'copy table'

5. Run SMART Analysis directly
    - copy and paste the table (ctrl+V) directly to the peak list section of https://smart3.ucsd.edu
    - Important: Apply one backslash to remove the additional space character that is imported with the NMR table'
    - Select experiment type (Normal HSQC or Edited HSQC)
    - (Optional) If you know the molecular weight of your compounds, please enter it.
    - Once the data and experimental condition are submitted, analysis is automatically processed.                                                                                                                            

![image](https://user-images.githubusercontent.com/51690359/118419405-bee1e600-b670-11eb-8d24-fd4170fd08c6.png)

    
**You will get the results from SMART 3.0 Analysis! :)**

**Please feel free to play around with the processing parameters such as including/excluding noise signals or signals from other minor compounds in case of mixtures or explore the differences of SMART results when referencing your spectra compared to tables without referencing. Overall SMART is designed to be very robust towards any of these changes as its training is not only based on the absolute position of the peaks, but the relative position of each peak towards every other peak (see also References 1 + 2).**

## Input Data Formatting using MestreNova

Please prepare your NMR peak lists of each compound using Excel or preferably notepad/wordpad and save the NMR table of each compound as a comma-separated value (.csv) or tab-separated value (.tsv) file. Furthermore, SMART does now support peaklists from TopSpin. Please always place 1H data in the first column and their corresponding 13C data in the second column. The first row will be left for strings “1H” and “13C” as table head. Please strictly keep 1H shifts in 2 decimals and 13C shifts in a single decimal only. You can name the .csv file with any name that your operating system accepts.

|     13C    |      1H     |  Intensity  |
|:----------:|:-----------:|:-----------:|
|   128.3    |     7.61    |   6861.17   |
|   106.8    |     6.16    |   7830.41   |
|   102.2    |     6.82    |   4734.83   |
|    45.8    |     3.69    | -14123.20   |
|    28.5    |     2.22    |  14767.24   |
|    20.5    |     2.34    |  21184.69   |

In the NMR table files, wherever there are diastereotopic protons on a methylene carbon (i.e., CH2 with two distinct proton shifts), please add a separate entry for both the carbon and proton:

|     13C    |      1H     |
|:----------:|:-----------:| 
| 43.2       | 3.12        |
| 43.2       | 4.30        | 

### Supported Formats for SMART Analysis (MestreNova)
SMART supports CSV and TSV formats for analysis. You can check your .csv or .tsv file by opening it with text editors such as wordpad or notepad. Your table should appear like this:

    1H,13C      
    1.09,14.3
    2.21,22.2
    3.41,56.9
    7.21,128.6
    7.29,123.4
    
    or for .tsv:
    
    1H	13C
    1.09	14.3
    2.21	22.2
    3.41	56.9
    7.21	128.6
    7.29	123.4


### Submit CSV or TSV file, or TopSpin peaklist for SMART Analysis

After the .csv or .tsv file or TopSpin peaklist is ready, simply upload it to the first window by drag&drop (analysis automatically starts) or copy and paste the table and click analyze. You can submit multiple jobs at one time. **Important: You have to allow pop-ups for https://smart.ucsd.edu/classic. If your analysis is still running after 20 seconds, pop-up blockers are the most likely reason why your analysis failed.** The result will show up in the next webpage as images with chemical structures, compound names, similarity scores, and molecular weights. The Top 100 hits are ranked by similarity scores. You can download the result by clicking the “download results” icon and open it with Excel.

## Troubleshooting

### Overview
This section addresses some common issues with the analysis workflows at SMART. If you run into any of these common issues, hopefully this page will give you actionable steps to address them. If this page cannot help you, please refer to the [forum](https://groups.google.com/forum/#!forum/smartnmr) to help answer your questions.

### SMART Analysis

**Failed Job**

1. If your analysis is still running after 20 seconds, pop-up blockers are the most likely reason why your analysis failed. You have to allow pop-ups for https://smart.ucsd.edu/classic.
2. The file format input is incorrect. Please make sure it is a supported file format for SMART (see Input Data Formatting)
When you enter the peak list as tab-separated or comma-separated table make sure that the first row is:
1H,13C or 
1H    13C, respectively. Do NOT include any additional spaces or " signs.
**Especially, if you copy&paste your table directly from MestreNova, apply one backslash to remove the additional space character that is imported with the NMR table**

**Downloaded Results file cannot be opened**

Yes you can :) You can open the downloaded results with any text editor (Excel, Notepad, Wordpad). Just make a right-click on the file, click on 'open with' and then choose one of the former mentioned programs.

**Results Incorrect**

The SMART is still at its early stage of development. It will become more and more accurate the more spectra you contribute to the Moliverse (see Contribute to SMART).
- If your results show strange suggestions and all of the have a cosine score of 1.0 you probably have switched the columns. Please be sure that your first column has proton shifts (1H) and your second column has carbon shifts (13C), where each 1H-13C pair is separated by a comma.
- If your CSV file won't show exactly two decimals for 1H, or exactly a single digit for 13C, you can mannually change the last digit of the chemical shift from 0 to 1. E.g., change 7.30 to 7.31. Slight variation of the last digit won't affect your test result.

## Visualize your result(s) in the 3D cluster space (Moliverse)

Please click the "Embed in the Moliverse" bar in the upper right side of the result page. Wait a while for the Embedding Projector to launch. To find your query molecule, please type "query" in the search bar in the upper right side of the Embedding Projector Webpage. Naturally, your query compound as well as its 100 closest analogues will light up in unique colors (decreasing cosine similarity will be displayed by purple-red-orange-yellow).

- To better observe the clustering effect, please click UMAP or T-SNE tabs in the lower left window. The clusters will organize themselves in real-time, clustering more similar compounds in close proximity and more dissimilar compounds apart from each other.
- Further you can access and enjoy the whole 'Moliverse' under https://smart.ucsd.edu/classic (click on 'Explore Moliverse')

## Visualize the 3D cluster space (Moliverse) with VR devices

Please click [Video_1](https://youtu.be/Pgfw-90t1t0), [Video_2](https://youtu.be/2qM_kFf-SRo) and [Video_3](https://www.youtube.com/watch?v=v_c-Z1a5-KU) to watch the VR movie. The VR software will be delivered soon!

[![VIDEO_1](http://img.youtube.com/vi/Pgfw-90t1t0/0.jpg)](http://www.youtube.com/watch?v=Pgfw-90t1t0 "Small-molecule Interactive Mapping Applications")

[![VIDEO_2](http://img.youtube.com/vi/2qM_kFf-SRo/0.jpg)](http://www.youtube.com/watch?v=2qM_kFf-SRo "Into the Moliverse")

[![VIDEO_3](http://img.youtube.com/vi/v_c-Z1a5-KU/0.jpg)](http://www.youtube.com/watch?v=v_c-Z1a5-KU "The Moliverse, A Digital Frontier to Reshape Natural Products Research")

## Contact

If you have any questions about how to use SMART please first ask the community at the SMART forum located here: [forum](https://groups.google.com/forum/#!forum/smartnmr). We will respond you as soon as we can. If you woud like to access our to-be-online SMART 3.0, please email Hyunwoo Kim (hwkim at-sign ucsd dot edu) for off-line trials.

## References

To reference the system please cite the papers listed below.

1. Zhang C\*, Idelbayev Y\*, Roberts N, Tao Y, Nannapaneni Y, Duggan BM, Min J, Lin EC, Gerwick EC, Cottrell GW, Gerwick WH. Small Molecule Accurate Recognition Technology (SMART) to Enhance Natural Products Research. *Scientific Reports.* 2017, **7(1)**, 14243.  [DOI: 10.1038/s41598-017-13923-x](https://doi.org/10.1038/s41598-017-13923-x) *\*These authors contributed equally to this work.*

2. Reher R\*, Kim H\*, Zhang C\*, Mao HH, Wang M, Nothias LF, Caraballo-Rodriguez AM, Glukhov E, Teke B, Leao T, Alexander KL, Duggan BM, Van Everbroeck EL, Dorrestein PC, Cottrell GW, Gerwick WH. A Convolutional Neural Network-based approach for the Rapid Characterization of Molecularly Diverse Natural Products. *Journal of the American Chemical Society.* 2020, **142(9)**, 4114-4120.  [DOI: 10.1021/jacs.9b13786](https://doi.org/10.1021/jacs.9b13786) *\*These authors contributed equally to this work.*  

3. Coming soon...

<a href="https://smart.ucsd.edu/classic" target="_blank"><img src="https://user-images.githubusercontent.com/20175888/70386594-ecd8dc00-194e-11ea-8378-ba1929e90ae4.png" alt="SMART_LOGO" title="USE SMART" align="right" width="250" height="100" border="10" /></a>
