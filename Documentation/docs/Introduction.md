# Welcome to SMART Documentation

SMART 3.0 (Small Molecule Accurate Recognition Technology) is the latest AI-based NMR analysis. 

SMART is a user-friendly, AI-based dereplication and analysis tool that uses 2D NMR data to rapidly associate newly isolated NPs with their known analogues. SMART has been designed to mimic the normal path of experiential learning in that additional 2D NMR spectral inputs can be used to enrich its database and improve its performance. In short, SMART aims to become an experienced associate to natural products researchers as well as other classes of organic chemists.

# How SMART3.0 is working?

![image](https://user-images.githubusercontent.com/51690359/118436403-a0421600-b695-11eb-85ed-98cca36bfa5f.png)

Figure 1. Overview of SMART 3.0. In the data construction phase, experimental 1H-13C HSQC data are converted to a constructed 1H-13C HSQC spectrum. In the structure annotation phase, the chemical fingerprints and molecular weights are predicted by SMART 3.0 and compared with structure databases to identify related NPs. In the class annotation phase, the compound class is predicted and provides support to the structure identification and annotation process.

# Overall framework

![image](https://user-images.githubusercontent.com/51690359/118436457-b2bc4f80-b695-11eb-9880-7e8ded31aeaf.png)

Figure 2. (A) The multi-task learning architecture of SMART 3.0. In the feature extraction step, the convolutional neural network extracts the features from HSQC spectra. Based on the extracted features, fully connected layers predict chemical fingerprints, molecular weights, and chemical classes. By using the predicted properties, structure annotation is performed. (B) The difference between HSQC and Edited HSQC spectra. On the Edited HSQC spectra, CH2 (methylene) correlations are shown separately from CH3/CH correlations. All of these correlations are binned and reconstructed on the image using 128 by 128 pixels.


# The performance of SMART 3.0 from the benchmark.

Testset from 3982 compounds is available on https://github.com/mwang87/SMART_NMR3/tree/master/Test_set

- Prediction performance

![image](https://user-images.githubusercontent.com/51690359/118436739-30805b00-b696-11eb-96bd-6eb583be8196.png)

Figure 3. Evaluation of the accuracy of SMART 3.0 to predict properties using a test set (n=3,982) not present in the training set. (A) Average (orange line) and median (blue line) of cosine scores between predicted and ground truth fingerprints for HSQC and Edited HSQC data input. (B) Linear regression between measured (x axis) and predicted molecular weights (y axis). (C) Confusion matrix of classification results using SMART 3.0 with HSQC data. (D) Confusion matrix of classification results using SMART 3.0 with Edited HSQC data.

- Structural identification/annotation performance

![image](https://user-images.githubusercontent.com/51690359/118436557-ded7d080-b695-11eb-9dfb-dc14823670be.png)

Figure 4. Performance evaluation of SMART 3.0 and other available tools with the same test set (n=3,982). Percentage of correctly identified structures (A) and annotated structures (B) found in the top k output of the different tools, for maximum rank k = 1, 2, …, 50. For the measurement of annotation rate, cosine score of 0.8 was set as the threshold.

- Performance in different solvent conditions

![image](https://user-images.githubusercontent.com/51690359/118436611-f44cfa80-b695-11eb-8b52-57c9ebeb888b.png)

Figure 5. Evaluation of SMART 3.0 analysis in different solvent conditions. (A) Venn diagram of the number of experimental NMR data obtained in chloroform-d and methanol-d4, respectively. (B) Experimental HSQC spectrum of neoline dissolved in chloroform-d¬ (blue) and methanol-d4 (red). (C) Identification (solid) and annotation (dashed) rates in total experimental data. (D) Identification (solid) and annotation (dashed) rates in compounds with NMR data recorded in both solvents.

# Frequently asked Questions

## Data Privacy

**Q:** Is the data I upload to analyze publicly visible?

**A:** Data you upload to analyze at SMART will remain private unless you explicitly make it public. The way to make your data public are to add individual annotated HSQC data to SMART database. See also the tab 'Contributing to SMART'.

## Analyzing Data

**Q:** Can I select more than one file at a time for each group in molecular networking?

**A:** Yes, you can drag&drop up to 8 CSV files in the analysis window. The analyses will run iteratively.

## NMR spectrometer Types (Bruker, JEOL, Varian...)

**Q:** Does SMART support data from different NMR vendors?

**A:** Yes, the only input you need for SMART analysis is the NMR table of 1H and 13C chemical shifts. These are vendor and processing software independent.

## Browser Support

We officially test on the latest Chrome browser. Other browsers, e.g. Firefox, Internet Explorer, Opera, Edge, are not officially supported but likely will not have any issues with SMART 3.0.


# Acknowledgements

We thank Advanced Chemical Design, Inc. for permission to utilize their Spectrus Processor 2017.2.1 software tool to predict HSQC spectra of various natural products. We further thank Dr. Kikuko Hayamizu for permission to utilize tabulated HSQC data for natural products from the CH-NMR-NP database. This work was supported by NIH grant GM107550 to G.W.C , P.C.D., and W.H.G. and by the Gordon and Betty Moore Foundation under grant GBMF7622 to G.W.C., P.C.D., and W.H.G. 

We sincerely thank the following people who have done great work to make the idea of web-based SMART become true. You guys are awesome. 
Names appear in alphabetic order.

| First       |  Middle  | Last                    |
|-------------|:--------:|-------------------------|
| Kelsey      |    L.    | Alexander               |
| Nuno        |          | Bandeira                |
| Wout        |          | Bittremieux             |
| Antoine     |    R.    | Blosse                  |
| Andrés      | Mauricio | Caraballo-Rodriguez     |
| Mitchell    |          | Christy                 |
| Garrison    |    W.    | Cottrell                |
| Alexander   |          | Dagman                  |
| Pieter      |    C.    | Dorrestein              |
| Brendan     |    M.    | Duggan                  |
| Joseph      |    M.    | Egan                    |
| Martha      |          | Gahl                    |
| Erik        |    C.    | Gerwick                 |
| Lena        |          | Gerwick                 |
| William     |    H.    | Gerwick                 |
| Michael     |    K.    | Gilson                  |
| Jenny       |          | Hamer                   |
| You Kyong   |          | Han                     |
| Yerlan      |          | Idelbayev               |
| Kyobin      |          | Kang                    |
| Hyunwoo     |          | Kim                     |
| Preston     |    B.    | Landon                  |
| Aaron       |          | Landon                  |
| Tiago       |    F.    | Leao                    |
| Ki Yong     |          | Lee                     |
| Bettina     |          | Lehman                  |
| Eugene      |    C.    | Lin                     |
| Roger       |    G.    | Linington               |
| Zheng       |          | Long                    |
| Huanru      |   Henry  | Mao                     |
| Jie         |          | Min                     |
| Anthony     |          | Mrse                    |
| Ben         |    C.    | Naman                   |
| Yashwanth   |          | Nannapaneni             |
| Louis-Félix |          | Nothias                 |
| Poornav     |    S.    | Purushothama            |
| Siddarth    |          | Ravichandran            |
| Raphael     |          | Reher                   |
| Nicholas    |    C.    | Roberts                 |
| Hyeji       |          | Shin                    |
| Yiwen       |          | Tao                     |
| Yoshinori   |          | Uekusa                  |
| Ezra        |    L.    | Van Everbroeck          |
| Vishal      |    T.    | Vasudevan               |
| Mingxun     |          | Wang                    |
| Yiran       |          | Xu                      |
| Chen        |          | Zhang                   |
| Jianping    |          | Zhao                    |


<a href="https://smart.ucsd.edu/classic" target="_blank"><img src="https://user-images.githubusercontent.com/20175888/70386594-ecd8dc00-194e-11ea-8378-ba1929e90ae4.png" alt="SMART_LOGO" title="USE SMART" align="right" width="250" height="100" border="10" /></a>
