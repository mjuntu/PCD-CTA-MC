# PCD-CTA-MC
Monte Carlo codes for assessing the cost-effectiveness of CT Angiography using a photon counting detector CT.

To run the codes, add the folder to MATLAB's search path (add to path --> folder and subfolders). The results are highly dependent on the used reclassification rates of CAD-RADS grades. For this work, the data from 10.1016/j.jcct.2024.10.011 is used (see supplementary material of this publication). To this end, the user may introduce their own reclassification data by modifying the table "Reclassification_data.xlsx". If one wants to include their own reclassification data, simply include your EID-CT, PCD-CT, and ICA CAD-RADS grades in columns E-G
