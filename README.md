# Dynamic prediction for multivariate longitudinal and time-to-event data using super learning and multivariate functional principal component analysis-based methods

This is a public repository containing the different documents related to Arnau Garcia's master's thesis for the [Master's degree in Statistics and Operations Research (MESIO) UPC]{https://mesioupcub.masters.upc.edu/en?set_language=en}. In this repository you will find both code files and pdf documents with the reports that were made in the course of the work.

## Rpository structure

* **Reports**: Contains the main reports produced during the master's thesis:

  * `SL_1stReport.pdf`, `SL_2ndReport.pdf`, `SL_3rdReport.pdf`
  * `multi_longitudinal_data_MMvsMFPCA.pdf`
  
    Presentation slides are included in:
    
  * `slides_SLinJM_meeting1.pdf`
  * `slides_SLinJM_meeting2.pdf`

* **PBC\_analysis folder**: Provides the code and materials used in the case studies based on the PBC dataset, corresponding to the analyses presented in Sections 3.3 and 4.3.4 of the master's thesis.

* **simul\_SL folder**: Contains R scripts for the simulation study evaluating super learner performance in joint modeling. The design, methodology, and results of this study are described in Section 3.2 of the work.

* **simul\_MFPCA\_missing folder**: Includes R code for the simulation study investigating the robustness of MFPCA-based predictions under varying missing data scenarios. This study is detailed in Section 4.2.1 of the master's thesis.

* **simul\_MFPCA folder**: Contains both R and Python scripts for the simulation study comparing MFPCA-based methods, the super learner approach, and multivariate joint models. The corresponding analysis is presented in Section 4.3.3 of the work.


