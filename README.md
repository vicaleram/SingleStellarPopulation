# SingleStellarPopulation

This repository contains the notebooks to analyse the our population of magnetic massive stars. The evolutionary models used for each star in the population are computed using Modules for Experiments in Stellar Astrophysics [(MESA)](https://docs.mesastar.org/en/release-r22.05.1/index.html), and were developed by [Z. Keszthelyi, et al. (2022)](https://academic.oup.com/mnras/article/517/2/2028/6701644). Recreating a population of stars would require too much computing time with MESA, therefore we used the poppulation synthesis code MOSS [(Götberg, Y. et al.)](https://ui.adsabs.harvard.edu/abs/2019A%26A...629A.134G/abstract) to handle both the creation of stars and the evolution of each of them in our simulated star cluster. 

The goal of our work is to explore the possible explanations for the observable 10% incidence of magnetic O, A, and B type stars [(Grunhut et al. 2017)](https://academic.oup.com/mnras/article/465/2/2432/2417464?login=false). Our experiments consist of exploring different initial magnetic field functions (IBF), and determine how they affect the observational incidence rates of our popoulation at different cluster ages. Lastly, we include the effects of current observational capabilities to detect the magnetic fields of massive stars through spectropolarimetry and determine what stellar charactertistics have greater influence in the detection of their surface magnetic field.

This repository presents the Jupyter notebooks used in our analysis and any other relevant code. In the examples folder you can find a series of notebooks that **attempt** to give an ordered approach to our analysis. 

Feel free to contact us if you have any comments, questions or feedback!

![alt-text](https://github.com/vicaleram/SingleStellarPopulation/blob/main//hr_strbrstBf25_01.gif)
