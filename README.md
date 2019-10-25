----------------------------------------Dataset for the article:----------------------------------------------------------------
--------Entrainment of mammalian motile cilia in the brain with hydrodynamic forces.--------------------------------------------
--by Nicola Pellicciotta--


We analysed the dynamics of ependymal ciliated cells under oscillatory flow. The complete dataset is composed of more than 2 Tb, so here we uploded only a complete dataset for one cell.
Each directory is named  accordingly to the flow treatment, i.e external frequency applied and the magnitude of the flow. For example:

60X_V2_O3_P40.27Nov2018_17.50.20

60X: stand for the objective that we used

V2: is the voltage applied at the solenoidal end, we applied from 2 to 5 Volts. For a convertion into the flow velocity see supplementaries of the article.
 
O3: this is not useful parameter, it is always the voltage plus 1, in this case V=2, then O=3.

P40:  this is the period of the external oscillatory flow. the frequency is 1000/P , in this case f_ex = 1000/40 = 25 Hz

So we now that in this directory we report cilia dynamics (as a series of tiffs) for a cell under an oscillatory flow with frequency 25Hz and voltage 2 (see paper for the corresponding flow velocity or check the matlab variable calibration_matrix.mat).

---------------------------- code for analysis ---------------------------------------------------------------------------------------

We report the Matlab code that I used to extract from each cell dataset an arnold plot (as in Figure 1 in the paper). 
The workflow is divided in 2 steps:

step1_fft.m:
for each flow treatment (for example 60X_V2_O3_P40.27Nov2018_17.50.20)  we find the CBF of the cell buy using FFT analysis of the pixel intensity over time on a region above the cell. All the CBF and fft spectra for each flow treatment are saved in a variable fft_results_fft.mat. Here we upload the results for this particular cell dataset as fft_resuts_fft_nicola.mat , you are welcome to use it.

step2_arnold.m:
From the fft spectra the code evaluate if the cilia were entrianed by the external flow. THe procedure is explained in the methods of the paper. For each entrianment event it plot a black filled circle, otherwise an empty circle. Eventually it save a variable Cell.mat that contain the entrainment information. (the one that I used is Cell_nicola.mat). 

This procedure was repeated for all the cell analysed (58 cells) and then statistical analysis is performed. The code used for the figures in the article are on github.   


---------------------------------------------------------------------------------------------------------------------------------------

