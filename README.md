# On Mouse Mating work 

## Project goals

### Central questions 

* Determine the Influence of parental weight for offspring weight gain.
* Should mothers loose weight before becoming pregnant?

## Formulated goals

* Test the following hypotheses:
  * _"Pre-pregnant maternal weight loss leads to heavier offspring especially when father gets HFD."_
  * _"Leanest children result when both parents taking the same diet, either HFD x HFD or chow x chow."_
  * _"Offspring resulting from HFD x chow mating are heavier than HFD x HFD"._
* Rebuild plots from slides in `/Users/paul/Documents/HM_MouseMating/communication/190910 Labmeeting Verpaarungsstudie.pdf`

## 23-11-2022

* created directory structure using `find . -type d -exec mkdir -p /Users/paul/Documents/HM_LOBB_BSC/{} \;`
* saved project background materials to `/Users/paul/Documents/HM_MouseMating/communication/190910 Labmeeting Verpaarungsstudie.pdf`
* saved project raw data to `/Users/paul/Documents/HM_MouseMating/analysis/raw_data/190719 Rohdaten Verpaarungsprojekt.xlsx`
* created `/scripts/inspect_data.R`
* got a mouse shape with CC0 License from [here](https://www.svgrepo.com/svg/142184/mouse-mammal-animal-shape), stored at `/manuscript/figure_drafts/mouse-mammal-animal-shape-svgrepo-com.svg`, also go many other shapes
* finished draft figure at `/manuscript/figure_drafts/fig_foo.svg`
* commit `a3c8514a6f62bcc73f20c4484e7bd154da5224b9`

## 12-12-2022

* finish draft figure in biorender
* see `/Users/paul/Documents/HM_MouseMating/manuscript/figure_drafts_biorender/221212_mouse_mating.pdf`
* see `/Users/paul/Documents/HM_MouseMating/manuscript/figure_drafts_biorender/221212_mouse_mating.png`
* got cleaned data from Sebasian
* see `/Users/paul/Documents/HM_MouseMating/analysis/raw_data/221212 Verpaarungsprojekt.xlsx`

## 14-12-2022

* received `/Users/paul/Documents/HM_MouseMating/analysis/raw_data/221212 Verpaarungsprojekt.xlsx`
* changed all dates to iso dates for import into R
* start import and data inspection in `/Users/paul/Documents/HM_MouseMating/analysis/scripts/inspect_data.R`
* commit `aa0f3bce64991cf869676d03be66eddf1693415d`

## 19-12-2022

* reshaping raw data in `221214_mice_pairing_resahped.xlsx`
* read data into `inspect_data.R`
* commit `dd2d532533ed7656314a1c010c76edbb5c8f724c`

## 20-12-2022

* in `inspect_data.R`
  * starting data inspection and plotting
  * commit prior to changing plotting code
  * commit `2f56668823fb3f6494c916343eca3d5a77cadb11`
  * restructured code and code comments to indicate how to proceed
  * minor edits to ANOVA tests 
  * commmit `bde3eb9b08451332e0b539fadc964cb4dd0af11a`

## 21-12-2022

* in `inspect_data.R`
  *  starting to explore polynomials
* commit `26f841e6b4bb024466ea8171cb9953ab97a85340`
* plotted multivariate linear model with polynomials in ggplot
* commit `c725d55d68521c737bf668f35153099fc1a145d3`

## 22-12-2022

* in `inspect_data.R`
  * extending model formula to inlcude diets
  * commit prior to chnaging plotting code to use `ggplot()` `c725d55d68521c737bf668f35153099fc1a145d3`
  * commit prior to removing `ggeffects()` plots and using only `ggplot()` `6b19c55a12fea4c3651811373fbfc1e66f5f9d93`
  * test model is finished, without random effects
  * commit `cb00b9c4b993dfc7396c8425dbe4774ca8da4644`
  * finished preliminary analysis
  * commit `b938ce5972986f9fedfa66bd842158f57775611e`
  * updated what to do next, commit `4d7dac255f96aab9c04e0c0671db21509fd7cc85`
  
## 02-01-2023

* in `inspect_data.R`
  * cleaned script
  * revised polynomial regression
  * implemented mixed effect modelling
  * implemented ANOVA-equivalent using OLS regression
* commit `34ea85a9cc6fd748bc8c1eee6d0b389771666570`

## 04-01-2023

* starting implemnting reporting of `inspect_data.R`
  * spun using `knitr::spin("/Users/paul/Documents/HM_MouseMating/analysis/scripts/inspect_data.R")`
  * after corceting errors using `rmarkdown::render("/Users/paul/Documents/HM_MouseMating/analysis/scripts/inspect_data.R", output_format = "html_notebook", output_dir = here("script_spun")`
* in `inspect_data.R` 
  * added `bartlett.test()` to test equal variance across groups
* commit `bcc98482c45ea662c6df800ecbe6be65d6a19594`

## 10-01-2023

 * read Appendix and read Chapter 5 of **Zuur AF, Ieno EN, Walker N, et al. 2009. Mixed effects models and extensions in ecology with R. New York, NY: Springer.**
 * found helpful hint on modeleling, in addition to Chapter 5 on [webpage (saved in Zotero)/](https://stats.stackexchange.com/questions/71087/analysis-of-a-time-series-with-a-fixed-and-random-factor-in-r)
 * updated todo queue
 * in `inspect_data.R`
   * added analysis remarks to all section
   * restructured data inspection of full data
   * addded code outline as per Appendix for data inspection and regression setup (gls)from Zuur et al. 
   * commit `6d5fb92e0bb80b1d4c34eb742ea2cc72e738cdde`
   * added Chapter 5 steps to GLS modelling section (last part of script.
   * started modelling using `nlme::lme()` - need to abort for AH
   * see comments in code (line 593) and also wehre to continue reading in book
   * commit `c3eddd7baa0cba9f45c99d5b291f07636244b0d1`   

## 23-01-2023

* [x] in `inspect_data.R` - verified modelling as per ISME with random effect modelling - see communications folder 09.01.2023 - and see the verious tutorial, and the `lme4` book.
  * [x] improve implement random effect modelling
  * [x] inspect model residuals
  * [ ] ~~select splines / polynomials based on AIC~~ **(splines won't work in `lme4`)**
  * [ ] ~~implemnet backwards selection using GLMER~~ **(used manual forward selection as in Bates (2015) and Zuur et al. (2005))**  
* commit `b8519802e991060688c2a88a6a88c5b3fbe2f9a8`

## 24-01-2023

* [x] in `inspect_data.R`
  * [x] confirm modelling with HI-MAG stats support group
  * [x] names of stats poeple and add to paper **Elmar Spiegel** and **Roman Le Glut**, **Maciej Rosolowski**
  * [x] updated next steps
  * commit 

## 25-01-2023

* [x] in `inspect_data.R`
  * [x] use conditional AICs to inpect models - lines `547` and `580`
  * [x] plot the best models using ggplot
  * [x] basic modelling script done - now needs to be cleaned and expanded to answer the question posed in the introduction

* commit `d5871657c9665c01aa421bb58015cb6dbb77635c`
* created `main_analysis.R` for refining lme4 solution and starts answering questions
* commit `d5871657c9665c01aa421bb58015cb6dbb77635c`
* [x] in `main_analysis.R` 
  * added more model fomulae with dietary interactions
  * re-sorted modelling table - about to delete models
  * commit `8fe8a92c37807244ed2729c3a437d9a80b11bf37`
* updated plotting 
* commit `39b035178f0277d89b0d3b568b1fdd9e0e2489aa`

## 26-01-2023

* uodated ToDo - added new ideas

## 30-01-2023

* updated ToDo - added new ideas

## 08-02-2023

* working in `main_analysis.R` 
  * restructured script 
  * commit `cd6a930f3a39199744609742bd7470c5a49240e7`

  

  
## Todo queue

* [ ] in `main_analysis.R `
  * [ ] find mothers with weight loss
  * [ ] add weight loss as factor to modelling data
  * [ ] sort models sytematically
  * [ ] plot fixed and random effects as with DINCH mice
  * [ ] consider usding splines or GAM
  * [ ] build and re-evaluate models with reagrds to maternal weight loss
  * [ ] build and re-evaluate models with reagrds to sex-sepecific weight gain - sex needs to be included or evaluated
  * [ ] to understand effect of polynomials and sex - plot random-conditioned fixed effects using the `effects` package
  * adress high VIFs - by recoding factor variable 
  * adress high VIFS - by ignoring factor variable (`https://stackoverflow.com/questions/33397689/multi-collinearity-for-categorical-variables`)
  * adress high VIFS - **by centering time variable** (?) 
  * [ ] answer last question
  * [ ] sort script for rendering
  * [ ] possibly rewrite script as done for DINCH mis, also using `GGpairs`

* in manuscript
  * add methods  
  * add IMISE consulting personell to acknowledgements: 
    * **Maciej Rosolowski** <maciej.rosolowski@imise.uni-leipzig.de>
    * Maryam Yahiaoui-Doktor <maryam.Yahiaoui-Doktor@imise.uni-leipzig.de>
  * add HM consulting personell to paper
    * **Elmar Spiegel**  
    * Roman Le Glut

## Ideas not yet pursued

* [ ] implement modelling steps as in Appendix **Zuur AF, Ieno EN, Walker N, et al. 2009. Mixed effects models and extensions in ecology*** 
* [ ] in `main_analysis.R ` - possibly play around further 
  * [ ] possibly add classification tree
  * [ ] list variance explained by factors - as per [web page](https://github.com/timnewbold/MResEcologicalModelling/blob/4adab861b41934e460f1a3f9bad17a4398acb068/1StatisticalModels/WorkshopExercises.md) in Zotero.
  * [ ] rewrite everything succeinctly as in DINCH mice?
 
## Considerations

* AH: _"Hier die Übersicht der Maiting study. Im Prinzip wäre es schön, wenn du die Plots aus den Folien in R erstellst (etwas hübscher machen) und eine geeignete Statistik drüber laufen lässt für die time series ANOVA."_

## Further reading

### Serial correlations and model verification
* **Zuur AF, Ieno EN, Walker N, et al. 2009. Mixed effects models and extensions in ecology with R. New York, NY: Springer.**

### Not relevant anymore 

* on multivariate polynomials
  * Kahle D. 2013. mpoly: Multivariate Polynomials in R. The R Journal 5: 162.

* On contrast coding: `https://marissabarlaz.github.io/portfolio/contrastcoding/#what-is-a-contrast`

* for modelling consider _"A guide to creating design matrices for gene expression experiments. https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html. Viewed 25 Nov 2022."_, section 6.1
  






