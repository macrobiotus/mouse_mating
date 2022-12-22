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

## 21-12-2022

* in `inspect_data.R`
  * extending model formula to inlcude diets
  * commit prior to chnaging plotting code to use `ggplot()` `c725d55d68521c737bf668f35153099fc1a145d3`
  * commit prior to removing `ggeffects()` plots and using only `ggplot()` `6b19c55a12fea4c3651811373fbfc1e66f5f9d93`
  * test model is finished, without random effects
  * commit `cb00b9c4b993dfc7396c8425dbe4774ca8da4644`
  * finished preliminary analysis
  * commit `b938ce5972986f9fedfa66bd842158f57775611e`


## Todo queue

* [ ] read on time series ANOVA
* [ ] do time series ANOVA in R using longtidinal data technique as described [here](https://stats.stackexchange.com/questions/197607/how-to-test-difference-between-times-series-does-time-series-anova-exist). Uses package `nlme::glm()`
* [ ] do time series ANOVA in R - possibly using `ggbetweenstats()` of [ggstatsplot](https://indrajeetpatil.github.io/ggstatsplot/) 


## Considerations

* AH: _"Hier die Übersicht der Maiting study. Im Prinzip wäre es schön, wenn du die Plots aus den Folien in R erstellst (etwas hübscher machen) und eine geeignete Statistik drüber laufen lässt für die time series ANOVA."_

* for modelling consider _"A guide to creating design matrices for gene expression experiments. https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html. Viewed 25 Nov 2022."_, section 6.1

## Further reading

* on multivariate polynomials
  * Kahle D. 2013. mpoly: Multivariate Polynomials in R. The R Journal 5: 162.

* On contrast coding: `https://marissabarlaz.github.io/portfolio/contrastcoding/#what-is-a-contrast`
  






