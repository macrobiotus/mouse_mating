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
  * starting to work on freshly defined Objective 2
  * finished Objective 2, as per `gam06` in section
  * updated `README.md` below
  * commit `84d6d5f422af4bd33b7a7babcda8bb81032eb750`

## 09-02-2023

* working in `main_analysis.R`
  * improved GAM and found better one - with interaction
  * after reading lower section this [website](https://marissabarlaz.github.io/portfolio/gam/#visualizing-gam-results).
  * commit `1a4172ed230ff5007a050afe746d5bfc9ce7f55f`
  * finished Objective 3, as per `gam10` in section
  * commit `137527f7afd1ed32581f8d87e6f38425b866e5d9`
  * starting on Objective 1
  * alamost deone define "mothers who have lost weight" - positing that thes have the curviest weight gain curve
  * commit `5064bf141548c4981a17bbd7596b98d8ac007f65`

## 10-02-2023

* working in `main_analysis.R`
  * started modelling for for Objective 1
     * done for mother's weight loss only
     * done for mother's weight loss and fathers diet
     * see `gam14`
     * reporting results in chat
     * emailling about results
 * all initial 3 Objecktives analysed.
 * commit `b9b23a3c341912b0f2ea478432ce004b6bfc1950`

## 13-02-2023

* after meeting with Nora, we agreed that writing can start regardeless of RNA-Seq data
* re-ran code in `main_analysis.R`to line 1027
* started drafting in word document (`????_main_text_short_version.docx`)
* saved everythimg beyond line 1027 of `main_analysis.R` to `scratch_code_main_analysis.R`
* saved work space
* commit `72efe63015822d83450c4b5a4cc0acfea946f2dd`
* minor chnages in `main_analysis.R` fot rendering

## 25-04-2023

* started working on manuscript outline, opening main analysis, likely for minor chnages and package installs

## 28-04-2023

* started outlining across multiple documents
* likely continuing coding next and scaffolding in parallel
* main document is  `[...] main_text_scaffolding`
* but see others
* starting to code and write in parallel
* main text document with growing manuscript is in `/Users/paul/Documents/HM_MouseMating/manuscript/230428_main_text_scaffolding.docx`
* starting to code in parallel, while adjusting existing code with script that correspond to text, namely:
  * `000_r_format_data.r`
  * `010_r_define_obesity.r`
  * `020_r_h1.r`
  * `030_r_h2.r`
  * ` 030_r_h3.r`             
* commmit `506381bbbeb6429dba2b4d9bb13aedfa645024be`

## 04-05-2023

* finished an approved manuscript introduction in `/Users/paul/Documents/HM_MouseMating/manuscript/230504_main_text_scaffolding.docx`
* finished initial version of `000_r_format_data.r`
* starting to work in `010_r_define_obesity.r`
* commit `18e41585b59a17220bf14eecbbbb0161f847e7c0`
* revised obesity definition script to include parents obesity status
* commit `ab82740ee5c3712dae66c66ebaf516d475f687aa`

## 05-05-2023

* finished growth curve analysis - updated all scripts and manuscript file accodingly
* commit `18e41585b59a17220bf14eecbbbb0161f847e7c0`
* finished logistic rgression for H1 in code and text
* inserted draft code for H2 and H3

## 10-05-2023

* inserted draft code for H2 and H3
* commit `f370392e69bdb3aeb5c6d177e8efb21f58557314`
* finished logistic regression for H1
* finished well fitting GAM model - updated main text. 
* commit `31fb7fdfb3f2a19da11322778f99166f5a663a32`
* some further inspections
* commit `d58abc6aeb4c1719417b4aea3d858c1920763c43`

## 12-05-2023

* handed off work for RNA seq analysis to AH
* analysis design follows
  * `/Users/paul/Documents/HM_MouseMating/communication/230512_RNAseq_analysis_draft_1.png`
  * `/Users/paul/Documents/HM_MouseMating/communication/230512_RNAseq_analysis_draft_2.png`
  * testing if possible after chatting with AH
  * tesing in script `040_r_h3.r`
  * asdivsing AH to test conditions of factor "ObeseParents" while neglegting onesity status of offspring "ObesityLogical" - not enough data available
  * see overview in "/Users/paul/Documents/HM_MouseMating/communication/230512_RNAseq_data_vs_h3_model_outcomes.png"
  
## 15-May-2023

  * emailled stats consulting unit
  * non-relavent chnages to 3rd hypothesis R script 
  * commit `ebc24c69c5658177b3253eb5c7d153ad4e50c88f`
  
## 16-May-2023

 * updated code for data analysis caffee (DAC)
 * DAC staff fine with model - interpret as linera - possibly check random effect
 * commit `07141cf9a13b35068165abbc7681dabd112f0de2`
 
## 24-May-2023
 
 * attempting to run AHs code 
 * starting to adjust AH's code in `/Users/paul/Documents/HM_MouseMating/analysis/scripts/050_r_array_analysis_ah.r`
 * started adjusting manuscript for males and for array analysis
 * commit `be2856715f6b5b10ba0d5a35dc022cffc2bf5fa8`
 
## 25-May-2023
 
 * updated manuscript with current analysis state
 * about to start PCA with added obesity factor highlight
 * started, and need to continue renaming of data - sample names - and tissue variable names - before continuing
 * commit `6f89b2156f9b3ce38c2590c7eea1388c80415929`
 
## 26-May-2023

 * implemneted overall PCA and updated manuscript 
 * commit `2bded37dc8b779ae7cc113e8a72aa6a61d211341`
 
## 30-May-2023

  * commit before continuing to step through code post initial PCA
  * commit `1454b2f64bff271cbf955243ab6c0b5c3827926f`
  * simplified PCA plotting code
  * ran all PCA plotting code started updating results in manuscript
  * commit `fd09fd425f02e38d05edefdfe58919517a8d50d3`
  
## 31-May-2023

  * started reading in AH's DGE results
  * prepared re-coding dietarey variables in line 175 of array analysis script 
  * commit `865e6c3ac30b286a40cc1c68b30711bfdcf936f4`

## 01-06-2023
  
  * updated todo queue
  * thinned clutter in thus file
  * commit prior to implemnteing DGE in script 50 for obesity variables, instead of AH dietary variables
  * commit `5fd8790e5024bce8e05f885d08f219b1c736ef58`
  * successfully read in and inspected (normalized) array data
  * added analysis ideas to script `050_r_array_analysis_ah.r`
  * commit `3192f2500adfade9f3e3a72a76dd5b1427ddf9ba`
  * tried modelling further but need to find other ideas - talk to TH
  * script 50 is marked at site where analysis ought to be continued - see this file below on what is to do further
  * in summary: experimented with DGE modelling - about to start correlation analysis as in Omentin paper
  * commit `5339189d3579e75926086fb8a81df7cf714bf840`

## 02-06-2023

 * updated todo list
 * isolated experimental model code
 * tested DGE of `ObesityLgl` - none found
 * started for DGE based on parents obesity
 * got plan on what to do
 * got package to annotate array targets
 * commit `0da0360e749065590b4656bce1753e3b85574544`

## 05-06-2023
 
 * finished getting DGE results
 * commit `297fd3d6c49c0f155191d4a1749e6786bb12b6bc`
 * improved DGE contarts definitions by implementing all possible contarstst
 * choose PCA to choose from meaningful results
 * commit `01ab2cb42c87f5eda63e240ca4653fd81fdb9c57`

## 06-06-2023
 
 * consolidated PCA code
 * analysed PC1 and appended results notes to manuscript
 * ready to work on DGE results - see Todo below
 * commit `f1eff498fb52fd33c2ebca18f5bab080a72031c5`
 * inspected coefficients in PCA analysis
 * got annotations for ExpressionData
 * commit `dc8a882bce0c9227c68097a06b80dff68b00e34e` 
 * commit after adding ToDo's to script `ebd8a371b8e586f7725f8f70f7b3bbc12a2baab3`

## 07-06-2023

 * added methods summary text to main text scaffold
 * implemented tissue specific DGE
 * **started** selected relavant DGEs sets based on PCA - select contarts in a way as to get signal in comparison to all others!
 * **next** do an export volcano plots
 * **next** implement KEGG and GO analysis
 * commit `897b6fc4594fd2d469e5ee1ec1362f2ba855d485`
 * updated todo and started selecting coefficients, saved
 * commit `50cb1c3dab4d090bed823bcbee9e56e5327df89d`
 
## Todo queue
 
  * [x] always keep in mind `/HM_MouseMating/manuscript/display_items/230512_RNAseq_data_vs_h3_model_outcomes.png`(or `/Users/paul/Documents/HM_MouseMating/manuscript/communication/190916 Probenliste Clariom S.xlsx`):
    * [x] parental diet conforms exactly with dietary variables considered by AH for array data
    * [x] offsprings obesity statuts does not conform with dietary variables considered by AH for array data
    * [x] the latter is needed but possibly the former data is the only one available
  * [ ] revise script `050_r_array_analysis_ah.r`:
    * [x] **Day 1** - **use to inform DGE and write this down** inspect PCA signals  - overall and tissue specific - for obesity variables -> ***use to choose from DGE results** 
    * [x] **Day 2** - implement DGE - tissue specific - get gene annotations from array targets -> use to inform KEGG and GO - look up how to check model assumptions / possibly email 
    * [x] **Day 3** - ~~finished DGE~~
    * [ ] **Day 4** - modify reference levels in modelling - structure models more - rewrite 
    * [ ] **Day 5** - subset DGE results based on coefficients - isolate each contrats against all others
    * [ ] **Day 6** - export volcano plots
    * [ ] **Day 7** - implement KEGG and GO for all contrast see Tamoxifen project 
    * [ ] **Day 8** - possibly implement upset plots 
    * [ ] **Day 9** - check for errors among expression data annotation get Isolate model code - email Elmar on how to model outcomes directly - on cluster from intensity data by tissue and gene
  * check if this is still valid
    * [ ] establish correlations between varaibles `ObesityLgcl` and possibly also `ObeseParents` and intensities - 
       * [ ] consult ITLN1 paper to find how AH does it - or check teams chat for package
       * [ ] possibly get numerical summaries of PCAs as noted in script 50
    * [ ] try to reimplemnet DGE with `ObesityLgcl` and possibly also `ObeseParents`
    * [ ] adjust AH's Vulcano plot code for new DGE results
    * [ ] adjust AH's UpsetPlot code for new DGE results 
    * [ ] **possibly** implement including other variables instead of body weight body fat
    * [ ] **possibly** use AH's dietary varibles but recode to obesity variable
    * [ ] **possibly** use AH's dietary variables neat
  * [ ] in manuscript
    * [ ] add IMISE consulting personell to acknowledgements:
      * **Maciej Rosolowski** <maciej.rosolowski@imise.uni-leipzig.de>
      * Maryam Yahiaoui-Doktor <maryam.Yahiaoui-Doktor@imise.uni-leipzig.de>
    * [ ] add HM consulting personell to paper
      * **Elmar Spiegel**
      * Roman Le Glut
      * Kainat Khowaja
  * [ ] archive old scripts and data files
    * `/Users/paul/Documents/HM_MouseMating/analysis/scripts/inspect_data.R`
    * `/Users/paul/Documents/HM_MouseMating/analysis/scripts/main_analysis.R`
    
## Ideas not yet pursued

 * foo 

## Further reading

### Modelling

 * **Zuur AF, Ieno EN, Walker N, et al. 2009. Mixed effects models and extensions in ecology with R. New York, NY: Springer.**
