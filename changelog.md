# On Mouse Mating work

## Project goals

### Central questions

* Determine the Influence of parental weight for offspring weight gain.
* Should mothers loose weight before becoming pregnant?

## Originally formulated goals

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
 * updated todo in README
 * commit `7e9d17ca02f1e54fba1afda157b5b85ccd091c13`
 * added testing all reference leveles in analysis of 1st components

## 08-06-2023

 * read PC-Regression results and designed contrasts
 * **next:** subset DGE results to noted contrasts 
 * **next:** compile one well-labelled DGE list
 * commit `28bdb73a2dece35ca9c1252ca9bf088a4355ec72` 
 * commit prior to extending parents' DGE list `75ea960ab4b05a7ad80a50e6d3421604bb4395cc`
 * finished getting a list with all DGE results - saved workspace
 * commit `58bf6dbc4989ca2e16a0a336860db5b534b891f0`

## 12-06-2023

 * started on Volcano plots
 * commit `58bf6dbc4989ca2e16a0a336860db5b534b891f0`

## 14-06-2023

 * adjusting and saving Volcano plots
 * continue by volcano updating plots in manuscript
 * from `/Users/paul/Documents/HM_MouseMating/analysis/plots`
 * commit `e41733678e8c6f30de6952ad7b5b290ed3c36d2d`

## 21-06-2023

 * continuing in `050_r_array_analysis_ah.r`
 * implemented KEGG and GO analysis
 * commit `26119d5a2c7718ff93a3d7441d92d15906f3ae63`

## 22-06-2023
 
 * thinned ou entire folder structure by erasing superfluous output files
 * did not erase any scripts
 * re-ran the following scripts consecutively 
   * `000_r_format_data.r`
   * `010_r_define_obesity.r`
   * `020_r_h1.r`
   * `030_r_h2.r`
   * `040_r_h3.r`
   * `050_r_array_analysis.r`
   * also see corresponding `.Rdata` files
 * commit `c8d1cd7e774b6a71f2b503049307676fe9fc4230`

## 28-06-2023

 * to get overview creating a results summary document with display items
 * name `/HM_MouseMating/manuscript/230628_results_summary.docx`
   * re-ran the following scripts consecutively 
   * `000_r_format_data.r` - minor changes to formatting
   * `010_r_define_obesity.r` - added DI export code - including table summaries - writes to `/Users/paul/Documents/HM_MouseMating/manuscript/display_items`
 * commit after reaching line 292 in `010_r_define_obesity.r` and prior to removing Mixed group in `000_r_format_data`
 * commit `b7bd0758c5c037169ed24dcae421bcad6c9cbcdf`  
 * re-ran `000_r_format_data.r`
   * removed mixed group among f1 and group information among f0 since it was confusing
 * commit `640149bdafaeb7e3f9ac7d61336a2a43edcce16b`

## 29-06-2023

 * adding Roxygen-style markdown to  `000_r_format_data.r` - see `https://bookdown.org/yihui/rmarkdown-cookbook/spin.html`
 * re-running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/000_r_format_data.r')`
 * re-running `010_r_define_obesity.r` from start
   * found f0 individuals with mixed diets `010_r_define_obesity.r` in line 362 
   * commit - filter in `000_r_format_data.r` re run to `010_r_define_obesity.r` line 362
   * commit `a54481c0fc4fd590d7a0f99e53f30fe6ce96123b`
 * adding some Markdown structure to `010_r_define_obesity.r`
 * running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/010_r_define_obesity.r')`
 * commit `be7b1fbed71ff2c4eeca66e0a46d631b85dcba39`
 * continued writing - finished main text methods - need to revise supplemental methods
 * commit `f6cc4a51b2021f9214638b322b90e3078d7013b3`

## 30-06-2023 - work day 44
  
 * in `000_r_format_data.r` checking filtering out  of f1 mice with mothers of mixed diets 
 * re-running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/000_r_format_data.r')`
 * re-running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/010_r_define_obesity.r')`
 * running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/020_r_h1.r')`
 * including model formula summary as per 
   * `https://datalorax.github.io/equatiomatic/`
   * `https://bookdown.org/yihui/rmarkdown-cookbook/equatiomatic.html`
   * `https://vincentarelbundock.github.io/modelsummary/articles/modelsummary.html`
 * running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/030_r_h2.r')`
 * running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/040_r_h3.r')`
 * started on results section 
 * before continuing writing - need refreshed array analysis results first
 * commit `97a250fa653f36b5a8452d1e56af4a3a8631a19e`
 * commit prior to altering script 50 due to missing mother obese in LIAT
 * commit `d7c10c120f5c516b362620484ab8e4e1d825b9e1`
 * after attempting to deal with missing data will ignore LIAT
 * running `git checkout d7c10c120f5c516b362620484ab8e4e1d825b9e1 -- 050_r_array_analysis.r`
 * including code to check availability of obesity variables for expression data
 * commit `bc46bba79b132b56c66a5be6c464f2f610dd47f0`

## 04-07-2023 - work day 45

 * re-ran `050_r_array_analysis.r` to break point in line 1010 and saved workspace image
 * updated results section in manuscript and SI figure
 * re-running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/050_r_array_analysis.r')`

## 05-07-2023 - work day 46

 * updated results section
 * commit prior to dealing with `050_r_array_analysis.r` and missing contrast in LIAT in script 
 * commit `cc0d696935d32eb267ea227572ea2c8971a67234`
 * getting DEG results falso fro LIAT with reduced set
 * commit `f8ca937c8142e37258185ec8d98cd5723b5a3572`
 * continuing to update manuscript
 * extending `050_r_array_analysis.r` - adding DGE list export to: 
   * `/Users/paul/Documents/HM_MouseMating/analysis/tables/050_r_array_analysis_BRAT_TTL_down.xlsx`
   * `/Users/paul/Documents/HM_MouseMating/analysis/tables/050_r_array_analysis_BRAT_TTL_sign.xlsx`
   * `/Users/paul/Documents/HM_MouseMating/analysis/tables/050_r_array_analysis_LIAT_TTL_down.xlsx`
   * `/Users/paul/Documents/HM_MouseMating/analysis/tables/050_r_array_analysis_LIAT_TTL_sign.xlsx`
 * commit `482459afbd07f16e6af1a7f94f95cd63103a75b7`
 * finished script `050_r_array_analysis`, let it run in full, and saved full image
 * updated results in manuscript - **only GO terms still need doing**
 * commit `224d11d65037d253652c9a131548fe68c10684f1`

## 06-07-2023 - work day 47 - at retreat 
 
 * extending `050_r_array_analysis.r`
 * drafted `md` formatting, by duplicating second-level-headings
 * rendered `md` version
 * commit `363b82da3bdb880b8ec2540aea6c495d96b103d2`
 
## 10-07-2023 - work day 48 

 * updating Figure 1 for manuscript, and script `010_r_define_obesity.r` for summary (see line 115)
 * commit `a55ace13bd77f28e3c7a2cb0a18a63f86e0076de`
 * started on Discussion in main text
 * corrected cases in captions of script `010_r_define_obesity.r`
 * commit `ef91e56bd34e14303fa7639b9fadfe8f6ba9feb4`
 * re-running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/050_r_array_analysis.r')`

## 11-07-2023 - work day 49

 * re-ran `050_r_array_analysis.r` after correcting table export function
 * commit `cce6f274fbb3f7f5de3d83d6993b20fd4b819ead`
 * added messy heat map code to `050_r_array_analysis.r`
 * re-running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/050_r_array_analysis.r')`
 * re-running `050_r_array_analysis.r`
 * commit `65bf8fb487d792018e1b0122f7860a85dabb1319`
 * appended messy heat map code to `050_r_array_analysis.r`
 * re-running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/050_r_array_analysis.r')`
 * commit `d1279a3497cad29886b0bc4e4fbebf34cd548b90`

## 12-07-2023 - work day 50

 * added script `060_r_arrange_display_items.R` 
 * which creates `060_r_arrange_display_items__fig_2.png` for inclusion into manuscript
 * re-ran `050_r_array_analysis.r` to get new input files for above script
 * commit `1781c5c70efd6e68b541b7382dcc8aa451540b2a`
 * re-running `rmarkdown::render(input = '/Users/paul/Documents/HM_MouseMating/analysis/scripts/050_r_array_analysis.r')`

## 14-07-2023 - work day 51

 * finished and mailed of manuscript - see `/Users/paul/Documents/HM_MouseMating/manuscript/230714_main_text_to_coauthors/230714_main_text.docx`
 * slightly rearranged folder structure in `/Users/paul/Documents/HM_MouseMating/manuscript`
 * updated todo below
 * commit `557e386cc0584dfc0e43fe114ace800703394b77`

## 21-08-2023 - work day 52

 * received edited manuscript version from NKB
 * stored:
   * mailed-off version (original): `/Users/paul/Documents/HM_MouseMating/manuscript/230714_main_text__reference.docx`
   * mailed-off version (revised):`/Users/paul/Documents/HM_MouseMating/manuscript/230714_main_text__edited_by_NKB.docx`
   * received version (to be revised further):`/Users/paul/Documents/HM_MouseMating/manuscript/230821_main_text.docx`
   * started accepting minor changes to above manuscript. 

## 22-08-2023 - work day 53
 
 * met with NKB - updated todo
 * commit `b96c627bd7e22f13357108cecf73d81ca4596f2d`
 * created and switched to analysis fork - `git checkout -b "more_bars"`
 * created additional version of SI Fig 1 b for main text
   * see text 22-Aug-2023 - newly inserted Fig 2
   * see fork `more_bars`, script `010_r_define_obesity.r`, line ~480
 * commit `3c9fbb5857428cd06a28060eecd779f411e77bf9`(on `more_bars`)

## 29-08-2023 - work day 54 - working in revision suggestions of NKB (see ToDo)

 * in `/Users/paul/Documents/HM_MouseMating/analysis/scripts/010_r_define_obesity.r`
   * slight re-arrangment of plot saving code
   * re-ran full script
 * commit `42db8c441a5018ee3bbeebcdf593e0ef507a78be`
 * created and switched to analysis fork - `git checkout -b "slim_mouse"`
 * in `/Users/paul/Documents/HM_MouseMating/analysis/scripts/010_r_define_obesity.r`
   * slight re-arrangment of plot saving code
   * updated plotting code - added diet to F0 
   * re-ran entire script
 * in `/Users/paul/Documents/HM_MouseMating/analysis/scripts/050_r_array_analysis.r`
   * adjusted heat map plotting code to plot z-scores, and updated figure manuscript
 * starting revising main text - proceeded to ened of Methods section
 * commit `13dc723cc44ae31212ffc5ce57e0a19b5cd9a435`
 * updated Todo queue 
 * commit `4d299b57cc06654b51b5b623e2608fbe9a967829`

## 30-08-2023 - work day 55 - working in revision suggestions of NKB (see ToDo)

 * corrected figure caption for heat maps in main text
 * in `/Users/paul/Documents/HM_MouseMating/analysis/scripts/010_r_define_obesity.r`
   * checking obesity definition in F0 - with CD never obese - with HFD depending on curve analysis - see line ~365 
   * checking obesity definition in F1 - always depending on curve analysis - see line ~390
 * worked through todo list below
 * commit before explorig odd obesity value assignment among F1 - see below
 * commit `e3af6b4daa580840f880b4b52b6f1cfd72dda613`
 * leaving branch `more_bars` open ended
 * **merging current branch** `slim_mouse` to `master` as per `https://www.w3docs.com/snippets/git/how-to-make-the-current-git-branch-a-master-branch.html`
 * **creating analysis fork** `adjust_cutoffs` to change curve cutoff of F1 to that of F0 `git checkout -b adjust_cutoffs`
 * working in `010_r_define_obesity.r` as listed in ToDo queue
 * re-running `010_r_define_obesity.r` entirely
 * re-running `020_r_h1.r` entirely
 * re-running `030_r_h2.r` entirely
 * re-running `040_r_h3.r` entirely - added exporting plot `gam6_plot` as `SI Figure2`
 * not re-running `050_r_array_analysis.r` - results shouldn't have changed
 * continuing to edit file `230830_main_text_b.docx`
 * commit `bf92b4919f9baaac0d892e3a116e06ea719f1364`
 * **merging current branch** `adjust_cutoffs` to `master`
 * commit `ce0887565d9dcaaa915b820b6fda5c72d49bc533`
 * copying and continuing to edit file `230830_main_text_b_cleaned.docx` - removing superflous tracked changes
 * updated ToDo
 * copying cleaned `230830_main_text_b_cleaned` to `230831_main_text`, mailing latter of to NKB
 * summary of current state in `/Users/paul/Documents/HM_MouseMating/communication/230830_on_revision_to_NKB.pdf`
 * commit `08c0d9f27dce5ef0d68ef7adf8426450c0dcadb2`
 * suspecting error in growth curve analysis 
 * in `010_r_define_obesity.r` adding inspection code at `Isolate measurement days ----`
   * using obesity definition by left and right of 0 result stays the same
   * corrected figures in manuscript file 31.08.2023
   * commit to new fork `new_curvatures` - `8bee3da04e795336a6a80a50beb6356428208519`

## 31-08-2023 - work day 56 - working on main text

 * expanded Acknowledgements, started funding section
 * received comments from NKB
 * updated ToDo 

## 04-09-2023 - work day 56 - working on main text

 * [x] working on submission version `/HM_MouseMating/manuscript/230904_manuscript_to_coauthors/230904_main_text.docx`
   * [x] revise text
   * [x] draft abstract
 * [x] updated ToDo

## 28-09-2023 - work day 57 - revising main text

 * working in `/Users/paul/Documents/HM_MouseMating/manuscript/230928_main_text.docx`
 * merged branch `new_curvatures`
 * opened branch `revise_array_analysis`
 * starting to adjust `050_r_array_analysis.r`
   * adjusted lfc and p-cutoffs in `get_one_volcanoplot`
   * adjusted lfc and p-cutoffs in heatmap code
   * commented out code that is irrelevant and starting to bit rot due to package updates
 * commit `61a1404cd9fb96dc0e2b9c3a208ce399a719af82`

## 09-10-2023 - work day 58 - revising main text

 * working in `/Users/paul/Documents/HM_MouseMating/manuscript/231009_main_text.docx`
 * continued to adjust `050_r_array_analysis.r`
   * updated packages, checking which packages aren't needed
   * corrected PCA plot labels and exported, break point set
   * corrected PCA plots in manuscript SI
   

## 10-10-2023 - inset - exporting data for other project

 * adjusted `050_r_array_analysis.r`
   * added export code to save prepared input data (around line 665)

## 13-10-2023 -  work day 59 -  improved boxplot
  
 * in `010_r_define_obesity.r`
   * improved boxplot
 * commit `1d7ad2fe18eff90e9ed62966dd4a31025aa30f95 `

## 17-10-2023 -  work day 60 - adjusted figures
 
  * in `050_r_array_analysis.r`
    * [x] re-run from first breakpoint in line 672 - environment not saved
    * [x] see break point in R script (environment saved up until there)
    * [x]see mark in manuscript
    * [x] ~~~for volcano plot plots possibly remove zero-count genes~~~
    * [x] adjusted GO and Kegg plots
    * [x] re-ran full script from start 
  * [x] adjusted and re-ran `060_r_arrange_display_items.R`
    * [x] now creates Figs 3, 4, 5
  * adjusted manuscript further - inserted new figures - modified text (significance definition for KEGG and GO analysis)
  * commit `88e9e299782a3e2a792ad535b692d0fa9de0df23`
 
## 18-10-2023 -  work day 61 - adjusted figures

  * in manuscript of todays date adjusted results and discussion
  * [x] in `050_r_array_analysis.r`
    * [x] adjusted plotting functions
    * [x] now alos exporting all required tables
    * [x] re-ran entirely
  * [x] adjusted and re-ran `060_r_arrange_display_items.R`
  * [x] commit `cfc340befcd6717c343985c347af1fae6712a610`

## 19-10-2023 -  work day 61 - revised results in manuscript

  * [x] in `050_r_array_analysis.r`
    * [x] adjusted go plotting functions - now exports table with gene annotations
  * [x] commit `0f26267c94138eee542162af37fb92c8e46acca5`

## 25-10-2023 -  work day 62 - preparing submission

  * in `/HM_MouseMating/manuscript/231025_submission_1`
    * finished manuscript
    * main tex figures
    * SOM tables
  * in `/HM_MouseMating/analysis/scripts`
    * updated `.gitingnore`
  * cleaned git repository using `git filter-repo --invert-paths --path path_to_file --force`
    * to retore file use backup before 25-Oct-2023 11:00 CET
    
## 03-11-2023 -  work day 63 - adjusting code for seminar talk

 * on branch `seminar_talk`, in `010_r_define_obesity.r`
   * adjusting code for `mice_f1_slct_xyplot_final`
   * exports males and females separately
   * to `/HM_miscellaneous/231103_hm_seminar/231103_hm_on_mated_obese_mice`
 * commit `a93ba018086272aebbe1cfb0fae4440356e4010a`
 * in `050_r_array_analysis.r`
   * added PCA plot code for seminar, exported to plots folder, exported as `.png` and annotated
 * commit `2fc00677eafd1c4824a1a60d981556c6ba3f18e3`

## 07-11-2023 -  work day 64 - gave seminar talk

 * gave seminar talk
   * `/231103_hm_seminar/231106_hm_on_mated_obese_mice_extended.pdf`
  * updated todo

## 13-11-2023 -  work day 65 - gave seminar talk

 * gave retreat talk
   * `/231103_hm_seminar/231103_hm_on_mated_obese_mice_compressed.pdf`
  * updated todo
  * commit `9de302375c7153ca6b1f8c1e0a972b4f282ccc6d`

## 14-11-2023 -  work day 66 - learned more about GAMS and PCA interperations across 2 statistics course
 
 * updated todo - revisions needs improvent for
   * actual diets
   * PCA interpretation - imprecise and possiby odd
   * GAM building of offsprings obesity, likely wrong

## 21-03-2024 -  work day 67 - starting on first revision

  * re-installed packages 
  * commit prior to changing code
  * started editing `/HM_MouseMating/manuscript/240321_submission_2_preparation/240321_revision_notes.docx`
  * adding studies for first revision to distinct sub-folder
  * updated Todo (below), for revision
  * created branch for first revision

## 22-03-2024 -  work day 68 - working on first revision

  * text editing with current date - no commit yet - only text work

## 25-03-2024 -  work day 69 - working on first revision

  * text editing with current date - no commit yet - only text work
  * see different versions with todays' date 
    * ended up outlining introduction first
    * saved introduction draft in sperate file 
    * **this latest file also needs to be used to build the introduction in the new, revised docoment**

## 10-04-2024 -  work day 70 - working on first revision

  * checking revised introduction in file `/HM_MouseMating/manuscript/240321_submission_2_preparation/240410_4_outline_introduction.docx`
  * adjusting previous manuscript file with new introduction and further notes throughout (`/Users/paul/Documents/HM_MouseMating/manuscript/240321_submission_2_preparation/240410_3_main_text.docx`)
  * working on code branch `mouse_mating/first_revision`
  * finished checking script `/000_r_format_data.r`
    * [x] found - check for dietary magnetic resonance imaging data and start carrying this through the analysis
    * [x] found - check for litter size information
  * for obesity definition and comparison saved reference weights table from `https://www.jax.org/jax-mice-and-services/strain-data-sheet-pages/body-weight-chart-005304`
    * file is at `/HM_MouseMating/analysis/raw_data/240410_reference_weights.xlsx`
    * caption was: "Groups of 80 males and 80 females were weighed weekly.  Mice were fed a diet containing 6% fat (LabDiet 5K52 formulation). Values represent mean and one standard deviation. Ages are ± 3 day"
  * previous file edited agin with information from orignal provider `https://www.taconic.com/products/mouse-rat/standard-strains-and-stocks/black-6-b6ntac#tabsmobiledropdown-155d2d1ed8-item-dd16099164-tab` 
    * also see info sheet downloaded to Zotero
  * plotted out weights with reference weight  in  `/010_r_define_obesity__mice_weights_references.pdf`
  
## 11-04-2024 -  work day 71 - working on first revision 

  * updated todo - multiple times
  * worked in script `000_r_test_saemix.R`
    * finished training for Nonlinear Mixed-Effects Growth Models
    * using `https://doi.org/10.5964/meth.7061`
    * to substitute hypothesis testing
      * defined in introduction
      * currently done in `020_r_h1.r` `030_r_h2.r` `040_r_h2.r` 
      * with saemix code
  * worked in script `010_r_define obesity.r`
    * checked wether or not Diet and Obesity variables overlap perfectly
      * line 467: **F0 8992 had HFD but no gain - is this important for the anaslysis?**
      * line 491: **all F1 have low weight gain, regardless of diet**
      * line 500: here is definition of "Obesity"
        * `F0 with any HFD and hi weight gain as classified "obese" and all others are "not obese"`
        * `F0 with any HFD and hi weight gain as classified "obese" and all others are "not obese"`
        * "Obesity" variable indicates HFD of parents, rightfully apart from **F0 8992**
    * during execution gathered notes for writing
      * "obesity" variable can be substituted with HFD throughout the text, for DEG analysis
      * remove mention of curve trajectories curve trajectories
      * where needed still: these were useful in confirming that parent individual 8992 didn't gain weight
      * effects of diet needs to be evaluated with saemix
      * lattice plots were re-exported without "obesity" and "weight gain" difinitions
      * tables were re-exported without "obesity" and "weight gain" difinitions
      * model descriptions can describe teh effcet instead of the lattice plot mentioned in the previous line

## 12-04-2024 -  work day 72 - working on first revision 

  * created script `015_r_check_accumulation_curves.R` to use with SAEMIX
  
  
## Todo queue (last updated 11-04-2024)
 
### **revision work** - after first submission
 
 * based on `/HM_MouseMating/manuscript/240321_submission_2_preparation/240321_revision_notes.docx`
   * [x] commit - branch repository - commit
   * [x] adress as many comments as reasonable without re-running code
   * [x] revise introduction
   * [ ] re-run analysis - see how "diets" could be used instead of "obesity" 
     * ~~[ ] consider using the the reference data beyond the initial plots~~
     * ~~[ ] starting point - fit Gompertz curves and compare significant changes in `alpha` parameter - see `https://www.ipb.pt/~vcadavez/websiteVC/tutorial/rcode/2019-04-28-gompertzmodel/`~~
     * ~~[ ] starting point - use `nlme` - see `https://www.r-bloggers.com/2019/09/fitting-complex-mixed-models-with-nlme-example-3/`~~
     * ~~[ ] starting point - compare different growth models - see `https://cran.r-project.org/web/packages/biogrowth/vignettes/v04_model_comparison.html`~~
     * [x] keep all current results in script `010_r_define obesity.r`
     * [ ] of script `010_r_define obesity.r` integrate re-exported items
     * [ ] learn how to correlate variables with trajectories
       * [x] look at `saemix` example
       * [ ] ~~correct GAM modeling of offspring obesity as in Gavin Simpsons rat hormone hgam example, Physalia course GAMs in R, day 4, 23.11.2023, possibly also pull GS repoitory from course and look at "chick example"~~
       * [ ] ~~possibly use GAM approach - see Winter and Wieling (2016), https://doi.org/10.1093/jole/lzv003~~
     * [x] using SAEMIX instead: 
     * [ ] see what to do with the weight trajectories
       * [ ] substitute `020_r_h1.r` `030_r_h2.r` `040_r_h2.r` with saemix code
       * [ ] both parents obese - offspring sex neglected - male  - female
       * [ ] one parent obese - offspring sex neglected - male - female
     * [ ] likely keep all old results in script `50_r_array_analysis.r`, 
       * [ ] but with relabel figures (in text or code)
       * [ ] check wether axis correct axis is being looked at in PCA
       * [ ] verify DEGs with human or web data
     * [ ] improve figure 1 - diet labelling
 * [ ] re -outline individual manuscript sections after analysis
   * [ ] consider the rebuttal document on or after 25.3.2025
   * [ ] revise methods
   * [ ] revise results
   * [ ] copy other sections from outline document first
 * [x] consider litter size 
   * [x] started to comment in revision file
   * [x] addressed in manuscript file
   * [x] to correct body weight after weaning consider (Zhang et al. 2012)
   * ~~[x] possibly report on the most profound effects of the aforementioned weight adjustments~~
   * [ ] consider NKBs comment on litter sizes - see communications folder 25-Mar-204
  

- note reference for GAM modeling, from Zotero
- in script 50 
  
  
### **revision work** - previously urgent

 * [ ] re-eavaluate all interpretation in the context of
   * [ ] correct diets - Western vs Western Control
   * [ ] correct interpretation of PCA
   * [ ] use `dimdesc()`as in course notes and code examples Multivariate Statistics 1 23.11.2023 (**see example sheet (R) for anova on treatment in PCA** as [here](/Users/paul/Documents/HM_miscellaneous/231116_cf-statcon_mv_statistics/231123_day_3__exercises.R))

### **revision work** - previously  additional notes

 * [ ] possibly revise description of contrast finsing and definition, refer to seminar talk 7-Nov- 2023
 * [ ] AJ: check diet definition and composition 
 * [ ] OB BS: Check challenge for ITT GTT to see if they have a metabolic syndrom
 * [ ] check communications folder for
 * [ ] consider diet descriptions and other comments NK
  * [ ] reviewer link AH
     
### **manuscript work** for first submission   

 * [x] circulate to MB 
 * [x] revise as per MB and AH
 * [x] rename heat map plot
 * [x] export KEGG and GO terms with gene names
 * [x] format plots
 * [x] implement GO term lookup
 * [x] circulate to coauthors 
 * [x] revise results and discussion further
 * [x] copy DEG and GSEA results tables
 * [x] upload code
 * [x] upload array data upon publication
 * [x] insert array dat link

### possible further anslysis work

 * [x] create analysis fork - `more_bars`
   * [x] for NKB create additional version of SI Fig 1 b for main text
   * [x] calculate weight deltas for each f1 mouse between day 0 and 100 ("A...")
   * [x] plot weight deltas by factor
   * [x] include in main text
   * ˜˜˜[ ] _possibly_ omit or use obesity factors instead, where available ˜˜˜ 
 
 * [x] create analysis fork - `slim_mouse`
   * [x] deal with mouse 8989
   * [x] check if it is removed form analysis as it should - if so add to main text
   * [x] if it isn't removed re-run analysis without it
 
  * ~~~[ ] _possibly_ create analysis fork - `weight_deltas`~~~
    * [ ] redefine obesity using weight deltas instead of curves and re-run everything
    * [ ] in script `010_r_define_obesity.r`
      * [ ] possibly define obesity using weight gain deltas instead of growth curves
      * [ ] re-check weight gain definitions among F0 and F1

 * [x] manuscript work
   * [x] see above 
   * [ ] improve presentation of weight variables
     * [x] add F0 weight deltas in `010_r_define_obesity.r`
     * [x] to F1 and F0 weight deltas attempt to add assigned obesity status, e.g. via differently-colored barplot-overlays, or sub figures, in `010_r_define_obesity.r`
       * [x] update Fig. 2 in main text, to Fig 2a and b. 
     * ~~~[ ] add diet to F1 curve plot in `010_r_define_obesity.r`, as already done for F0~~~ **not needed, all CD, added to caption**
     * [x] update SI Figure 1b in manuscript
     * [x] check caption of SI Figure 1b
     * [x] add curvature summaries to SI Fig 1
     * [x] explore odd obesity value assignment among F1
       * [x] commit to save state `81033ebad16be523b7096e8a841697e7fd037b33`
       * [x] backing up `010_r_define_obesity__mice_derivatives_densities.pdf` to `010_r_define_obesity__mice_derivatives_densities_backup.pdf` 
       * [x] create analysis fork `adjust_cutoffs` to change curve cutoff of F1 to that of F0 
       * [x] in  `010_r_define_obesity.r`  
         * [x] check and remove missing data among F1 - week 16
         * [x] re-run script and look for changes, decide what to do next
         * [x] save new version of manuscript file, named `230830_main_text_b.docx`
         * [x] update manuscript file 
   * [x] add count values to manuscript
     * [x] see `010_r_define_obesity__mice_f1_slct__obesity.xlsx` and  line ~415 of `010_r_define_obesity.r`
     * [x] add parents and their obesity status to current SI table 1 which is `040_r_h3__rna_seq_sample.xlsx`
     * [x] update count values in text from `040_r_h3__rna_seq_sample.xlsx`
   * [x] validate caption of Fig. 3  
   * [x] update all analysis and text
   * [x] revise draft with coauthor's suggestions
   * [x] circulate to NKB
     * [x] see `/Users/paul/Documents/HM_MouseMating/communication/230831_on_finishing.pdf`
     * [x] meet with NKB
     * [x] work in NKB edits
     * [x] use commnets from `/Users/paul/Documents/HM_MouseMating/manuscript/230831_main_text_NK.docx`
     * [x] apply commnets in file `/Users/paul/Documents/HM_MouseMating/manuscript/230831_main_text.docx`
     * [x] add abstract
     * [x] add acknowledgements, funding
   
   * [ ] communicate and **possibly** 
     * [x] in fork `new_curvatures` 
       * [x] implement old analysis with corrected figures
       * [ ] explore correction of growth curve analysis
       * [ ] and re-run analysis with corrected obesity definition
     * [ ] in fork `weight_delta` fork explore alternative to growth curve analysis and re-run analysis, use weight deltas
   * [ ] possibly - duplicate and load duplicated environment of `050_r_array_analysis.r`, in heatmap remove "inferred obesity status of the sequenced individual"
   * [ ] correct figure labels - in `pdfs`'s and code - use 
     * [ ] epigonal visceral (EVAT), 
     * [ ] subcutaneous (SCAT), 
     * [ ] liver (L) and 
     * [ ] brown adipose tissue (BAT)
   * [ ] expand SI table 1 to female and litter information 
   * [ ] prepare submissions **or**     
   
 * [x] always keep in mind `/HM_MouseMating/manuscript/display_items/230512_RNAseq_data_vs_h3_model_outcomes.png`(or `/Users/paul/Documents/HM_MouseMating/manuscript/communication/190916 Probenliste Clariom S.xlsx`):
   * [x] parental diet conforms exactly with dietary variables considered by AH for array data
   * [x] offsprings obesity statuts does not conform with dietary variables considered by AH for array data
   * [x] the latter is needed but possibly the former data is the only one available

## Additional ideas, if desired, or required

 * [ ] in script `050_r_array_analysis.r`
   * [ ] check if DGE, KEGG, GO results can be interpreted automatically with some sort of comparison data
   * [ ]  possibly use running score and preranked list of GSEA result [see here](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html).
   * [ ]  clean out heat map code - create function, name variables properly
   * [ ] **someday** - possibly Upset plots
   * [ ] **someday** - possibly network plots
   * [ ] **someday** - possibly implement SWAMP instad of PCA
 
 * [ ] check if maternal weight loss can be adressed with statistical analyses using mixed group
 
 * [ ] archive old scripts and data files
   * `/Users/paul/Documents/HM_MouseMating/analysis/scripts/inspect_data.R`
   * `/Users/paul/Documents/HM_MouseMating/analysis/scripts/main_analysis.R`
