# MSB1015-Assignment-2

[![GitHub License](https://img.shields.io/github/license/Rrtk2/MSB1015-Assignment-2)](https://github.com/Rrtk2/MSB1015-Assignment-2/blob/master/LICENSE.md) ![](https://img.shields.io/badge/Status-Wrapping_up-green) [![GitHub Watches](https://img.shields.io/github/watchers/Rrtk2/MSB1015-Assignment-2.svg?style=social&label=Watch&maxAge=2592000)](https://github.com/Rrtk2/MSB1015-Assignment-2/watchers) 


#### What is this project about
This repository is the final product of assignment 2, requested by the course MSB1015 (Scientific Programming). The goal is to create a repository which contains all information, documentation and files needed to run a wikidata query using html/javascript as a new user.
The setup of a [GitHub page](https://rrtk2.github.io/MSB1015-Assignment-2/) needs the following documents: _site.yml and index.Rmd. _site.yml contains the architecture of the github page. index.Rmd contains the actual R markdow script, which will be transferred into the website. When these files are set up, the rmarkdown::render_site() function is ran in the main (_site.yml) folder, this assembles the site and a _site folder will be created. The contents of the folder need to be transferred into the ~/docs/ folder of the repo. Github pages need to be enabled and set to the ~/docs/ folder of the repo. 


#### Usage
Visit the [GitHub page](https://rrtk2.github.io/MSB1015-Assignment-2/) to see the result of this project. To run the code on your local computer, go to the [GitHub raw code page](https://rrtk2.github.io/MSB1015-Assignment-2/raw.html) and insert this code into your local R version. This will install and load the required libraries, do the query, use machine learning and create 2 models which can be used to predict the boiling temperatures. In block 2 (SETTINGS) several parameters can be adjusted, such as the training size or query call. 

#### Expected output
Using the default settings the output should be very similar as observed in the [GitHub page](https://rrtk2.github.io/MSB1015-Assignment-2/). Any updates which alters the result of the query will alter the results when the script is re-ran.

#### Project structure
##### Where does data come from?
The query asks information from [Wikidata](http://wikidata.org) in a similar fashion the dedicated [Wikidata database query](https://query.wikidata.org/) works using the SPARQL language. Data on wikidata is published under the [Creative Commons Zero](https://creativecommons.org/share-your-work/public-domain/cc0) license, stating 'others may freely build upon, enhance and reuse the works for any purposes without restriction under copyright or database law'.

##### How is data shared, in what format, with what protocols?
Using the tool developed in this project, data is shared using the [wikidata-sdk](https://www.wikidata.org/w/api.php). 

##### Workflow
The following workflow is applied, coherent to the script blocks:

1) Install packages

2) Define user settings

3) Set user functions

4) Define WikiData query call & extract data

5) Convert units to Kelvin

6) (Optional) filtering and data cleanup

7) Visualize raw data 

8) Data enrichment (using rcdk)

9) Latent variable filtering

10) Data subsetting for machine learning

11) Machine learning (PLS model)

12) Machine learning (RandomForest model)


#### Contact
ra.reijnders@student.maastrichtuniversity.nl


#### License and contributing guidelines
[License](/LICENSE.md) 

[Contributing guidelines](/CONTRIBUTING.md) 


#### Who is involved, and what are their roles.
RRtK2 (owner and contributor)


#### Status of project
Wrapping up. Minor adjustments to the [GitHub page](https://rrtk2.github.io/MSB1015-Assignment-2/) might occur, as well as to the scripts. Bugfixes.


#### Copyright and authors
All code and documents in the MSB1015-Assignment-2 folder was created by [these author(s)](/AUTHORS.md).
