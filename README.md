# stomata-independence

This repository contains source code associated with the manuscript:

[Developmental integration cannot explain major features of stomatal anatomical evolution in seed plants](https://doi.org/10.1101/XXXXX). *bioRxiv*.

## Author contributions

[Christopher D. Muir](https://cdmuir.netlify.app) designed the study, compiled data, analyzed data, and wrote the manuscript with input from all authors. 

These authors contributed data:

* Miquel Àngel Conesa
* Jeroni Galmés
* Varsha S. Pathare
* PatriciaRivera
* Rosana López Rodríguez
* Teresa Terrazas
* Dongliang Xiong

## Contents

This repository has the following file folders:

- `alignments`: output (sequence alignments and partitions) from [PyPHLAWD](https://github.com/FePhyFoFum/PyPHLAWD) for tree inference using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)
- `figures`: figures generated from *R* code
- `objects`: saved objects generated from *R* code
- `processed-data`: processed data generated from *R* code
- `ms`: manuscript input (e.g. `ms.Rmd` and `stomata-independence.bib`) and output (e.g. `ms.pdf`) files
- `r`: *R* scripts for all data processing and analysis
- `raxml`: out from tree inference using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)
- `template`: manuscript style files

## Prerequisites (copied from https://github.com/daijiang/workflow_demo):

To run code and render manuscript:

- [*R*](https://cran.r-project.org/) version >4.1.0 and [*RStudio*](https://www.rstudio.com/)
- [LaTeX](https://www.latex-project.org/): you can install the full version or try [**tinytex**](https://yihui.org/tinytex/)
- [GNU Make](https://www.gnu.org/software/make/): In terminal, you can just type `make paper` to render the manuscript. You can also use it to re-run all scripts.

Before running scripts, you'll need to install the following *R* packages:

```
source("r/install-packages.R")
```

To fit **brms** model, set up [**cmdstanr**](https://mc-stan.org/cmdstanr/).

To infer phylogenies in RAxML using code in script `raxml/run-raxml.sh`, you will need to follow instructions here: [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)

**Not recommended:** To get sequence clusters and alignments using [PyPHLAWD](https://github.com/FePhyFoFum/PyPHLAWD), you will need to follow directions to set up PyPHLAWD, blast, and download the GenBank plant and fungal sequence data base. Scripts for this are not in this repository because they require a bunch of dependencies to work and the database is too large. However, the output is deposited here for inspection.

## Downloading data and code 

1. Download or clone this repository to your machine.

```
git clone git@github.com:cdmuir/stomatal-independence.git
```

2. Open `stomatal-independence.Rproj` in [RStudio](https://www.rstudio.com/)

## Rendering manuscript

### Software requirements

At minimum, you will need [R](https://cran.r-project.org/) installed on your machine. Install additional packages by running `r/install-packages.R`.

### Rendering manuscript with pre-saved outout

Open `ms/ms.Rmd` and knit using [RStudio](https://www.rstudio.com/).

You can also run the following code from the R console:

```{r}
rmarkdown::render(
  input = "ms/ms.Rmd",
  output_dir = "ms"
)
```

or use `make`

```
make paper
```

### Generating all results

You can re-run all analyses, figures, etc. using [GNU make](https://www.gnu.org/software/make/). Type `make --version` in a terminal/shell to see if it is already installed.

```
# Clear out previously saved output
make cleanall
# This will take a long time to run
make
```
