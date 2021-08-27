# A good simple tutorial about Make can be found at http://kbroman.org/minimal_make/ 
R_OPTS=--no-save --no-restore --no-init-file --no-site-file
all: data model paper
data: objects/pair_div.rds 
model: objects/fit.rds
paper: ms/ms.pdf figures/concepts.pdf figures/h1-raw.pdf figures/h1.pdf figures/h2.pdf

ms/ms.pdf: ms/ms.Rmd ms/stomata-independence.bib figures/concepts.pdf figures/h1-raw.pdf figures/h2-raw.pdf figures/h1.pdf figures/h2.pdf
	Rscript -e 'rmarkdown::render("ms/ms.Rmd", output_format = "bookdown::pdf_document2", output_file = "ms.pdf")'

figures/concepts.pdf: r/16_plot-concepts.R
	Rscript -e 'source("r/16_plot-concepts.R")'

figures/h1-raw.pdf: r/08_plot-h1-raw.R
	Rscript -e 'source("r/08_plot-h1-raw.R")'

figures/h2-raw.pdf: r/09_plot-h2-raw.R
	Rscript -e 'source("r/09_plot-h2-raw.R")'

figures/h1.pdf: r/11_prepare-plotting.R r/12_plot-h1-pairs.R r/14_plot-h1.R
	Rscript -e 'source("r/11_prepare-plotting.R")'
	Rscript -e 'source("r/12_plot-h1-pairs.R")'
	Rscript -e 'source("r/14_plot-h1.R")'

figures/h2.pdf: r/11_prepare-plotting.R r/13_plot-h2-pairs.R r/15_plot-h2.R
	Rscript -e 'source("r/11_prepare-plotting.R")'
	Rscript -e 'source("r/13_plot-h2-pairs.R")'
	Rscript -e 'source("r/15_plot-h2.R")'

objects/fit.rds: r/10_fit-pairs.R processed-data/trimmed-data.rds processed-data/trimmed-phylogeny.rds objects/pair_div.rds
	Rscript -e 'source("r/10_fit-pairs.R")'

processed-data/trimmed-data.rds: r/06_connect-data.R processed-data/resolved_names.rds processed-data/phy2.rds
	Rscript -e 'source("r/06_connect-data.R")'
	
processed-data/trimmed-phylogeny.rds: r/06_connect-data.R processed-data/resolved_names.rds processed-data/phy2.rds
	Rscript -e 'source("r/06_connect-data.R")'

objects/pair_div.rds: r/07_make-pairs.R	processed-data/trimmed-data.rds processed-data/trimmed-phylogeny.rds
	Rscript -e 'source("r/07_make-pairs.R")'

processed-data/resolved_names.rds: r/02_resolve-names.R processed-data/taxize_output.rds
	Rscript -e 'source("r/02_resolve-names.R")'

processed-data/phy2.rds: r/05_make-subtrees.R processed-data/phy1.rds objects/polytomies.rds
	Rscript -e 'source("r/05_make-subtrees.R")'
	
processed-data/taxize_output.rds: r/01_taxize-data.R
	Rscript -e 'source("r/01_taxize-data.R")'

processed-data/phy1.rds: r/04_prepare-subtrees.R processed-data/phy.rds
	Rscript -e 'source("r/04_prepare-subtrees.R")'
	source raxml/run-raxml.sh

objects/polytomies.rds: r/04_prepare-subtrees.R processed-data/phy.rds
	Rscript -e 'source("r/04_prepare-subtrees.R")'

processed-data/phy.rds: r/03_make-megatree.R processed-data/resolved_names.rds
	Rscript -e 'source("r/03_make-megatree.R")'

clean: 
	\rm -f *~ *.Rout */*~ */*.Rout .RData Rplots.pdf
	
cleanall: 
	\rm -f *.aux *.bbl *.blg *.log *.pdf *~ *.Rout */*~ */*.Rout figures/*.png figures/*.pdf ms/ms.pdf ms/ms.tex objects/*.rds processed-data/*.rds */*.aux */*.log 

