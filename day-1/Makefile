all: slides purl

slides: index.Rmd slides.css
	Rscript -e 'library("rmarkdown"); render("index.Rmd")'

purl: index.Rmd
	Rscript -e "knitr::purl(\"index.Rmd\")"

