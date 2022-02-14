all: slides purl copy

slides: index.Rmd slides.css
	Rscript -e 'library("rmarkdown"); render("index.Rmd")'

purl: index.Rmd
	Rscript -e "knitr::purl(\"index.Rmd\")"

copy: index.html slides.css macros.js
	cp -R -u index_files index.html macros.js slides.css libs resources ~/work/web/jekyll/blog/teaching/courses/physalia/2022/gam/day-1/
