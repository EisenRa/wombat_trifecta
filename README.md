# Wombat trifecta
Genome-resolved metagenomics of faecal samples from all three wombat species.


To render the bookdown document, run these lines of code:
```
library(bookdown)
library(htmlwidgets)
library(webshot)

render_book(input = "code/", output_format = "bookdown::gitbook", output_dir = "../docs")

``` 
