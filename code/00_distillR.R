## Running distillR
# Raphael Eisenhofer 

#library(devtools)
#install_github("anttonalberdi/distillR")
library(rairtable)
library(tidyverse) 
library(distillR)
library(vroom)

##BASH one liners to download files and wrangle (was having issues downloading using read_tsv)

```
cd data

while read file; do wget $file; done < anno_urls.tsv

cat *anno.tsv.gz > annotations.tsv.gz

gzcat annotations.tsv.gz | grep -v 'scaffold' >> annos.tsv
gzcat annotations.tsv.gz | head -1 > header.tsv
cat header.tsv annos.tsv > all_annotations.tsv
gzip all_annotations.tsv

rm header.tsv && rm annos.tsv && rm annotations.tsv.gz
```


genome_annotations <- vroom("../annos/merged_annotations.tsv.gz")

genome_gifts_raw <- distill(genome_annotations,
                            GIFT_db,genomecol=2,
                            annotcol=c(9,10,19),
                            verbosity=F)

save(genome_gifts_raw, file = "data/genome_gifts_raw.Rdata")

# Identifiers in the annotation table: 460 
# Identifiers in the database: 1547 
# Identifiers in both: 28 
# Percentage of annotation table identifiers used for distillation: 6.09%
# Percentage of database identifiers used for distillation: 1.81%

