## Running distillR
# Raphael Eisenhofer 

#library(devtools)
#install_github("anttonalberdi/distillR")
library(tidyverse) 
library(rairtable)
library(distillR)
library(vroom)

sample_metadata <- read_tsv("data/DMB0173_metadata.tsv.gz") %>%
  rename(sample=1)

read_counts <- read_tsv("data/DMB0173_counts.tsv.gz") %>%
  rename(genome=1) %>% 
  select(one_of(c("genome",sample_metadata$sample)))

genome_metadata <- read_tsv("data/DMB0173_mag_info.tsv.gz") %>%
  rename(length=mag_size)%>%
  semi_join(., read_counts, by = "genome") %>% 
  arrange(match(genome,read_counts$genome))

##Only needs to be done once, outputs saved in 'data'
#Pull individual annotation files from AirTable,
anno_urls <- airtable("MAGs",
         base = "appWbHBNLE6iAsMRV",
         view = "viwMD6lnM4YFxkAmi") %>% #points to specific view (faster)
        read_airtable(., fields = c("ID","mag_name","number_genes","anno_url"), id_to_col = TRUE) %>% #get 3 columns from MAGs table
        filter(mag_name %in% paste0(genome_metadata$genome,".fa")) %>% #filter by MAG name
        filter(number_genes > 0) %>% #genes need to exist
        select(anno_url)

write_tsv(anno_urls, "data/anno_urls.tsv", col_names = F)

##BASH one liners to download files and wrangle (was having issues downloading using read_tsv)

``
cd data

while read file; do wget $file; done < anno_urls.tsv

cat *anno.tsv.gz > annotations.tsv.gz

gzcat annotations.tsv.gz | grep -v 'scaffold' >> annos.tsv
gzcat annotations.tsv.gz | head -1 > header.tsv
cat header.tsv annos.tsv > all_annotations.tsv
gzip all_annotations.tsv

rm header.tsv && rm annos.tsv && rm annotations.tsv.gz
``

genome_annotations <- vroom("data/all_annotations.tsv.gz")

genome_gifts_raw <- distill(genome_annotations,
                            GIFT_db,genomecol=2,
                            annotcol=c(9,10,19),
                            verbosity=F)

save(genome_gifts_raw, file = "data/genome_gifts_raw.Rdata")

# Identifiers in the annotation table: 1616
# Identifiers in the database: 1819
# Identifiers in both: 259
# Percentage of annotation table identifiers used for distillation: 16.03%
# Percentage of database identifiers used for distillation: 14.24%

