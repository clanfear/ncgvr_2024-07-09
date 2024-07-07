library(tidyverse)



read_cdc <- function(input_file){
  input_text <- readLines(input_file)
  df <- read_tsv(input_file, n_max =  str_which(input_text, "Total")[1]-2)
  df <- df |>
    mutate(type = str_to_sentence(str_extract(input_file, "(all|firearm)"))) |>
    select(year = Year, deaths = Deaths, population = Population, rate = `Crude Rate`, type)
  return(df)
}

cdc_mortality <- map(list.files("./data/raw/cdc", full.names = TRUE), read_cdc) |>
  list_rbind() |>
  distinct()
save(cdc_mortality, file = "./data/cdc_mortality.RData")
