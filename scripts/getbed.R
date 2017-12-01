library(tidyverse)

df = read_tsv(snakemake@input[[1]], col_names=FALSE) %>% 
        select(chrom=X1, start=X2, end=X3, summit=X6, expression=X7, strand=X4, tfiib=X9) %>% 
        group_by(chrom, start, end, summit, expression, strand) %>% slice(1) %>% ungroup() %>% 
        mutate_at(vars(tfiib), funs(if_else(.<0, FALSE, TRUE))) %>% 
        mutate_at(vars(start), funs(if_else(strand=="+", .+summit, .))) %>% 
        mutate_at(vars(end), funs(if_else(strand=="-", start+summit, .))) %>% 
        mutate(name=".")

with_tfiib = df %>% filter(tfiib) %>%
                select(chrom,start,end,name,expression,strand) %>% 
                arrange(desc(expression)) %>% 
                write_tsv(snakemake@output[["yes"]], col_names=FALSE)
no_tfiib = df %>% filter(!tfiib) %>%
            select(chrom,start,end,name,expression,strand) %>% 
            arrange(desc(expression)) %>% 
            write_tsv(snakemake@output[["no"]], col_names=FALSE)
