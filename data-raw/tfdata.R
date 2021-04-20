library(clusterProfiler)
c5 <- read.gmt("c3.tft.v6.2.symbols.gmt")
library(stringr)

ttt <- unique(c5$term)[is.na(str_match(unique(c5$term), "UNKNOWN"))]

ttt <- as.character(ttt)

c6 <- c5[c5$term %in% ttt, ]

c6$term <- as.character(c6$term)

unique(c6$term)

for (i in 1:nrow(c6)) {
  a <- str_count(c6$term[i], "_")
  if (a == 1) {
    c6$term[i] <- str_split(c6$term[i], "_")[[1]][1]
  }
  if (a == 2) {
    c6$term[i] <- str_split(c6$term[i], "_")[[1]][2]
  }
}

names(c6) <- c("tf", "gene")

tfdata <- c6

usethis::use_data(tfdata, overwrite = TRUE)
