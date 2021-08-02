
x <- load_example_dataset()

head(x$expr)
head(x$pdata)
dim(x$expr)
dim(x$pdata)

del <- c('GSM2286248','GSM2286249',
         'GSM2286316','GSM2286317',
         'GSM2286396','GSM2286397',
         'GSM2286398','GSM2286399')

TCM_expr <- x$expr[,!names(x$expr)%in%del]

pdata<- x$pdata[!x$pdata$geo_accession%in%del,]

# Handle weird words
pdata$title <- sub("渭", "μ", pdata$title)
pdata$title <- gsub("灏戼㹥", "β", pdata$title)
pdata$perturbagen <- sub("渭", "μ", pdata$perturbagen)
pdata$title <- sub("尾", "β", pdata$title)
pdata$perturbagen <- sub("尾", "β", pdata$perturbagen)

TCM_pdata <- pdata

usethis::use_data(TCM_pdata)

usethis::use_data(TCM_expr)
