
suppressMessages(library("ChemmineR"))
suppressMessages(library("ChemmineOB"))

data(smisample) # Loads the same SMIset provided by the library
smiset <- smisample

smi <- data.frame(SMILES = as.character(smiset[1:100]))

# smiles_matrix <- read.csv(args[1], header = FALSE)
smiles_matrix <- smi
print("### Lipinski Rule of Five - Druglikeness Violations ###")
sdf <- apply(smiles_matrix, 1, smiles2sdf) # 1 = rows, convert smiles to sdfs

for (i in 1:length(sdf)) {
  cid(sdf[[i]]) <- names(sdf)[i]
}

cid(sdf[[2]])

plot(sdf[[2]], regenCoords = TRUE, print = FALSE) # 'print=TRUE' returns SDF summaries
