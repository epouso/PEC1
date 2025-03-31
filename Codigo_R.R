#### PEC1: CÓDIGO R UTILIZADO ####

# Paquetes:
library(SummarizedExperiment)

## IMPORTACIÓN DE DATOS Y CREACIÓN DE SUMMARIZEDEXPERIMENT ##

# Importar la tabla de datos. El archivo de origen es en formato .txt.

datos <- read.table("ST001618_AN002653_Results.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Eliminación de todos los puntos (.) en la tabla. Los valores de intensidad contienen (.) 
# de forma aleatoria. En un incio indicaban la separación de miles, pero en el documento .txt
# habian perdido el sentido.

datos[-1] <- lapply(datos[-1], function(x) gsub("\\.", "", x))

# Crear el colData (SummarizedExperiment) a partir de los grupos a los que
# pertenece cada muestra (Fila 1 del df datos).

col_data <- datos[1, , drop = F] # Separar los datos de la fila 1 en un nuevo df.
col_data <- as.data.frame(t(col_data)) # Transponer el df.
colnames(col_data) <- c("Grupo") # Modificar el header

col_data <- col_data[-1, , drop = F] # Eliminar la primera fila. Eliminar????

col_data$Grupo <- gsub("Phenotype:", "", col_data$Grupo)

# Crear el rowData (SummarizedExperiment) a partir de los metabolitos que se
# encuentran en la fila 1. Como el experimento no dispone de datos para los metabolitos
# (más alla de su identificador), únicamente los numeraremos.

row_data <- datos[, 1, drop = F] # Separa los datos de la columna 1.
row_data$numeracion <- seq(from = 1, to = dim(row_data)[1]) # Crea la columna "numeracion".

# Crear los metadatos del experimento.

metadatos <- list(ExperimentName = "Metabolomics Analysis: Opioid Addiction Project (Golestan Cohort Study)", ResearcherName = "Susan Sumner", Institution = "University of North Carolina at Chapel Hill", Department = "Nutrition", Contact = "susan_sumner@unc.edu", Species = "Homo Sapiens", AgeRange = "40-75", Gender = "Male & Female", SampleType = "Urine")
                
# Crear el SummarizedExperiment.

sumexp <- SummarizedExperiment(assays = list(Opioids = datos), rowData = row_data, colData = col_data, metadata = metadatos)

# Guardar el SummarizedExperiment en formato binario (.Rda).

save(sumexp, file = "sumexp.Rda")

## ANÁLISIS EXPLORATORIO: DESCRIPTIVO ##

# Dimensiones.

dim(datos)
dim(col_data)
dim(row_data)

# Valores NA.

sum(is.na(assay(sumexp)))


