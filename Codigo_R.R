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
col_data$Grupo <- gsub("Phenotype:", "", col_data$Grupo)

# Crear el rowData (SummarizedExperiment) a partir de los metabolitos que se
# encuentran en la fila 1. Como el experimento no dispone de datos para los metabolitos
# (más alla de su identificador), únicamente los numeraremos.

row_data <- datos[, 1, drop = F] # Separa los datos de la columna 1.
row_data <- as.data.frame(row_data[-1,])
row_data$numeracion <- seq(from = 1, to = dim(row_data)[1]) # Crea la columna "numeracion".

# Crear los metadatos del experimento.

metadatos <- list(ExperimentName = "Metabolomics Analysis: Opioid Addiction Project (Golestan Cohort Study)", ResearcherName = "Susan Sumner", Institution = "University of North Carolina at Chapel Hill", Department = "Nutrition", Contact = "susan_sumner@unc.edu", Species = "Homo Sapiens", AgeRange = "40-75", Gender = "Male & Female", SampleType = "Urine")
                
# Eliminar la primera fila de los datos.

datos <- datos[-1,]

# Transformación de los datos de carateres a numéricos.

datos[,-1] <- as.data.frame(lapply(datos[,-1], as.numeric))
View(datos)
str(datos)

# Crear el SummarizedExperiment.

sumexp <- SummarizedExperiment(assays = list(Opioids = datos), rowData = row_data, colData = col_data, metadata = metadatos)

# Guardar el SummarizedExperiment en formato binario (.Rda).

save(sumexp, file = "sumexp.Rda")

## ANÁLISIS EXPLORATORIO: DESCRIPTIVO ##

# Dimensiones y tipo de variables.

dim(datos)
dim(col_data)
dim(row_data)

str(datos)

# Valores NA.

sum(is.na(assay(sumexp)))

## Correlación y clustering de variables para agrupar información similar
## PCA

matriz <- assay(sumexp) # Creamos una variable con los datos del assay.
matriz_num <- matriz[, sapply(matriz, is.numeric)] # Creamos una variables con los datos numéricos del assay.
pca <- prcomp(matriz_num, scale = TRUE)  # Realizamos el analisis de componentes principales, y normalizamos las variables.
summary(pca) # Resultado pca (Desviación estandar y % de varianza sobre el total).
names(pca)
pca$rotation[,1] # Carga de cada variable al PC1
pca$rotation[,2] # Carga de cada variable al PC2
var_pc1 <- pca$sdev[1]^2/sum(pca$sdev^2) # % de variabilidad explicado por el PC1
var_pc2 <- pca$sdev[2]^2/sum(pca$sdev^2) # % de variabilidad explicado por el PC2

col_grupo <- colData(sumexp)[-1, "Grupo"]

plot(pca$x[,1:2], col = as.factor(col_grupo), main = "PCA: Proyecciones PC1 & PC2", xlab = paste("PC1 ", round(var_pc1, 2),"%"), ylab = paste("PC2 ", round(var_pc2,2), "%"))  # Visualizar las primeras dos componentes
legend("topright", legend = levels(as.factor(col_grupo)), 
       fill = rainbow(length(levels(as.factor(col_grupo)))))

## Clustering

matriz_num_scal <- scale(matriz_num) # Escalar los datos
sd_matriz <- apply(matriz_num_scal, MARGIN=1, FUN="sd") # Calcular la desviación estandar de cada fila.
sel_matriz <- (sd_matriz > quantile(sd_matriz, 0.975)) # Determinar las variables con una desviación estandar mayor al 1%.
matriz_clust <- matriz_num_scal[sel_matriz,] # Seleccionar las variables con una desviación estandar mayor al 1%.
dim(matriz_clust)

t_matriz_clust <- t(matriz_clust) # Transponer la matriz para hacer el clustering (Observaciones en filas y variables en columnas).
dim(t_matriz_clust)

dist_matriz <- dist(t_matriz_clust) # Determinar la distancia euclidea de las observaciones.
clust <- hclust(dist_matriz, method = "average") # Realizar el clustering.
plot(clust, main = "Dendograma", xlab = "Observaciones", ylab = "Distancia", labels = F) # Crear el dendograma. Sin etiquetas.
plot(clust, main = "Dendograma", xlab = "Observaciones", ylab = "Distancia") # Crear el dendograma. Con etiquetas.
