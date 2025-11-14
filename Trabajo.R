# PaulaPortilla_Trabajo2.R
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento

# 1. Cargar librerías (si necesarias) y datos del archivo "datos_biomed.csv". (0.5 pts)

# ggplot2 es una librería para crear gráficos y los construye por capas
# reshape2 es una librería diseñada para cambiar la estructura o forma de los datos
if (!require(ggplot2)) { install.packages("ggplot2"); library(ggplot2) } else { library(ggplot2) }
if (!require(reshape2)) { install.packages("reshape2"); library(reshape2) } else { library(reshape2) }

datos <- read.csv("datos_biomed.csv", sep = ";", header = TRUE)
# El archivo lo separé por ; y esa es la razón para el uso de sep= ";"

head(datos) # Lo usamos para imprimir las primeras filas y comprobar que se leyó correctamente

# 2. Exploración inicial con las funciones head(), summary(), dim() y str(). ¿Cuántas variables hay? ¿Cuántos tratamientos? (0.5 pts)

str(datos)       # Muestra los tipos de datos de cada columna: numérico, factor, etc
summary(datos)   # Muestra estadísticas básicas: mínimo, media, mediana, máximo, cuartiles
dim(datos)       # Muestra dimensiones del dataset: número de filas y columnas

num_variables <- ncol(datos)
cat("Número total de variables:", num_variables, "\n")
# Usamos este código para contar cuántas variables hay

num_tratamientos <- length(unique(datos$Tratamiento))
cat("Número de tratamientos diferentes:", num_tratamientos, "\n")
# Usamos este código para contar cuántos tratamientos hay

# 3. Una gráfica que incluya todos los boxplots por tratamiento. (1 pt)
# Un boxplot muestra la distribución de una variable numérica, incluyendo mediana, cuartiles y valores atípicos

# Primero convertimos el tratamiento en factor para que R lo interprete como categorías
datos$Tratamiento <- as.factor(datos$Tratamiento)

# Convertimos las columnas en numéricas
cols_num <- c("Glucosa", "Presion", "Colesterol")
for(col in cols_num){
  if(col %in% names(datos)){
    # sustituir comas decimales por punto y convertir
    datos[[col]] <- as.numeric(gsub(",", ".", datos[[col]]))
  } else {
    stop(paste0("No se encontró la columna '", col, "'. Revisa nombres en el CSV."))
  }
}

#Pasamos a formato largo con melt (reshape2)
datos_largo <- reshape2::melt(datos,
                              id.vars = "Tratamiento",
                              measure.vars = cols_num,
                              variable.name = "Parametro",
                              value.name = "Valor")

# Para el boxplot único, usamos facet_wrap para mostrar en un solo panel dividido por parámetro
p <- ggplot(datos_largo, aes(x = Tratamiento, y = Valor, fill = Tratamiento)) +
  geom_boxplot() +
  facet_wrap(~Parametro, scales = "free_y") +
  labs(title = "Boxplots de Glucosa, Presión y Colesterol por Tratamiento",
       x = "Tratamiento", y = "Valor") +
  theme_minimal()
ggsave("03_boxplots_unico_por_tratamiento.png", p, width = 10, height = 5, dpi = 300)
print(p)

# 4. Realiza un violin plot (investiga qué es). (1 pt)
# Un violin plot muestra la densidad de la distribución además del IQR/mediana
# Combina la información de densidad con la mediana/dispersión

p_violin <- ggplot(datos_largo, aes(x = Tratamiento, y = Valor, fill = Tratamiento)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ Parametro, scales = "free_y") +
  labs(title = "Violin plots por tratamiento",
       subtitle = "Muestra la densidad y la distribución de los valores") +
  theme_minimal()
print(p_violin)

# 5. Realiza un gráfico de dispersión "Glucosa vs Presión". Emplea legend() para incluir una leyenda en la parte inferior derecha. (1 pt)

p_dispersion <- ggplot(datos, aes(x = Presion, y = Glucosa, color = Tratamiento)) +
  geom_point(alpha = 0.8, size = 2) + 
  labs(title = "Gráfico de Dispersión: Glucosa vs Presión",
       x = "Presión",
       y = "Glucosa") +
  theme_minimal() +
  # Ajustamos la posición de la leyenda a la esquina inferior derecha
  theme(legend.position = c(0.95, 0.05), 
        legend.justification = c("right", "bottom")) 
print(p_dispersion)

# 6. Realiza un facet Grid (investiga qué es): Colesterol vs Presión por tratamiento. (1 pt)
#La función facet grid en ggplot2 forma una matriz de paneles definidos por variables de faceting de filas y columnas

p_facet <- ggplot(datos, aes(x = Presion, y = Colesterol)) +
  geom_point(aes(color = Tratamiento)) +
  facet_wrap(~ Tratamiento, scales = "free") +
  labs(title = "Colesterol vs Presión por Tratamiento",
       x = "Presión", y = "Colesterol") +
  theme_minimal()
print(p_facet)

# 7. Realiza un histogramas para cada variable. (0.5 pts)

p_hist <- ggplot(datos_largo, aes(x = Valor, fill = Parametro)) +
    geom_histogram(bins = 15, color = "black") +
# 'scales = "free"' permite que los ejes X e Y se ajusten a cada gráfico
  facet_wrap(~ Parametro, scales = "free") +
  labs(title = "Histogramas de Parámetros Biomédicos",
       x = "Valor",
       y = "Frecuencia",
       fill = "Parámetro") +  # Esto le da un título claro a la leyenda
  theme_minimal()
print(p_hist)

# 8. Crea un factor a partir del tratamiento. Investifa factor(). (1 pt)
# Factor() convierte caracteres en categorías y es útil para análisis por grupo.

trat_factor <- factor(datos$Tratamiento)
cat("Clase de trat_factor:", class(trat_factor), "\n")
cat("Niveles del factor:\n"); print(levels(trat_factor))
cat("Tabla de conteo por nivel:\n"); print(table(trat_factor))

# 9. Obtén la media y desviación estándar de los niveles de glucosa por tratamiento. Emplea aggregate() o apply(). (0.5 pts)

# Media de la glucosa por tratamiento
media_glucosa <- aggregate(Glucosa ~ Tratamiento, data = datos, FUN = mean, na.action = na.omit)
print(media_glucosa)

#Desviación estándar de la glucosa por tratamiento
sd_glucosa <- aggregate(Glucosa ~ Tratamiento, data = datos, FUN = sd, na.action = na.omit)
print(sd_glucosa)

# 10. Extrae los datos para cada tratamiento y almacenalos en una variable. Ejemplo todos los datos de Placebo en una variable llamada placebo. (1 pt)

# Extraemos los datos de 'Placebo', 'FármacoA, 'FármacoB'
datos_placebo <- datos[datos$Tratamiento == "Placebo", ]
datos_farmacoA <- datos[datos$Tratamiento == "FarmacoA", ]
datos_farmacoB <- datos[datos$Tratamiento == "FarmacoB", ]

# Verificamos la creación de las variables e imprimimos las primeras filas y el conteo
cat("\n--- Resumen de datos 'Placebo' ---\n")
print(head(datos_placebo))
cat(paste("Total filas Placebo:", nrow(datos_placebo), "\n"))

cat("\n--- Resumen de datos 'FarmacoA' ---\n")
print(head(datos_farmacoA))
cat(paste("Total filas FarmacoA:", nrow(datos_farmacoA), "\n"))

cat("\n--- Resumen de datos 'FarmacoB' ---\n")
print(head(datos_farmacoB))
cat(paste("Total filas FarmacoB:", nrow(datos_farmacoB), "\n"))

# 11. Evalúa si los datos siguen una distribución normal y realiza una comparativa de medias acorde. (1 pt)

# Análisis para glucosa
cat("\n--- Medias y Desviación Estándar por Tratamiento ---\n")
desc_glucosa <- aggregate(Glucosa ~ Tratamiento, data = datos, FUN = function(x) c(Media = mean(x), DE = sd(x)))
print(desc_glucosa)
shapiro_glucosa <- sapply(split(datos$Glucosa, datos$Tratamiento), function(x) shapiro.test(x)$p.value)
bartlett_glucosa <- bartlett.test(Glucosa ~ Tratamiento, data = datos)$p.value
cat("\n--- Resumen de Supuestos ---\n")
cat("Normalidad (p-values por grupo):", paste(round(shapiro_glucosa, 4), collapse = ", "), "\n")
cat("Homocedasticidad (Bartlett p-value):", round(bartlett_glucosa, 4), "\n\n")
if (all(shapiro_glucosa >= 0.05) && bartlett_glucosa >= 0.05) {
  cat("--> TEST RECOMENDADO: ANOVA (Paramétrico)\n")
  resultado_glucosa <- aov(Glucosa ~ Tratamiento, data = datos)
  print(summary(resultado_glucosa))
} else if (all(shapiro_glucosa >= 0.05) && bartlett_glucosa < 0.05) {
  cat("--> TEST RECOMENDADO: ANOVA de Welch\n")
  resultado_glucosa <- oneway.test(Glucosa ~ Tratamiento, data = datos, var.equal = FALSE)
  print(resultado_glucosa)
} else {
  cat("--> TEST RECOMENDADO: Kruskal-Wallis (No Paramétrico)\n")
  resultado_glucosa <- kruskal.test(Glucosa ~ Tratamiento, data = datos)
  print(resultado_glucosa)
}

# Análisis para presión 

cat("\n--- Medias y Desviación Estándar por Tratamiento ---\n")
desc_presion <- aggregate(Presion ~ Tratamiento, data = datos, FUN = function(x) c(Media = mean(x), DE = sd(x)))
print(desc_presion)
shapiro_presion <- sapply(split(datos$Presion, datos$Tratamiento), function(x) shapiro.test(x)$p.value)
bartlett_presion <- bartlett.test(Presion ~ Tratamiento, data = datos)$p.value
cat("\n--- Resumen de Supuestos ---\n")
cat("Normalidad (p-values por grupo):", paste(round(shapiro_presion, 4), collapse = ", "), "\n")
cat("Homocedasticidad (Bartlett p-value):", round(bartlett_presion, 4), "\n\n")
if (all(shapiro_presion >= 0.05) && bartlett_presion >= 0.05) {
  cat("--> TEST RECOMENDADO: ANOVA (Paramétrico)\n")
  resultado_presion <- aov(Presion ~ Tratamiento, data = datos)
  print(summary(resultado_presion))
} else if (all(shapiro_presion >= 0.05) && bartlett_presion < 0.05) {
  cat("--> TEST RECOMENDADO: ANOVA de Welch\n")
  resultado_presion <- oneway.test(Presion ~ Tratamiento, data = datos, var.equal = FALSE)
  print(resultado_presion)
} else {
  cat("--> TEST RECOMENDADO: Kruskal-Wallis (No Paramétrico)\n")
  resultado_presion <- kruskal.test(Presion ~ Tratamiento, data = datos)
  print(resultado_presion)
}

# Análisis para colesterol

cat("\n--- Medias y Desviación Estándar por Tratamiento ---\n")
desc_colesterol <- aggregate(Colesterol ~ Tratamiento, data = datos, FUN = function(x) c(Media = mean(x), DE = sd(x)))
print(desc_colesterol)
shapiro_colesterol <- sapply(split(datos$Colesterol, datos$Tratamiento), function(x) shapiro.test(x)$p.value)
bartlett_colesterol <- bartlett.test(Colesterol ~ Tratamiento, data = datos)$p.value
cat("\n--- Resumen de Supuestos ---\n")
cat("Normalidad (p-values por grupo):", paste(round(shapiro_colesterol, 4), collapse = ", "), "\n")
cat("Homocedasticidad (Bartlett p-value):", round(bartlett_colesterol, 4), "\n\n")
if (all(shapiro_colesterol >= 0.05) && bartlett_colesterol >= 0.05) {
  cat("--> TEST RECOMENDADO: ANOVA (Paramétrico)\n")
  resultado_colesterol <- aov(Colesterol ~ Tratamiento, data = datos)
  print(summary(resultado_colesterol))
} else if (all(shapiro_colesterol >= 0.05) && bartlett_colesterol < 0.05) {
  cat("--> TEST RECOMENDADO: ANOVA de Welch\n")
  resultado_colesterol <- oneway.test(Colesterol ~ Tratamiento, data = datos, var.equal = FALSE)
  print(resultado_colesterol)
} else {
  cat("--> TEST RECOMENDADO: Kruskal-Wallis (No Paramétrico)\n")
  resultado_colesterol <- kruskal.test(Colesterol ~ Tratamiento, data = datos)
  print(resultado_colesterol)
}

# 12. Realiza un ANOVA sobre la glucosa para cada tratamiento. (1 pt)

# Primero ejecutamos el ANOVA
anova_glucosa_final <- aov(Glucosa ~ Tratamiento, data = datos)
cat("\n--- 1. Resultado del ANOVA ---\n")
print(summary(anova_glucosa_final))

# Realizamos el test Post-Hoc que compara pares de tratamientos si el ANOVA fue significativo
cat("\n--- 2. Test Post-Hoc Tukey HSD (Comparación de pares) ---\n")
tukey_resultado <- TukeyHSD(anova_glucosa_final)
print(tukey_resultado)


