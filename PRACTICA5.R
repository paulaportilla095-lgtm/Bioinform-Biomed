#############################################################################
#
# PRACTICA R
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")
a
# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) #filas y columnas conteo
head(data) #seis primeras filas y columnas
tail(data) #seis ultimas filas y colunmas

# Hacemos un primer histograma para explorar los datos
hist(data, col = "blue", main="GSE5583 - Histogram")

# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
#La transformación logarítmica nos permite organizar los datos de otra manera, tanto de los datos como de lo que estamos visualizando
#Sirve para mejorar la visualización de los datos
data2 = log2(data)

hist(data2, col = "purple", main="GSE5583 (log2) - Histogram")


# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
#Data2 es el logaritmo visual que hemos usado, col=c es para aplicar los colores a cada columna, main es para agregar el nombre del título 
# y las=2 sirve para poner los nombres de cada caja en forma vertical
# ¿Qué es un boxplot?
#Un boxplot es una representación visual de una distribución de datos
#Un boxplot cuenta de cajas que presenta cada dato, las líneas dentro de las cajas representan la mediana
#Los bigotes muestran la dispersion y los puntos que están fuera de los bigotes y de las cajas, representan valores atípicos
boxplot(data2, col=c("blue", "blue", "blue",
	"orange", "orange", "orange"),
	main="GSE5583 - boxplots", las=2)


boxplot(data, col=c("blue", "blue", "blue",
	"orange", "orange", "orange"),
	main="GSE5583 - boxplots", las=2)

	

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación ç
# de los valores de expresión. ¿Es correcta la separación?
# Si es correcta la separación porque se pueden ver las muestras agrupadas correctamente segun su tipo biológico: WT y KO
hc = hclust(as.dist(1-cor(data2)))
plot(hc, main="GSE5583 - Hierarchical Clustering")


#######################################
# Análisis de Expresión Diferencial 
#######################################

head(data)

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
# Hemos usado los datos brutos del wildtype y del knockout, el matrix significa que hemos usado una tabla que está
# formada por filas y columnas y el array son los vectores usados para la tabla
wt <- data[,1:3]
ko <- data[,4:6]
class(wt)

# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)
head(wt.mean) 
head(ko.mean)

# ¿Cuál es la media más alta?
# La media mas alta es 37460.5
limit = max(wt.mean, ko.mean)
limit

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO",
	main = "GSE5583 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Añadir una línea diagonal con abline
abline(0, 1, col = "red")

# ¿Eres capaz de añadirle un grid?
# Si
grid()
#abline(a, b): línea de pendiente b y ordenada en el origen a
#abline(h=y): línea horizontal
#abline(v=x): línea vertical
abline(1, 2, col = "red")      #línea y = 2x + 1
abline(h = 2, col = "green")   #línea y = 2
abline(v = 3, col = "violet")  #línea x = 3

# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean

# Hacemos un histograma de las diferencias de medias
hist(diff.mean, col = "gray")

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# Porque los datos que tienen una distribución normal no necesitan transformación, dado a que el t-test compara directamente las medias
# ¿Cuántas valores tiene cada muestra?
# Cada muestra tiene 3 valores
pvalue = NULL 
tstat = NULL 
for(i in 1 : nrow(data)) { #Para cada gen
	x = wt[i,] # gene wt número i
	y = ko[i,] # gene ko número i
	
	# Hacemos el test
	t = t.test(x, y)
	
	# Añadimos el p-value a la lista
	pvalue[i] = t$p.value
	# Añadimos las estadísticas a la lista
	tstat[i] = t$statistic
}

head(pvalue)

# Ahora comprobamos que hemos hecho TODOS los cálculos
length(pvalue)

# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10?
# Aplicar el -log10 es una forma de transformar los datos
hist(pvalue,col="gray")
hist(-log10(pvalue), col = "gray")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
# Si se puede representar en el gráfico
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline(v = diff.mean_cutoff, col = "blue", lwd = 3)
#abline(v = -diff.mean_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
# Hay 426 genes que cumplen los criterios de diferencia de medias ≥ 2 y p-value ≤ 0.01
# Y cada uno de ellos tiene 6 valores de expresión: 3 de WT y 3 de KO
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]),col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
# Porque en el código los genes con valores negativos de diferencia de medias se pintarán de rojo y los positivos en azul
# Pero en la práctica los genes sobreexpresados suelen estar en rojo y los reprimidos en azul
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0],
	-log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0],
	-log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")


# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# Filtered sirve para mostrar la matriz de datos que queremos ver, en este caso, genes y muestras
# Rowv = rowv sirve para definir el orden de las filas
# Colv = colv sirve para definir el orden de las columnas
# cexCol= 0.7 es el tamaño de las columnas
# labRow=FALSE sirve para no mostrar el nombre de las filas
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE, col=hcl.colors(100, "Viridis"))

heatmap(filtered)


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

# Hacemos nuestro heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row")

# Lo guardamos en un archivo PDF
pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row",labRow=FALSE)
dev.off()
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,col = redgreen(75), scale = "row",labRow=FALSE)

# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t",quote = FALSE)
