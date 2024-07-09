#!/usr/bin/env Rscript



# SCRIPT: Sumstat_Converter
# AUTOR: Francisco Javier Suarez Lopez (fjaviersuarez@correo.ugr.es), Universidad de Granada.
# VERSION: 1.1 (Nota: anadida opcion --nohla para eliminar el HLA de los analisis)
# ULTIMA MODIFICACION: 9/07/2024 (DD/MM/YYYY)
# DESCRIPCION: Sumstat_Converter es un programa para la conversion automatica de archivos harmonizados de NHGRI-EBI GWAS Catalogue para poder ser utilizados en PRSice-2: Polygenic Risk Score Software for Biobank-Scale Data.





# Unicamente para el help y sus colores
colorear <- function(texto, color) {
  lista_colores <- c(
    "rojo" = "\033[31m",
    "rosa" = "\033[95m",
    "azul" = "\033[34m",
    "reset" = "\033[0m"
  )
  
  if (color %in% names(lista_colores)) {
    cat(paste0(lista_colores[color], texto, lista_colores["reset"]))
  } else {
    cat(texto)
  }
}
    nombre <- "
╔═╗┬ ┬┌┬┐┌─┐┌┬┐┌─┐┌┬┐   ┌─┐┌─┐┌┐┌┬  ┬┌─┐┬─┐┌┬┐┌─┐┬─┐
╚═╗│ ││││└─┐ │ ├─┤ │    │  │ ││││└┐┌┘├┤ ├┬┘ │ ├┤ ├┬┘
╚═╝└─┘┴ ┴└─┘ ┴ ┴ ┴ ┴────└─┘└─┘┘└┘ └┘ └─┘┴└─ ┴ └─┘┴└─
"
help_color <- function() {
  cat("\n")


  colorear(nombre, "rojo")

  colorear("\n\n * Proposito:\n", "azul")
  cat("Convertir datos de cualquier enfermedad de GWAS Catalog al mismo formato que los .assoc.logistic de la afeccion de testeo para PRS\n")
  colorear(" * Sintaxis:\n", "azul") 
  cat("Rscript Sumstat_Converter.r --base ENFERMEDAD.assoc.logistic --target ENFERMEDAD_GWASCATALOG.txt --out SALIDA\n")
  colorear(" * Argumentos:\n", "azul") 

  cat("--base: corresponde al fichero .assoc.logistic de la enfermedad que va a ser usado por el script como base para reescribir los archivos\n")
  cat("--target corresponde al fichero de GWAS Catalog ya descomprimido (con gunzip)\n")
  cat("--out corresponde a la salida SIN EXTENSIoN, el script ya le coloca la salida necesaria a cada archivo (.log, .prob y .meta para los resultados del meta analisis, y .txt para el archivo final de salida generado con los SNP y fichero, ambos corregidos)\n")
  cat("--nohla es una opcion adicional que permite eliminar el HLA (chr6, bp >20M <40M), en caso de especificarse la opcion. Si no se hace, lo incluye en los analisis\n\n")

  colorear("Dudas, sugerencias o errores, por favor --> fjaviersuarez@correo.ugr.es\n", "rosa") 
}



# Este comando es para coger argumentos del comando
args <- commandArgs(trailingOnly = TRUE)


# Inicializaciones
base <- NULL
out <- NULL
target <- NULL
HLA_switch <- FALSE

for (i in 1:length(args)){
    # Verifica cada comando y que justo despues del comando haya algo (que se le da un argumento)
    if(args[i] == "--base" && !is.null(args[i+1])){
    base <- args[i+1]
    }else if(args[i] == "--out" && !is.null(args[i+1])){
    out <- args[i+1]
    }else if(args[i] == "--target" && !is.null(args[i+1])){
    target <- args[i+1]
    }else if(args[i] == "--h" || args[i] == "--help"){
    help_color()    
    q(save = "no")
    }else if(args[i] == "--nohla"){
    HLA_switch <- TRUE
    }

}

# Comprobaciones
if(is.null(base)){
    stop("! Advertencia: No se ha proporcionado el archivo base a usar de guia")
}
if(is.null(out)){
    stop("! Advertencia: No se ha proporcionado nombre para la salida")
}
if(is.null(target)){
    stop("! Advertencia: No se ha proporcionado el archivo target/enfermedad a transformar/adaptar")
}





# Esta funcion es para seleccionar unicamente las columnas necesarias de todas las descargadas del GWAS Catalog
editar <- function(datos, out){
    # Lectura de datos desde el archivo
    columnas_seleccionadas <- unlist(strsplit("hm_chrom hm_variant_id hm_pos hm_effect_allele hm_odds_ratio standard_error p_value hm_other_allele", " ")) # Esto sirve para separar el vector con las columnas que se proporcionan y, ademas, es lo que establece el orden. Si quieres usar otras columnas u otro orden cambialo aqui (ten en cuenta que debes tambien ajustar las funciones siguientes de renombre)
    if (!all(columnas_seleccionadas %in% colnames(datos))) { # Si no estan todas las columnas necesarias devuelve error 
        stop("!! ERROR: No todas las columnas a seleccionar estan presentes en los datos. Revisa los datos manualmente, faltan columnas")
    }
    # Se filtran todas las filas y las columnas seleccionadas
    datos_filtrados <- datos[, columnas_seleccionadas, drop = FALSE] # el formato es set_de_datos[filas,columnas,opciones], si en filas lo dejamos vacio las coge todas, y si en columnas ponemos el vector que tenemos solamente va a guardar las columnas que coincidan con este vector
    # Se actualiza
    write.table(datos_filtrados, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    cat("\n\n\n► Columnas editadas\n")
}

# Esta funcion es para renombrar las columnas editadas anteriormente
renombrar <- function(datos, out){
    nuevos_nombres <- strsplit("CHR SNP BP A1 OR SE P A2", " ")[[1]]
    # Comprobaciones de columnas
    if (length(nuevos_nombres) != ncol(datos)) {
        stop(paste("!! ERROR: No coinciden en cantidad las columnas a introducir (", length(nuevos_nombres), " columnas) con las del archivo (", ncol(datos), " columnas), revisa posibles errores en el edit", sep = ""))
    }
    # Nueva asignacion de los headers
    colnames(datos) <- nuevos_nombres

    write.table(datos, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    cat("► Columnas renombradas\n")
}

# Esta funcion es para crear la columna SNP concatenando las demas columnas (chrN:BP:A2:A1)
SNP <- function(datos, out){ 
    columnas_necesarias <- c("CHR", "BP", "A1", "A2")
    # Comprobaciones
    if (!all(columnas_necesarias %in% colnames(datos))) { 
        columnas_faltan <- columnas_necesarias[!columnas_necesarias %in% colnames(datos)]
        stop(paste("!! ERROR: Faltan las siguientes columnas necesarias para SNP:", paste(columnas_faltan, collapse = ", ")))} # Para que sea mas legible
    # Se actualiza la columna SNP o se crea si no existe, a partir de concatenar lo necesario
    datos$SNP <- paste("chr",datos$CHR, ":", datos$BP, ":",datos$A2,":",datos$A1, sep = "") 

   
    write.table(datos, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("► Columna SNP creada\n")
}


# Eliminar valores NA
limpiar <- function(datos, out) {
    datos_sin_na <- na.omit(datos)
  write.table(datos_sin_na, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}



meta_analizar <- function(base, out){
    system(paste("plink1.9 --meta-analysis", paste0(out,".txt"),base, "--out", out)) # System permite que el programa "escriba" por pantalla y si le damos la misma sintaxis que para plink1.9 o cualquier otro script, va a funcionar igual. Vease EasyPRS.r que funciona igual para otro ejemplo
    cat("► Meta-analisis finalizado\n")
}


arreglar <- function(datos, datosprob, out, datosbase) {
  # Como datosprob no tiene header se le anade para poder filtrar
  colnames(datosprob) <- c("archivo", "SNP", "error")
  
  # Se filtran los datosprob que tengan el error ALLELE_MISMATCH, descartando asi los BAD_ES/BAD_SE
  datosprob_filtrado <- datosprob[datosprob$error == "ALLELE_MISMATCH", ]
  
  # Se extrae el CHR y BP de cada SNP en datosprob_filtrado
  partes_prob <- strsplit(as.character(datosprob_filtrado$SNP), ":") # strsplit separa un string en base a un delimitador (: porque los SNP tienen ese formato, CHR:BP:A1:A2, si es _ se sustituye a este nivel de codigo)
  chr_prob <- sapply(partes_prob, "[", 1) # sapply sirve para aplicar una funcion a todos los elementos de un vector y devuelve otro vector, el corchete y el 1 indican que acceda al primer elemento, o sea, al cromosoma para guardarlo en el vector chr_prob
  bp_prob <- sapply(partes_prob, "[", 2)
  
  # Se crea el indice para cada coincidencia
  indice_datos <- paste0("chr", datos$CHR) %in% chr_prob & datos$BP %in% bp_prob # indice_datos es booleano! sera TRUE si el cromosoma y el BP de datos se encuentra en datosprob: asi aseguramos que ese SNP en concreto se encuentre alli. No lo hacemos buscando directamente por SNP porque el formato puede estar intercambado
  
  # Se intercambian A1 y A2 en esos indices
  temp <- datos$A2[indice_datos] # Los datos que cumplan la condicion de arriba (TRUE) van a reemplazarse guardando el A2 en una variable temporal para no perderlo (si cambiamos A1 --> A2 sin haberlo guardado habremos perdido esa informacion y ahora tendremos dos columnas A1)
  datos$A2[indice_datos] <- datos$A1[indice_datos]
  datos$A1[indice_datos] <- temp
  
  # Se invierte la OR. Se convierte a numerico porque en el set de datos no lo es y devuelve error al invertirla (1/string = error). Otra posible solucion era anadir la columna odds_ratio no harmonizada porque es la inversa a la hm_odds_ratio, y que en las coincidencias la OR se reemplazase por la odds_ratio no harmonizada
  indice_or <- indice_datos & !is.na(datos$OR)
  datos$OR[indice_or] <- 1/as.numeric(datos$OR[indice_or])
  
  write.table(datos, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("► Fichero corregido\n")

}


noHLA <- function(datos, out){
  partes_snp <- strsplit(as.character(datos$SNP), ":") 
  chr_snp <- sapply(partes_snp, "[", 1)
  bp_snp <- sapply(partes_snp, "[", 2)

  bp_snp <- as.numeric(bp_snp)


  indice_datos <-  (chr_snp != "chr6") | (chr_snp == "chr6" & (bp_snp < 20000000 | bp_snp > 40000000)) # FIltro para los SNPs que o bien no sean el chr6, o bien lo sean pero esten fuera del HLA

  datos_sin_HLA <- datos[indice_datos,]

  write.table(datos_sin_HLA, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  cat("► HLA eliminado\n")


}


mover <- function(datos, out){
    # El vector de las subopciones se separa
    columnas_seleccionadas <- unlist(strsplit("SNP CHR BP A1 A2 OR P SE", " "))

    # Se filtran todas las filas y las columnas seleccionadas
    datos_filtrados <- datos[, columnas_seleccionadas, drop = FALSE]
    # Se actualiza

    write.table(datos_filtrados, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    cat("► Columnas reordenadas\n")


    # Esta linea facilita encontrar el fichero generado entre todos los que se generan y evita obligar a hacer ls ni escribir para encontrarlo: se copia y se pega
    cat("► Programa finalizado. Se ha generado el siguiente fichero .txt: ",paste0(out,".txt"),"\n\n")

}





# Llamada a funciones


# Nota: datos se carga tantas veces como se necesite despues, si solo se cargase una vez se seguiria trabajando con el mismo set para cada funcion
colorear(nombre, "rojo")


datos <- read.table(target, header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

editar(datos, out)

datos <- read.table(paste0(out,".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

renombrar(datos,out)

datos <- read.table(paste0(out,".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

SNP(datos, out)

datos <- read.table(paste0(out,".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

limpiar(datos,out)

datos <- read.table(paste0(out,".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

meta_analizar(base, out)

datosprob <- read.table(paste0(out, ".prob"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

datosbase <- read.table(base, header = TRUE, sep = "", stringsAsFactors = FALSE, colClasses = "character", check.names=FALSE)

datos <- read.table(paste0(out,".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

datosprob <- read.table(paste0(out, ".prob"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

arreglar(datos, datosprob, out, datosbase)

datos <- read.table(paste0(out,".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

if(HLA_switch) noHLA(datos,out)

datos <- read.table(paste0(out,".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, colClasses = "character")

mover(datos,out)
