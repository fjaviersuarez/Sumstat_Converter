#!/usr/bin/env Rscript


# Unicamente para el help y sus colores
color_text <- function(text, color) {
  colors <- c(
    "rojo" = "\033[31m",
    "rosa" = "\033[95m",
    "azul" = "\033[34m",
    "reset" = "\033[0m"
  )
  
  if (color %in% names(colors)) {
    cat(paste0(colors[color], text, colors["reset"]))
  } else {
    cat(text)
  }
}

mostrar_ayuda_color <- function() {
    cat("\n")
    nombre <- "
╔═╗┬ ┬┌┬┐┌─┐┌┬┐┌─┐┌┬┐   ┌─┐┌─┐┌┐┌┬  ┬┌─┐┬─┐┌┬┐┌─┐┬─┐
╚═╗│ ││││└─┐ │ ├─┤ │    │  │ ││││└┐┌┘├┤ ├┬┘ │ ├┤ ├┬┘
╚═╝└─┘┴ ┴└─┘ ┴ ┴ ┴ ┴────└─┘└─┘┘└┘ └┘ └─┘┴└─ ┴ └─┘┴└─
"

    color_text(nombre, "rojo")

    color_text("\n\n * Propósito:\n", "azul")
    cat("Convertir datos de cualquier enfermedad de GWAS Catalog al mismo formato que los .assoc.logistic de la afeccion de testeo para PRS\n")

    color_text(" * Sintaxis:\n", "azul") 
    cat("Rscript Sumstat_Converter.r --base ENFERMEDAD.assoc.logistic --target ENFERMEDAD_GWASCATALOG.txt --out SALIDA\n")
    color_text(" * Argumentos:\n", "azul") 

    cat("--base: corresponde al fichero .assoc.logistic de la enfermedad que va a ser usado por el script como base para reescribir los archivos\n")
    cat("--target: corresponde al fichero de GWAS Catalog ya descomprimido (con gunzip)\n")
    cat("--out: corresponde a la salida SIN EXTENSIÓN, el script ya le coloca la salida necesaria a cada archivo (.log, .prob y .meta para los resultados del meta análisis, y .txt para el archivo final de salida generado con los SNP y fichero, ambos corregidos)\n\n")
    color_text("Dudas, sugerencias o errores, por favor --> fjaviersuarez@correo.ugr.es\n", "rosa") 
}



# Este comando es para coger argumentos del comando
args <- commandArgs(trailingOnly = TRUE)


# Inicializaciones
base <- NULL
out <- NULL
target <- NULL


for (i in 1:length(args)){
    # Verifica cada comando y que justo después del comando haya algo (que se le da un argumento)
    if(args[i] == "--base" && !is.null(args[i+1])){
    base <- args[i+1]
    }else if(args[i] == "--out" && !is.null(args[i+1])){
    out <- args[i+1]
    }else if(args[i] == "--target" && !is.null(args[i+1])){
    target <- args[i+1]
    }else if(args[i] == "--h" || args[i] == "--help"){
    mostrar_ayuda_color()    
    q(save = "no")

    }
}

# Comprobaciones
if(is.null(base)){
    stop("Advertencia: No se ha proporcionado el archivo base (NOA) a usar de guía")
}
if(is.null(out)){
    stop("Advertencia: No se ha proporcionado nombre para la salida")
}
if(is.null(target)){
    stop("Advertencia: No se ha proporcionado el archivo target/enfermedad a transformar/adaptar")
}





# Esta función es para seleccionar unicamente las columnas necesarias de todas las descargadas del GWAS Catalog
editar <- function(datos, out){
    # Lectura de datos desde el archivo
    # El vector de las subopciones se separa
    columnas_seleccionadas <- unlist(strsplit("hm_chrom hm_variant_id hm_pos hm_effect_allele odds_ratio standard_error p_value hm_other_allele", " "))
    if (!all(columnas_seleccionadas %in% colnames(datos))) {
        stop("ERROR: No todas las columnas a seleccionar están presentes en los datos. Revisa los datos manualmente, faltan columnas")
    }
    # Se filtran todas las filas y las columnas seleccionadas
    datos_filtrados <- datos[, columnas_seleccionadas, drop = FALSE]
    # Se actualiza

    write.table(datos_filtrados, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    cat("Columnas editadas\n")
}

# Esta función es para renombrar las columnas editadas anteriormente
renombrar <- function(datos, out){
    nuevos_nombres <- strsplit("CHR SNP BP A1 OR SE P A2", " ")[[1]]
    # Comprobaciones de columnas
    if (length(nuevos_nombres) != ncol(datos)) {
        stop(paste("ERROR: No coinciden en cantidad las columnas a introducir (", length(nuevos_nombres), " columnas) con las del archivo (", ncol(datos), " columnas), revisa posibles errores en el edit", sep = ""))
    }
    # Nueva asignación de los headers
    colnames(datos) <- nuevos_nombres

    write.table(datos, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    cat("Columnas renombradas\n")
}

# Esta función es para crear la columna SNP concatenando las demás columnas (chrN:BP:A2:A1)
SNP <- function(datos, out){ 
    columnas_necesarias <- c("CHR", "BP", "A1", "A2")
    # Comprobaciones
    if (!all(columnas_necesarias %in% colnames(datos))) {
        columnas_faltan <- columnas_necesarias[!columnas_necesarias %in% colnames(datos)]
        stop(paste("ERROR: Faltan las siguientes columnas necesarias para SNP:", paste(columnas_faltan, collapse = ", ")))} # Para que sea más legible
    # Se actualiza la columna SNP o se crea si no existe, a partir de concatenar lo necesario
    datos$SNP <- paste("chr",datos$CHR, ":", datos$BP, ":",datos$A2,":",datos$A1, sep = "") 

   
    write.table(datos, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Columna SNP creada\n")
}

# Eliminar AT/TA/CG/GC
limpiar <- function(datos, out) {

    ORmiss <- sum(is.na(datos$OR))
    SEmiss <- sum(is.na(datos$SE))
    cat("OR NA: ", ORmiss)
    cat("\nSE NA: ", SEmiss)
    datos_limpios <- datos[!is.na(datos$SE), ]
    datos_limpios <- datos_limpios[!is.na(datos$OR), ]


  write.table(datos_limpios, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nEliminados valores NA de SE/OR\n")

}

# Para generar el archivo .prob que se va a utilizar después
meta_analizar <- function(base, out){
    system(paste("plink1.9 --meta-analysis", paste0(out,".txt"),base, "--out", out))

}

# Esta función es para arreglar los SNP generados anteriormente para abordar los cambios de hebra, los alelos intercambiados, los allele mismatch, y cualquier problema
arreglar <- function(datos, datosprob, out, datosbase) {
contador <- 0
    # Le pone header al .prob para poder buscar en función de ellos
    colnames(datosprob)<-c("archivo", "SNP", "error")

    for(snp in datosprob$SNP){
        # Partes:               CHR : BP : A2 : A1
        partes_prob <- unlist(strsplit(as.character(snp), ":"))
        indice <- which((paste0("chr",datos$CHR) == partes_prob[1])&(datos$BP == partes_prob[2]))
        partes_snp <- unlist(strsplit(as.character(datos[indice,"SNP"]), ":"))
        indice_prob <- which(snp==datosprob$SNP)

        if (length(indice) > 0) {

            # Si el snp está bien escrito en todos lados y aún así está aquí + su código de error es allele_mismatch
            if((snp%in%datosbase$SNP)&&(snp%in%datos$SNP)&&(snp%in%datosprob$SNP)&&((datosprob[indice_prob,"error"]=="ALLELE_MISMATCH"))){
                # A2 se almacena en temp
                temp <- datos[indice,"A2"]
                # A2 se sobreescribe con A1
                datos[indice,"A2"] <- datos[indice,"A1"]
                # A1 se sobreescribe con temp (A2)
                datos[indice,"A1"] <- temp
                if(!is.na(datos[indice,"OR"])){
                    OR<-as.numeric(datos[indice,"OR"])
                    OR <- 1/OR
                } 

            } 
                
        
        }
            contador<-contador+1 # A cada vuelta de bucle = cada linea en prob, se suma 1
            porcentaje <- (contador/nrow(datosprob))*100

            cat("Porcentaje completado: ", sprintf("%.2f", porcentaje), "% . Líneas procesadas: ",contador," / ", nrow(datosprob),"\r")
            flush.console() # Y eso solo fuerza para que se muestre
        
        
        }
            
        write.table(datos, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

}


mover <- function(datos, out){
    # Lectura de datos desde el archivo
    # El vector de las subopciones se separa
    columnas_seleccionadas <- unlist(strsplit("SNP CHR BP A1 A2 OR P SE", " "))

    # Se filtran todas las filas y las columnas seleccionadas
    datos_filtrados <- datos[, columnas_seleccionadas, drop = FALSE]
    # Se actualiza

    write.table(datos_filtrados, file = paste0(out, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    cat("Columnas editadas\n")
}


####
# Llamada a funciones
####

# Nota: datos se carga tantas veces como se necesite después, si solo se cargase una vez se seguiría trabajando con el mismo set para cada función
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
mover(datos,out)
