#' @title Combina datos de DDCT y de Desviacion Estandar para Analisis de qRT PCR.
#' 
#' @description Trabaja sobre objetos de R previamente creados.
#' Inicia analisis y transformacion de datos por corrida
#' Combina datos de desviacion estandar conjunta (control y tratamiento/condicion)
#' Hace grafico de barras de tratamiento vs control, con desviacion estandar conjunta, guarda resultados en winsys
#' 
#'
#' @param ruta Ruta de carpeta origen en winsys
#' @param target Gen problema (por ejemplo cornichon CN1)
#' @param normalizador Gen normalizador (por ejemplo tubulina, Tb, o factor de elongacion, Ef)
#' @return Datos combinados por las diferentes corridas, listos para graficacion y analisis 
#'
#' @export


# segunda funcion para analisis de qRT PCR
combinar_promediar_qRT_PCR <- function(ruta,
                                       target,
                                       normalizador){
  
  # guardar objeto de nombres de data frames de desviaciones
  todos_objetos <- ls(envir = parent.frame())
  
  ##################### combinar desviaciones estandar #######################
  
  # clave sd_conjunta
  sd_clave1 <- "sd_conjunta_"
  
  # filtrar por palabra clave "3h"
  filtrados1 <- todos_objetos[grep(sd_clave1, todos_objetos)]
  
  # crear listas con objetos por clave
  sd_combinados <- mget(filtrados1, envir = parent.frame())
  
  # unirlos como data frame
  sd_combinados_func <- as.data.frame(do.call(rbind, sd_combinados))
  
  # volver global
  sd_combinados <<- sd_combinados_func
  
  # sacar promedios de datos_promedios_final
  sd_combinados_promedio <<- aggregate(x = sd ~ id,
                                       data = sd_combinados_func,
                                       FUN = mean)
  
  ########################### combinar datos Cqs #############################
  
  # clave datos
  datos_clave2 = "datos_"
  
  # filtrar por palabra clave "sd_conjunta"
  filtrados2 <- todos_objetos[grep(datos_clave2, todos_objetos)]
  
  # crear listas de objetos por clave
  datos_combinados <- mget(filtrados2, envir = parent.frame())
  
  # unirlos como data frame
  datos_combinados_func <- as.data.frame(do.call(rbind, datos_combinados))
  
  # volver global
  datos_combinados <<- datos_combinados_func
  
  # obtener promedios de datos combinados
  datos_combinados_promedio <<- aggregate(x = Cq ~ Id,
                                          data = datos_combinados_func,
                                          FUN = mean)

  ######################## promediar/combinar datos DDCT ########################
  
  # clave datos
  datos_clave2 = "DDCT_"
  
  # filtrar por palabra clave "sd_conjunta"
  filtrados2 <- todos_objetos[grep(datos_clave2, todos_objetos)]
  
  ##############################################################################
  
  # corrida1
  corrida1 <- get(filtrados2[1], envir = parent.frame())
  
  # corrida2
  corrida2 <- get(filtrados2[2], envir = parent.frame())
  
  # eliminar columna de ID de cada df
  corrida_1 <- corrida1[, "exp2DDCT"]
  corrida_2 <- corrida2[, "exp2DDCT"]
    
  # obtener promedios de ambas corridas en expr2^DDCT
  DDCT_combinados <- data.frame(
    exp2DDCT_promedio = rowMeans(cbind(corrida_1, corrida_2))
    )
  
  # agregar ID
  DDCT_combinados$ID <- str_remove(corrida1$ID[1], "_corrida1") 
  
  # objeto para trabajar en la funcion
  DDCT_combinados_func <- DDCT_combinados
  
  # volver global
  DDCT_combinados <<- DDCT_combinados 
   
  ############################ agregar controles ###############################
  
  # crear data frame de controles
  controles <- data.frame(matrix(nrow = nrow(DDCT_combinados_func),
                                 ncol = ncol(DDCT_combinados_func)))
  
  # nombre de controles
  colnames(controles) <- colnames(DDCT_combinados_func)
  
  # target = "CN1"
  # normalizador = "Ea"
  # sustituir gen con Ctl en columna ID
  controles$ID <- str_replace(DDCT_combinados_func$ID, target, "Ctl")
  
  # poblar columna de datos
  controles$exp2DDCT_promedio <- 1
  
  # unir data frames
  DDCT_combinados_finales <<- rbind(DDCT_combinados_func,
                                    controles)
  
}
