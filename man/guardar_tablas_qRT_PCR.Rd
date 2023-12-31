% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/guardar_tablas_qRT_PCR.R
\name{guardar_tablas_qRT_PCR}
\alias{guardar_tablas_qRT_PCR}
\title{Guarda en el Sistema de Archivos de Windows los Archivos Principales Generados en R Durante (y para) el Analisis de Expresion Relativa por qRT-PCR}
\usage{
guardar_tablas_qRT_PCR(ruta, tratamiento_condicion, target, normalizador)
}
\arguments{
\item{ruta}{Ruta de carpeta origen en winsys}

\item{tratamiento_condicion}{Tratamiento y condicion (por ejemplo: Tratamiento: NaCl. Condicion: 3h)}

\item{target}{Gen problema (por ejemplo cornichon CN1)}

\item{normalizador}{Gen normalizador (por ejemplo tubulina, Tb o factor de elongacion, Ef)}
}
\value{
Graficos y resumen de resultados de anova
}
\description{
Guarda los objetos de R principales, en winsys
Incluyen:
Tablas de desviacion estandar de todas las corridas y promediadas
Tablas usadas para graficas
Tablas usadas para analisis estadistico de ANOVA y Tukey
Tablas de desviaciones estandar (ds) por corrida, combinadas (por tratamiento y contol) y ds promedio
}
