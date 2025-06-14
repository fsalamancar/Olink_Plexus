#Limpieza de datos del PLEXUS, Preestadisticos y Pre-visuales####

#Carga de paquetes####
packages <- function(requirements,quiet=FALSE){
  has   <- requirements %in% rownames(installed.packages())
  if(any(!has)){
    message("Installing packages...")
    setRepositories(ind=c(1:7))
    r <- getOption("repos")
    r["CRAN"] <- "https://cran.uni-muenster.de/"
    #options(install.packages.check.source = "no")
    install.packages(requirements[!has],repos=r)
  }
  if(quiet){
    for(r in requirements){suppressMessages(require(r,character.only=TRUE))}
  }else for(r in requirements){message(paste(r,suppressMessages(require(r,character.only=TRUE)),sep=': '))}
}
packages(c('tidyverse',"tidyr","readxl","readr","tibble","forcats","purrr","stringr","readxl",'data.table',"ggrepel","dplyr","openxlsx","fs","ggplot2","snpStats","data.table","qqman","CMplot"))




#Definir funciones utiles####
crear_carpeta <- function(nombre_carpeta){
  # Crear la carpeta si no existe
  if (!dir.exists(nombre_carpeta)) {
    dir.create(nombre_carpeta)
    cat("La carpeta se ha creado correctamente en:", nombre_carpeta, "\n")
  } else {
    cat("La carpeta ya existe en:", nombre_carpeta, "\n")
  }
}


#Carga y limpieza de datos####

##Tabla de Covariables####
plexus_covariables_Raw <- read_excel(file.path("1.Datos_Raw", "Metadata_plexus.xlsx")) %>%as_tibble()

##Tabla para comparacion entre UCvsCD####

#Crear una subcarpeta para los datos
crear_carpeta("2.Datos_Limpios_UCvsCD")

#Carga de tabla expresiones del Plexus 
plexus_raw <- read.csv(file.path("1.Datos_Raw", "Plexus_rawdata.csv"), sep = ";") %>%as_tibble()
plexus_raw_inflammation <-  plexus_raw %>% filter(str_detect(Panel,"Inflammation") & QC_Warning == "PASS") ##QUE PANEL ES EL CORRECTO
plexus_raw_inflammation <- plexus_raw_inflammation %>% dplyr::select(SampleID,Index,OlinkID,Assay,NPX) %>% rename(ID=SampleID)

#Carga de tabla de metadatos del Plexus
plexus_metadata <- read.delim(file.path("1.Datos_Raw","Plexus_metadata.tsv"))%>% dplyr::rename(ID=RAW_DATA_FILE_NAME)
plexus_metadata[,1] <- str_replace(plexus_metadata[,1], "\\.0$", "") 

#Join para dejar solo una tabla, ordenar y remover "IBD Unclassified"
Plexus_InflammatoryPanel <- plexus_raw_inflammation %>% left_join(plexus_metadata, by = "ID") %>% drop_na()
Plexus_InflammatoryPanel <- Plexus_InflammatoryPanel %>% dplyr::select(ID,Assay,NPX,SEX,DIAGNOSIS,DISEASE_LOCATION,ETHNICITY) %>% pivot_wider(names_from = Assay, values_from = NPX) %>% drop_na()
Plexus_InflammatoryPanel <- Plexus_InflammatoryPanel %>%  filter(str_detect(DIAGNOSIS,"Crohn|Ulcerative|Ileal|Colonic"))

#Escribir la tabla para la comparacion UCvsCD
write.csv(Plexus_InflammatoryPanel, file.path("2.Datos_Limpios_UCvsCD", "Plexus_InflammatoryPanel.csv"))

#Transformar la tabla Plexus_InflammatoryPanel a formato largo 
Plexus_InflammatoryPanel_Long <- Plexus_InflammatoryPanel%>%pivot_longer(cols = -c(1:5), names_to = "PROTEIN", values_to = "VALUE") %>% mutate(VALUE=as.numeric(str_replace(VALUE,",","."))) %>% drop_na()%>% select(ID,SEX,DIAGNOSIS,ETHNICITY,PROTEIN,VALUE)

#Guardar la tabla larga con texto
write.csv(Plexus_InflammatoryPanel_Long, file.path("2.Datos_Limpios_UCvsCD", "Plexus_InflammatoryPanel_Long.csv"))

#Convertir variables de la columna SEX a MALE=0, FEMALE=1
Plexus_InflammatoryPanel_Long_Binary <- Plexus_InflammatoryPanel_Long
Plexus_InflammatoryPanel_Long_Binary$SEX <-  ifelse(Plexus_InflammatoryPanel_Long_Binary$SEX == "Female", 1, 0)

#Convertir variables de la columna DIAGNOSIS A Crohn Disease =0 Ulcerative Colitis=1
Plexus_InflammatoryPanel_Long_Binary$DIAGNOSIS <- ifelse(Plexus_InflammatoryPanel_Long_Binary$DIAGNOSIS == "Ulcerative Colitis", 1, 0)

#Escribir la tabla en formato largo y binario
write.csv(Plexus_InflammatoryPanel_Long_Binary, file.path("2.Datos_Limpios_UCvsCD", "Plexus_InflammatoryPanel_Long_Binary.csv"))


##Tabla para comparacion entre casos de CD####

#Filtrar tabla para las comparaciones Zonales para CD
crear_carpeta("2.1.Datos_Limpios_CD_Comparissons")

#Filtrar por unicamente Crohn Disease: Ileocolonic, Colonic, Ileal
Plexus_InflammatoryPanel_ONLY_CD <- Plexus_InflammatoryPanel %>% filter(str_detect(DISEASE_LOCATION,"Ileocolonic|Colonic|Ileal"))

#Escribir la tabla para las comparaciones Zonales para CD (Todos son CD)
write.csv(Plexus_InflammatoryPanel_ONLY_CD, file.path("2.1.Datos_Limpios_CD_Comparissons", "Plexus_InflammatoryPanel_ONLY_CD.csv"))

#Transformar la tabla Plexus_InflammatoryPanel_ONLY_CD a formato largo 
Plexus_InflammatoryPanel_ONLY_CD_Long <- Plexus_InflammatoryPanel_ONLY_CD%>%pivot_longer(cols = -c(1:5), names_to = "PROTEIN", values_to = "VALUE") %>% mutate(VALUE=as.numeric(str_replace(VALUE,",","."))) %>% drop_na()%>% select(ID,SEX,DIAGNOSIS,DISEASE_LOCATION,ETHNICITY,PROTEIN,VALUE)

#Guardar la tabla larga con texto
write.csv(Plexus_InflammatoryPanel_ONLY_CD_Long , file.path("2.1.Datos_Limpios_CD_Comparissons", "Plexus_InflammatoryPanel_ONLY_CD_Long.csv"))

#Convertir variables de la columna SEX a MALE=0, FEMALE=1
Plexus_InflammatoryPanel_ONLY_CD_Long_Binary <- Plexus_InflammatoryPanel_ONLY_CD_Long
Plexus_InflammatoryPanel_ONLY_CD_Long_Binary$SEX <-  ifelse(Plexus_InflammatoryPanel_ONLY_CD_Long_Binary$SEX == "Female", 1, 0)

#Convertir variables de la columna DISEASE LOCATION A "Ileocolonic"=0 "Colonic"=1  "Ileal"=2     
Plexus_InflammatoryPanel_ONLY_CD_Long_Binary <- Plexus_InflammatoryPanel_ONLY_CD_Long_Binary  %>%mutate(DISEASE_LOCATION = case_when(
    DISEASE_LOCATION == "Ileocolonic" ~ 0,
    DISEASE_LOCATION == "Colonic" ~ 1,
    DISEASE_LOCATION == "Ileal" ~ 2 
  ))

#Escribir la tabla en formato largo y binario(categoricas como numericas)
write.csv(Plexus_InflammatoryPanel_ONLY_CD_Long_Binary, file.path("2.1.Datos_Limpios_CD_Comparissons", "Plexus_InflammatoryPanel_ONLY_CD_Long_Binary.csv"))

#Pre-Visualizaciones####
crear_carpeta("3.Previsualizaciones")

#Boxplot
Boxplot_uc_vs_cd <- function(){
  Plexus_InflammatoryPanel_Long$DIAGNOSIS <- as.factor(Plexus_InflammatoryPanel_Long$DIAGNOSIS)
  ggplot(Plexus_InflammatoryPanel_Long, aes(x = DIAGNOSIS, y = value, fill = DIAGNOSIS)) +
    geom_boxplot() +
    labs(title = "Distribución de valores de proteínas por diagnóstico",
         x = "Diagnóstico", y = "Valor de la proteína") +
    theme_minimal()
}
plot <- Boxplot_uc_vs_cd()

ggsave(filename = file.path("3.Previsualizaciones", "Boxplots_por_diagnostic.jpg"), plot = plot, width = 8, height = 6, dpi = 300)

#Histograma
Histogram_uc_vs_cd <- function(){ 
  Plexus_InflammatoryPanel_Long$DIAGNOSIS <- as.factor(Plexus_InflammatoryPanel_Long$DIAGNOSIS)
  ggplot(Plexus_InflammatoryPanel_Long, aes(x = value, fill = SEX)) +
    geom_histogram(binwidth = 0.5, position = "dodge") +
    facet_wrap(~ DIAGNOSIS) +
    labs(title = "Distribución de valores de proteínas por diagnóstico y sexo",
         x = "Valor de la proteína", y = "Frecuencia") +
    theme_minimal()
}

plot2 <- Histogram_uc_vs_cd()
ggsave(filename = file.path("3.Previsualizaciones", "Histogram.jpg"), plot = plot2, width = 8, height = 6, dpi = 300)







