#Modelamiento de Biomarcadores proteicos UCvsCD y entre CD (Ileal, ileocolonic, colonic) usando datos del Plexus####

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
    cat("La carpeta ya existe en:", nombre_carpeta ,". Por lo tanto, se sobreescribira", "\n")
  }
}


#Modelamiento de biomarcadores proteicos UCvsCD####

#Cargar Data. | Crohn Disease =0 Ulcerative Colitis=1 | MALE=0, FEMALE=1. 
Plexus_InflammatoryPanel_Long_Binary <- read.csv("2.Datos_Limpios_UCvsCD/Plexus_InflammatoryPanel_Long_Binary.csv") %>% select(ID,SEX,DIAGNOSIS,PROTEIN,VALUE) 
Plexus_InflammatoryPanel_Long_Binary$DIAGNOSIS <- as.character(Plexus_InflammatoryPanel_Long_Binary$DIAGNOSIS) 
Plexus_InflammatoryPanel_Long_Binary$SEX <- as.character(Plexus_InflammatoryPanel_Long_Binary$SEX) 

#Modelamiento con Regresion Logistica
getAssocs_UCvsCD_Plexus <-  function(outliers){
  
  #Remover outliers deseados por el usuario
  #Calcula los cuantiles:
  #Calcula el cuantil Inferior
  qci <- Plexus_InflammatoryPanel_Long_Binary%>%group_by(PROTEIN)%>%summarise(qci=quantile(VALUE,(0+outliers),na.rm=T))
  
  #Calcula el cuantil Superior
  qcf <- Plexus_InflammatoryPanel_Long_Binary%>%group_by(PROTEIN)%>%summarise(qcf=quantile(VALUE,(1-outliers),na.rm=T))
  

  #Actualizar los valores atipicos encontrados como "NA"
  Plexus_InflammatoryPanel_Long_Binary_qc <- Plexus_InflammatoryPanel_Long_Binary%>%
    left_join(qci,by=c("PROTEIN"))%>%
    left_join(qcf,by=c("PROTEIN"))%>%
    mutate(value_qc=case_when(VALUE<qci~NA_real_,VALUE>qcf~NA_real_,TRUE~VALUE))
  
  #Escribir El Dataset generado en una nueva carpeta.
  nombre_carpeta <- paste0("Test_UCvsCD_outliers_",outliers)
  crear_carpeta(nombre_carpeta)
  write.csv(Plexus_InflammatoryPanel_Long_Binary_qc, file = paste0(nombre_carpeta,"/Plexus_InflammatoryPanel_Long_Binary_qc",outliers,"_.csv"), row.names=FALSE)
  
  #Histograma de los datos filtrados sin outliers 
  p2 <- ggplot(Plexus_InflammatoryPanel_Long_Binary_qc,aes(x=value_qc,fill=DIAGNOSIS))+
    geom_histogram()+
    facet_wrap(~PROTEIN,scales='free')+
    theme_minimal()
  
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_histograms_rawdata_qc.pdf"),p2,width=25,height=25)
  
  #Grafico de outliers eliminados
  p3 <- Plexus_InflammatoryPanel_Long_Binary_qc%>%group_by(PROTEIN)%>%
    filter(is.na(value_qc))%>%ggplot(aes(x=PROTEIN,fill=DIAGNOSIS))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_Number_of_outliers_removed_by_protein.pdf"),p3,width=15,height=7)

  
  #Realizar el modelo, Modelo Logistico. Se calcula el assoc de cada valor de proteina con el diagnosis, se corrige por SEX.
  #Realizar un bucle para cada una de las proteinas.
  assocs <- tibble()
  #Pasar La variable a predecir a numerica
  Plexus_InflammatoryPanel_Long_Binary_qc$DIAGNOSIS <- as.numeric(Plexus_InflammatoryPanel_Long_Binary_qc$DIAGNOSIS)
  
  for (p in unique(Plexus_InflammatoryPanel_Long_Binary_qc$PROTEIN)){
    #Cuadrar el modelo logistico
    assoc <- glm(DIAGNOSIS ~ value_qc + SEX , data = Plexus_InflammatoryPanel_Long_Binary_qc%>%filter(PROTEIN==p), family = binomial(link = "logit"))
    #Recibir el resumen del modelo
    summary(assoc)
    #Realizar el dataframe que sera el output
    assoc_sum <- summary(assoc)$coefficients%>%as.data.frame()%>%
      rownames_to_column("variable")%>%
      mutate(PROTEIN=p,.before=1)%>%
      rename(Pvalue="Pr(>|z|)")%>%
      mutate(variable=str_replace_all(variable,"\\(|\\)",""))%>%
      as_tibble() 
    colnames(assoc_sum) <- c("protein","variable","Estimate","StdError","Zvalue","Pvalue")
    assocs <- bind_rows(assocs,assoc_sum)
  }
  
  assocs%>%write_tsv(paste0(nombre_carpeta,"/ol_",outliers,"_assocs.txt"))
  
  
  ##Ajustar P valor mediante FDR
  assocdf <- assocs%>%filter(variable=="value_qc")%>%mutate(PvalueFDR=p.adjust(Pvalue,method="fdr"))
  
  #Grafico de puntos del estimado de las proteinas, con color basado en el pvalor (pvalor< 0.05 en rojo)
  sortedplot <- assocdf%>%arrange((Estimate))
  p4 <- ggplot(assocdf%>%mutate(protein=factor(protein,levels=sortedplot$protein)),aes(x=Estimate,y=protein,color=PvalueFDR<0.05))+geom_point(cex=4)+theme(axis.text.y = element_text(size=13))+scale_color_manual(values=c("#53545b","red"))+theme(legend.position = "none")+theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_dotplot_associations.pdf"),p4,width=8,height=15)
  
  # volcano plot # put the name of the proteins with pvalue < 0.05 in the plot with lines and do not overlap
  p5 <- ggplot(assocdf,aes(x=Estimate,y=-log10(Pvalue),label=ifelse(Pvalue<0.05,protein,""),color=PvalueFDR<0.05))+geom_point()+geom_text_repel()+theme(axis.text.y = element_text(size=8))+scale_color_manual(values=c("#53545b","red"))+theme(legend.position = "none")+theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_volcanoplot_associations.pdf"),p5,width=8,height=8)#
  
  #Obtener datasets con proteinas significativas.
  result <- subset(assocdf, PvalueFDR < 0.05) #Cambiar a PValue o PValueFDR segun si se quiere las significativas corregidas. 
  write.csv(result, file = paste0 (nombre_carpeta,"/Proteinas_significativas_CONFDR",outliers,"_.csv"), row.names = FALSE)
}

##llamada de la funcion
getAssocs_UCvsCD_Plexus(outliers=0.03)























#Modelamiento de biomarcadores proteicos  CD (ileal vs ileocolonic)####

#Cargar Data. | Ileocolonic"=0. | "Colonic"=1 | "Ileal"=2     
Plexus_InflammatoryPanel_ONLY_CD_Long_Binary <- read.csv("2.1.Datos_Limpios_CD_Comparissons/Plexus_InflammatoryPanel_ONLY_CD_Long_Binary.csv") %>% select(ID,SEX,DISEASE_LOCATION,PROTEIN,VALUE) 

#Filtrar unicamente por Ileal(2), y Ileocolonic(0)
Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo <- Plexus_InflammatoryPanel_ONLY_CD_Long_Binary%>% filter(DISEASE_LOCATION == 2 | DISEASE_LOCATION == 0 )

#El modelo unicamente recibe valores 0 y 1, por lo que Ileal sera (1)
Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo %>% mutate(DISEASE_LOCATION = case_when(DISEASE_LOCATION == 2 ~ 1,TRUE ~ DISEASE_LOCATION))


#Modelamiento con Regresion Logistica
getAssocs_IlealvsIleo_Plexus <-  function(outliers){
  
  #Remover outliers deseados por el usuario
  #Calcula los cuantiles:
  #Calcula el cuantil Inferior
  qci <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo%>%group_by(PROTEIN)%>%summarise(qci=quantile(VALUE,(0+outliers),na.rm=T))
  
  #Calcula el cuantil Superior
  qcf <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo%>%group_by(PROTEIN)%>%summarise(qcf=quantile(VALUE,(1-outliers),na.rm=T))
  
  
  #Actualizar los valores atipicos encontrados como "NA"
  Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo%>%
    left_join(qci,by=c("PROTEIN"))%>%
    left_join(qcf,by=c("PROTEIN"))%>%
    mutate(value_qc=case_when(VALUE<qci~NA_real_,VALUE>qcf~NA_real_,TRUE~VALUE))
  
  #Escribir El Dataset generado en una nueva carpeta.
  nombre_carpeta <- paste0("Test_CD_IlealvsIleo_outliers_",outliers)
  crear_carpeta(nombre_carpeta)
  write.csv(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc, file = paste0(nombre_carpeta,"/Plexus_InflammatoryPanel_Long_Binary_qc",outliers,"_.csv"), row.names=FALSE)
  
  #Pasar a Caracter para que sea posible la realizacion de graficas
  Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc$DISEASE_LOCATION <- as.character(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc$DISEASE_LOCATION) 
  Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc$SEX <- as.character(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc$SEX) 
  
  
  #Histograma de los datos filtrados sin outliers 
  p2 <- ggplot(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc,aes(x=value_qc,fill=DISEASE_LOCATION))+
    geom_histogram()+
    facet_wrap(~PROTEIN,scales='free')+
    theme_minimal()
  
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_histograms_rawdata_qc.pdf"),p2,width=25,height=25)
  
  #Grafico de outliers eliminados
  p3 <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc%>%group_by(PROTEIN)%>%
    filter(is.na(value_qc))%>%ggplot(aes(x=PROTEIN,fill=DISEASE_LOCATION))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_Number_of_outliers_removed_by_protein.pdf"),p3,width=15,height=7)
  
  
  #Realizar el modelo, Modelo Logistico. Se calcula el assoc de cada valor de proteina con el diagnosis, se corrige por SEX.
  #Realizar un bucle para cada una de las proteinas.
  assocs <- tibble()
  #Pasar La variable a predecir a numerica
  Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc$DISEASE_LOCATION <- as.numeric(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc$DISEASE_LOCATION)
  
  for (p in unique(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc$PROTEIN)){
    #Cuadrar el modelo logistico
    assoc <- glm(DISEASE_LOCATION ~ value_qc + SEX , data = Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsIleo_qc%>%filter(PROTEIN==p), family = binomial(link = "logit"))
    #Recibir el resumen del modelo
    summary(assoc)
    #Realizar el dataframe que sera el output
    assoc_sum <- summary(assoc)$coefficients%>%as.data.frame()%>%
      rownames_to_column("variable")%>%
      mutate(PROTEIN=p,.before=1)%>%
      rename(Pvalue="Pr(>|z|)")%>%
      mutate(variable=str_replace_all(variable,"\\(|\\)",""))%>%
      as_tibble() 
    colnames(assoc_sum) <- c("protein","variable","Estimate","StdError","Zvalue","Pvalue")
    assocs <- bind_rows(assocs,assoc_sum)
  }
  
  assocs%>%write_tsv(paste0(nombre_carpeta,"/ol_",outliers,"_assocs.txt"))
  
  
  ##Ajustar P valor mediante FDR
  assocdf <- assocs%>%filter(variable=="value_qc")%>%mutate(PvalueFDR=p.adjust(Pvalue,method="fdr"))
  
  #Grafico de puntos del estimado de las proteinas, con color basado en el pvalor (pvalor< 0.05 en rojo)
  sortedplot <- assocdf%>%arrange((Estimate))
  p4 <- ggplot(assocdf%>%mutate(protein=factor(protein,levels=sortedplot$protein)),aes(x=Estimate,y=protein,color=PvalueFDR<0.05))+geom_point(cex=4)+theme(axis.text.y = element_text(size=13))+scale_color_manual(values=c("#53545b","red"))+theme(legend.position = "none")+theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_dotplot_associations.pdf"),p4,width=8,height=15)
  
  # volcano plot # put the name of the proteins with pvalue < 0.05 in the plot with lines and do not overlap
  p5 <- ggplot(assocdf,aes(x=Estimate,y=-log10(Pvalue),label=ifelse(Pvalue<0.05,protein,""),color=PvalueFDR<0.05))+geom_point()+geom_text_repel()+theme(axis.text.y = element_text(size=8))+scale_color_manual(values=c("#53545b","red"))+theme(legend.position = "none")+theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_volcanoplot_associations.pdf"),p5,width=8,height=8)#
  
  #Obtener datasets con proteinas significativas.
  result <- subset(assocdf, PvalueFDR < 0.05) #Cambiar a PValue o PValueFDR segun si se quiere las significativas corregidas. 
  write.csv(result, file = paste0 (nombre_carpeta,"/Proteinas_significativas_CONFDR",outliers,"_.csv"), row.names = FALSE)
}

##llamada de la funcion
getAssocs_IlealvsIleo_Plexus(outliers=0.03)












#Modelamiento de biomarcadores proteicos  CD (ileal vs colonic)####

#Cargar Data. | Ileocolonic"=0. | "Colonic"=1 | "Ileal"=2     
Plexus_InflammatoryPanel_ONLY_CD_Long_Binary <- read.csv("2.1.Datos_Limpios_CD_Comparissons/Plexus_InflammatoryPanel_ONLY_CD_Long_Binary.csv") %>% select(ID,SEX,DISEASE_LOCATION,PROTEIN,VALUE) 

#Filtrar unicamente por Ileal(2), y Colonic(1)
Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic <- Plexus_InflammatoryPanel_ONLY_CD_Long_Binary%>% filter(DISEASE_LOCATION == 2 | DISEASE_LOCATION == 1 )

#El modelo unicamente recibe valores 0 y 1, por lo que Ileal sera (0)
Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic %>% mutate(DISEASE_LOCATION = case_when(DISEASE_LOCATION == 2 ~ 0,TRUE ~ DISEASE_LOCATION))


#Modelamiento con Regresion Logistica
getAssocs_IlealvsColonic_Plexus <-  function(outliers){
  
  #Remover outliers deseados por el usuario
  #Calcula los cuantiles:
  #Calcula el cuantil Inferior
  qci <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic%>%group_by(PROTEIN)%>%summarise(qci=quantile(VALUE,(0+outliers),na.rm=T))
  
  #Calcula el cuantil Superior
  qcf <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic%>%group_by(PROTEIN)%>%summarise(qcf=quantile(VALUE,(1-outliers),na.rm=T))
  
  
  #Actualizar los valores atipicos encontrados como "NA"
  Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic%>%
    left_join(qci,by=c("PROTEIN"))%>%
    left_join(qcf,by=c("PROTEIN"))%>%
    mutate(value_qc=case_when(VALUE<qci~NA_real_,VALUE>qcf~NA_real_,TRUE~VALUE))
  
  #Escribir El Dataset generado en una nueva carpeta.
  nombre_carpeta <- paste0("Test_CD_IlealvsColonic_outliers_",outliers)
  crear_carpeta(nombre_carpeta)
  write.csv(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc, file = paste0(nombre_carpeta,"/Plexus_InflammatoryPanel_Long_Binary_qc",outliers,"_.csv"), row.names=FALSE)
  
  #Pasar a Caracter para que sea posible la realizacion de graficas
  Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc$DISEASE_LOCATION <- as.character(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc$DISEASE_LOCATION) 
  Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc$SEX <- as.character(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc$SEX) 
  
  
  #Histograma de los datos filtrados sin outliers 
  p2 <- ggplot(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc,aes(x=value_qc,fill=DISEASE_LOCATION))+
    geom_histogram()+
    facet_wrap(~PROTEIN,scales='free')+
    theme_minimal()
  
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_histograms_rawdata_qc.pdf"),p2,width=25,height=25)
  
  #Grafico de outliers eliminados
  p3 <- Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc%>%group_by(PROTEIN)%>%
    filter(is.na(value_qc))%>%ggplot(aes(x=PROTEIN,fill=DISEASE_LOCATION))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_Number_of_outliers_removed_by_protein.pdf"),p3,width=15,height=7)
  
  
  #Realizar el modelo, Modelo Logistico. Se calcula el assoc de cada valor de proteina con el diagnosis, se corrige por SEX.
  #Realizar un bucle para cada una de las proteinas.
  assocs <- tibble()
  #Pasar La variable a predecir a numerica
  Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc$DISEASE_LOCATION <- as.numeric(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc$DISEASE_LOCATION)
  
  for (p in unique(Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc$PROTEIN)){
    #Cuadrar el modelo logistico
    assoc <- glm(DISEASE_LOCATION ~ value_qc + SEX , data = Plexus_InflammatoryPanel_ONLY_CD_Long_IlealvsColonic_qc%>%filter(PROTEIN==p), family = binomial(link = "logit"))
    #Recibir el resumen del modelo
    summary(assoc)
    #Realizar el dataframe que sera el output
    assoc_sum <- summary(assoc)$coefficients%>%as.data.frame()%>%
      rownames_to_column("variable")%>%
      mutate(PROTEIN=p,.before=1)%>%
      rename(Pvalue="Pr(>|z|)")%>%
      mutate(variable=str_replace_all(variable,"\\(|\\)",""))%>%
      as_tibble() 
    colnames(assoc_sum) <- c("protein","variable","Estimate","StdError","Zvalue","Pvalue")
    assocs <- bind_rows(assocs,assoc_sum)
  }
  
  assocs%>%write_tsv(paste0(nombre_carpeta,"/ol_",outliers,"_assocs.txt"))
  
  
  ##Ajustar P valor mediante FDR
  assocdf <- assocs%>%filter(variable=="value_qc")%>%mutate(PvalueFDR=p.adjust(Pvalue,method="fdr"))
  
  #Grafico de puntos del estimado de las proteinas, con color basado en el pvalor (pvalor< 0.05 en rojo)
  sortedplot <- assocdf%>%arrange((Estimate))
  p4 <- ggplot(assocdf%>%mutate(protein=factor(protein,levels=sortedplot$protein)),aes(x=Estimate,y=protein,color=PvalueFDR<0.05))+geom_point(cex=4)+theme(axis.text.y = element_text(size=13))+scale_color_manual(values=c("#53545b","red"))+theme(legend.position = "none")+theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_dotplot_associations.pdf"),p4,width=8,height=15)
  
  # volcano plot # put the name of the proteins with pvalue < 0.05 in the plot with lines and do not overlap
  p5 <- ggplot(assocdf,aes(x=Estimate,y=-log10(Pvalue),label=ifelse(Pvalue<0.05,protein,""),color=PvalueFDR<0.05))+geom_point()+geom_text_repel()+theme(axis.text.y = element_text(size=8))+scale_color_manual(values=c("#53545b","red"))+theme(legend.position = "none")+theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_volcanoplot_associations.pdf"),p5,width=8,height=8)#
  
  #Obtener datasets con proteinas significativas.
  result <- subset(assocdf, PvalueFDR < 0.05) #Cambiar a PValue o PValueFDR segun si se quiere las significativas corregidas. 
  write.csv(result, file = paste0 (nombre_carpeta,"/Proteinas_significativas_CONFDR",outliers,"_.csv"), row.names = FALSE)
}

##llamada de la funcion
getAssocs_IlealvsColonic_Plexus(outliers=0.03)










#Modelamiento de biomarcadores proteicos  CD (Ileocolonic vs colonic)####

#Cargar Data. | Ileocolonic"=0. | "Colonic"=1 | "Ileal"=2     
Plexus_InflammatoryPanel_ONLY_CD_Long_Binary <- read.csv("2.1.Datos_Limpios_CD_Comparissons/Plexus_InflammatoryPanel_ONLY_CD_Long_Binary.csv") %>% select(ID,SEX,DISEASE_LOCATION,PROTEIN,VALUE) 

#Filtrar unicamente por Ileocolonic(0), y Colonic(1)
Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic <- Plexus_InflammatoryPanel_ONLY_CD_Long_Binary%>% filter(DISEASE_LOCATION == 0 | DISEASE_LOCATION == 1 )

#Modelamiento con Regresion Logistica
getAssocs_IleovsColonic_Plexus <-  function(outliers){
  
  #Remover outliers deseados por el usuario
  #Calcula los cuantiles:
  #Calcula el cuantil Inferior
  qci <- Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic%>%group_by(PROTEIN)%>%summarise(qci=quantile(VALUE,(0+outliers),na.rm=T))
  
  #Calcula el cuantil Superior
  qcf <- Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic%>%group_by(PROTEIN)%>%summarise(qcf=quantile(VALUE,(1-outliers),na.rm=T))
  
  
  #Actualizar los valores atipicos encontrados como "NA"
  Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc <- Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic%>%
    left_join(qci,by=c("PROTEIN"))%>%
    left_join(qcf,by=c("PROTEIN"))%>%
    mutate(value_qc=case_when(VALUE<qci~NA_real_,VALUE>qcf~NA_real_,TRUE~VALUE))
  
  #Escribir El Dataset generado en una nueva carpeta.
  nombre_carpeta <- paste0("Test_CD_IleovsColonic_outliers_",outliers)
  crear_carpeta(nombre_carpeta)
  write.csv(Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc, file = paste0(nombre_carpeta,"/Plexus_InflammatoryPanel_Long_Binary_qc",outliers,"_.csv"), row.names=FALSE)
  
  #Pasar a Caracter para que sea posible la realizacion de graficas
  Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc$DISEASE_LOCATION <- as.character(Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc$DISEASE_LOCATION) 
  Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc$SEX <- as.character(Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc$SEX) 
  
  
  #Histograma de los datos filtrados sin outliers 
  p2 <- ggplot(Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc,aes(x=value_qc,fill=DISEASE_LOCATION))+
    geom_histogram()+
    facet_wrap(~PROTEIN,scales='free')+
    theme_minimal()
  
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_histograms_rawdata_qc.pdf"),p2,width=25,height=25)
  
  #Grafico de outliers eliminados
  p3 <- Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc%>%group_by(PROTEIN)%>%
    filter(is.na(value_qc))%>%ggplot(aes(x=PROTEIN,fill=DISEASE_LOCATION))+
    geom_bar()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_Number_of_outliers_removed_by_protein.pdf"),p3,width=15,height=7)
  
  
  #Realizar el modelo, Modelo Logistico. Se calcula el assoc de cada valor de proteina con el diagnosis, se corrige por SEX.
  #Realizar un bucle para cada una de las proteinas.
  assocs <- tibble()
  #Pasar La variable a predecir a numerica
  Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc$DISEASE_LOCATION <- as.numeric(Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc$DISEASE_LOCATION)
  
  for (p in unique(Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc$PROTEIN)){
    #Cuadrar el modelo logistico
    assoc <- glm(DISEASE_LOCATION ~ value_qc + SEX , data = Plexus_InflammatoryPanel_ONLY_CD_Long_IleovsColonic_qc%>%filter(PROTEIN==p), family = binomial(link = "logit"))
    #Recibir el resumen del modelo
    summary(assoc)
    #Realizar el dataframe que sera el output
    assoc_sum <- summary(assoc)$coefficients%>%as.data.frame()%>%
      rownames_to_column("variable")%>%
      mutate(PROTEIN=p,.before=1)%>%
      rename(Pvalue="Pr(>|z|)")%>%
      mutate(variable=str_replace_all(variable,"\\(|\\)",""))%>%
      as_tibble() 
    colnames(assoc_sum) <- c("protein","variable","Estimate","StdError","Zvalue","Pvalue")
    assocs <- bind_rows(assocs,assoc_sum)
  }
  
  assocs%>%write_tsv(paste0(nombre_carpeta,"/ol_",outliers,"_assocs.txt"))
  
  
  ##Ajustar P valor mediante FDR
  assocdf <- assocs%>%filter(variable=="value_qc")%>%mutate(PvalueFDR=p.adjust(Pvalue,method="fdr"))
  
  #Grafico de puntos del estimado de las proteinas, con color basado en el pvalor (pvalor< 0.05 en rojo)
  sortedplot <- assocdf%>%arrange((Estimate))
  p4 <- ggplot(assocdf%>%mutate(protein=factor(protein,levels=sortedplot$protein)),aes(x=Estimate,y=protein,color=PvalueFDR<0.05))+geom_point(cex=4)+theme(axis.text.y = element_text(size=13))+scale_color_manual(values=c("#53545b","red"))+theme(legend.position = "none")+theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_dotplot_associations.pdf"),p4,width=8,height=15)
  
  # volcano plot # put the name of the proteins with pvalue < 0.05 in the plot with lines and do not overlap
  p5 <- ggplot(assocdf,aes(x=Estimate,y=-log10(Pvalue),label=ifelse(Pvalue<0.05,protein,""),color=PvalueFDR<0.05))+geom_point()+geom_text_repel()+theme(axis.text.y = element_text(size=8))+scale_color_manual(values=c("#53545b","red"))+theme(legend.position = "none")+theme_minimal()
  ggsave(paste0(nombre_carpeta,"/ol_",outliers,"_volcanoplot_associations.pdf"),p5,width=8,height=8)#
  
  #Obtener datasets con proteinas significativas.
  result <- subset(assocdf, PvalueFDR < 0.05) #Cambiar a PValue o PValueFDR segun si se quiere las significativas corregidas. 
  write.csv(result, file = paste0 (nombre_carpeta,"/Proteinas_significativas_CONFDR",outliers,"_.csv"), row.names = FALSE)
}

##llamada de la funcion
getAssocs_IleovsColonic_Plexus(outliers=0.03)
















