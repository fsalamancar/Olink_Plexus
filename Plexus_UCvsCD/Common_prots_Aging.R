#Resultados de los analisis
Proteinas_significativas_UCvsCD_Plexus <- read.csv("Test_UCvsCD_outliers_0.03/Proteinas_significativas_CONFDR0.03_.csv") %>% select(protein) %>% as.data.frame()
Proteinas_significativas_IleovsColonic_Plexus <- read.csv("Test_CD_IleovsColonic_outliers_0.03/Proteinas_significativas_CONFDR0.03_.csv") %>% select(protein)
Proteinas_significativas_IlealvsIleo_Plexus <- read.csv("Test_CD_IlealvsIleo_outliers_0.03/Proteinas_significativas_CONFDR0.03_.csv") %>% select(protein)
Proteinas_significativas_IlealvsColonic_Plexus <- read.csv("Test_CD_IlealvsColonic_outliers_0.03/Proteinas_significativas_CONFDR0.03_.csv") %>% select(protein)


#Proteinas Aging Clock,
Data_result <- read_excel("1.Datos_Raw/Angie_Proteins/Data_result.xlsx") %>% drop_na()%>% as.data.frame()
Data_Plasma_protein_based <- read_excel("1.Datos_Raw/Angie_Proteins/Data_Plasma protein-based.xlsx")%>% drop_na()


#Buscar Proteinas del analisis en proteinas Aging Clock
#UCvsCD
UCvsCD_Result<- left_join(Proteinas_significativas_UCvsCD_Plexus,Data_result, by = "protein")%>% drop_na()
UCvsCD_Result_Plasma_protein<- left_join(Proteinas_significativas_UCvsCD_Plexus,Data_Plasma_protein_based, by = "protein")%>% drop_na()


#IleovsColonic
IleovsColonic_Result <- left_join(Proteinas_significativas_IleovsColonic_Plexus,Data_result, by = "protein")%>% drop_na()
IleovsColonic_Result_Plasma_protein<- left_join(Proteinas_significativas_IleovsColonic_Plexus,Data_Plasma_protein_based, by = "protein")%>% drop_na()



#IlealvsIleo
IlealvsIleo_Result<- left_join(Proteinas_significativas_IlealvsIleo_Plexus ,Data_result, by = "protein")%>% drop_na()
IlealvsIleo_Result_Plasma_protein<- left_join(Proteinas_significativas_IlealvsIleo_Plexus ,Data_Plasma_protein_based, by = "protein")%>% drop_na()


#IlealvsColonic
IlealvsColonic_Result<- left_join(Proteinas_significativas_IlealvsColonic_Plexus,Data_result, by = "protein")%>% drop_na()
IlealvsColonic_Result_Plasma_protein<- left_join(Proteinas_significativas_IlealvsColonic_Plexus,Data_Plasma_protein_based, by = "protein")%>% drop_na()



dos_papers <- left_join(Data_result,Data_Plasma_protein_based, by = "protein")%>% drop_na()
