library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

if(!require(nationalparkcolors)){
  devtools::install_github("katiejolly/nationalparkcolors")
  library(nationalparkcolors)
}

pal<-c(park_palette('DeathValley'),park_palette("GrandTeton"))

##read in imaging features 
imgFeat <-readxl::read_xlsx('data/cptac_hnc_pyrad_features_v2.xlsx')
##convert to long format
imgFeat <- imgFeat %>%
  pivot_longer(3:ncol(imgFeat),names_to='patient',values_to='value')

##read in proteomics
impProts <- read.table('data/S054_HNSCC_imputed_0920.tsv',sep='\t',header=T,check.names=F)

##convert to long format
tumProts <- impProts%>%
  dplyr::select(Gene,ends_with('-T')) 

##convert to long format
normProts <- impProts%>%
  dplyr::select(Gene,ends_with('-N'))

normProts <- normProts%>%
  pivot_longer(2:ncol(normProts),names_to='patient',values_to='inferredNormAbund')%>%
  mutate(patient=stringr::str_remove(patient,'-N'))

tumProts <- tumProts%>%
  pivot_longer(2:ncol(tumProts),names_to='patient',values_to='inferredTumAbund')%>%
  mutate(patient=stringr::str_remove(patient,'-T'))


##can we get lcinical too?
clin.dat<-read.table('data/HNSCC_Feb2021_clinical_data.csv',sep=',',header=T,check.names=F,quote='"')%>%
  dplyr::select(case_id,contains("lymph_node"))%>%### this gets us 8 different categories
  rename(extranodal_extension="baseline/lymph_nodes_extranodal_extension",
         staging="baseline/pathologic_staging_regional_lymph_nodes_pn",
         dissection_performed="baseline/lymph_node_neck_dissection_performed",
         dissection_left="baseline/method_of_lymph_node_neck_dissection_left",
         dissection_right="baseline/method_of_lymph_node_neck_dissection_right",
         ln_examined="baseline/number_of_lymph_nodes_examined",
         positive_by_he="baseline/number_of_lymph_nodes_positive_for_tumor_by_he",
         postiive_by_ihc="baseline/number_of_lymph_nodes_positive_for_tumor_by_ihc_staining_only",
         patient='case_id')%>%
  mutate(frac_ln=positive_by_he/ln_examined)%>%
  subset(staging!='')%>%
  mutate(staging=toupper(staging))



##reduced clinical data
red.dat<-clin.dat%>%
  select(patient,staging,positive_by_he)

red.dat$lnCounts<-sapply(red.dat$positive_by_he,function(x) ifelse(x==0,'zeroLN',ifelse(x>2,'moreThanTwo','oneOrTwo')))


##now combine the abundances to the same table
fullProts <- normProts%>%inner_join(tumProts,by=c('Gene','patient'))%>%
  mutate(tumNormDiff=inferredTumAbund/inferredNormAbund)%>%
  left_join(clin.dat,by='patient')


##helper functions
##quick function to do PCA on a column value
pcaFromDf<-function(df,column){
  fn=c(mean)
  names(fn)<-column
  mat <-df%>%
    dplyr::select('Gene','patient',column)%>%
    pivot_wider(names_from='patient',values_from=column,values_fn=fn)%>%
    tibble::column_to_rownames('Gene')%>%
    as.matrix()
  return(data.frame(prcomp(t(mat))$x)%>%tibble::rownames_to_column('patient'))
}


## read in gene list
invList <-read.table('data/pouliquenEtAlProteins.tsv',sep='\t',header=T,check.names=F)