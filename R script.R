##########################################
## R script related to Hermanson et al. ##
##########################################

set.seed(1)

librar(ape)
library(boot)
library(cluster)

#Read character-taxon matrix
turtle_mat <- read.table('SupplementaryFile1_CharacterTaxonMatrix.csv',header = T, sep=',',
                         row.names = 1)
head(turtle_mat)


#Replace '-' and '?' with NAs before analysis

turtle_mat_full <- t(apply ( turtle_mat , 1, function(x) ifelse(x == '-' | x == '?',NA,x)  ))
turtle_mat_full <- data.frame(apply(turtle_mat_full[,1:15],2,as.integer), 'Clade'=turtle_mat_full[,"Clade"])
head(turtle_mat_full)
  
#Principal Coordinate Analysis (PCoA)
dist_mat_full <- daisy(turtle_mat_full[,1:15],metric = 'gower')

pcoa_full <- pcoa(dist_mat_full,correction='cailliez')


#Define major groups to analyse (some were not included in our final analyses...)

groups_full <- lapply ( levels(as.factor(turtle_mat$Clade)) , 
                        function(x) rownames(subset(turtle_mat,Clade==x)) )
names(groups_full) <- levels(as.factor(turtle_mat$Clade)) 
  groups_full$Trionychidae <- groups_full$Trionychia[-1]
  
  groups_full$Testudinata <- rownames(turtle_mat_full)[turtle_mat_full$Clade!='non_shelled']
  
  groups_full$Testudinata_excl_marine <- rownames(turtle_mat_full)[turtle_mat_full$Clade!='non_shelled' & turtle_mat_full$Clade!='Chelonioidea']
  groups_full$Testudinata_excl_marine <- groups_full$Testudinata_excl_marine[-13]
  
  groups_full$Testudinata_excl_terrestrial <- rownames(turtle_mat_full)[turtle_mat_full$Clade!='non_shelled' & turtle_mat_full$Clade!='Testudinidae']
  
  groups_full$Pan_Testudines <- rownames(turtle_mat_full)
  
  groups_full$Testudinata_excl_Chelydr <- rownames(turtle_mat_full)[turtle_mat_full$Clade!='non_shelled' & turtle_mat_full$Clade!='Chelydroidea']
  
  groups_full$crown <- rownames(turtle_mat_full)[turtle_mat_full$Clade!='non_shelled' & turtle_mat_full$Clade!='Triassic_stem']
  groups_full$crown_minus_marine <- groups_full$crown[-c(10,18,19)] #remove chelonioids and Carettochelys
  groups_full$crown_minus_terr <- groups_full$crown[-c(25:30)] #remove testudinids
  
  order_groups <- c('Pan_Testudines','non_shelled','Testudinata','Testudinata_excl_marine',
                    'Testudinata_excl_terrestrial','crown',
                    'crown_minus_marine','crown_minus_terr',
                    'Testudinata_excl_Chelydr',
                    'Triassic_stem','Chelidae','Pelomedusoides','Trionychia','Trionychidae',
                  'Chelydroidea','Chelonioidea',
                  'Emysternia','Testudinidae','Geoemydidae')

groups_full <- groups_full[order_groups]

#Disparity analysis

#Function to calculate sum of variances (as the disparity metric). Will be used to run the bootstrapped calculations of the metric
sum_var <- function(data,indices) {
  d <- data[indices,]
  sov <- sum(apply(d,2,var))
  return(sov)
}

#Run 1000 bootstrap replicates of calculations of sum of variances using our PCoA vectors as the input data
sov_boot <- lapply ( groups_full, function(x) boot(data= pcoa_full$vectors[x,],statistic = sum_var, R=1000)$t  )

#Create individual objects to store the sum of variances calculations replicates for key groups

sov_Pan <-  sov_boot$Pan_Testudines #for Pan-Testudines (whole group)
sov_non_shelled <- sov_boot$non_shelled #for non-shelled stem turtles
sov_Testudinata <- sov_boot$Testudinata #for shelled turtles (Testudinata)
sov_crown <- sov_boot$crown #for crown group turtles (Testudines)
sov_crown_minus_marine <- sov_boot$crown_minus_marine #for crown group turtles, but excluding the highly aquatic species (chelonioids + Carettochelys)
sov_crown_minus_terr <- sov_boot$crown_minus_terr #for crown group turtles, but excluding the highly terrestrial species (testudinids)

#Run ANOVA comparisons

comparisons <- list(

#Pan-Testudines vs. Testudinata
'Pan_vs_Testudinata'= summary(aov(c(sov_Pan,sov_Testudinata) ~ rep(c('Pan','Testudinata'),each=1000)  )),

#Pan-Testudines vs. crown
'Pan_vs_crown' = summary(aov(c(sov_Pan,sov_crown) ~ rep(c('Pan','crown'),each=1000)  )),

#Testudinata vs. crown
'Testudinata_vs_crown'= summary(aov(c(sov_Testudinata,sov_crown) ~ rep(c('Testudinata','crown'),each=1000)  )),

#Crown vs. crown minus marine
'Crown_vs_crown_minus_marine'= summary(aov(c(sov_crown,sov_crown_minus_marine) ~ rep(c('crown','crown_minus_marine'),each=1000)  )),

#Crown vs. crown minus terrestrial
'Crown_vs_crown_minus_terr' = summary(aov(c(sov_crown,sov_crown_minus_terr) ~ rep(c('crown','crown_minus_terr'),each=1000)  ))

)

comparisons

#Adjusted P-values using Bonferroni correction

p.adjust( unlist( lapply(comparisons,function(x) x[[1]][['Pr(>F)']][1] )) , method = 'bonferroni')


#Code to produce Figure 29 of the main text (Disparity analyses). Later edited in Adobe Illustrator...

pdf('.pdf',width = 7,height = 7,useDingbats = F)

par(mfrow=c(2,2))

my_pch_full <- c(21,21,21,21,21,25,21,21,22,21)
cols_pco <- c(rep('#bb5566',5),'white','#bb5566','#bb5566','grey60','#bb5566')

scree.data <- apply ( pcoa_full$vectors , 2 , var ) / sum ( apply ( pcoa_full$vectors , 2 , var ) ) * 100

plot(pcoa_full$vectors[,1:2], pch=my_pch_full[as.factor(turtle_mat_full$Clade)],
     col='grey50',
     bg = adjustcolor( cols_pco[as.numeric(as.factor(turtle_mat_full$Clade))],alpha.f = 0.8), 
     cex=1.5,cex.axis=0.8,
     xlab = paste0('PCo1 [',round(scree.data,1)[1],'%]'), ylab=paste0('PCo2 [',round(scree.data,1)[2],'%]') )

legend('bottomright',legend=c('non-shelled turtles',
                             'non-Testudines\ntestudinatans',
                             'Testudines\n(crown turtles)'),
       cex=0.8, pch=c(25,22,21), col='grey50', pt.bg=c('white','grey20','#bb5566'),bty='n',
       pt.cex=1.25)

#Boxplots for disparity bootstraps

to_plot <- names(sov_boot)[c(1,2,10,3,6,7,8)] #only key groups discussed in the main text
cols_viol <- c('grey90','grey90','grey40', 'grey40','#bb5566','#bb5566','#bb5566')

plot('n',ylim=range(unlist(sov_boot)),
     xlim = c(0,7), xaxt='n', ylab = 'Disparity (sum of variances)',xlab='',
     cex.axis=0.8)
#axis(1,at=mean(c(0,7)),lab='',line = NA,tick = F)

for ( i in 1:length(to_plot)  ){
  
  stripchart(sov_boot[[to_plot[i]]][,1] , vertical=T , at = i-0.25,
             pch=16, 
             col= adjustcolor(cols_viol[i] , alpha.f = 0.3) , add=T , method='jitter' )
  
  boxplot(sov_boot[[to_plot[i]]][,1] , at = i-0.25,
          xlab=NULL, names=NA,add=T, 
          col = adjustcolor(cols_viol[i],alpha.f = 0.5),ylab=NULL,yaxt='n',pch=NA)
  
}

dev.off()
