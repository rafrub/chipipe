library("ChIPseeker")
library("org.At.tair.db")
library("TxDb.Athaliana.BioMart.plantsmart28")
library("clusterProfiler")
library("pathview")

#Definimos todos los parámetros
args<- commandArgs(trailingOnly = T)
peak.file <- args[[1]]
chromosomes <- as.character(args[[2]])
tssup<-as.numeric(args[[3]])
tssdown<-as.numeric(args[[4]])
peak.type<-as.numeric(args[[5]])

#Definimos las regiones promotoras, tiene que leer el fichero de picos

#Defining promoter regions
peaks <- readPeakFile(peakfile= peak.file, header=F)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
promoter <- getPromoters(TxDb= txdb, upstream = tssup, downstream = tssdown)

#Calculating the peak distribution along the genome
covplot(peaks, weightCol = "V5")
#la columna v5 es la que nos interesa ya que es la que nos dice la posición


#Annotating peaks according to the types of DNA regions they bind to
peakAnno<- annotatePeak(peak=peaks, tssRegion=c(-(tssdown), tssup), 
                        TxDb = txdb, annoDb = "org.At.tair.db")

plotAnnoPie(peakAnno,main="Distribution of TF binding regions",
            outer=T,line=-2)

plotDistToTSS(peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

## Saving peaks that bind to promoter regions.
df.annotation <- as.data.frame(peakAnno)
#escribimos el peakanno como data frame, y una de las columnas de data frame
#es annotation
{
  if (peak.type == 1)
    #este if dice que si el tipo de pico es 1 (narrow, FT), coge solo los que se unan
    #al promotor, porque los demás picos pueden ser ruido.
    
  {
    promoter.df.annotation <- subset(df.annotation,annotation=="Promoter")
    #subset hace que coja solo los que en anotación tenga promoter
    #subset se utiliza para filtrar tablas, metes los datos que quieres y la condición
    #que quieres que se dé. Ve a la columna annotation y todo lo que tenga el nombre 
    #promoter coger esa fila.
  }
  else if (peak.type == 2)
  {
    annotation.to.exclude <- c("Distal Intergenic", "Downstream (1-2kb)","Downstream (2-3kb)","Downstream (<1kb)")
    #en un espacio intergénico no nos interesa si hay una marca de histonas porque no 
    #podemos decir qué está regulando, por tanto excluimos ese tipo de picos
    annotation.to.include <- setdiff(unique(df.annotation$annotation),annotation.to.exclude)
    #setdiff: metes dos conjuntos, devuelve los elementos del primer conjunto que no están en el segundo
    #si lo primero es nuestra anotación completa, nos quedamos solo con la parte de la tabña que no tenga
    #lo que hemos definido antes
    promoter.df.annotation <- subset(df.annotation,df.annotation$annotation %in% annotation.to.include)
    #el in busca si un elemento está en una det lista o no
  }
}


## Listing genes affected by the TF (its regulome).
regulome <- promoter.df.annotation$geneId
write(regulome,file = "regulome.txt")

## Defining universe for GO & kegg terms enrichment.
genes.atha <- as.data.frame(genes(txdb))
{
  if (chromosomes == "all")
  {
    my.universe<-genes.atha$gene_id
  }
  else
  {
    chromosomes.for.universe <- as.vector(strsplit(chromosomes,",")[[1]])
    #str split sirve para cortar una cadena mediante comas o espacios
    genes.of.my.universe <- subset(genes.atha,seqnames==chromosomes.for.universe)
    my.universe <- genes.of.my.universe$gene_id
  }
}



## GO terms BP enrichment.
ego <- enrichGO(gene          = regulome,
                universe      = my.universe,
                OrgDb         = org.At.tair.db,
                ont           = "BP",
                keyType = "TAIR", pvalueCutoff = 0.05)
#pvaluecutoff en le tenemos que poner el pvalue que nosotros queramos
GO.enrichment <- as.data.frame(ego)
#guardamos el enriquecimiento como data frame
write.table(GO.enrichment,file = "go_terms.tsv",sep="\t",row.names = F)

#si no sale ningún término de enriquecimiento, da error al crear el data frame
#y entonces ya no ejecuta nada de lo demás, por tanto ponemos el siguiente if
{
  
  if(nrow(GO.enrichment) == 0)
  {
    print("No enrichment of GO terms for biological processes detected.")
  }
  
  else
  {
    pdf(file = "plots_go_bp.pdf",width = 14, height = 14)
    #a partir de aquí todo lo que ponga lo va a escribir en este pdf,
    #por tanto ya no lo devuelve dentro de R
    
    if (nrow(GO.enrichment) > 5)
    {
      b <- goplot(ego)
      plot(b)
    }
    c <- barplot(ego,showCategory = 30,title = "Barplot of detected GO terms for biological processes")
    plot(c)
    
    d <- dotplot(ego, showCategory=30,title = "Dotplot of detected GO terms for biological processes")
    plot(d)
    
    e <- cnetplot(ego,colorEdge=T,layout = "circle",circular=T,cex_category=2,cex_label_category=2)
    plot(e)
    
    dev.off() #esto es para cerrar el pdf y volver a devolver todo a través de R
  }
}



#Volvemos a hacer lo mismo para enriquecimiento en molecular function y cellular component

#GO terms MF enrichment.
ego.mf <- enrichGO(gene       = regulome,
                   universe      = my.universe,
                   OrgDb         = org.At.tair.db,
                   ont           = "MF",
                   keyType = "TAIR", pvalueCutoff = 0.05)
GO.enrichment.mf <- as.data.frame(ego.mf)
write.table(GO.enrichment.mf,file = "go_terms_mf.tsv",sep="\t",row.names = F)

{
  
  if(nrow(GO.enrichment.mf) == 0)
  {
    print("No enrichment of GO terms for molecular functions detected.")
  }
  
  else
  {
    pdf(file = "plots_go_mf.pdf",width = 14, height = 14)
    
    if (nrow(GO.enrichment.mf) > 5)
    {
      b <- goplot(ego.mf)
      plot(b)
    }
    c <- barplot(ego.mf,showCategory = 30,title = "Barplot of detected GO terms for molecular functions")
    plot(c)
    
    d <- dotplot(ego.mf, showCategory=30,title = "Dotplot of detected GO terms for molecular functions")
    plot(d)
    
    e <- cnetplot(ego.mf,colorEdge=T,layout = "circle",circular=T,cex_category=2,cex_label_category=2)
    plot(e)
    
    dev.off()
  }
}

#GO terms CC enrichment.
ego.cc <- enrichGO(gene          = regulome,
                   universe      = my.universe,
                   OrgDb         = org.At.tair.db,
                   ont           = "CC",
                   keyType = "TAIR", pvalueCutoff = 0.05)
GO.enrichment.cc <- as.data.frame(ego.cc)
write.table(GO.enrichment.cc,file = "go_terms_cc.tsv",sep="\t",row.names = F)

{
  
  if(nrow(GO.enrichment.cc) == 0)
  {
    print("No enrichment of GO terms for celular components detected.")
  }
  
  else
  {
    pdf(file = "plots_go_cc.pdf",width = 14, height = 14)
    
    if (nrow(GO.enrichment.cc) > 5)
    {
      b <- goplot(ego.cc)
      plot(b)
    }
    c <- barplot(ego.cc,showCategory = 30,title = "Barplot of detected GO terms for cellular components")
    plot(c)
    
    d <- dotplot(ego.cc, showCategory=30,title = "Dotplot of detected GO terms for cellular components")
    plot(d)
    
    e <- cnetplot(ego.cc,colorEdge=T,layout = "circle",circular=T,cex_category=2,cex_label_category=2)
    plot(e)
    
    dev.off()
  }
}


## KEGG terms enrichment.
pathway.enrich <- as.data.frame(enrichKEGG(gene=regulome, keyType = "kegg", 
                                           organism="ath", pvalueCutoff = 0.1, pAdjustMethod = "BH"))
write.table(pathway.enrich,file = "kegg_terms.tsv",sep="\t",row.names = F)
{
  if (nrow(pathway.enrich) == 0)
  {
    print("No enrichment of KEGG pathways detected.")
  }
  else
  {
    for (i in 1:nrow(pathway.enrich))
    {
      pathway.current.id <- pathway.enrich[i,1]
      pathview(gene.data = regulome, pathway.id = pathway.current.id,species = "ath",
               gene.idtype = "TAIR", kegg.dir = "kegg_images/")
    }
  }
}
