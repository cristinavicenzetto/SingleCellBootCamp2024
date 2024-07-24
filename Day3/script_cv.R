setwd("/home/pkf/gendata2/bioinformatica/rcarriero/corso")

#Upload libraries

library(Seurat)
library(dplyr)
library(clustree)
library(ggraph)

#Upload Seurat object

tcells<- readRDS('Treg.rds')
#per quante cellule
ncol(tcells)
#per quanti gene
nrow(tcells)

tcells@assays$RNA # c'è data e fetures
#si può assegnare all'oggetto con la funzione ident per dire come suddividerli, qui settato per responder e no

#si può può cambiare l'ident, in questo caso per campione
Idents(tcells) <- tcells$orig.ident

#Find variable features sulla matrice normalizzata già
tcells <- FindVariableFeatures(tcells, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tcells), 10)
#VariableFeatures(tcells) funztion to extraxct the 10 highetst
# plot variable features with and without labels
pdf("plot_variable_features.pdf", width = 15, height = 10)
plot1 <- VariableFeaturePlot(tcells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

#Scale data
#in this case si scalano tutti per evitare di perdere dati che possono essere di interesse. Si può usre solo the most variable perchè son quelli che si usano per PCA e UMAP. Meglio fare prima opzione perchè poi si può voler vedere quelli che non sono ipervariable (che magari non aiutano nella clusterizzazione)
all.genes <- rownames(tcells)
tcells <- ScaleData(tcells, features = all.genes)

#Run PCA

tcells <- RunPCA(tcells, features = VariableFeatures(object = tcells))

# Examine and visualize PCA results a few different ways
#vediamo le 5 prime feateure dei primi 5 pc, si vedono positive correlation e negative correlation
print(tcells[["pca"]], dims = 1:5, nfeatures = 5)

#vedere il pca con
pdf("vizdim_pca_plot.pdf")
VizDimLoadings(tcells, dims = 1:2, reduction = "pca")#vedere i geni che più sono rilevanti nelle pc selezionati, si può colorare la pca in base all'active.ident, che va settato nel dataset (di default non c'è)
DimPlot(tcells, reduction = "pca")
DimHeatmap(tcells, dims = 1, cells = 500, balanced = TRUE)# si vedono le heatmap per dimensioni per capire quali sono i geni p
DimHeatmap(tcells, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#compute the pvalue per each dimension and to choose how many pc
tcells <- JackStraw(tcells, num.replicate = 100)
tcells <- ScoreJackStraw(tcells, dims = 1:20)

pdf("JackStrawPlot.pdf")
JackStrawPlot(tcells, dims = 1:20)
ElbowPlot(tcells)
dev.off()

#Select significant components
tcells <- FindNeighbors(tcells, dims = 1:14)#per trovare le distanze fra due cellule in base al pattern di espressione, si può aggiungere le dimensione significative che si vogliono analizzare

cells_clust_tree <- FindClusters(tcells, resolution = seq(from=0, to=2, by=0.1), print.output = 0, save.SNN = T)#si raggruppano le cellule in base alle similarities dei vicini; parameter importante: resolution, c'è uno standard, ma è meglio scelierlo bene per avere cluster troppo grandi o cluster stroppo frammentati, quindi con poche celllule. Si fa una seq da 0 a risoluzione 2 (molto alto) e suddividere per quanti liveli di resolution servono
#per scegliere la risuoluzione si usa clustree che è la versione grafica

pdf("CL7_tregs_Resolution_tree_up_to_2.0_pca_1_17.pdf",  width=17, height=17)
clustree(cells_clust_tree) #vedi tutti i colori
clustree(cells_clust_tree, node_colour = sc3_stability)#colorare per stability score
dev.off()
#il cluster 0 ha tutte e cell in una popolazione, si deve guardare come le cellule si muovono in base al numero dicluster che aumentano. The edges are considered base on the cells in the cluater. purple arrow: cells in this cluster can be good in another cluster by increasing resolution. Ad un certo punto ci sono molte frecce perchè non c'è più stabilità e le cellule fittano in tutti i cluster.
#criteri:
#lecel of stable resolution in base al numero delle frecce
#in base al colore: più è intenso più cellule possono muoversi
#si devono considerare il numero di cellule di partenz

#scelto 0.4 perchè tirando unariga  si vede che è il più stabile divendendo 
#controlla bene la vignetta, si può separare per un metadata di una serie di geni che identificano una data popolazione

#si possono scegliere più risoluzione (2 o3) e poi fare le analisi per capire la separazione

#identificazione cluster dopo aver scelto la risoluzione si fa la umap 
tcells <- FindClusters(tcells, resolution = 0.4)
tcells <- RunUMAP(tcells, dims = 1:14) #scegleindo solo le pc significative

pdf("umap.pdf")
plot_1 <- DimPlot(tcells, reduction = "umap", label = T)
#per aggiungere il filtro per condizione per poi capie se alcuni cluster sono più rielvanti per tipo, si crea la variabile con la clusterizzazionw
tcells$Condition <- tcells$`tcells@active.ident`
plot_2 <- DimPlot(tcells, reduction = "umap", label = T, group.by = 'Condition')
plot_1+plot_2
dev.off()

#per vedere un particolare gene nella umap
FeaturePlot(tcells, features = 'CXCL13', label = T)


#cercare i marker per ogni cluster per arrivare 
tcells.markers <- FindAllMarkers(tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#si possono usare i valori di default, ma sipossono cambiare:
#onmly.pos =true only genes that are upregulated to in one cluster in respect to all the other, is easier to explain. It is possible to use the opposite, but it says that in that cluster the gean is down respect to other
#min.pct= see expression of that marker at leat in 25% of samples in the cluster, if lower it is not a marker of the cluster; come minino si può settare fino a 10%. questo valore cambia in base alla tecnologia perchè cambia come l'RNA è catturato
#logf-threr= si può mettere 0 e poi cambiare per tenere anche geni non differenzialmente espressi e poi si filtra, ma conviene farlo prima

#per vedere i cluster in base ad una divisione per un metadata di interesse, LR0linear regression perchè è una condizione pèer comparare i marcatori trovati senza correzione e per la variabile di interesse per confrontarli con latent.vat
tcells.markers_cond <- FindAllMarkers(tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use= "LR", latent.vars = "Condition")

#per vedere i risultati basati
tcells.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
tcells.markers_cond %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

#fare la tabella per i top10 features per cluster
tcells.markers %>% group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
#fare heat map
pdf("heatmap.pdf")
DoHeatmap(tcells, features = top10$gene) + NoLegend()
dev.off()
#per esportare la tablla con i geni. occhio, i marker gene could be marker of multiple cluster. The list exported rinomina i geni con unfercore il cluster in base al cluster se sono in più cluster con _numero di cluster nella prima colonna ,meglio tenere l'ultima colonna perchè c'è il nome originale del gene. Occhio che i nomi delle colonne sono shiftati di uno. Tabella avg_foldchanfùge come il cluster cambia rispetto gli altri, pct1 la percentuale delle cellule che lo esprimono su pct.1, pct.2, pvalue occhio che ci sono anche dei geni non significativi, se si vogliono solo i sign si può mettere un filtro in findallmarkers)
write.table(tcells.markers, "tcells.markers.txt", sep="\t")

#compare two cluster e poi usa gli altri metodi per vedere la tabella
tcells.markers_clust <- FindMarkers(tcells, ident.1 = "0", ident.2 = c("1","2"), only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Plot of marker genes

pdf("VlnPlot_FOXP3.pdf")
VlnPlot(tcells, features = 'FOXP3')
dev.off()

pdf("FeaturePlot_CXCL13.pdf")
FeaturePlot(tcells, features = 'CXCL13', label = T)
dev.off()


pdf("DotPlot_IL1R2.pdf")
DotPlot(tcells, features = 'IL1R2')#per cluster la percentuale di dell che esprimono il marker e il level of expresspion
dev.off()

#Responder vs Non responder --> all cells, without considering clustering analysis, si assegnano i livelli di ident
Idents(tcells)

levels(tcells@meta.data$`tcells@active.ident`)[1] <-"Non_responder"
levels(tcells@meta.data$`tcells@active.ident`)[2] <-"Responder"
table(tcells@meta.data$`tcells@active.ident`)

Idents(tcells)<- tcells@meta.data$`tcells@active.ident`#change the slot in the object of the identity per cambiare il livello
R_vs_NR<- FindMarkers(tcells, ident.1 = 'Responder', ident.2 = 'Non_responder')
head(R_vs_NR)
write.table(R_vs_NR, file='R_vs_NR.txt', sep="\t")

#Responder vs Non responder--> cluster 0
Idents(tcells)<- tcells$RNA_snn_res.0.4#cambio ident sui cluster
tcells_c0<- subset(tcells, idents = '0')#subsect l'object per un cluster#vi set the identity of only cluster 0, ad esempio se si vuole suddividere ulteriormente, es è il cluster CD8 e si può rianalizzare, si può fare un cluster alla votla
Idents(tcells_c0)<-tcells_c0@meta.data$`tcells@active.ident`#per rifare la suddivisione su un cluster in un particolare cluater devo riportare l'ident sulla suddivisione per gruppi di pts
Idents(tcells_c0)
c0_R_vs_NR<- FindMarkers(tcells_c0, ident.1 = 'Responder', ident.2 = 'Non_responder')#creo la variabile con i geni diff expr fra resp e nn respond incluster 0
head(c0_R_vs_NR)
write.table(c0_R_vs_NR, file='c0_R_vs_NR.txt', sep="\t")
DoHeatmap(tcells_c0, features = top10$gene) + NoLegend()
