BiocManager::install(c("limma", "edgeR", "EnhancedVolcano", "ComplexHeatmap", "tximport", "readxl", "dplyr", "EnsDb.Hsapiens.v86"), force = TRUE)

library(tximport)
quant_files <- list.files(
  "/Users/andrea/rnaseq_analysis/transcript_quant_MCF-7", 
  pattern = "quant.sf", 
  recursive = TRUE, 
  full.names = TRUE
)
sample_names <- basename(dirname(quant_files))

txi <- tximport(
  files = setNames(quant_files, sample_names),
  type = "salmon",
  txOut = TRUE
)

head(txi$counts)

library(readxl)
sample_annotation <- read_excel("/Users/andrea/rnaseq_analysis/sample_annotation_MCF-7.xlsx")
class(sample_annotation) <- "data.frame"
rownames(sample_annotation) <- sample_annotation$Sample

raw_counts <- txi$counts

library(limma)
library(edgeR)
# --- Analisi Limma-Voom per miRNA e combinazione risultati ---

# Lista di miRNA unici
miRNAs <- unique(sample_annotation$miRNA)

degs_matrix_gene_list <- list()
degs_limma_list <- list()
degs_matrix_list <- list()

for (current_mirna in miRNAs) {
  # Filtra annotazioni e conteggi per il miRNA corrente
  subset_annotation <- sample_annotation[sample_annotation$miRNA == current_mirna, ]
  subset_counts <- raw_counts[, rownames(subset_annotation)]
  
  # Ricostruisci l'oggetto DGEList e la matrice di disegno
  dge <- DGEList(counts = subset_counts)
  experiment <- model.matrix(~ 0 + Condition, data = subset_annotation)
  
  # Filtraggio ed normalizzazione (come prima, ma con i subset)
  keep <- filterByExpr(dge, design = experiment)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  norm_dge <- voom(dge, experiment)
  
  # Analisi Limma-Voom (come prima, ma con i subset)
  fit <- lmFit(norm_dge, experiment)
  contrast_matrix <- makeContrasts(ConditionOE - ConditionCTL, levels = experiment)
  cont_fit <- contrasts.fit(fit, contrasts = contrast_matrix)
  ebayes_fit <- eBayes(cont_fit)
  
  # Estrazione risultati
  degs <- topTable(ebayes_fit, coef=1, n=Inf)
  degs_limma <- degs[degs$adj.P.Val < 0.05 & abs(degs$logFC) >= log2(1.5),]
  
  library(ComplexHeatmap)
  library(org.Hs.eg.db)
  library(EnsDb.Hsapiens.v86)
  library(dplyr)
  degs_matrix <- subset_counts[rownames(degs_limma),]
  degs_matrix <- t(scale(t(degs_matrix)))
  
  # Raggruppamento delle righe per ID inserendo i valori medi:
  
  # Converti degs_matrix in un data.frame
  degs_matrix <- as.data.frame(degs_matrix) 
  
  # Estrai gli ID Ensembl dai nomi delle righe di degs_matrix
  ensembl_ids <- sapply(strsplit(rownames(degs_matrix), "\\."), function(x) x[1])
  
  # Aggiungi una colonna ENSEMBL a degs_matrix
  degs_matrix$ENSEMBL <- ensembl_ids
  
  # Raggruppa per ENSEMBL e calcola la media per tutte le altre colonne
  degs_matrix_gene <- degs_matrix %>%
    group_by(ENSEMBL) %>%
    summarize(across(everything(), mean))
  degs_matrix_gene <- data.frame(degs_matrix_gene)
  
  # Estrazione della colonna ENSEMBL per usarla come nomi delle righe
  gene_names <- degs_matrix_gene$ENSEMBL
  
  # Conversione di degs_matrix_gene in un data.frame standard
  degs_matrix_gene <- as.data.frame(degs_matrix_gene)
  
  # Impostazione dei nomi delle righe di degs_matrix_gene
  rownames(degs_matrix_gene) <- gene_names
  
  # Rimozione della colonna ENSEMBL (ora che i nomi delle righe sono impostati)
  degs_matrix_gene$ENSEMBL <- NULL
  
  # Gestione dei geni duplicati causati dal mapping id-gene:
  
  # Salvamento dei row names originali (ENSEMBLTRANS IDs)
  original_rownames <- rownames(degs_matrix_gene)
  
  symbols <- mapIds(EnsDb.Hsapiens.v86, keys = original_rownames, column = c('SYMBOL'), keytype = 'TXID')
  symbols <- symbols[!is.na(symbols)]
  
  # Creazione di un vettore per i nuovi nomi di riga, inizializzato con i simboli mappati
  new_rownames <- symbols[original_rownames]
  
  # Gestione dei nomi di riga duplicati
  duplicati <- duplicated(new_rownames) | duplicated(new_rownames, fromLast = TRUE) # Trova tutti i duplicati, la prima e l'ultima occorrenza
  new_rownames_duplicati <- new_rownames # Crea una copia per modificarla
  
  # Aggiunta di ENSEMBLTRANS ID ai nomi di riga duplicati per renderli unici
  for(i in which(duplicati)) {
    new_rownames_duplicati[i] <- paste0(new_rownames[i], "_", original_rownames[i])
  }
  
  # Impostazione deii nuovi nomi di riga (unici) nella matrice
  rownames(degs_matrix_gene) <- new_rownames_duplicati
  
  # Rimozione dei geni non riconosciuti
  righe_da_rimuovere <- startsWith(rownames(degs_matrix_gene), "NA_")
  degs_matrix_gene <- degs_matrix_gene[!righe_da_rimuovere, ]
  degs_matrix_gene <- as.matrix(degs_matrix_gene)
  
  degs_limma_list[[current_mirna]] <- degs_limma
  degs_matrix_list[[current_mirna]] <- degs_matrix
  degs_matrix_gene_list[[current_mirna]] <- degs_matrix_gene
  
}

degs_matrix_gene_miR_455_3p <- degs_matrix_gene_list[["miR-455-3p"]]
degs_matrix_gene_miR_874 <- degs_matrix_gene_list[["miR-874"]]
degs_matrix_gene_miR_9 <- degs_matrix_gene_list[["miR-9"]]

degs_limma_miR_455_3p <- degs_limma_list[["miR-455-3p"]]
degs_limma_miR_874 <- degs_limma_list[["miR-874"]]
degs_limma_miR_9 <- degs_limma_list[["miR-9"]]

degs_matrix_miR_455_3p <- degs_matrix_list[["miR-455-3p"]]
degs_matrix_miR_874 <- degs_matrix_list[["miR-874"]]
degs_matrix_miR_9 <- degs_matrix_list[["miR-9"]]

# EnhancedVolcano(degs_limma, lab = rownames(degs_limma), x = 'logFC', y = 'adj.P.Val', title = "Overexpression vs Control", pCutoff = 0.05, FCcutoff = log2(1.5))

# Heatmap(degs_matrix_gene, name = "Expression", row_names_gp = gpar(fontsize = 4))