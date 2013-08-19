#########################################
# gene_expression_funnet.R - analysis of genes expression based on classied subjects with Funnet
# Author: Hoai Tuong Nguyen
# Created: 05/07/2013
# Modified: 12/07/2013
# CMD: R --no-save --no-restore --slave -f gene_expression_funnet.R "--args input.dir='?' output.dir='?'"
#########################################

#First read and parse the arguments for on the fly run
args=(commandArgs(TRUE))

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}


library(mc)
library(FunNet)


#default working path for inner run
#setwd("MetaCARDIS/")
#setwd("C:/Users/TUONG/Documents/metacardis/7.scripts/")

#default input and output paths
if (!exists("input.dir"))
  input.dir="../1.data/metacardis"
if (!exists("output.dir"))
  output.dir="../3.results/CLASSIFICATION/Subjects/2013-06-17"

dir.create(output.dir, showWarnings = FALSE)


#Modified version of Funnet
fun<-function (wd = "", org = "hsa", two.lists = TRUE, up.frame = NULL, 
               down.frame = NULL, genes.frame = NULL, restrict = FALSE, 
               ref.list = NULL, logged = FALSE, discriminant = FALSE, go.bp = TRUE, 
               go.cc = TRUE, go.mf = TRUE, kegg = TRUE, annot.method = "specificity", 
               annot.details = TRUE, direct = FALSE, enriched = TRUE, fdr = NA, 
               build.annot.net = TRUE, coexp.matrix = NULL, coexp.method = "spearman", 
               estimate.th = FALSE, hard.th = NA, soft.th = NA, topological = FALSE, 
               keep.sign = FALSE, level = NA, annot.clust.method = "umilds", 
               annot.prox.measure = "unilat.pond.norm.mean", test.recovery = FALSE, 
               test.robust = FALSE, replace.annot = NA, random.annot = FALSE, 
               build.gene.net = FALSE, gene.clust.method = "hclust", gene.net.details = FALSE, 
               gene.clusters = NA, alpha = 0.05, RV = 0.9, sigma = NA, keep.rdata = FALSE, 
               zip = TRUE,fileprefix) 
{
  if (org == "HS") {
    org <- "hsa"
  }
  if (org == "MM") {
    org <- "mmu"
  }
  if (org == "RN") {
    org <- "rno"
  }
  if (org == "SC") {
    org <- "sce"
  }
  parameter.list <- list(analysis.date = date(), package.version = .funnet.version, 
                         org = org, annot.date = annot.date, annot.method = annot.method, 
                         annot.clust.method = annot.clust.method, annot.prox.measure = annot.prox.measure, 
                         direct = direct, enriched = enriched, fdr = fdr, two.lists = two.lists, 
                         restrict = restrict, go.bp = go.bp, go.cc = go.cc, go.mf = go.mf, 
                         kegg = kegg, discriminant = discriminant, logged = logged, 
                         annot.details = annot.details, estimate.th = estimate.th, 
                         hard.th = hard.th, soft.th = soft.th, coexp.method = coexp.method, 
                         topological = topological, keep.sign = keep.sign, build.annot.net = build.annot.net, 
                         level = level, test.recovery = test.recovery, test.robust = test.robust, 
                         replace.annot = replace.annot, random.annot = random.annot, 
                         build.gene.net = build.gene.net, gene.clust.method = gene.clust.method, 
                         gene.clusters = gene.clusters, gene.net.details = gene.net.details, 
                         alpha = alpha, RV = RV, sigma = sigma, keep.rdata = keep.rdata, 
                         zip = zip)
  .check.parameters(parameter.list, coexp.matrix, up.frame, 
                    down.frame, ref.list, genes.frame)
  if (discriminant) {
    two.lists <- TRUE
    restrict <- TRUE
  }
  if (!is.null(up.frame) & !is.null(down.frame) & is.null(genes.frame)) {
    two.lists <- TRUE
  }
  if (!is.null(genes.frame) & (is.null(up.frame) | is.null(down.frame))) {
    two.lists <- FALSE
    discriminant <- FALSE
  }
  cat(paste("\n\tFunNet started at: ", date(), sep = ""))
  cat(paste("\n\t\tUsing annotations updated on: ", annot.date, 
            sep = ""))
  if (wd != "") {
    setwd(wd)
  }
  results.dir <- paste(fileprefix,"_", format(Sys.time(), "%Y_%b_%d_%H-%M-%S"), 
                       sep = "")
  dir.create(paste(getwd(), "/", results.dir, sep = ""))
  dir.create(paste(getwd(), "/", results.dir, "/html", sep = ""))
  dir.create(paste(getwd(), "/", results.dir, "/images", sep = ""))
  try(write.table(as.matrix(print(parameter.list)), col.names = F, 
                  file = paste(getwd(), "/", results.dir, "/", "parameters_list.txt", 
                               sep = ""), sep = "\t"))
  wd <- getwd()
  locus.name <- annot.base[[org]]$locus.name[, 1:2]
  locus.symbol <- annot.base[[org]]$locus.name[, c(1, 3)]
  rownames(locus.name) <- locus.name[, 1]
  rownames(locus.symbol) <- locus.symbol[, 1]
  if (two.lists) {
    up.down <- .filter.genes(up.frame = up.frame, down.frame = down.frame, 
                             two.lists = TRUE, locus.name = locus.name, logged = logged)
    up.frame <- up.down$up.frame
    down.frame <- up.down$down.frame
    rm(up.down)
    if (discriminant) {
      ref.list <- c(as.character(up.frame[, 1]), as.character(down.frame[,1]))
    }
    else if (restrict & !discriminant) {
      ref.list <- .filter.genes(restrict = TRUE, ref.list = ref.list, 
                                locus.name = locus.name)
    }
    else if (!restrict & !discriminant) {
      ref.list <- NULL
    }
  }
  else {
    genes.frame <- .filter.genes(genes.frame = genes.frame, 
                                 two.lists = FALSE, locus.name = locus.name)
    if (restrict) {
      ref.list <- .filter.genes(restrict = TRUE, ref.list = ref.list, 
                                locus.name = locus.name)
    }
    else {
      ref.list <- NULL
    }
  }
  cat(paste("\n\tSaving start-up environment... ", format(Sys.time(), 
                                                          "%X"), sep = ""))
  save(up.frame, down.frame, ref.list, genes.frame, parameter.list, 
       coexp.matrix, locus.name, file = paste(getwd(), "/", 
                                              results.dir, "/", "start-up_environment.RData", sep = ""), 
       compress = T)
  save(up.frame, down.frame, ref.list, genes.frame, parameter.list, 
       coexp.matrix, locus.name, file = paste(getwd(), "/", 
                                              "start-up_environment.RData", sep = ""), compress = T)
  if (estimate.th) {
    datas <- NULL
    if (!is.null(up.frame) & !is.null(down.frame)) {
      datas <- rbind(up.frame, down.frame)
    }
    if (!is.null(genes.frame) & (is.null(up.frame) | is.null(down.frame))) {
      datas <- genes.frame
    }
    rownames(datas) <- datas[, 1]
    datas <- datas[, 2:ncol(datas)]
    cat("\n\tHard thresholding...\n")
    try(hard.th <- .PickHardThreshold(datExpr1 = t(datas), 
                                      coexp.method = coexp.method))
    try(write.table(hard.th$tablou, file = paste(getwd(), 
                                                 "/", results.dir, "/", coexp.method, "_hard_threshold.txt", 
                                                 sep = ""), append = FALSE, col.names = TRUE, , row.names = F, 
                    sep = "\t"))
    try(write(paste("\n\nHard threshold estimate: ", hard.th$estimate, 
                    "\n", sep = ""), file = paste(getwd(), "/", results.dir, 
                                                  "/", coexp.method, "_hard_threshold.txt", sep = ""), 
              append = TRUE))
    try(write(paste("Do not trust this automated estimation without checking it!\n", 
                    "Do not hesitate to select another threshold depending on the associated connectivity values.\n", 
                    "Then please restart FunNet interaction analysis with your selected threshold.\n", 
                    sep = ""), file = paste(getwd(), "/", results.dir, 
                                            "/", coexp.method, "_hard_threshold.txt", sep = ""), 
              append = TRUE))
    cat("\n\tSoft thresholding...\n")
    try(soft.th <- .PickSoftThreshold(datExpr1 = t(datas), 
                                      coexp.method = coexp.method))
    try(write.table(soft.th$tablou, file = paste(getwd(), 
                                                 "/", results.dir, "/", coexp.method, "_soft_threshold.txt", 
                                                 sep = ""), append = FALSE, col.names = TRUE, , row.names = F, 
                    sep = "\t"))
    try(write(paste("\n\nSoft threshold estimate: ", soft.th$estimate, 
                    "\n", sep = ""), file = paste(getwd(), "/", results.dir, 
                                                  "/", coexp.method, "_soft_threshold.txt", sep = ""), 
              append = TRUE))
    try(write(paste("Do not trust this automated estimation without checking it!\n", 
                    "Do not hesitate to select another threshold depending on the associated connectivity values.\n", 
                    "Then please restart FunNet interaction analysis with your selected threshold.\n", 
                    sep = ""), file = paste(getwd(), "/", results.dir, 
                                            "/", coexp.method, "_soft_threshold.txt", sep = ""), 
              append = TRUE))
    print("Estimation of the co-expression threshold finished!")
    print("Please restart FunNet with your chosen threshold.")
    if (!keep.rdata) {
      try(unlink(paste(getwd(), "/", results.dir, "/", 
                       list.files(path = paste(getwd(), "/", results.dir, 
                                               "/", sep = ""), pattern = "[:print:]*.RData"), 
                       sep = ""), recursive = TRUE))
    }
    if (zip) {
      try(unlink(paste(getwd(), "/", results.dir, "/html", 
                       sep = ""), recursive = TRUE))
      try(unlink(paste(getwd(), "/", results.dir, "/images", 
                       sep = ""), recursive = TRUE))
      try(system(command = paste("zip -r9q ", results.dir, 
                                 ".zip ", "./", results.dir, "/*", sep = "")))
      try(unlink(paste(getwd(), "/", results.dir, sep = ""), 
                 recursive = TRUE))
    }
    options(show.error.messages = FALSE)
    stop()
  }
  if (is.null(coexp.matrix) & (build.annot.net | build.gene.net)) {
    cat(paste("\n\tComputing co-expression matrix... ", format(Sys.time(), 
                                                               "%X"), sep = ""))
    datas <- NULL
    if (!is.null(up.frame) & !is.null(down.frame)) {
      datas <- rbind(up.frame, down.frame)
    }
    if (!is.null(genes.frame) & (is.null(up.frame) | is.null(down.frame))) {
      datas <- genes.frame
    }
    rownames(datas) <- datas[, 1]
    datas <- datas[, 2:ncol(datas)]
    
    if (coexp.method %in% c("spearman", "pearson", "kendall")) {
      coexp.matrix <- rcorr(t(datas), type = coexp.method)$r
      
    }
    else if (coexp.method == "euclid") {
      coexp.matrix <- 1 - (as.matrix(dist(datas))/max(dist(datas), 
                                                      na.rm = TRUE))
    }
    try(save(coexp.matrix, file = paste(getwd(), "/", results.dir, 
                                        "/", "brut_coexp_matrix.RData", sep = ""), compress = T))
  }
  if (build.annot.net | build.gene.net) {
    sign.matrix <- coexp.matrix/abs(coexp.matrix)
    coexp.matrix <- abs(coexp.matrix)
    if (annot.clust.method %in% c("umilds", "spectral")) {
      if (!is.na(hard.th)) {
        coexp.matrix[coexp.matrix >= hard.th] <- 1
        coexp.matrix[coexp.matrix < hard.th] <- 0
      }
      else if (!is.na(soft.th)) {
        coexp.matrix <- coexp.matrix^soft.th
      }
      if (topological) {
        coexp.matrix <- 1 - .TOMdist(adjmat1 = coexp.matrix)
      }
      else if (keep.sign) {
        coexp.matrix <- coexp.matrix * sign.matrix
      }
    }
    else if (annot.clust.method == "ucknn" & !is.na(hard.th)) {
      coexp.matrix[coexp.matrix < hard.th] <- 0
    }
    try(save(coexp.matrix, sign.matrix, parameter.list, file = paste(getwd(), 
                                                                     "/", results.dir, "/", "gene_adj_matrix.RData", sep = "")))
  }
  if (kegg) {
    terms.name <- KEGG.terms.name
    rownames(terms.name) <- terms.name[, 1]
    file.annot <- annot.base[[org]]$KEGG.file.annot
    taxoname <- "KEGG"
    if (two.lists) {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = FALSE, results.dir = results.dir, 
                 alpha = alpha, locus.name = locus.name, annot.clust.method = annot.clust.method, 
                 up.frame = up.frame, down.frame = down.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = NA, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
    else {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = FALSE, results.dir = results.dir, 
                 annot.clust.method = annot.clust.method, alpha = alpha, 
                 locus.name = locus.name, genes.frame = genes.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = NA, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
  }
  go.name <- c("GO Biological Process", "GO Cellular Component", 
               "GO Molecular Function")
  terms.name <- GO.terms.name
  rownames(terms.name) <- terms.name[, 1]
  if (go.bp) {
    file.annot <- annot.base[[org]]$GO.DIR.BP.file.annot
    taxoname <- go.name[1]
    if (two.lists == TRUE) {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 alpha = alpha, locus.name = locus.name, annot.clust.method = annot.clust.method, 
                 up.frame = up.frame, down.frame = down.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
    else {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 annot.clust.method = annot.clust.method, alpha = alpha, 
                 locus.name = locus.name, genes.frame = genes.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
  }
  if (go.cc) {
    file.annot <- annot.base[[org]]$GO.DIR.CC.file.annot
    taxoname <- go.name[2]
    if (two.lists == TRUE) {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 alpha = alpha, locus.name = locus.name, annot.clust.method = annot.clust.method, 
                 up.frame = up.frame, down.frame = down.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
    else {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 annot.clust.method = annot.clust.method, alpha = alpha, 
                 locus.name = locus.name, genes.frame = genes.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
  }
  if (go.mf) {
    file.annot <- annot.base[[org]]$GO.DIR.MF.file.annot
    taxoname <- go.name[3]
    if (two.lists == TRUE) {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 alpha = alpha, locus.name = locus.name, annot.clust.method = annot.clust.method, 
                 up.frame = up.frame, down.frame = down.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
    else {
      .main.loop(file.annot = file.annot, taxoname = taxoname, 
                 annot.method = annot.method, terms.name = terms.name, 
                 direct = direct, fdr = fdr, go = TRUE, results.dir = results.dir, 
                 annot.clust.method = annot.clust.method, alpha = alpha, 
                 locus.name = locus.name, genes.frame = genes.frame, 
                 restrict = restrict, ref.list = ref.list, annot.details = annot.details, 
                 level = level, build.annot.net = build.annot.net, 
                 test.recovery = test.recovery, test.robust = test.robust, 
                 replace.annot = replace.annot, locus.symbol = locus.symbol, 
                 annot.prox.measure = annot.prox.measure, coexp.matrix = coexp.matrix, 
                 parameter.list = parameter.list, org = org, gene.net.details = gene.net.details, 
                 RV = RV, sigma = sigma, random.annot = random.annot)
    }
  }
  if (build.gene.net & !is.null(coexp.matrix)) {
    try(clusters <- .build.coexp.net(coexp.matrix = coexp.matrix, 
                                     locus.name = locus.name, locus.symbol = locus.symbol, 
                                     gene.clust.method = gene.clust.method, gene.clusters = gene.clusters))
    try(save(clusters, parameter.list, file = paste(getwd(), 
                                                    "/", results.dir, "/", "co-expression_clusters.RData", 
                                                    sep = ""), compress = T))
    net.matrix <- coexp.matrix
    rownames(net.matrix) <- locus.symbol[rownames(net.matrix), 
                                         2]
    colnames(net.matrix) <- locus.symbol[colnames(net.matrix), 
                                         2]
    try(.cyto.sym(net.matrix = net.matrix, file.net = paste(getwd(), 
                                                            "/", results.dir, "/", "co-expression_net.txt", sep = ""), 
                  diagonal = FALSE, thresh = NULL))
    rm(net.matrix)
    try(centrality <- .genes.centrality(adj.matrix = coexp.matrix, 
                                        clusters = clusters, taxoname = taxoname, locus.symbol = locus.symbol, 
                                        results.dir = results.dir, coexp = TRUE))
    if (two.lists) {
      up.down <- rbind(matrix(1, nrow(up.frame), 1), matrix(0, 
                                                            nrow(down.frame), 1))
      rownames(up.down) <- c(as.character(up.frame[, 1]), 
                             as.character(down.frame[, 1]))
      try(write.table(cbind(rownames(clusters$gene.connect), 
                            as.vector(locus.symbol[rownames(clusters$gene.connect), 
                                                   2]), as.vector(locus.name[rownames(clusters$gene.connect), 
                                                                             2]), up.down, clusters$gene.connect, centrality[rownames(clusters$gene.connect), 
                                                                                                                             ]), file = paste(getwd(), "/", results.dir, 
                                                                                                                                              "/", "co-expression_net_info.txt", sep = ""), 
                      sep = "\t", col.names = c("geneid", "symbol", 
                                                "name", "up(1)_down(0)", colnames(clusters$gene.connect), 
                                                colnames(centrality)), row.names = F))
      try(rm(clusters, up.down, centrality))
    }
    else {
      try(write.table(cbind(rownames(clusters$gene.connect), 
                            as.vector(locus.symbol[rownames(clusters$gene.connect), 
                                                   2]), as.vector(locus.name[rownames(clusters$gene.connect), 
                                                                             2]), clusters$gene.connect, centrality[rownames(clusters$gene.connect), 
                                                                                                                    ]), file = paste(getwd(), "/", results.dir, 
                                                                                                                                     "/", "co-expression_net_info.txt", sep = ""), 
                      sep = "\t", col.names = c("geneid", "symbol", 
                                                "name", colnames(clusters$gene.connect), colnames(centrality)), 
                      row.names = F))
      try(rm(clusters, centrality))
    }
    cat(paste("\n\tCo-expression net building finished... ", 
              date(), sep = ""))
    rm()
  }
  if (!keep.rdata) {
    try(unlink(paste(getwd(), "/", results.dir, "/", list.files(path = paste(getwd(), 
                                                                             "/", results.dir, "/", sep = ""), pattern = "[:print:]*.RData"), 
                     sep = ""), recursive = TRUE))
  }
  if (zip) {
    try(system(command = paste("zip -r9q ", results.dir, 
                               ".zip ", "./", results.dir, "/*", sep = "")))
    try(unlink(paste(getwd(), "/", results.dir, sep = ""), 
               recursive = TRUE))
  }
  cat(paste("\n\tEnd  of treatment at: ", date(), "\n", sep = ""))
  rm()
}

#function for analysis of genes expression based on classied subjects with Funnet
genes.updown.funnet<-function(genes,prode.entrez,high,low,stars,fileprefix,output.dir){
  m<-mean(as.numeric(as.matrix(genes)))
  res<-as.data.frame(rbind(c("Probe_ID","High_Risk_Group","Low_Risk_Group"),cbind(colnames(genes),ifelse(colMeans(genes[high,])>m,1,0),ifelse(colMeans(genes[low,])>m,1,0))))
  #high.up
  ge.high.up<-res[rownames(res)%in%stars & res[,2]==1,]
  high.up<-t(genes[,colnames(genes)%in%as.character(ge.high.up[,1])])
  rownames(high.up)<-prode.entrez[prode.entrez[,1]%in%rownames(high.up),2]
  high.up<-high.up[na.omit(rownames(high.up)),]
  
  #high.up
  ge.high.up<-res[rownames(res)%in%stars & res[,2]==1,]
  high.up<-t(genes[,colnames(genes)%in%as.character(ge.high.up[,1])])
  rownames(high.up)<-prode.entrez[prode.entrez[,1]%in%rownames(high.up),2]
  high.up<-high.up[na.omit(rownames(high.up)),]
  
  #high.up
  ge.high.down<-res[rownames(res)%in%stars & res[,2]==0,]
  high.down<-t(genes[,colnames(genes)%in%as.character(ge.high.down[,1])])
  rownames(high.down)<-prode.entrez[prode.entrez[,1]%in%rownames(high.down),2]
  high.down<-high.down[na.omit(rownames(high.down)),]
  
  #low.up
  ge.low.up<-res[rownames(res)%in%stars & res[,3]==1,]
  low.up<-t(genes[,colnames(genes)%in%as.character(ge.low.up[,1])])
  rownames(low.up)<-prode.entrez[prode.entrez[,1]%in%rownames(low.up),2]
  low.up<-low.up[na.omit(rownames(low.up)),]  
  #low.down
  ge.low.down<-res[rownames(res)%in%stars & res[,3]==0,]
  low.down<-t(genes[,colnames(genes)%in%as.character(ge.low.down[,1])])
  rownames(low.down)<-prode.entrez[prode.entrez[,1]%in%rownames(low.down),2]
  low.down<-low.down[na.omit(rownames(low.down)),]
  
  high.refs<-c(rownames(high.up),rownames(high.down))
  low.refs<-c(rownames(low.up),rownames(low.down))

  high.up.frame<-as.data.frame(high.up)
  colnames(high.up.frame)[1]<-"GeneID"
  high.down.frame<-as.data.frame(high.down)
  colnames(high.down.frame)[1]<-"GeneID"
  
  low.up.frame<-as.data.frame(low.up)
  colnames(low.up.frame)[1]<-"GeneID"  
  low.down.frame<-as.data.frame(low.down)
  colnames(low.down.frame)[1]<-"GeneID"
  

  write.table(high.up,file=sprintf("%s/%s.high.up.genes.txt",output.dir,fileprefix),col.names=F,row.names=T,sep="\t",quote=F)
  write.table(high.down,file=sprintf("%s/%s.high.down.genes.txt",output.dir,fileprefix),col.names=F,row.names=T,sep="\t",quote=F)
  write.table(low.up,file=sprintf("%s/%s.low.up.genes.txt",output.dir,fileprefix),col.names=F,row.names=T,sep="\t",quote=F)
  write.table(low.down,file=sprintf("%s/%s.low.down.genes.txt",output.dir,fileprefix),col.names=F,row.names=T,sep="\t",quote=F)
  write.table(high.refs,file=sprintf("%s/%s.high.refs.genes.txt",output.dir,fileprefix),col.names=F,row.names=F,sep="\t",quote=F)
  write.table(low.refs,file=sprintf("%s/%s.low.refs.genes.txt",output.dir,fileprefix),col.names=F,row.names=F,sep="\t",quote=F)
  write.table(res,file=sprintf("%s/updown.genes.%s.txt",output.dir,fileprefix),col.names=F,row.names=F,sep="\t",quote=F)
  
  
  high.up.frame<-as.data.frame(read.table(sprintf("%s/%s.high.up.genes.txt",output.dir,fileprefix),sep="\t",header=F))
  high.down.frame<-as.data.frame(read.table(sprintf("%s/%s.high.down.genes.txt",output.dir,fileprefix),sep="\t",header=F))
  low.up.frame<-as.data.frame(read.table(sprintf("%s/%s.low.up.genes.txt",output.dir,fileprefix),sep="\t",header=F))
  low.down.frame<-as.data.frame(read.table(sprintf("%s/%s.low.down.genes.txt",output.dir,fileprefix),sep="\t",header=F))
    
 
  fun(org="HS", two.lists=TRUE, up.frame=high.up.frame, down.frame=high.down.frame,
      genes.frame=NULL, restrict=TRUE, ref.list=data.frame(unique(ge.prode.entrez[,2])), logged=FALSE,
      discriminant=TRUE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,
      annot.method="specificity", annot.details=TRUE,
      direct=FALSE, enriched=TRUE, fdr=NA, build.annot.net=TRUE,
      coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE,
      hard.th=0.8, soft.th=NA, topological = FALSE, keep.sign=FALSE, level=1,
      annot.clust.method="umilds", annot.prox.measure="dynamical",
      test.recovery=FALSE, test.robust=FALSE, replace.annot=NA,
      build.gene.net=TRUE, gene.clust.method="hclust", gene.net.details=TRUE,
      gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE,fileprefix=sprintf("%s.high.genes",fileprefix))


  fun(org="HS", two.lists=TRUE, up.frame=low.up.frame, down.frame=low.down.frame,
      genes.frame=NULL, restrict=TRUE, ref.list=data.frame(unique(ge.prode.entrez[,2])), logged=FALSE,
      discriminant=TRUE, go.bp=TRUE, go.cc=TRUE, go.mf=TRUE, kegg=TRUE,
      annot.method="specificity", annot.details=TRUE,
      direct=FALSE, enriched=TRUE, fdr=NA, build.annot.net=TRUE,
      coexp.matrix=NULL, coexp.method="spearman", estimate.th=FALSE,
      hard.th=0.8, soft.th=NA, topological = FALSE, keep.sign=FALSE, level=1,
      annot.clust.method="umilds", annot.prox.measure="dynamical",
      test.recovery=FALSE, test.robust=FALSE, replace.annot=NA,
      build.gene.net=TRUE, gene.clust.method="hclust", gene.net.details=TRUE,
      gene.clusters=NA, alpha=0.05, RV=0.90, sigma=NA, keep.rdata=FALSE, zip=TRUE,fileprefix=sprintf("%s.low.genes",fileprefix))

  
  return(res)
}



#all gene expression matrix
ge.data<-read.table(sprintf("%s/data.d.all.txt",input.dir),sep="\t",header=T)
ge.prode.entrez<-ge.data[,2:3]

#T0 gene expression matrix
ge.data.t0<-get(load(sprintf("%s/GE_t0.RData",input.dir)))

#No NA - T0 gene expression matrix
#ge.data.t0.nona<-na.omit(ge.data.t0)
n.na<-sapply(1:ncol(ge.data.t0),function(x) sum(is.na(ge.data.t0[,x]))) 
ge.data.t0.nona<-ge.data.t0[,which(n.na==0)]



##############
# LINEAR MODELS
##############

#genes by groups
ge.lm.classified.adfm.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
ge.lm.classified.wfmh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lm.classified.adwh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
ge.lm.classified.adfmh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lm.classified.adh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir))
ge.lm.classified.adbmifma.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir))

#genes by groups after testing
ge.lm.classified.adfm.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_auto-test.csv",output.dir),sep=";")
ge.lm.classified.wfmh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_auto-test.csv",output.dir),sep=";")
ge.lm.classified.adwh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test.csv",output.dir),sep=";")
ge.lm.classified.adfmh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test.csv",output.dir),sep=";")
ge.lm.classified.adh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_Height2_auto-test.csv",output.dir),sep=";")
ge.lm.classified.adbmifma.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto-test.csv",output.dir),sep=";")

#genes significantly different by groups
ge.lm.classified.adfm.t0.autotest.star<-ge.lm.classified.adfm.t0.autotest[ge.lm.classified.adfm.t0.autotest[,4]!="",]
ge.lm.classified.wfmh2.t0.autotest.star<-ge.lm.classified.wfmh2.t0.autotest[ge.lm.classified.wfmh2.t0.autotest[,4]!="",]
ge.lm.classified.adwh2.t0.autotest.star<-ge.lm.classified.adwh2.t0.autotest[ge.lm.classified.adwh2.t0.autotest[,4]!="",]
ge.lm.classified.adfmh2.t0.autotest.star<-ge.lm.classified.adfmh2.t0.autotest[ge.lm.classified.adfmh2.t0.autotest[,4]!="",]
ge.lm.classified.adh2.t0.autotest.star<-ge.lm.classified.adh2.t0.autotest[ge.lm.classified.adh2.t0.autotest[,4]!="",]
ge.lm.classified.adbmifma.t0.autotest.star<-ge.lm.classified.adbmifma.t0.autotest[ge.lm.classified.adbmifma.t0.autotest[,4]!="",]






#analysis of genes expression based on classied subjects with Funnet
updown.genes.ge.lm.classified.wfmh2.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                            prode.entrez=ge.prode.entrez,
                                                            high=which(ge.lm.classified.adfm.t0[,4]==0),
                                                            low=which(ge.lm.classified.adfm.t0[,4]==1),
                                                            stars=ge.lm.classified.adfm.t0.autotest.star[,1],
                                                            fileprefix="ge.lm.classified.adfm.t0",
                                                            output.dir=output.dir)
updown.genes.ge.lm.classified.adwh2.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                            prode.entrez=ge.prode.entrez,
                                                            high=which(ge.lm.classified.adwh2.t0[,4]==0),
                                                            low=which(ge.lm.classified.adwh2.t0[,4]==1),
                                                            stars=ge.lm.classified.adwh2.t0.autotest.star[,1],
                                                            fileprefix="ge.lm.classified.adwh2.t0",
                                                            output.dir=output.dir)
updown.genes.ge.lm.classified.adfmh2.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                             prode.entrez=ge.prode.entrez,
                                                             high=which(ge.lm.classified.adfmh2.t0[,4]==0),
                                                             low=which(ge.lm.classified.adfmh2.t0[,4]==1),
                                                             stars=ge.lm.classified.adfmh2.t0.autotest.star[,1],
                                                             fileprefix="ge.lm.classified.adfmh2.t0",
                                                             output.dir=output.dir)
updown.genes.ge.lm.classified.adh2.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                           prode.entrez=ge.prode.entrez,
                                                           high=which(ge.lm.classified.adh2.t0[,4]==0),
                                                           low=which(ge.lm.classified.adh2.t0[,4]==1),
                                                           stars=ge.lm.classified.adh2.t0.autotest.star[,1],
                                                           fileprefix="ge.lm.classified.adh2.t0",
                                                           output.dir=output.dir)
updown.genes.ge.lm.classified.adbmifma.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                               prode.entrez=ge.prode.entrez,
                                                               high=which(ge.lm.classified.adbmifma.t0[,4]==0),
                                                               low=which(ge.lm.classified.adbmifma.t0[,4]==1),
                                                               stars=ge.lm.classified.adbmifma.t0.autotest.star[,1],
                                                               fileprefix="ge.lm.classified.adbmifma.t0",
                                                               output.dir=output.dir)








##############
# LOWESS
##############


#genes by groups
ge.lw.classified.adfm.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_points.csv",output.dir))
ge.lw.classified.adwh2.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_points.csv",output.dir))
ge.lw.classified.adfmh2.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lw.classified.adh2.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Height2_points.csv",output.dir))
ge.lw.classified.adbmifma.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_points.csv",output.dir))

#genes by groups after testing
ge.lw.classified.adfm.t0.autotest<-read.table(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_auto-test.csv",output.dir),sep=";")
ge.lw.classified.adwh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Waist_Height2_auto-test.csv",output.dir),sep=";")
ge.lw.classified.adfmh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Fat-mass_Height2_auto-test.csv",output.dir),sep=";")
ge.lw.classified.adh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_Height2_auto-test.csv",output.dir),sep=";")
ge.lw.classified.adbmifma.t0.autotest<-read.table(sprintf("%s/ge-nona_lw_t0_Adipocyte-diameter_vs_BMI-Fat-mass-Age_auto.csv",output.dir),sep=";")

#genes significantly different by groups
ge.lw.classified.adfm.t0.autotest.star<-ge.lw.classified.adfm.t0.autotest[ge.lw.classified.adfm.t0.autotest[,4]!="",]
ge.lw.classified.adwh2.t0.autotest.star<-ge.lw.classified.adwh2.t0.autotest[ge.lw.classified.adwh2.t0.autotest[,4]!="",]
ge.lw.classified.adfmh2.t0.autotest.star<-ge.lw.classified.adfmh2.t0.autotest[ge.lw.classified.adfmh2.t0.autotest[,4]!="",]
ge.lw.classified.adh2.t0.autotest.star<-ge.lw.classified.adh2.t0.autotest[ge.lw.classified.adh2.t0.autotest[,4]!="",]
ge.lw.classified.adbmifma.t0.autotest.star<-ge.lw.classified.adbmifma.t0.autotest[ge.lw.classified.adbmifma.t0.autotest[,4]!="",]





#analysis of genes expression based on classied subjects with Funnet
updown.genes.ge.lw.classified.adfm.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                     prode.entrez=ge.prode.entrez,
                                                     high=which(ge.lw.classified.adfm.t0[,4]==0),
                                                     low=which(ge.lw.classified.adfm.t0[,4]==1),
                                                     stars=ge.lw.classified.adfm.t0.autotest.star[,1],
                                                     fileprefix="ge.lw.classified.adfm.t0",
                                                     output.dir=output.dir)
updown.genes.ge.lw.classified.adwh2.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                     prode.entrez=ge.prode.entrez,
                                                     high=which(ge.lw.classified.adwh2.t0[,4]==0),
                                                     low=which(ge.lw.classified.adwh2.t0[,4]==1),
                                                     stars=ge.lw.classified.adwh2.t0.autotest.star[,1],
                                                     fileprefix="ge.lw.classified.adwh2.t0",
                                                     output.dir=output.dir)
updown.genes.ge.lw.classified.adfmh2.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                      prode.entrez=ge.prode.entrez,
                                                      high=which(ge.lw.classified.adfmh2.t0[,4]==0),
                                                      low=which(ge.lw.classified.adfmh2.t0[,4]==1),
                                                      stars=ge.lw.classified.adfmh2.t0.autotest.star[,1],
                                                      fileprefix="ge.lw.classified.adfmh2.t0",
                                                      output.dir=output.dir)
updown.genes.ge.lw.classified.adh2.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                    prode.entrez=ge.prode.entrez,
                                                    high=which(ge.lw.classified.adh2.t0[,4]==0),
                                                    low=which(ge.lw.classified.adh2.t0[,4]==1),
                                                    stars=ge.lw.classified.adh2.t0.autotest.star[,1],
                                                    fileprefix="ge.lw.classified.adh2.t0",
                                                    output.dir=output.dir)
updown.genes.ge.lw.classified.adbmifma.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                        prode.entrez=ge.prode.entrez,
                                                        high=which(ge.lw.classified.adbmifma.t0[,4]==0),
                                                        low=which(ge.lw.classified.adbmifma.t0[,4]==1),
                                                        stars=ge.lw.classified.adbmifma.t0.autotest.star[,1],
                                                        fileprefix="ge.lw.classified.adbmifma.t0",
                                                        output.dir=output.dir)


ge.lm.classified.wfmh2.t0<-read.csv(sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lm.classified.wfmh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lm_t0_Waist_vs_Fat-mass_Height2_auto-test.csv",output.dir),sep=";")
ge.lm.classified.wfmh2.t0.autotest.star<-ge.lm.classified.wfmh2.t0.autotest[ge.lm.classified.wfmh2.t0.autotest[,4]!="",]

updown.genes.ge.lm.classified.wfmh2.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                            prode.entrez=ge.prode.entrez,
                                                            high=which(ge.lm.classified.wfmh2.t0[,4]==0),
                                                            low=which(ge.lm.classified.wfmh2.t0[,4]==1),
                                                            stars=ge.lm.classified.wfmh2.t0.autotest.star[,1],
                                                            fileprefix="ge.lm.classified.wfmh2.t0",
                                                            output.dir=output.dir)



ge.lw.classified.wfmh2.t0<-read.csv(sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_points.csv",output.dir))
ge.lw.classified.wfmh2.t0.autotest<-read.table(sprintf("%s/ge-nona_lw_t0_Waist_vs_Fat-mass_Height2_auto-test.csv",output.dir),sep=";")
ge.lw.classified.wfmh2.t0.autotest.star<-ge.lw.classified.wfmh2.t0.autotest[ge.lw.classified.wfmh2.t0.autotest[,4]!="",]

updown.genes.ge.lw.classified.wfmh2.t0<-genes.updown.funnet(genes=ge.data.t0.nona,
                                                            prode.entrez=ge.prode.entrez,
                                                            high=which(ge.lw.classified.wfmh2.t0[,4]==0),
                                                            low=which(ge.lw.classified.wfmh2.t0[,4]==1),
                                                            stars=ge.lw.classified.wfmh2.t0.autotest.star[,1],
                                                            fileprefix="ge.lw.classified.wfmh2.t0",
                                                            output.dir=output.dir)

#Checking UP-DOWN
m<-mean(as.numeric(as.matrix(ge.data.t0.nona)))
updown<-cbind(colnames(ge.data.t0.nona),colMeans(ge.data.t0.nona[which(ge.lm.classified.adbmifma.t0[,4]==1),]),colMeans(ge.data.t0.nona[which(ge.lm.classified.adbmifma.t0[,4]==0),]),ifelse(colMeans(ge.data.t0.nona[which(ge.lm.classified.adbmifma.t0[,4]==1),])>m,1,0),ifelse(colMeans(ge.data.t0.nona[which(ge.lm.classified.adbmifma.t0[,4]==0),])>m,1,0))

