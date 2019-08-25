##############################################################################
## Copyright (C) Xbiome Biotech Co., Ltd. - All Rights Reserved
## Unauthorized copying of this file, via any medium is strictly prohibited
## Proprietary and confidential
##
## File description: text
##
## Written by Wei Jiang <w4356y@163.com>, March 2019
## Updated by Wei Jiang <w4356y@163.com>, 2019-03-21
##############################################################################



############################################################# 16S Analysis Processing #######################################################################

# ROXYGEN_START

## Process as single end

#' Process Single Reads
#'
#' This function is used to process output results from dada2,
#' which will trim the sequence with length trim_len from the end of
#' reverse reads as features.
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param seqRes_double_exp output from dada2 results, list object containing
#' @param trim_len the length of features that will be trimmed from, the end of reverse reads. Technically, this value should be no more
#' than 217(or the length of single reads)
#' @return An otu table with rows being samples and columns being features.
#' @export
#' @examples

process_se <- function(seqRes_double_exp, trim_len){
  ## Note: This function is used to process output results from dada2,
  ## which will trim the sequence with length trim_len from the end of
  ## reverse reads as features.

  ## Inputs: seqRes_double_exp, output from dada2 results, list object containing
  ## clean_reads (any more info: https://gitlab.com/xbiome/data-analysis/bj_cancer/pipeline-tuning/blob/develop/main.R).

  ## Inputs: trim_len, the length of features that will be trimmed from
  ## the end of reverse reads. Technically, this value should be no more
  ## than 217(or the length of single reads)

  ## Output: An otu table with rows being samples and columns being features.

  library(stringr)
  otu <- seqRes_double_exp$clean_reads
  colnames(otu) <- str_sub(colnames(otu), -trim_len, -1)
  otu_derep <- t(rowsum(t(otu), group = colnames(otu)))
  row.names(otu_derep) <- str_replace(row.names(otu_derep), '(.*)_1.fastq.gz', '\\1')
  return (otu_derep)
}


## Process Results from dada2

#' Process as Merged Paired End
#'
#' This function is used to to process results from dada2 results,
#' which will trim the sequence with length trim_len from the end of
#' reverse reads as features. Mostly the same with process_se, but this
#' function can trim more than single reads.
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param seqRes_double_exp output from dada2 results
#' @param trim_len reads whose length is less than this value will be dumped
#' @return A list contain Abd and Rel_Abd
#' @export
#' @examples

pre_process <- function(seqRes_double_exp, trim_len){

  ## Note: This function is used to to process results from dada2 results,
  ## which will trim the sequence with length trim_len from the end of
  ## reverse reads as features. Mostly the same with process_se, but this function
  ## can trim more than single reads.

  ## Inputs: seqRes_double_exp, output from dada2 results;

  ## Inputs: trim_len, reads whose length is less than this value will be dumped.

  ## Output: A list contain Abd and Rel_Abd.
  library(stringr)
  name <- c()
  id <- c()
  seq <- c()
  for (i in 1:length(seqRes_double_exp$sequence_track)){
    name <- c(name, names(seqRes_double_exp$sequence_track[[i]]$nomerge_sequence))
    id <- c(id, Reduce(rbind, seqRes_double_exp$sequence_track[[i]]$nomerge_sequence))
    seq <- c(seq, seqRes_double_exp$sequence_track[[i]]$sequence)
  }
  df_mt <- unique(data.frame(id, name, seq))

  with_id <- c()
  for(k in 1:length(colnames(seqRes_double_exp$clean_reads))){
    if (nchar(as.character(filter(df_mt,id == colnames(seqRes_double_exp$clean_reads)[k])$seq))>217){
      with_id[length(with_id)+1] = k
      seq_name = as.character(filter(df_mt,id == colnames(seqRes_double_exp$clean_reads)[k])$seq)
      colnames(seqRes_double_exp$clean_reads)[k]=str_sub(seq_name,length(seq_name)-218,-1)
  }

}

  otu_rel <- sweep(seqRes_double_exp$clean_reads,1,rowSums(seqRes_double_exp$clean_reads),'/')[,with_id]
  otu_count <- seqRes_double_exp$clean_reads[,with_id]
  res <- list("Abd"=otu_count,"Rel_Abd"=otu_rel)
  return (res)

}


### Merge Feature Sequence to Genus Level

#' Group Sequence to Genus Level
#'
#' This function is used to group feature sequence into Genus level,
#' Sequence with same Genus taxonomy will be add together.
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param taxa, a data frame/matrix with row names being sequece consistent
#' with the column name of otu table.
#' @param otu, otu table with rows being samples and columns being sequence
#' @return An otu table with rows being samples and columns being Genus names
#' @export
#' @examples


read_to_genus <- function(taxa, otu){
  ## Note: This function is used to group feature sequence into Genus level,
  ## Sequence with same Genus taxonomy will be add together.

  ## Input: taxa, a data frame/matrix with row names being sequece consistent
  ## with the column name of otu table.

  ## Input: otu, otu table with rows being samples and columns being feature
  ## sequence.

  ## Output: An otu table with rows being samples and columns being Genus names.

  map_data <- taxa[colnames(otu), ] %>% as.data.frame() %>%
    mutate(id=paste(.$Family, .$Genus, sep = " "), names = colnames(otu))
  otu_genus <- NULL
  for (i in 1:length(unique(map_data$id))){
    if(is.na(str_match(unique(map_data$id)[i],"NA"))){
      #otu_merge %>% mutate(unique(map_data$id))[i]=sum(.[,filter(map_data,id==unique(map_data$id))[i])]))
      x=unique(filter(map_data,id == unique(map_data$id)[i])$names)
      if(length(x) == 1){otu_genus[[unique(map_data$id)[i]]] = otu[,x]}
      else{
        otu_genus[[unique(map_data$id)[i]]] = unique(filter(map_data, id == unique(map_data$id)[i])$names) %>% otu[,.] %>% apply(.,1,sum) %>% as.numeric()
      }
    }
  }
  otu_genus <- as.data.frame(otu_genus)

}

### Rarefaction of Otu Table

#' Rarefaction of Otu Table
#'
#' Rarefy the original Otu table with minimal sample reads
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param otu Absolute count
#' @return rarefied otu
#' @export
#' @examples


rarefy_otu <- function(otu){
  ## Note: Rarefy the original Otu table with minimal sample reads.

  ## Input: otu, Absolute count

  ## Output: otu, rarefied otu.

  library("GUniFrac")
  otu_rared <- Rarefy(otu,depth = min(rowSums(otu)))
  return (otu_rared$otu.tab.rff)
}

## Write an otu table to Biom format.

#' Write Biom File
#'
#' write otu table into a biom file
#'
#' @author Wei Jiang <jiangwei@xbiome.cn>
#' @param otu otu table for write
#' @param rep_fp file path to write representative sequence to
#' @return
#' @export
#' @examples

write_biom_from_otu_table <- function(otu, rep_fp, biom_fp){
  ## Note: write otu table into a biom file.

  ## Input: otu, otu table for write.

  ## Input: rep_fp, file path to write representative sequence to.

  ## Input: biom_fp, file path to write biom file to.

  library(biomformat)
  seq <- colnames(otu)
  for(i in 1:length(seq)){
    write.table(paste(">seq", i, sep = ""), file = rep_fp, append = T, row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(as.character(seq[i]), file = rep_fp, append = T, quote = FALSE, col.names = FALSE, row.names = FALSE)
  }
  colnames(otu) = paste("seq",c(1:length(colnames(otu))), sep = "")
  biom_data = make_biom(t(otu))
  write_biom(biom_data,biom_file = biom_fp)
}

### Remove Batch Effect

#' Remove Batch Effect
#'
#' tu table with column being feature and rows being sample;
#' replace zeros with half the minimal frequency;
#' method with removeBatchEffect from limma;
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param otu_count otu table with absolute count
#' @param batch factor or vector indicating batches
#' @return otu table de-batch-effected
#' @export
#' @examples

remove_batch=function(otu_count,batch){
  ## Note: otu table with column being feature and rows being sample;
  ## replace zeros with half the minimal frequency;
  ## method with removeBatchEffect from limma;

  ## Input: otu_count, otu table with absolute count

  ## Input: batch, factor or vector indicating batches.

  ## Output: otu table de-batch-effected.

  library(limma)
  min_v=otu_count[otu_count>0] %>% min()
  otu_count[otu_count==0]=min_v/2
  batched_otu=removeBatchEffect(otu_count %>% log() %>% t(),batch) %>% exp() %>% t() %>% as.data.frame()
  return(batched_otu)
}


###  Distance Metric Function

#' Cosine Distance of Samples
#'
#' Give the OTU table, return the cosine distance matrix;
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param x OTU table, with rows being samples and columns being features.
#' @return distance matrix
#' @export
#' @examples

cosine <- function(x) {
  ## Note: Give the OTU table, return the cosine distance matrix;

  ## Input: x, OTU table, with rows being samples and columns being features.

  ## Output: distance matrix;

  x <- t(x)
  x <- as.matrix(x)
  y <- t(x) %*% x
  res <- 1 - y / (sqrt(diag(y)) %*% t(sqrt(diag(y))))
  res <- as.dist(res)
  attr(res, "method") <- "cosine"
  return(res)
}

#' Bray-Curtis Distance of Samples
#'
#' Give the OTU table, return the bray-curtis distance matrix
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param x OTU table, with rows being samples and columns being features
#' @return distance matrix
#' @export
#' @examples

bray <- function(x){
  ## Note: Give the OTU table, return the bray-curtis distance matrix;

  ## Input: x, OTU table, with rows being samples and columns being features.

  ## Output: distance matrix;

  x <- t(x)
  res = vegan::vegdist(t(x), method = "bray")
  attr(res, "method") <- "bray"
  return(res)
}

#' Jaccard Distance of Samples
#'
#' Give the OTU table, return the jaccard distance matrix
#'
#' @author Wei Jiang <jiangwei@xbiome.cn>
#' @param x OTU table, with rows being samples and columns being features.
#' @return distance matrix
#' @export
#' @examples

jaccard <- function(x){
  ## Note: Give the OTU table, return the jaccard distance matrix;

  ## Input: x, OTU table, with rows being samples and columns being features.

  ## Output: distance matrix;

  x <- t(x)
  res = vegan::vegdist(t(x), method = "jaccard")
  attr(res, "method") <- "jaccard"
  return(res)
}


############################################################# Ploting #######################################################################

pacman::p_load("pvclust","stats","ggplot2","vegan","ggbiplot","ade4","phyloseq")

## Principal Component Analysis

#' Principal Component Analysis
#'
#' Give otu table and meta table to generate PCA plot
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param data OTU Table with rows being samples and columns being feature;
#' @param meta meta table
#' @param column a column in meta table used as reference
#' @return
#' @export
#' @examples

plot_pca <- function(data,meta,column){
  ## Give otu table and meta table to generate PCA plot

  ## Input: data, OTU Table,

  ## Input: meta, meta table

  ## Input: column, a column in meta table used as reference.

  library(ggbiplot)
  pca <- prcomp(data, scale. = TRUE)
  gp = ggbiplot(pca,
                var.axes = FALSE,
                obs.scale = 1,
                groups = meta[row.names(pca$x),column] %>%
                  unlist() %>%
                  as.character(), ellipse = TRUE) +
    theme_bw()
  res <- list("plot" = gp,
              "data" = pca)
  return (res)
}

## Permanova Test

#' Permanova Test
#'
#' Permanova test used to test the dispersion difference
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param d OTU TABLE with rows being sample and columns being features
#' @param meta metadata
#' @param feature a column name in meta
#' @param metric a distance metric used to compute the distance between samples
#' @return a statistac res, p value is get by res$aov.tab$`Pr(>F)`[1]
#' @export
#' @examples

permanova_test <- function(d, meta, feature, metric){
  ## Note: Permanova test used to test the dispersion difference

  ## Input: d, OTU TABLE with rows being sample and columns being features.

  ## Input: meta, metadata;

  ## Input: feature, a column name in meta;

  ## Input: metric, a distance metric used to compute the distance between samples.

  ## Output: a statistac res, p value is get by res$aov.tab$`Pr(>F)`[1]

  id = which(!is.na(meta[,feature]) & meta[,feature] != "Others")
  res = adonis(d[id, ] ~ as.character(unlist(meta[id, feature])),
               data = data.frame(meta),
               permutations = 999,
               method = metric)
  return(res)
}

## Permanova Test with TWO variable

#' Permanova Test with TWO variabl
#'
#' Permanova test used to test the dispersion difference FOR two variables
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param d OTU TABLE with rows being sample and columns being features
#' @param meta metadata
#' @param feature1 a column name in meta
#' @param feature2 a column name in meta
#' @return a statistac res, p value is get by res$aov.tab$`Pr(>F)`[1]
#' @export
#' @examples

advanced_permanova_test = function(d, meta, feature1, feature2, metric){
  ## Note: Permanova test used to test the dispersion difference FOR two variables.

  ## Input: d, OTU TABLE with rows being sample and columns being features.

  ## Input: meta, metadata;

  ## Input: feature1, a column name in meta;

  ## Input: feature2, a column name in meta

  ## Input: metric, a distance metric used to compute the distance between samples.

  ## Output: a statistac res, p value is get by res$aov.tab$`Pr(>F)`[1]

  f = paste("d",
            " ~ ",
            feature1,
            "*",
            feature2,
            sep = "")
  res = do.call("adonis2",
                list(as.formula(f),
                     data = data.frame(meta),
                     permutations = 999,
                     method = "bray"))
  #res=adonis2(f,data=meta,permutations = 999,method=metric)
  return(res)
}


## Principal Cordinate Analysis

#' Principal Cordinate Analysis
#'
#' This function is used to plot PCoA graph of otu table
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param data OTU TABLE with rows being sample and columns being features
#' @param meta metadata
#' @param column a column in metadata used to add color
#' @param metric the distance metric used to calculate the distance between samples,
#' including "bray", "jaccard", "unifrac", "wunifrac"
#' @param trefile tree file, .nwk or .tre file showing the phylogenetic tree of
#' features, this has to be specified if the metric is "unifrac" or "wunifrac".
#' @return A list contain pcoa results(pcoa), data frame merged by pcoa and meta(data_ord), pcoa plot
#' @export
#' @examples

plot_pcoa <- function(data, meta, column, metric, trefile=""){
  ## Note: This function is used to plot PCoA graph of otu table.

  ## Input: data, OTU TABLE with rows being sample and columns being features.

  ## Input: meta, metadata;

  ## Input: column, a column in metadata used to add color.

  ## Input: metric, the distance metric used to calculate the distance between samples.

  ## Input: trefile, tree file, .nwk or .tre file showing the phylogenetic tree of
  ## features, this has to be specified if the metric is "unifrac" or "wunifrac".

  ## Output: A list contain pcoa results(pcoa), data frame merged by pcoa and meta(data_ord), pcoa plot.

  ## Dependency: permanova_test function has to be pre-defined to calculate the dispersion.

  perm_res = permanova_test(data, meta, column, metric)
  p_value = perm_res$aov.tab$`Pr(>F)`[1]
  if (!metric %in% c("bray", "jaccard", "unifrac", "wunifrac")){
    stop("distance metric should be one of (Bray-Curtis, Jaccard, Unifrac or Weighted-Unifrac)")
  }else{
    if (metric %in% c("unifrac", "wunifrac")){
      if (is.null(tree)){
        stop("Phylogenic tree should be specified for Unifrac or Weighted-Unifrac distance")
      }else{
        library(phyloseq)
        sample_data = sample_data(meta)
        otu_table = otu_table(data, taxa_are_rows = TRUE)
        tree=phy_tree(trefile)
        phylo1 = phyloseq(sample_data, otu_table, tree)
        pcoa = cmdscale(distance(phylo1, method = metric),
                        k = 2, eig = T)
      }

    }else{
      pcoa = cmdscale(vegdist(data, method = metric),
                      k = 2, eig = T)
    }
    mds1 = as.data.frame(pcoa$points)
    colnames(mds1) = c("pc1", "pc2")
    data_ord=merge(mds1, meta, by = "row.names")
    row.names(data_ord) = data_ord$Row.names
    p= ggplot(data=data_ord, aes(x = pc1, y = pc2)) +
      geom_point(aes(color = as.character(unlist(meta[row.names(data_ord), column])), size = 3)) +
      labs(color = column) +
      guides(size = FALSE) +
      xlab(paste("PC1 ", round(100*as.numeric(pcoa$eig[1] / sum(pcoa$eig)), 2), "%", sep = "")) +
      ylab(paste("PC2 ", round(100 * as.numeric(pcoa$eig[2] / sum(pcoa$eig)), 2), "%", sep = "")) +
      ggtitle(metric) +
      theme(plot.title = element_text(hjust = 1,
                                      size = 10)) +
      guides(color = guide_legend(override.aes = list(size = 4))) +
      annotate("text",
               x = max(data_ord$pc1),
               y = max(data_ord$pc2),
               label = paste("p = ", round(p_value, 4), sep = ""),
               hjust = 1) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.text = element_text(size = 12)) +
      ggsci::scale_color_jco() +
      ggsci::scale_fill_jco()
  }
  res= list("pcoa"=pcoa,
            "data"=data_ord,
            "plot"=p)
  return (res)
}


## Principal Cordinate Analysis with Shape

#' Principal Cordinate Analysis with 2 Confounding Factors
#'
#' This function is used to plot PCoA with 2 variables, display in graph as color and shape
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param data OTU TABLE with rows being sample and columns being features
#' @param meta metadata
#' @param column column in meta table used to plot color
#' @param col_control column in meta table used to plot shape
#' @param metric distance measure to compute the distance between samples
#' @param trefile tree file, .nwk or .tre file showing the phylogenetic tree of
#' features, this has to be specified if the metric is "unifrac" or "wunifrac"
#' @return
#' @export
#' @examples

plot_pcoa_with_shape = function(data, meta, column, col_control, metric, trefile = ""){
  ## Note: This function is used to plot PCoA with 2 variables, display in graph as color
  ## shape;

  ## Input: data, OTU TABLE with rows being sample and columns being features.

  ## Input: meta, metadata;

  ## Input: column, column in meta table used to plot color

  ## Input: col_control, column in meta table used to plot shape

  ## Input: metric, distance measure to compute the distance between samples.

  ## Input: trefile, tree file, .nwk or .tre file showing the phylogenetic tree of
  ## features, this has to be specified if the metric is "unifrac" or "wunifrac".

  ## Dependency: advanced_permanova_test function has to be pre-defined to calculate the dispersion.

  meta <- data.frame(meta)
  if(meta[,col_control] %>% is.na() %>% sum() > 0){
    meta[,col_control] = addNA(meta[,col_control])
  }
  if(is.numeric(meta[,col_control])){
    meta[, col_control] = factor(meta[, col_control])
  }
  if(F){
    perm_res = advanced_permanova_test(data, meta, column, col_control, metric)
    p_value=perm_res$`Pr(>F)`[3]
  }
  if(T){
    perm_res=permanova_test(data, meta, col_control, metric)
    p_value=perm_res$aov.tab$`Pr(>F)`[1]
  }
  if (!metric %in% c("bray", "jaccard", "unifrac", "wunifrac")){
    stop("distance metric should be one of (Bray-Curtis, Jaccard, Unifrac or Weighted-Unifrac)")
  }else{
    if (metric %in% c("unifrac", "wunifrac")){
      if (is.null(tree)){
        stop("Phylogenic tree should be specified for Unifrac or Weighted-Unifrac distance")
      }else{
        library(phyloseq)
        sample_data = sample_data(meta)
        otu_table = otu_table(data, taxa_are_rows = TRUE)
        tree = phy_tree(trefile)
        phylo1 = phyloseq(sample_data, otu_table, tree)
        pcoa = cmdscale(distance(phylo1, method = metric),
                        k = 2, eig = T)
      }

    }else{
      pcoa=cmdscale(vegdist(data, method = metric),
                    k = 2, eig = T)
    }
    mds1 = as.data.frame(pcoa$points)
    colnames(mds1) = c("pc1", "pc2")
    data_ord=merge(mds1, meta, by = "row.names")
    row.names(data_ord) = data_ord$Row.names
    shape_values=c(16, 17, 0, 3, 4, 7, 8, 9)
    ggplot(data=data_ord, aes(x = pc1, y = pc2)) +
      geom_point(aes(color = meta[row.names(data_ord), column],
                     shape = meta[row.names(data_ord), col_control],
                     size=3),
                 na.rm = FALSE) +
      labs(color = column,
           shape = col_control) +
      guides(size = FALSE) +
      xlab(paste("PC1 ", round(100 * as.numeric(pcoa$eig[1] / sum(pcoa$eig)), 2), "%", sep = "")) +
      ylab(paste("PC2 ", round(100 * as.numeric(pcoa$eig[2] / sum(pcoa$eig)), 2), "%", sep = "")) +
      scale_shape_manual(values = c(16, 17, 0, 3, 4, 7, 8, 9),
                         na.value = shape_values[levels(meta[, col_control]) %>% length()]) +
      ggtitle(metric) +
      theme(plot.title = element_text(hjust = 1,
                                      size = 10)) +
      guides(color = guide_legend(override.aes = list(size = 5)),
             shape = guide_legend(override.aes = list(size = 4))) +
      annotate("text",
               x = max(data_ord$pc1),
               y = max(data_ord$pc2),
               label = paste("p = ", round(p_value, 4), sep = ""),
               hjust = 1) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.text = element_text(size = 12)) +
      ggsci::scale_color_jco() +
      ggsci::scale_fill_jco()
  }
}


## Plot PCoA and combine the PCoA plot with sequence Distribution

#' PCoA plot with sequence Distribution
#'
#' This function is used to combine PCoA graph and the distribution plot of
#' highly related sequence in PC1 axis
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param pcoa a pcoa object generate from plot_pcoa function
#' @param otu_df, the otu table used to calculate the correlation between PC1 component
#' and sequence feature across all samples. Relative abundant otu table is advised.
#' Must be a data frame.
#' @param Response environment variable used as color reference in PCoA graph. No quote!!!!
#' @param is_sequence, whether the colnames of otu table is sequence or not, if TRUE, taxonomy
#' is conducted using assignTaxonomy function from dada2 package
#' @return Taxonoy Results is displayed
#' @export
#' @examples

plot_pcoa_plus_seq = function(pcoa,otu_df,Response,is_sequence=TRUE){
  ## Note: This function is used to combine PCoA graph and the distribution plot of
  ## highly related sequence in PC1 axis;

  ## Input: pcoa, a pcoa object generate from plot_pcoa function;

  ## Input: otu_df, the otu table used to calculate the correlation between PC1 component
  ## and sequence feature across all samples. Relative abundant otu table is advised.
  ## Must be a data frame;

  ## Input: Response, environment variable used as color reference in PCoA graph. No quote!!!!

  ## Input: is_sequence, whether the colnames of otu table is sequence or not, if TRUE, taxonomy
  ## is conducted using assignTaxonomy function from dada2 package.

  ## Output: Taxonoy Results is displayed.

  library(dada2)
  library(dplyr)
  library(ggplot2)
  library(grid)
  var_col = enquo(Response)
  cor_list = NULL
  for (i in 1:ncol(otu_df)){
    cor_list[i] = cor(otu_df[row.names(pcoa$data), i], pcoa$data$pc1, method = "spearman")
  }
  df_cor = data.frame(cor = cor_list, row.names = colnames(otu_df))
  df_cor$abs_cor = abs(df_cor$cor)
  if(is_sequence){
    taxa1 = df_cor[order(df_cor$abs_cor,decreasing = T),][1:10,] %>%
      round() %>%
      as.matrix() %>%
      t() %>%
      assignTaxonomy(.,"/home/weijiang/xbiome/tumor/pipeline-tuning/data/silva_nr_v132_train_set.fa.gz")
  }else{
    taxa1 = df_cor[order(df_cor$abs_cor,decreasing = T),][1:10,] %>%
      as.matrix()
  }
  df = pcoa$data
  df$Type = "PCoA"
  seq_name = c()
  for(i in 1:4){
    s_name = paste("seq", i, sep = "")
    seq_name = c(seq_name, rep(s_name, nrow(df)))
  }
  df1 = data.frame(pc1 = rep(df$pc1, 4),
                 pc2 = rep(df$pc2, 4),
                 abd = stack(otu_df[row.names(pcoa$data), row.names(taxa1)[1:4]]) %>%
                   .$values,
                 #   abd=c(otu_df[row.names(pcoa$data),row.names(taxa1)[1]],otu_df[row.names(pcoa$data),row.names(taxa1)[2]],otu_df[row.names(pcoa$data),row.names(taxa1)[3]],otu_df[row.names(pcoa$data),row.names(taxa1)[4]]),
                 seq = factor(seq_name,levels=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8","seq9","seq10"))
  )
  #   seq=c(rep("seq1",nrow(df)),rep("seq2",nrow(df)),rep("seq3",nrow(df)),rep("seq4",nrow(df)),rep("seq5",nrow(df))))
  gp1 = ggplot(df, aes(x = pc1, y = pc2)) +
    geom_point(aes(color = !!var_col, size = 3)) +
    theme_bw() +
    ylab(paste("PC2   ", round(100 * as.numeric(pcoa$pcoa$eig[2] / sum(pcoa$pcoa$eig)), 2), "%", sep = "")) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_text(size=15),
      axis.text.y = element_text(size=12),
      legend.title = element_text(size=13),
      legend.text = element_text(size=10))+
    facet_grid(Type ~ .) +
    guides(size=FALSE, color = guide_legend(override.aes = list(size = 4))) +
    ggsci::scale_color_jco()

  gp2 = ggplot(data=df1,aes(x = pc1, y = abd, color = seq)) +
    geom_line() +
    facet_grid(rows = vars(seq)) +
    theme_bw() +
    xlab(paste("PC1   ", round(100 * as.numeric(pcoa$pcoa$eig[1] / sum(pcoa$pcoa$eig)), 2), "%", sep = "")) +
    ylab("Rel_abd") +
    theme(
      axis.title = element_text(size = 15),
      axis.text=element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 10)) +
    ggsci::scale_color_d3()
  g = rbind(ggplotGrob(gp1), ggplotGrob(gp2),size = "last")
  id <- g$layout$t[grep(pattern = "panel", g$layout$name)]
  g$heights[7] = unit(7, units = "null")
  grid.newpage()
  grid.draw(g)
  return (taxa1)
}



### Non-metric Multi-Dimensional Analysis

#' Non-metric Multi-Dimensional Analysis(NMDS)
#'
#' plot Non-metric Multi-Dimensional Analysis results
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param data OTU table with rows being sample names, and columns being features
#' @param meta meta table, with rows being sample names
#' @param column a scalar, quoted, denoted which environment factor should be plot with color
#' @return
#' @export
#' @examples
#' fun(1, 1)
#' fun(10, 1)

plot_nmds = function(data,
                     meta,
                     column){
  ## Note: plot Non-metric Multi-Dimensional Analysis results;

  ## Input: data, OTU table with rows being sample names, and columns being features.

  ## Input: meta, meta table, with rows being sample names;

  ## column: a scalar, quoted, denoted which environment factor should be plot with color.

  ndms1 = metaMDS(data)
  data_mds = merge(as.data.frame(ndms1$points), meta, by = "row.names")
  row.names(data_mds) = data_mds$Row.names
  ggplot(data = data_mds, aes(x = MDS1, y = MDS2)) +
    geom_point(aes(color = as.character(unlist(meta[row.names(data_mds), column])),
                                                            size=3)) +
    labs(color = column) +
    guides(size = FALSE) +
    theme_bw() +
    theme(axis.title = element_text(size=15),
          axis.text = element_text(size=12)) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}

############################################################# Modelling #######################################################################

### plot ROC curve

#' ROC curve
#'
#' plot ROC curve in train set
#'
#' @author Wei Jiang <jiangwei@xbiome.cn>
#' @param model a model after trained
#' @param rain_selected training dataset
#' @return
#' @export
#' @examples

plot_roc = function(model, train_selected){
  ## Note: plot ROC curve in train set

  ## Input: model, a model after trained;

  ## Input: train_selected, training dataset;

  library(pROC)
  preds = predict(model, type = "prob")$R
  roc1 = roc(train_selected$label,
             preds,
           plot = TRUE,
           print.thres = TRUE,
           print.auc = TRUE,
           percent = T,
           asp = NA,
           grid = TRUE)
}

#' Test ROC
#'
#' plot ROC curve in test set
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param model a model after trained
#' @param test_selected test dataset
#' @return the AUC value in test set
#' @export
#' @examples

plot_test_roc=function(model, test_data){
  ## Note: plot ROC curve in test set

  ## Input: model, a model after trained;

  ## Input: test_selected, test dataset;

  ## Output: Return the AUC value in test set.
  library(pROC)
  preds = predict(model, type = "prob", newdata = test_data)$R
  r1 = roc(test_data$label,
         preds,
         plot=TRUE,
         print.thres = TRUE,
         print.auc = TRUE,
         percent = T,
         asp = NA,
         grid = TRUE)
  return (r1$auc[1])
}


### Random Foret Model

#' Random Foret Model
#'
#' apply Random Forest Model for trainning in train set and test in test set
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param train train SET
#' @param test test SET
#' @param n seed
#' @return a list containing model, confusion matrix, true predicting case and test AUC
#' @export
#' @examples

rf_model=function(train, test, n){
  ## Note: apply Random Forest Model for trainning in train set and test in test set;

  ## Input: train, train SET;

  ## Input: test, test SET;

  ## Input: n, seed;

  ## Output: a list containing model, confusion matrix, true predicting case and test AUC;

  set.seed(n)
  control <- trainControl(method = "LGOCV",
                          #number = 8,
                          #repeats = 50,
                          search = "grid",
                          classProbs = TRUE,
                          summaryFunction = twoClassSummary,
                          savePredictions = TRUE,
                          allowParallel = FALSE)
  rf_grid <- expand.grid(.mtry = seq(2, 30, 2))
  res_list=c()
  modellist <- list()
  tree_list=c('100','500','800','1000','1200','1500','2000','2500','3000')
  for(ntree in c(100, 500, 800, 1000, 1200, 1500, 2000, 2500, 3000)) {
    model <- train(label ~. ,
                   data = train,
                   #data=d1,
                   method = "rf",
                   metric = "ROC",
                   tuneGrid = rf_grid,
                   trControl = control,
                   ntree = ntree
                   #num.trees = ntree
    )
    pred=predict(model, newdata = test)
    res_list[length(res_list) + 1] = sum(pred == test$label)
    key <- toString(ntree)
    modellist[[key]] <- model

  }
  a = tree_list[which(res_list == max(res_list))[1]]
  plot_roc(modellist[[a]], train)
  test_auc = plot_test_roc(modellist[[a]], test)
  test_reaction = caret::confusionMatrix(predict(modellist[[a]], test), as.factor(test$label))
  res = list("model" = modellist, "confM" = test_reaction,"true_n" = max(res_list), "test_auc" = test_auc)
  return(res)
}


###Svm Model Trainning

#' SVM Model
#'
#' Apply Support Vector Machine Linear for trainning in train set and test in test set
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param train train SET
#' @param test test SET
#' @param  n seed
#' @param method method for tuning hyperprameters
#' @param fold the number of folds used for cross validation
#' @return a list containing model, confusion matrix, and test AUC
#' @export
#' @examples

svm_model=function(train, test, n=100, method="cv", fold=3){
  ## Note: apply Support Vector Machine Linear for trainning in train set and test in test set;

  ## Input: train, train SET;

  ## Input: test, test SET;

  ## Input: n, seed;

  ## Input: method, method for tuning hyperprameters

  ## Input: fold, the number of folds used for cross validation

  ## Output: a list containing model, confusion matrix, and test AUC;

  set.seed(n)
  control <- trainControl(method = method,
                          number = fold,
                          #repeats = 50,
                          search = "grid",
                          classProbs = TRUE,
                          summaryFunction = twoClassSummary,
                          savePredictions = TRUE)
  grid<- expand.grid(C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,5))
  svm_model <- train(label ~. ,
                     data = train,
                     method = "svmLinear",
                     metric = "ROC",
                     tuneGrid = grid,
                     trControl = control,
                     preProc=c("center","scale")
                     #ntree = ntree
                     #num.trees = ntree
  )
  plot_roc(svm_model,train)
  test_reaction=confusionMatrix(predict(svm_model,test),as.factor(test$label))
  test_auc=plot_test_roc(svm_model,test)
  res=list("model"=svm_model,"conF"=test_reaction,"auc"=test_auc)
  return(res)
}

###Pernalized Logistic Regression

#' Pernalized Logistic Regression
#'
#' Apply Logistic Regression for trainning in train set and test in test set
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param train train SET
#' @param test test SET
#' @param n seed
#' @param method method for tuning hyperprameters
#' @param fold the number of folds used for cross validation
#' @return a list containing model, confusion matrix, and test AUC
#' @export
#' @examples

lr_model=function(train, test, n = 100, method = "cv", fold = 3){
  ## Note: apply Logistic Regression for trainning in train set and test in test set;

  ## Input: train, train SET;

  ## Input: test, test SET;

  ## Input: n, seed;

  ## Input: method, method for tuning hyperprameters

  ## Input: fold, the number of folds used for cross validation

  ## Output: a list containing model, confusion matrix, and test AUC;

  set.seed(n)
  cctrl1 <- trainControl(method = method, number = fold,returnResamp = "all",
                         classProbs = TRUE,
                         summaryFunction = twoClassSummary
  )
  lrGrid <- expand.grid(lambda = c(0, 0.0001, 0.001, 0.01, 0.1, 0.5, 1),cp = c("bic"))
  lr_model=train(train[, -ncol(train)], train$label, method = "plr", tuneGrid = lrGrid, trControl = cctrl1, metric = "ROC")
  plot_roc(lr_model, train)
  test_auc=plot_test_roc(lr_model, test)
  test_reaction=confusionMatrix(predict(lr_model, test), as.factor(test$label))
  res=list("model" = lr_model, "conF" = test_reaction, "auc"=test_auc)
  return (res)

}


###XGB MODEL

#' XGBoost Model
#'
#' Apply XGB Model for trainning in train set and test in test set
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param train train SET
#' @param test test SET
#' @return a list containing model, confusion matrix
#' @export
#' @examples


xgb_model = function(train, test, n){
  ## Note: apply XGB Model for trainning in train set and test in test set;

  ## Input: train, train SET;

  ## Input: test, test SET;

  ## Input: n, seed;

  ## Output: a list containing model, confusion matrix

  set.seed(n)
  cctrl1 <- trainControl(method = "cv", number = 3, returnResamp = "all",
                         classProbs = TRUE,
                         summaryFunction = twoClassSummary
  )

  xgbGrid <- expand.grid(nrounds = c(20, 50, 100, 200, 500),
                         max_depth = c(2, 4, 6, 8, 10, 12, 14),
                         eta = 0.05,
                         colsample_bytree = 0.90,
                         rate_drop=0.1,skip_drop=0.1,
                         min_child_weight = 2,
                         subsample = 0.75,
                         gamma = 0.01)
  xgb_model=train(train[,-ncol(train)],train$label,
                  method = "xgbDART",
                  trControl = cctrl1,
                  tuneGrid = xgbGrid,
  )
  plot_roc(xgb_model,train)
  test_reaction=confusionMatrix(predict(xgb_model,test),as.factor(test$label))
  res=list("model"=xgb_model,"conF"=test_reaction)
  return(res)
}



####Feature selection by rfe

#' Feature Selection by rfe
#'
#' Apply rfe in random forest for feature selection
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param train train SET
#' @param n seed
#' @return Feature Selection Results
#' @export
#' @examples

rfe_select=function(dataset,n){
  ## Note: apply rfe in random forest for feature selection;

  ## Input: train, train SET;

  ## Input: n, seed;

  ## Output: Feature Selection Results

  set.seed(n)
  control <- rfeControl(functions = rfFuncs, method = "LGOCV", repeats = 500)
  results <- rfe(dataset[, -ncol(dataset)], factor(dataset$label), sizes = c(1:10, seq(15, 40, 5)), rfeControl = control)
  res=list("results" = results)
  return(res)
}



###feature selection by importance

#' Feature Selection by Importance
#'
#' Select feature with importance
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param dataset Train Dataset
#' @param n seed
#' @param top_n top n features selected
#' @return Feature Selection Results
#' @export
#' @examples

imp_select=function(dataset, n, top_n = 10){
  ## Note:  Select feature with importance;

  ## Input: dataset, Train Dataset;

  ## Input: n, seed;

  ## Input: top_n, top n features selected;

  ## Output: Feature Selection Results

  set.seed(n)
  control <- trainControl(method = "repeatedcv", number = 5, repeats = 20)
  #set.seed(123)
  model <- train(label ~. ,
                 data = dataset,
                 method = "rf",
                 trControl = control,
                 ntree = 1000
  )

  res = varImp(model)[1]
  species = row.names(res$importance)[order(res$importance,decreasing = TRUE)]
  contributor = res$importance[order(res$importance,decreasing = TRUE),]
  feature=species[1:top_n]
  #train_selected=dataset[row.names(matson_data),c(str_replace_all(feature,'`',''),"label")]
  #test_selected=dataset[row.names(BJ_data),c(str_replace_all(feature,'`',''),"label")]
  #res=list("model"=model,"feature"=feature,"train"=train_selected,"test"=test_selected)
  res=list("feature" = feature)
  return(res)
}


###Plot heatmap

#' Plot Heatmap of OTU table
#'
#' Plot Heatmap of OTU table
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param data an otu table
#' @param feature feature used in otu table
#' @param lable a column in data used to plot color bar
#' @return
#' @export
#' @examples

plot_heatmap=function(data, feature, label){
  ## Plot Heatmap of OTU table

  ## Input: feature, feature used in otu table;

  ## Input: lable, a column in data used to plot color bar;

  ## Input: data, an otu table;

  anno_df = data[,"label"]
  anno_df = as.data.frame(anno_df,row.names(data))
  colnames(anno_df) = c("Response")
  col=list(Response = c("R" = "lightgreen", "NR" = "red"))
  ha <- HeatmapAnnotation(anno_df, col = col)
  log(t(data[, feature])+1)%>%
    Heatmap(top_annotation = ha,
            column_names_gp = gpar(fontsize = 7),
            row_names_gp = gpar(fontsize = 10), heatmap_legend_param = list(title= "log(1000*rel_abd+1)"),
    )

}




###Glmnet Model

#' Glmnet Model
#'
#' Apply Generalized Linear Model for trainning in train set and test in test set
#'
#' @author Wei Jiang <w4356y@163.com>
#' @param train train SET
#' @param test test SET
#' @param method method for tuning hyperprameters
#' @param n seed
#' @param method method for tuning hyperprameters
#' @param fold the number of folds used for cross validation
#' @return A list containing model, confusion matrix, test AUC and variables used
#' @export
#' @examples

glmnet_model=function(train,test,n=100,method="cv",fold=3){
  ## Note: apply Generalized Linear Model for trainning in train set and test in test set;

  ## Input: train, train SET;

  ## Input: test, test SET;

  ## Input: n, seed;

  ## Input: method, method for tuning hyperprameters

  ## Input: fold, the number of folds used for cross validation

  ## Output: a list containing model, confusion matrix, test AUC and variables used;

  set.seed(n)
  cctrl1 <- trainControl(method = method, number = fold, returnResamp = "all",
                         classProbs = TRUE, summaryFunction = twoClassSummary)
  model=train(label ~ .,
              data=train,
              method = "glmnet",
              trControl = cctrl1,
              metric = "ROC",
              preProc = c("center", "scale"),
              tuneGrid = expand.grid(.alpha = seq(.05, 1, length = 15),
                                     .lambda = c((1:5)/10)))
  plot_roc(model,train)
  test_reaction=caret::confusionMatrix(predict(model,test),as.factor(test$label))
  #res=list("model"=model,"conF"=test_reaction)
  test_auc=plot_test_roc(model,test)
  coefficients <- coef(model$finalModel, model$bestTune$lambda)
  variables <- names(coefficients[which(coefficients[,1] != 0),])
  res=list("model"=model,"conF"=test_reaction,"Var"=variables,"auc"=test_auc)
  return(res)

}





plot_bj_sample_description = function(phylo){
  otu = otu_table(phylo1)
  otu_rel= sweep(otu, 1, rowSums(otu),'/')
  sample_meta= sample_data(phylo)
  merged_otu = merge(otu_rel, sample_meta, by="row.names")
  merged_otu$differential_day = difftime(merged_otu$CollectionDate,merged_otu$StartDate,units = "days") %>% as.numeric()
  merged_otu$differential_week = round(merged_otu$differential_day/7)
  id1 = merged_otu %>% group_by(SubjectID) %>% arrange(CollectionTime, .by_group = T) %>% slice(1) %>% select(SubjectID,differential_day,TTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGACGCTCAACGTCTGCACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGCGCAACTGACGCTGAAGCTCGAAAGTGCGGGTATCGAACAGG)  %>% filter(differential_day < 0.5) %>% arrange(TTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGACGCTCAACGTCTGCACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGCGCAACTGACGCTGAAGCTCGAAAGTGCGGGTATCGAACAGG) %>% select(SubjectID) %>% unlist() %>% as.character()
  id2 = merged_otu %>% group_by(SubjectID) %>% arrange(CollectionTime, .by_group = T) %>% slice(1) %>% select(SubjectID,differential_day,TTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGACGCTCAACGTCTGCACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGCGCAACTGACGCTGAAGCTCGAAAGTGCGGGTATCGAACAGG)  %>% filter(differential_day >= 0.5) %>% arrange(TTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGACGCTCAACGTCTGCACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGCGCAACTGACGCTGAAGCTCGAAAGTGCGGGTATCGAACAGG) %>% select(SubjectID) %>% unlist() %>% as.character()
  p_level = c(id2,id1)
  merged_otu$SubjectID=factor(merged_otu$SubjectID,levels=p_level)
  merged_otu$type= ifelse(merged_otu$SubjectID %in% id1,"pre","post")

  metadata_plot2 <- merged_otu %>%
    select(SubjectID, starts_with("L"), CollectionTime) %>%
    gather(key = "Assessment", value = "Time", -SubjectID) %>%
    unique() %>%
    filter(!is.na(Time)) %>%
    mutate(AssessmentTime = Time) %>%
    mutate(AssessmentTime = str_replace(AssessmentTime, "W(\\d+).+$", "\\1")) %>%
    mutate(AssessmentTime = as.numeric(AssessmentTime)) %>%
    mutate(Event = str_replace(Time, "W\\d+ (.+)$", "\\1")) %>%
    mutate(Event = ifelse(Assessment == 'CollectionTime', 'Sample', Event)) %>%
    mutate(Event = factor(Event, levels = c('PR', 'SD', 'PD', 'Sample'))) %>%
    filter(Event != "Sample")
  p = merged_otu %>% group_by(type) %>% arrange(TTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGACGCTCAACGTCTGCACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGCGCAACTGACGCTGAAGCTCGAAAGTGCGGGTATCGAACAGG, .by_group = T) %>%
    ggplot(.,aes(x=differential_week, y=SubjectID)) + geom_point(aes(size=TTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGACGCTCAACGTCTGCACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGAGCGCAACTGACGCTGAAGCTCGAAAGTGCGGGTATCGAACAGG)) +
    guides(size=FALSE) + theme_bw() + geom_hline(yintercept = 60.5,color="red") +ggsci::scale_color_jco()
  p + geom_point(data=metadata_plot2,aes(x= AssessmentTime, y = SubjectID,shape=Event),size=3) +
    scale_shape_manual(values = c(0, 1, 2, 20)) + geom_line(data=metadata_plot2,aes(x= AssessmentTime, y = SubjectID),color="blue")

}


# ROXYGEN_STOP



