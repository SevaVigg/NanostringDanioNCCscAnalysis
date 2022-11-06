# This scripts performes integration of several RNA-Seq datasets. It runs Seurat 3.0 (please not that Nanostring scripts run Seurat 2.x!!!)
#
# Written by Leonid A. Uroshlev, 26.10.2022


library(batchelor)
library(scran)
library(Seurat)
library(patchwork)
library(biomaRt)
library(ggplot2)
library(cowplot)


FeaturePlotN <- function (object, features, genenames = NULL, dims = c(1, 2), cells = NULL, cols = if (blend) {
    c("lightgrey", "#ff0000", "#00ff00")
} else {
    c("lightgrey", "blue")
}, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA, 
    reduction = NULL, split.by = NULL, keep.scale = "feature", 
    shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5, 
    label = FALSE, label.size = 4, repel = FALSE, ncol = NULL, 
    coord.fixed = FALSE, by.col = TRUE, sort.cell = NULL, interactive = FALSE, 
    combine = TRUE, raster = NULL) 
{
    if (!is.null(x = sort.cell)) {
        warning("The sort.cell parameter is being deprecated. Please use the order ", 
            "parameter instead for equivalent functionality.", 
            call. = FALSE, immediate. = TRUE)
        if (isTRUE(x = sort.cell)) {
            order <- sort.cell
        }
    }
    if (interactive) {
        return(IFeaturePlot(object = object, feature = features[1], 
            dims = dims, reduction = reduction, slot = slot))
    }
    if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature", 
        "all"))) {
        stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
    }
    no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", 
            size = 14, margin = margin(r = 7)))
    reduction <- reduction %||% DefaultDimReduc(object = object)
    if (length(x = dims) != 2 || !is.numeric(x = dims)) {
        stop("'dims' must be a two-length integer vector")
    }
    if (blend && length(x = features) != 2) {
        stop("Blending feature plots only works with two features")
    }
    if (blend) {
        default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
        cols <- switch(EXPR = as.character(x = length(x = cols)), 
            `0` = {
                warning("No colors provided, using default colors", 
                  call. = FALSE, immediate. = TRUE)
                default.colors
            }, `1` = {
                warning("Only one color provided, assuming specified is double-negative and augmenting with default colors", 
                  call. = FALSE, immediate. = TRUE)
                c(cols, default.colors[2:3])
            }, `2` = {
                warning("Only two colors provided, assuming specified are for features and agumenting with '", 
                  default.colors[1], "' for double-negatives", 
                  call. = FALSE, immediate. = TRUE)
                c(default.colors[1], cols)
            }, `3` = cols, {
                warning("More than three colors provided, using only first three", 
                  call. = FALSE, immediate. = TRUE)
                cols[1:3]
            })
    }
    if (blend && length(x = cols) != 3) {
        stop("Blending feature plots only works with three colors; first one for negative cells")
    }
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% colnames(x = object)
    data <- FetchData(object = object, vars = c(dims, "ident", 
        features), cells = cells, slot = slot)
    if (ncol(x = data) < 4) {
        stop("None of the requested features were found: ", paste(features, 
            collapse = ", "), " in slot ", slot, call. = FALSE)
    }
    else if (!all(dims %in% colnames(x = data))) {
        stop("The dimensions requested were not found", call. = FALSE)
    }
    features <- colnames(x = data)[4:ncol(x = data)]
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
            feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
            feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
        max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
        ]$maxcolors, no = length(x = cols))
    data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), 
        FUN = function(index) {
            data.feature <- as.vector(x = data[, index])
            min.use <- SetQuantile(cutoff = min.cutoff[index - 
                3], data.feature)
            max.use <- SetQuantile(cutoff = max.cutoff[index - 
                3], data.feature)
            data.feature[data.feature < min.use] <- min.use
            data.feature[data.feature > max.use] <- max.use
            if (brewer.gran == 2) {
                return(data.feature)
            }
            data.cut <- if (all(data.feature == 0)) {
                0
            }
            else {
                as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                  breaks = brewer.gran)))
            }
            return(data.cut)
        })
    colnames(x = data)[4:ncol(x = data)] <- features
    rownames(x = data) <- cells
    data$split <- if (is.null(x = split.by)) {
        RandomName()
    }
    else {
        switch(EXPR = split.by, ident = Idents(object = object)[cells, 
            drop = TRUE], object[[split.by, drop = TRUE]][cells, 
            drop = TRUE])
    }
    if (!is.factor(x = data$split)) {
        data$split <- factor(x = data$split)
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    plots <- vector(mode = "list", length = ifelse(test = blend, 
        yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
    xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, 
        dims[1]])))
    ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, 
        dims[2]])))
    if (blend) {
        ncol <- 4
        color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold, 
            negative.color = cols[1])
        cols <- cols[2:3]
        colors <- list(color.matrix[, 1], color.matrix[1, ], 
            as.vector(x = color.matrix))
    }
    for (i in 1:length(x = levels(x = data$split))) {
        ident <- levels(x = data$split)[i]
        data.plot <- data[as.character(x = data$split) == ident, 
            , drop = FALSE]
        if (blend) {
            features <- features[1:2]
            no.expression <- features[colMeans(x = data.plot[, 
                features]) == 0]
            if (length(x = no.expression) != 0) {
                stop("The following features have no value: ", 
                  paste(no.expression, collapse = ", "), call. = FALSE)
            }
            data.plot <- cbind(data.plot[, c(dims, "ident")], 
                BlendExpression(data = data.plot[, features[1:2]]))
            features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
        }
        for (j in 1:length(x = features)) {
            feature <- features[j]
            if (blend) {
                cols.use <- as.numeric(x = as.character(x = data.plot[, 
                  feature])) + 1
                cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
            }
            else {
                cols.use <- NULL
            }
            data.single <- data.plot[, c(dims, "ident", feature, 
                shape.by)]
            plot <- SingleDimPlot(data = data.single, dims = dims, 
                col.by = feature, order = order, pt.size = pt.size, 
                cols = cols.use, shape.by = shape.by, label = FALSE, 
                raster = raster) + scale_x_continuous(limits = xlims) + 
                scale_y_continuous(limits = ylims) + theme_cowplot() + 
                CenterTitle()
            if (label) {
                plot <- LabelClusters(plot = plot, id = "ident", 
                  repel = repel, size = label.size)
            }
            if (length(x = levels(x = data$split)) > 1) {
                plot <- plot + theme(panel.border = element_rect(fill = NA, 
                  colour = "black"))
                plot <- plot + if (i == 1) {
                  labs(title = feature)
                }
                else {
                  labs(title = NULL)
                }
                if (j == length(x = features) && !blend) {
                  suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident), 
                    limits = ylims) + no.right)
                }
                if (j != 1) {
                  plot <- plot + theme(axis.line.y = element_blank(), 
                    axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                    axis.title.y.left = element_blank())
                }
                if (i != length(x = levels(x = data$split))) {
                  plot <- plot + theme(axis.line.x = element_blank(), 
                    axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                    axis.title.x = element_blank())
                }
            }
            else {
print(j)
                plot <- plot + labs(title = genenames[genenames[,2] == features[[j]], 1])+theme(axis.text=element_text(size=36), axis.title=element_text(size=36), plot.title = element_text(size = 48), legend.text=element_text(size=24), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1, 'cm'))
            }
            if (!blend) {
                plot <- plot + guides(color = NULL)
                cols.grad <- cols
                if (length(x = cols) == 1) {
                  plot <- plot + scale_color_brewer(palette = cols)
                }
                else if (length(x = cols) > 1) {
                  unique.feature.exp <- unique(data.plot[, feature])
                  if (length(unique.feature.exp) == 1) {
                    warning("All cells have the same value (", 
                      unique.feature.exp, ") of ", feature, ".")
                    if (unique.feature.exp == 0) {
                      cols.grad <- cols[1]
                    }
                    else {
                      cols.grad <- cols
                    }
                  }
                  plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad, 
                    guide = "colorbar"))
                }
            }
            if (!(is.null(x = keep.scale)) && keep.scale == "feature" && 
                !blend) {
                max.feature.value <- max(data.single[, feature])
                min.feature.value <- min(data.single[, feature])
                plot <- suppressMessages(plot & scale_color_gradientn(colors = cols, 
                  limits = c(min.feature.value, max.feature.value)))
            }
            if (coord.fixed) {
                plot <- plot + coord_fixed()
            }
            plot <- plot
            plots[[(length(x = features) * (i - 1)) + j]] <- plot
        }
    }
    if (blend) {
        blend.legend <- BlendMap(color.matrix = color.matrix)
        for (ii in 1:length(x = levels(x = data$split))) {
            suppressMessages(expr = plots <- append(x = plots, 
                values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) > 
                  1, yes = levels(x = data$split)[ii], no = "")), 
                  expand = c(0, 0)) + labs(x = features[1], y = features[2], 
                  title = if (ii == 1) {
                    paste("Color threshold:", blend.threshold)
                  } else {
                    NULL
                  }) + no.right), after = 4 * ii - 1))
        }
    }
    plots <- Filter(f = Negate(f = is.null), x = plots)
    if (is.null(x = ncol)) {
        ncol <- 2
        if (length(x = features) == 1) {
            ncol <- 1
        }
        if (length(x = features) > 6) {
            ncol <- 3
        }
        if (length(x = features) > 9) {
            ncol <- 4
        }
    }
    ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol, 
        no = length(x = features))
    legend <- if (blend) {
        "none"
    }
    else {
        split.by %iff% "none"
    }
    if (combine) {
        if (by.col && !is.null(x = split.by) && !blend) {
            plots <- lapply(X = plots, FUN = function(x) {
                return(suppressMessages(expr = x + theme_cowplot() + 
                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""), 
                  limits = ylims) + no.right))
            })
            nsplits <- length(x = levels(x = data$split))
            idx <- 1
            for (i in (length(x = features) * (nsplits - 1) + 
                1):(length(x = features) * nsplits)) {
                plots[[i]] <- suppressMessages(expr = plots[[i]] + 
                  scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]), 
                    limits = ylims) + no.right)
                idx <- idx + 1
            }
            idx <- 1
            for (i in which(x = 1:length(x = plots)%%length(x = features) == 
                1)) {
                plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
                  theme(plot.title = element_text(hjust = 0.5))
                idx <- idx + 1
            }
            idx <- 1 
            if (length(x = features) == 1) {
                for (i in 1:length(x = plots)) {
                  plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
                    theme(plot.title = element_text(hjust = 0.5))
                  idx <- idx + 1
                }
                ncol <- 1
                nrow <- nsplits
            }
            else {
                nrow <- split.by %iff% length(x = levels(x = data$split))
            }
            plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots), 
                f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
            plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
            if (!is.null(x = legend) && legend == "none") {
                plots <- plots & NoLegend()
            }
        }
        else {
            plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% 
                length(x = levels(x = data$split)))
        }
        if (!is.null(x = legend) && legend == "none") {
            plots <- plots & NoLegend()
        }
        if (!(is.null(x = keep.scale)) && keep.scale == "all" && 
            !blend) {
            max.feature.value <- max(data.plot[, features])
            min.feature.value <- min(data.plot[, features])
            plots <- suppressMessages(plots & scale_color_gradientn(colors = cols, 
                limits = c(min.feature.value, max.feature.value)))
        }
    }
    return(plots)
}

uplot <- function(object, cols) {
	p <- ggplot(as.data.frame(object@reductions$umap@cell.embeddings),aes(x=UMAP_1, y=UMAP_2, color=object@active.ident)) + theme_bw() + geom_point(size=0.001)
	return(p)
}

FeaturePlotL <- function (object, features, col, genenames = NULL, dims = c(1, 2), cells = NULL, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA, 
    reduction = NULL, keep.scale = "feature", 
    shape.by = NULL, slot = "data", blend = FALSE, blend.threshold = 0.5, 
    label = FALSE, label.size = 4, repel = FALSE, ncol = NULL, 
    coord.fixed = FALSE, by.col = TRUE, sort.cell = NULL,
    combine = TRUE, raster = NULL, split.by = NULL) 
{
	cols = c("lightgrey", "blue")
    no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
        axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", 
            size = 14, margin = margin(r = 7)))
    reduction <- reduction %||% DefaultDimReduc(object = object)
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% colnames(x = object)
    data <- FetchData(object = object, vars = c(dims, "ident", features), cells = cells, slot = slot)
    features <- colnames(x = data)[4:ncol(x = data)]
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
            feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
            feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
        max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
        ]$maxcolors, no = length(x = cols))
    data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), 
        FUN = function(index) {
            data.feature <- as.vector(x = data[, index])
            min.use <- SetQuantile(cutoff = min.cutoff[index - 
                3], data.feature)
            max.use <- SetQuantile(cutoff = max.cutoff[index - 
                3], data.feature)
            data.feature[data.feature < min.use] <- min.use
            data.feature[data.feature > max.use] <- max.use
            if (brewer.gran == 2) {
                return(data.feature)
            }
            data.cut <- if (all(data.feature == 0)) {
                0
            }
            else {
                as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), breaks = brewer.gran)))
            }
            return(data.cut)
        })
    colnames(x = data)[4:ncol(x = data)] <- features
    rownames(x = data) <- cells
    plots <- vector(mode = "list", length = ifelse(test = blend, yes = 4, no = length(x = features)))
    xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
    ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
    for (j in 1:length(x = features)) {
		data.plot <- data
		data2 <- object
		feature <- features[j]
		data.single <- data.plot[, c(dims, "ident", feature)]
##		data.plot[,"ident"] = as.factor(col)
		aa = col
		aa[which(object@assays$RNA@counts[feature,] == 0)] <- "Not Expressed"
		data2@active.ident <- factor(aa, levels=c(sort(unique(col)), "Not Expressed"))
		data.plot[,"ident"] = factor(aa, levels=c(sort(unique(col)), "Not Expressed"))
        data.single <- data.plot[, c(dims, "ident", feature)]
        plot <- uplot(data2, aa)
        plot <- plot + CenterTitle() + labs(title = genenames[genenames[,2] == features[[j]], 1])  +scale_color_discrete(name="") +
        scale_color_manual(values = c(
"Brombin A, Chandra T, Patton EE (GSE178364)" = "#f8766d",
"Howard AA, et al. (GSE152906)" = "#b79f00",
"Kenny C, Cornell RA (GSE198791)" = "#00ba38",
"Lencer E, Prekeris R, Artinger K (GSE163907)"  = "#00bfc4",
"Saunders LM, Parichy DM, Trapnell C (GSE131136)" = "#619cff",
"Wagner DE, et al. (GSE112294)" = "#f564e3",
"Not Expressed" = "grey"
))
        plots[[j]] <- plot+NoLegend()+ guides(colour = guide_legend(override.aes = list(size=25)))+ theme(axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25))+theme(plot.title = element_text(size=50))
#theme(legend.title=element_blank(), legend.text = element_text(size=40),legend.key.size = unit(4, 'cm'), plot.title = element_text(size=40))+ guides(colour = guide_legend(override.aes = list(size=10)))
    }
    plots <- Filter(f = Negate(f = is.null), x = plots)
    if (is.null(x = ncol)) {
        ncol <- 2
        if (length(x = features) == 1) {
            ncol <- 1
        }
        if (length(x = features) > 6) {
            ncol <- 3
        }
        if (length(x = features) > 9) {
            ncol <- 4
        }
    }
    ncol <- ifelse(test = is.null(x = split.by) || blend, yes = ncol, no = length(x = features))
    legend <- split.by %iff% "none"
    if (combine) {
        if (by.col && !is.null(x = split.by) && !blend) {
            plots <- lapply(X = plots, FUN = function(x) {
                return(suppressMessages(expr = x + theme_cowplot() + 
                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = ""), 
                  limits = ylims) + no.right))
            })
            nsplits <- length(x = levels(x = data$split))
            idx <- 1
            for (i in (length(x = features) * (nsplits - 1) + 
                1):(length(x = features) * nsplits)) {
                plots[[i]] <- suppressMessages(expr = plots[[i]] + 
                  scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]), 
                    limits = ylims) + no.right)
                idx <- idx + 1
            }
            idx <- 1
            for (i in which(x = 1:length(x = plots)%%length(x = features) == 
                1)) {
                plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
                  theme(plot.title = element_text(hjust = 0.5))
                idx <- idx + 1
            }
            idx <- 1 
            if (length(x = features) == 1) {
                for (i in 1:length(x = plots)) {
                  plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
                    theme(plot.title = element_text(hjust = 0.5))
                  idx <- idx + 1
                }
                ncol <- 1
                nrow <- nsplits
            }
            else {
                nrow <- split.by %iff% length(x = levels(x = data$split))
            }
            plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots), 
                f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
            plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
            if (!is.null(x = legend) && legend == "none") {
                plots <- plots & NoLegend()
            }
        }
        else {
            plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% 
                length(x = levels(x = data$split)))
        }
        if (!is.null(x = legend) && legend == "none") {
            plots <- plots & NoLegend()
        }
        if (!(is.null(x = keep.scale)) && keep.scale == "all" && 
            !blend) {
            max.feature.value <- max(data.plot[, features])
            min.feature.value <- min(data.plot[, features])
            plots <- suppressMessages(plots & scale_color_gradientn(colors = cols, 
                limits = c(min.feature.value, max.feature.value)))
        }
    }
    return(plots)
}
#Loading data
z1 <- read.csv("GSE178364_zf_WT_24hpf_Unt-vs-ErbBi.csv", row.names=1)
z2<-Matrix::readMM("GSM4990514_matrix.mtx")
z2 <- as.matrix(z2)
genes <- read.table("GSM4990514_features.tsv")
rownames(z2) <- genes$V1
barcode <- read.table("GSM4990514_barcodes.tsv")
colnames(z2) <- barcode$V1
z3<-Matrix::readMM("GSE131136_all_NC_matrix.mtx")
genes <- read.table("GSE131136_all_NC_genes.tsv", stringsAsFactors=F, sep="\t")
rownames(z3) <- genes$V1
features <- read.table("GSE131136_all_NC_barcodes.tsv", stringsAsFactors=F)
colnames(z3) <- features$V1
z4<-Matrix::readMM("GSM4629198_48h_matrix.mtx")
z4 <- as.matrix(z4)
genes <- read.table("GSM4629198_48h_genes.tsv", stringsAsFactors=F, sep="\t")
rownames(z4) <- genes$V1
features <- read.table("GSM4629198_48h_barcodes.tsv", stringsAsFactors=F)
colnames(z4) <- features$V1
z5<-Matrix::readMM("GSM4629199_68h_matrix.mtx")
z5 <- as.matrix(z5)
genes <- read.table("GSM4629199_68h_genes.tsv", stringsAsFactors=F, sep="\t")
rownames(z5) <- genes$V1
features <- read.table("GSM4629199_68h_barcodes.tsv", stringsAsFactors=F)
colnames(z5) <- features$V1
z7 <- Matrix::readMM("GSM5956969_matrix.mtx")
z7 <- as.matrix(z7)
genes <- read.table("GSM5956969_features.tsv")
rownames(z7) <- genes$V1
barcode <- read.table("GSM5956969_barcodes.tsv")
colnames(z7) <- barcode$V1
z6 <- read.table("Wagner.csv", row.names=1, header=T)
rownames(z6) <- unlist(lapply(X=rownames(z6), FUN=function(x) { return(strsplit(x, "[.]")[[1]][1]) }))
#Aggregate datasets
inters <- Reduce(intersect,list(rownames(z1),rownames(z2),rownames(z3),rownames(z4),rownames(z5), rownames(z6), rownames(z7)))
z1_1 <- z1[rownames(z1) %in% inters,]
z2_1 <- z2[rownames(z2) %in% inters,]
z3_1 <- z3[rownames(z3) %in% inters,]
z4_1 <- z4[rownames(z4) %in% inters,]
z5_1 <- z5[rownames(z5) %in% inters,]
z6_1 <- z6[rownames(z6) %in% inters,]
z7_1 <- z7[rownames(z7) %in% inters,]
s1 <- SingleCellExperiment(list(counts=z1_1))
s2 <- SingleCellExperiment(list(counts=z2_1))
s3 <- SingleCellExperiment(list(counts=z3_1))
s4 <- SingleCellExperiment(list(counts=z4_1))
s5 <- SingleCellExperiment(list(counts=z5_1))
s6 <- SingleCellExperiment(list(counts=z6_1))
s7 <- SingleCellExperiment(list(counts=z7_1))
sizeFactors(s1) <- runif(ncol(z1_1))
sizeFactors(s2) <- runif(ncol(z2_1))
sizeFactors(s3) <- runif(ncol(z3_1))
sizeFactors(s4) <- runif(ncol(z4_1))
sizeFactors(s5) <- runif(ncol(z5_1))
sizeFactors(s6) <- runif(ncol(z6_1))
sizeFactors(s7) <- runif(ncol(z7_1))
#Batch corrections
quick.corrected <- quickCorrect(s1, s2, s3, s4, s5, s6, s7, PARAM=FastMnnParam(BSPARAM=BiocSingular::RandomParam()))
#Seurat processing
pbmc <- CreateSeuratObject(as.matrix(assay(quick.corrected$corrected)), project = "pbmc3k", min.cells = 1)
all.genes <- rownames(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 20000)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) 
pbmc <- RunUMAP(pbmc, dims=1:35)
png("UMAP.png", width=1000, height=1000)
UMAPPlot(pbmc)
dev.off()

reads <- cbind(z1_1, z2_1, z3_1, z4_1, z5_1, z6_1, z7_1)
colnames(reads) <- colnames(pbmc@assays$RNA@counts)
pbmc@assays$RNA@counts <- as.matrix(reads)
pbmc@assays$RNA@data <- as.matrix(reads)


#Loading genes for picture
s_g <- read.table("new_list.txt")
ensembl <- useMart("ensembl")
ensembl = useDataset("drerio_gene_ensembl",mart=ensembl)
sg <- getBM(attributes=c('external_gene_name','ensembl_gene_id'), filters = 'external_gene_name', values = s_g$V1, mart = ensembl)
sg <- sg[order(sg$external_gene_name),]
sg <- sg[!duplicated(sg$external_gene_name),]


png("FeaturePlot.png", height=20000, width=10000)
FeaturePlotN(pbmc, sg, features = sg$ensembl_gene_id, max.cutoff = 12)
dev.off()

col <- c(rep("Brombin A, Chandra T, Patton EE (GSE178364)", ncol(z1_1)), rep("Lencer E, Prekeris R, Artinger K (GSE163907)", ncol(z2_1)),
rep("Saunders LM, Parichy DM, Trapnell C (GSE131136)", ncol(z3_1)),
rep("Howard AA, Baker PA, Ibarra-Garcia-Padilla R, Moore JA, Rivas LJ, Singleton EW, Westheimer JL, Corteguera JA, Tallman JJ, Uribe RA (GSE152906)", ncol(z4_1)),
rep("Howard AA, Baker PA, Ibarra-Garcia-Padilla R, Moore JA, Rivas LJ, Singleton EW, Westheimer JL, Corteguera JA, Tallman JJ, Uribe RA (GSE152906)", ncol(z5_1)),
rep("Wagner DE, Weinreb C, Collins ZM, Megason SG, Klein AM (GSE112294)", ncol(z6_1)),
rep("Kenny C, Cornell RA (GSE198791)", ncol(z7_1)))


png("AuthorsPlot.png", height=20000, width=10000)
FeaturePlotL(pbmc, col, sg, features = sg$ensembl_gene_id, max.cutoff = 12)
dev.off()

