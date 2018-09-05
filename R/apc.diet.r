apc.diet <- function (x, preyfile, palette = NULL, check = FALSE) 
{

  if (!inherits(x, "diet"))
    stop("Not an object of diet class.")
 
  if(is.character(preyfile))  
     preynms <- read.csv(preyfile)
  else
    preynms <- preyfile
  
  

  # Sort the prey if there is a sort column labelled "Sort"
  if("Sort" %in% names(preynms))
     preynms <- preynms[order(preynms$Sort),] 
  levls <- as.vector(preynms$PreyTaxonSort)
  preyID <- match(levels(x$Group), preynms$Group)

  if(any(is.na(preyID)) | length(levls) != length(levels(x$Group)))
    stop("Number of prey groups does not match the number in the dataset. Check your preyfile.")
  levels(x$Group) <- preynms$PreyTaxonSort[preyID]
  x$Group <- factor(as.vector(x$Group), levels = levls)  
  plen <- length(unique(preynms$PreyTaxBroad))

  # Palettes
  if (plen > 4 & is.null(palette)) {
    rgb.palette <- colorRampPalette(c("red", "orange",
                                      "blue"), space = "rgb")
  
    node.cols <- rgb.palette(nrow(preynms))
    names(node.cols) <- levels(preynms$PreyTaxonSort)
    node.cols <- node.cols[match(as.vector(preynms$PreyTaxonSort), names(node.cols))]
  }
  else if (plen <= 4 & is.null(palette)) {
    if (plen > 1) {      
      grps <- table(preynms$PreyTaxBroad)
      if (plen == 2) {
        assigncols <- Quote(c(topo.colors(2 * n1)[1:n1], 
                              terrain.colors(2 * n2)[1:n2]))
        this.assigncols <- do.call("substitute", list(assigncols, 
                                                      list(n1 = grps[1], n2 = grps[2])))
      }
      else if (plen == 3) {
        assigncols <- Quote(c(topo.colors(2 * n1)[1:n1], 
                              terrain.colors(2 * n2)[1:n2], heat.colors(2 * 
                                n3)[1:n3]))
        this.assigncols <- do.call("substitute", list(assigncols, 
                                                      list(n1 = grps[1], n2 = grps[2], n3 = grps[3])))
      }
      else {
        assigncols <- Quote(c(topo.colors(2 * n1)[1:n1], 
                              terrain.colors(2 * n2)[1:n2], heat.colors(2 * 
                                n3)[1:n3], cm.colors(2 * n4)[1:n4]))
        this.assigncols <- do.call("substitute", list(assigncols, 
                                                      list(n1 = grps[1], n2 = grps[2], n3 = grps[3], 
                                                           n4 = grps[4])))
      }
      node.cols <- eval(this.assigncols)
      orderprey <- preynms$Sort[order(preynms$PreyTaxBroad)]
      names(node.cols) <- as.vector(preynms$PreyTaxonSort)[orderprey]
      node.cols <- node.cols[match(as.vector(preynms$PreyTaxonSort), names(node.cols))]
    }
    else {
      node.cols <- topo.colors(nrow(preynms))
      names(node.cols) <- as.vector(preynms$PreyTaxonSort)
     
    }
   

  }
  else {
    assigncols <- Quote(p)
    node.cols <- do.call("substitute", list(assigncols, list(p = palette)))
    if (length(node.cols) != length(as.vector(preynms$PreyTaxonSort))) 
      stop("Your palette does not reflect the number of taxa present in the data.")
    names(node.cols) <- as.vector(preynms$PreyTaxonSort)
  }

  if(check){
   
    g_legend <- function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      legend
    }
    
    z <- rep(1, length(node.cols))
    y <- 1:length(node.cols)
    cols <- rev(node.cols)
    df <- data.frame(x = z, y = y, Species = names(cols))
    p <- ggplot(df, aes(x, y, color = Species)) + geom_point(shape = 15, size = 10)  +
      scale_color_manual(values = cols) + xlab("") + ylab("") +
      theme_void() + ggtitle("Species Colours")  + 
      theme(plot.title = element_text(hjust = 0.5), plot.margin = unit(c(0,0,0,0), "mm"))
    pleg <- g_legend(p)
    grid.arrange(pleg)
    
}
  
  list(cols = node.cols, x = x)
}