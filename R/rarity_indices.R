######################################
######    Irr calculation      #######
######################################

setGeneric("Irr", 
           def = function(assemblages,
                          W,
                          abundance = F,
                          Wmin = min(W),
                          Wmax = max(W)
           )
           {
             standardGeneric( "Irr" )
           }
)

setMethod("Irr",
          signature(assemblages = "vector", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), Wmax = max(W))
          {
            if(length(assemblages) != length(W))
            {
              cat("Remark: Number of species different between assemblages and W, combining with species names...\n")
              if(any(!(names(assemblages) %in% names(W))))
              {
                stop("Species names of assemblages do not correspond to species names of W")
              }
              W <- W[match(names(assemblages), names(W))]
            }
            if(sum(assemblages) == 0)
            {
              IrrValue <- 0
            } else
            {
              PA.assemblages <- assemblages
              PA.assemblages[PA.assemblages > 0] <- 1
              
              if(abundance)
              {
                IrrValue <- ((sum(assemblages * W) / sum(assemblages)) - Wmin) / (Wmax - Wmin) # Formula with abundance
              } else
              {
                IrrValue <- ((sum(PA.assemblages * W) / sum(PA.assemblages)) - Wmin) / (Wmax - Wmin) # Formula with presence-absence data       
              }      
            }
            return(c(Irr = IrrValue, Richness = sum(PA.assemblages)))
          }
)

setMethod("Irr",
          signature(assemblages = "vector", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            if(length(assemblages) != nrow(W))
            {
              cat("Remark: Number of species different between assemblages and W, combining with species names...\n")
              if(any(!(names(assemblages) %in% rownames(W))))
              {
                stop("Species names of assemblages do not correspond to species names of W")
              }
              W <- W[match(names(assemblages), rownames(W)), ]
            }
            if(any(colnames(W) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              W <- W[, -which(colnames(W) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))), drop = F]
            }
            
            IrrValue <- NULL
            for (x1 in 1:ncol(W))
            {
              if(sum(assemblages) == 0)
              {
                IrrValue <- c(IrrValue,
                              0)
              } else
              {
                PA.assemblages <- assemblages
                PA.assemblages[PA.assemblages > 0] <- 1
                
                if(abundance)
                {
                  IrrValue <-  c(IrrValue,
                                 ((sum(assemblages * W[, x1]) / sum(assemblages)) - Wmin[x1]) / (Wmax[x1] - Wmin[x1])) # Formula with abundance
                } else
                {
                  IrrValue <-  c(IrrValue,
                                 ((sum(PA.assemblages * W[, x1]) / sum(PA.assemblages)) - Wmin[x1]) / (Wmax[x1] - Wmin[x1])) # Formula with presence-absence data       
                }   
              }
              names(IrrValue)[length(IrrValue)] <- paste("Irr_", colnames(W)[x1], sep ="")
            }
            IrrValue <- c(IrrValue,
                          Richness = sum(PA.assemblages))
            return(IrrValue)
          }
)

setMethod("Irr",
          signature(assemblages = "vector", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            W <- as.matrix(W)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)

setMethod("Irr",
          signature(assemblages = "matrix", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), Wmax = max(W))
          {
            if(nrow(assemblages) != length(W))
            {
              cat("Remark: Number of species different between assemblages and W, combining with species names...\n")
              if(any(!(rownames(assemblages) %in% names(W))))
              {
                stop("Species names of assemblages do not correspond to species names of W")
              }
              W <- W[match(rownames(assemblages), names(W))]
            }
            PA.assemblages <- assemblages
            PA.assemblages[PA.assemblages > 0] <- 1
            
            if(abundance)
            {
              IrrValue <- apply(assemblages, 2, function(x, Weights, wmin, wmax)
              {
                if(sum(x) == 0)
                {
                  0
                } else
                {
                  (((sum(x * Weights) / sum(x)) - wmin) / (wmax - wmin))
                }
              }, Weights = W, wmin = Wmin, wmax = Wmax) # Formula with abundance
            } else
            {
              IrrValue <- apply(PA.assemblages, 2, function(x, Weights, wmin, wmax)
              {
                if(sum(x) == 0)
                {
                  0
                } else
                {
                  (((sum(x * Weights) / sum(x)) - wmin) / (wmax - wmin))
                }
              }, Weights = W, wmin = Wmin, wmax = Wmax) # Formula with presence-absence data       
            }
            
            IrrValue <- cbind(IrrValue,
                              Richness = apply(PA.assemblages, 2, sum))
            return(IrrValue)
          }
)

setMethod("Irr",
          signature(assemblages = "matrix", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            if(nrow(assemblages) != nrow(W))
            {
              cat("Remark: Number of species different between assemblages and W, combining with species names...\n")
              if(any(!(rownames(assemblages) %in% rownames(W))))
              {
                stop("Species names of assemblages do not correspond to species names of W")
              }
              W <- W[match(rownames(assemblages), rownames(W)), ]
            }
            if(any(colnames(W) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))))
            {
              W <- W[, -which(colnames(W) %in% c("Q", "R", "cut.off", paste("Q", 1:1000, sep = ""), paste("R", 1:1000, sep = ""),  paste("cut.off", 1:1000, sep = ""))), drop = F]
            }
            
            IrrValue <- NULL
            PA.assemblages <- assemblages
            PA.assemblages[PA.assemblages > 0] <- 1
            
            for (x1 in 1:ncol(W))
            {
              if(abundance)
              {
                IrrValue <- cbind(IrrValue,
                                  apply(assemblages, 2, function(x, Weights, wmin, wmax)
                                  {
                                    if(sum(x) == 0)
                                    {
                                      0
                                    } else
                                    {
                                      (((sum(x * Weights) / sum(x)) - wmin) / (wmax - wmin))
                                    }
                                  }, Weights = W[, x1], wmin = Wmin[x1], wmax = Wmax[x1]) # Formula with abundance
                )
              } else
              {
                IrrValue <- cbind(IrrValue,
                                  apply(PA.assemblages, 2, function(x,  Weights, wmin, wmax)
                                  {
                                    if(sum(x) == 0)
                                    {
                                      0
                                    } else
                                    {
                                      (((sum(x * Weights) / sum(x)) - wmin) / (wmax - wmin))
                                    }
                                  }, Weights = W[, x1], wmin = Wmin[x1], wmax = Wmax[x1]) # Formula with presence-absence data   
                )
              }
              colnames(IrrValue)[ncol(IrrValue)] <- paste("Irr_", colnames(W)[x1], sep ="")
              
            }
            IrrValue <- cbind(IrrValue,
                              Richness = apply(PA.assemblages, 2, sum))
            return(IrrValue)
          }
)

setMethod("Irr",
          signature(assemblages = "matrix", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            W <- as.matrix(W)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)

setMethod("Irr",
          signature(assemblages = "data.frame", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), Wmax = max(W))
          {
            assemblages <- as.matrix(assemblages)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)

setMethod("Irr",
          signature(assemblages = "data.frame", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            assemblages <- as.matrix(assemblages)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)

setMethod("Irr",
          signature(assemblages = "data.frame", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), Wmax = apply(W, 2, max))
          {
            assemblages <- as.matrix(assemblages)
            W <- as.matrix(W)
            Irr(assemblages, W, abundance = abundance, Wmin = Wmin, Wmax = Wmax)
          }
)