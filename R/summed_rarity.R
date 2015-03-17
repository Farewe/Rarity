######################################
######    Isr calculation      #######
######################################

setGeneric("Isr", 
           def = function(assemblages,
                          W,
                          abundance = FALSE,
                          Wmin = min(W),
                          normalise = FALSE
           )
           {
             standardGeneric( "Isr" )
           }
)

setMethod("Isr",
          signature(assemblages = "vector", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), normalise = F)
          {
            if(length(assemblages) != length(W))
            {
              cat("Remark: Number of species different between assemblages and W, combining with species names...\n")
              if(any(!(names(assemblages) %in% names(W))))
              {
                stop("Species names of assemblages do not correspond to species names of W")
              }
              full.W <- W
              W <- W[match(names(assemblages), names(W))]
            }
            if(sum(assemblages) == 0)
            {
              IsrValue <- 0
            } else
            {
              PA.assemblages <- assemblages
              PA.assemblages[PA.assemblages > 0] <- 1
                            
              
              if(abundance)
              {
                if(standardise)
                {
                  stop("The index cannot be standardised between 0 and 1 when abundance = TRUE. Set abundance = FALSE or normalise = FALSE")
                }
                IsrValue <- sum(assemblages * W) # Formula with abundance
              } else
              {
                if(normalise)
                {
                  IsrValue <- sum(PA.assemblages * W - Wmin) / (sum(full.W - Wmin)) # Formula with presence-absence data 
                } else
                {
                  IsrValue <- sum(PA.assemblages * W) # Formula with presence-absence data 
                }     
              }      
            }
            return(c(Isr = IsrValue, Richness = sum(PA.assemblages)))
          }
)

setMethod("Isr",
          signature(assemblages = "vector", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            if(length(assemblages) != nrow(W))
            {
              cat("Remark: Number of species different between assemblages and W, combining with species names...\n")
              if(any(!(names(assemblages) %in% rownames(W))))
              {
                stop("Species names of assemblages do not correspond to species names of W")
              }
              full.W <- W
              W <- W[match(names(assemblages), rownames(W)), ]
            }
            if(any(colnames(W) %in% c("Q", "R", paste("Q", 1:100, sep = ""), paste("R", 1:100, sep = ""))))
            {
              full.W <- full.W[, -which(colnames(full.W) %in% c("Q", "R", paste("Q", 1:100, sep = ""), paste("R", 1:100, sep = ""))), drop = F]
              W <- W[, -which(colnames(W) %in% c("Q", "R", paste("Q", 1:100, sep = ""), paste("R", 1:100, sep = ""))), drop = F]
            }
            
            IsrValue <- NULL
            for (x1 in 1:ncol(W))
            {
              if(sum(assemblages) == 0)
              {
                IsrValue <- c(IsrValue,
                              0)
              } else
              {
                PA.assemblages <- assemblages
                PA.assemblages[PA.assemblages > 0] <- 1
                
                if(abundance)
                {
                  if(standardise)
                  {
                    stop("The index cannot be standardised between 0 and 1 when abundance = TRUE. Set abundance = FALSE or normalise = FALSE")
                  }
                  IsrValue <-  c(IsrValue,
                                 sum(assemblages * W[, x1])) # Formula with abundance
                } else
                {
                  if(normalise)
                  {
                    IsrValue <-  c(IsrValue,
                                   sum(PA.assemblages * W[, x1] - Wmin[x1]) / (sum(full.W[, x1]) - Wmin[x1])) # Formula with presence-absence data       
                  } else
                  {
                    IsrValue <-  c(IsrValue,
                                   sum(PA.assemblages * W[, x1])) # Formula with presence-absence data       
                  }
                }
              }
              names(IsrValue)[length(IsrValue)] <- paste("Isr_", colnames(W)[x1], sep ="")
            }
            IsrValue <- c(IsrValue,
                          Richness = sum(PA.assemblages))
            return(IsrValue)
          }
)

setMethod("Isr",
          signature(assemblages = "vector", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            W <- as.matrix(W)
            if(length(assemblages) != nrow(W))
            {
              cat("Remark: Number of species different between assemblages and W, combining with species names...\n")
              if(any(!(names(assemblages) %in% rownames(W))))
              {
                stop("Species names of assemblages do not correspond to species names of W")
              }
              full.W <- W
              W <- W[match(names(assemblages), rownames(W)), ]
            }
            if(any(colnames(W) %in% c("Q", "R", paste("Q", 1:100, sep = ""), paste("R", 1:100, sep = ""))))
            {
              full.W <- full.W[, -which(colnames(full.W) %in% c("Q", "R", paste("Q", 1:100, sep = ""), paste("R", 1:100, sep = ""))), drop = F]
              W <- W[, -which(colnames(W) %in% c("Q", "R", paste("Q", 1:100, sep = ""), paste("R", 1:100, sep = ""))), drop = F]
            }
            
            IsrValue <- NULL
            for (x1 in 1:ncol(W))
            {
              if(sum(assemblages) == 0)
              {
                IsrValue <- c(IsrValue,
                              0)
              } else
              {
                PA.assemblages <- assemblages
                PA.assemblages[PA.assemblages > 0] <- 1
                
                if(abundance)
                {
                  if(standardise)
                  {
                    stop("The index cannot be standardised between 0 and 1 when abundance = TRUE. Set abundance = FALSE or normalise = FALSE")
                  }
                  IsrValue <-  c(IsrValue,
                                 sum(assemblages * W[, x1])) # Formula with abundance
                } else
                {
                  if(normalise)
                  {
                    IsrValue <-  c(IsrValue,
                                   sum(PA.assemblages * W[, x1] - Wmin[x1]) / (sum(full.W[, x1]) - Wmin[x1])) # Formula with presence-absence data       
                  } else
                  {
                    IsrValue <-  c(IsrValue,
                                   sum(PA.assemblages * W[, x1])) # Formula with presence-absence data       
                  }
                }                   
              }
              names(IsrValue)[length(IsrValue)] <- paste("Isr_", colnames(W)[x1], sep ="")
            }
            IsrValue <- c(IsrValue,
                          Richness = sum(assemblages))
            return(IsrValue)
          }
)

setMethod("Isr",
          signature(assemblages = "matrix", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), normalise = F)
          {
            if(nrow(assemblages) != length(W))
            {
              cat("Remark: Number of species different between assemblages and W, combining with species names...\n")
              if(any(!(rownames(assemblages) %in% names(W))))
              {
                stop("Species names of assemblages do not correspond to species names of W")
              }
              full.W <- W
              W <- W[match(rownames(assemblages), names(W))]
            }
            PA.assemblages <- assemblages
            PA.assemblages[PA.assemblages > 0] <- 1
            
            
            if(abundance)
            {
              if(normalise)
              {
                stop("The index cannot be standardised between 0 and 1 when abundance = TRUE. Set abundance = FALSE or normalise = FALSE")
              }
              IsrValue <- apply(assemblages, 2, function(x, Weights)
              {
                sum(x * Weights)
              }, Weights = W) # Formula with abundance
            } else
            {
              if(normalise)
              {
                IsrValue <- apply(PA.assemblages, 2, function(x, Weights, full.W, wmin)
                {
                  sum(x * Weights - wmin) / (sum(full.W - wmin))
                }, Weights = W, full.W = full.W, wmin = Wmin) # Formula with presence-absence data       
              } else {
                IsrValue <- apply(PA.assemblages, 2, function(x, Weights)
                {
                  sum(x * Weights)
                }, Weights = W)
              }
              
            }
            
            IsrValue <- cbind(IsrValue,
                              Richness = apply(PA.assemblages, 2, sum))
            return(IsrValue)
          }
)

setMethod("Isr",
          signature(assemblages = "matrix", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            if(nrow(assemblages) != nrow(W))
            {
              cat("Remark: Number of species different between assemblages and W, combining with species names...\n")
              if(any(!(rownames(assemblages) %in% rownames(W))))
              {
                stop("Species names of assemblages do not correspond to species names of W")
              }
              full.W <- W
              W <- W[match(rownames(assemblages), rownames(W)), ]
            }
            if(any(colnames(W) %in% c("Q", "R", paste("Q", 1:100, sep = ""), paste("R", 1:100, sep = ""))))
            {
              full.W <- full.W[, -which(colnames(full.W) %in% c("Q", "R", paste("Q", 1:100, sep = ""), paste("R", 1:100, sep = ""))), drop = F]
              W <- W[, -which(colnames(W) %in% c("Q", "R", paste("Q", 1:100, sep = ""), paste("R", 1:100, sep = ""))), drop = F]
            }
            
            IsrValue <- NULL
            PA.assemblages <- assemblages
            PA.assemblages[PA.assemblages > 0] <- 1
            
            for (x1 in 1:ncol(W))
            {
              if(abundance)
              {
                if(normalise)
                {
                  stop("The index cannot be standardised between 0 and 1 when abundance = TRUE. Set abundance = FALSE or normalise = FALSE")
                }
                IsrValue <- cbind(IsrValue,
                                  apply(assemblages, 2, function(x, Weights)
                                  {
                                    sum(x * Weights)
                                  }, Weights = W[, x1]) # Formula with abundance
                )
              } else
              {
                if(normalise)
                {
                IsrValue <- cbind(IsrValue,
                                  apply(PA.assemblages, 2, function(x, Weights, full.W, wmin)
                                  {
                                    sum(x * Weights - wmin) / (sum(full.W - wmin))
                                  }, Weights = W[, x1], full.W = full.W[, x1], wmin = Wmin[x1]) # Formula with presence-absence data   
                )
                } else {
                  IsrValue <- cbind(IsrValue,
                                    apply(PA.assemblages, 2, function(x, Weights)
                                    {
                                      sum(x * Weights)
                                    }, Weights = W[, x1]) # Formula with presence-absence data   
                  )
                }
              }
              colnames(IsrValue)[ncol(IsrValue)] <- paste("Isr_", colnames(W)[x1], sep ="")
            }
            IsrValue <- cbind(IsrValue,
                              Richness = apply(PA.assemblages, 2, sum))
            return(IsrValue)
          }
)

setMethod("Isr",
          signature(assemblages = "matrix", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            W <- as.matrix(W)
            Isr(assemblages, W, abundance = abundance, Wmin = Wmin, normalise = normalise)
          }
)

setMethod("Isr",
          signature(assemblages = "data.frame", W = "vector"),
          function(assemblages, W, abundance = F, Wmin = min(W), normalise = F)
          {
            assemblages <- as.matrix(assemblages)
            Isr(assemblages, W = W, abundance = abundance, Wmin = Wmin, normalise = normalise)
          }
)

setMethod("Isr",
          signature(assemblages = "data.frame", W = "matrix"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            assemblages <- as.matrix(assemblages)
            Isr(assemblages, W, abundance = abundance, Wmin = Wmin, normalise = normalise)
          }
)

setMethod("Isr",
          signature(assemblages = "data.frame", W = "data.frame"),
          function(assemblages, W, abundance = F, Wmin = apply(W, 2, min), normalise = F)
          {
            assemblages <- as.matrix(assemblages)
            W <- as.matrix(W)
            Isr(assemblages, W, abundance = abundance, Wmin = Wmin, normalise = normalise)
          }
)