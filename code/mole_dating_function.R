mole_dating <- function(prefix){
  
  main_dir_file <- paste0(main_dir,dir_file)
  
  metadata <- data
  
  model_to_test <- model
  
  iter <- no_iter
  
  chr <- chrom

  #loading the results from gubbins in 
  cluster_name <- prefix
  
  clust_no <- sub("^clust","", cluster_name)
  
  t <- loadGubbins(prefix = paste0(main_dir_file,prefix,".",chr))
  
  tip_names <- t$tip.label
  
  #find the corresponding cluster no of each tip
  if (any(grepl("#", tip_names, fixed = TRUE))) {
    
    cluster <- metadata$cluster_poppunk_updated[match(tip_names, paste0(metadata$name_Sanger, "_1"))]
    
    #find the index of the cluster that is not the same as cluster of interest (i.e outgroup)
    outgroup_ind <- which(cluster != clust_no)
    
    #get the name of the outgroup
    outgroup <- tip_names[outgroup_ind]
    
    #drop the outgroup from the input tree
    tree <- drop.tip(t, outgroup)
    
    #find the collection date of each tip
    d <- metadata$year[match(tree$tip.label, paste0(metadata$name_Sanger, "_1"))]
    
  } else {
    
    #get the name of the outgroup
    outgroup <- c("Reference")
    
    #drop the outgroup from the input tree
    tree <- drop.tip(t, outgroup)
    
    #find the collection date of each tip
    d <- metadata$year[match(tree$tip.label, metadata$name)]
    
  }
  
  d <- as.numeric(d)
  
  #try to root the tree at the most optimised position
  rooted <- initRoot(tree, d)
  
  dated <- bactdate(rooted, d, useRec = T, updateRoot = F, model = model_to_test, nbIts = iter)
  
  return(dated)
}