###0-1. rotation 
coord_rotate <- function(seurat.object,theta,clockwise=TRUE) {
  seurat.object$xcoord <- seurat.object@reductions$spatial@cell.embeddings[,1]
  seurat.object$ycoord <- seurat.object@reductions$spatial@cell.embeddings[,2]
  
  coordinates <- seurat.object@meta.data[c("xcoord","ycoord")]
  
  theta=theta*pi/180
  if(clockwise){
    rotation_matrix <- matrix(c(cos(theta), sin(theta),-sin(theta), cos(theta)),nrow = 2)
  }else{
    rotation_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)),nrow = 2)
  }
  
  
  # rotated coordinates
  rotated_coordinates <- as.matrix(coordinates) %*% rotation_matrix
  # new x and y
  seurat.object$x_rotated <- rotated_coordinates[, 1]
  seurat.object$y_rotated <- rotated_coordinates[, 2]
  #set min[x,y] as [0,0]
  if (min(seurat.object$x_rotated)<0){
    seurat.object$x_rotated <- abs(min(seurat.object$x_rotated)) + seurat.object$x_rotated
  }
  if (min(seurat.object$y_rotated)<0){
    seurat.object$y_rotated <- abs(min(seurat.object$y_rotated)) + seurat.object$y_rotated
  }
  if (min(seurat.object$x_rotated)>0){
    seurat.object$x_rotated <-  seurat.object$x_rotated-min(seurat.object$x_rotated)
  }
  if (min(seurat.object$y_rotated)>0){
    seurat.object$y_rotated <-  seurat.object$y_rotated-min(seurat.object$y_rotated)
  }
  
  coordinates_new <- data.frame(xcoord=seurat.object$x_rotated,ycoord=seurat.object$y_rotated)
  
  rownames(coordinates_new) <- colnames(seurat.object)
  colnames(coordinates_new) <- paste0("Spatial_",1:2)
  assay.name=DefaultAssay(seurat.object)
  rotated.object<-CreateSeuratObject(counts=GetAssayData(seurat.object,assay = assay.name,layer = 'counts'),assay = 'Spatial')
  rotated.object[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(coordinates_new), key = "Spatial",assay = "Spatial")
  rotated.object<-AddMetaData(rotated.object,seurat.object@meta.data)
  return(rotated.object)
}

###0-2.mirror (if necessary)
coord_mirrow <- function(seurat.object,by=c('x','y')){
  coordinates <- data.frame(x=seurat.object@reductions$spatial@cell.embeddings[,1],
                            y=seurat.object@reductions$spatial@cell.embeddings[,2])
  if(length(by)==2){
    coordinates$x <- abs(coordinates$x-max(coordinates$x))
    coordinates$y <- abs(coordinates$y-max(coordinates$y))
  }else{
    if(by=='x'){
      coordinates$x <- abs(coordinates$x-max(coordinates$x))
    }
    if(by=='y'){
      coordinates$y <- abs(coordinates$y-max(coordinates$y))
    }
  }
  
  rownames(coordinates) <- colnames(seurat.object)
  colnames(coordinates) = paste0("Spatial_",1:2)
  assay.name=DefaultAssay(seurat.object)
  mirror.object<-CreateSeuratObject(counts=GetAssayData(seurat.object,assay = assay.name,layer = 'counts'),assay = 'Spatial')
  mirror.object[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(coordinates), key = "Spatial",assay = "Spatial")
  AddMetaData(mirror.object,seurat.object@meta.data) -> mirror.object
  return(mirror.object)
}


###1. set origin (if necessary)
coord_set0 <- function(seurat.object) {
  seurat.object$xcoord <- seurat.object@reductions$spatial@cell.embeddings[,1]
  seurat.object$ycoord <- seurat.object@reductions$spatial@cell.embeddings[,2]
  
  #set min[x,y] as [0,0]
  if (min(seurat.object$xcoord)<0){
    seurat.object$xcoord_new <- abs(min(seurat.object$xcoord)) + seurat.object$xcoord
  }
  if (min(seurat.object$ycoord)<0){
    seurat.object$ycoord_new <- abs(min(seurat.object$ycoord)) + seurat.object$ycoord
  }
  if (min(seurat.object$xcoord)>0){
    seurat.object$xcoord_new <-  seurat.object$xcoord-min(seurat.object$xcoord)
  }
  if (min(seurat.object$ycoord)>0){
    seurat.object$ycoord_new <-  seurat.object$ycoord-min(seurat.object$ycoord)
  }
  
  coordinates_new <- data.frame(xcoord=seurat.object$xcoord_new,ycoord=seurat.object$ycoord_new)
  
  rownames(coordinates_new) <- colnames(seurat.object)
  colnames(coordinates_new) <- paste0("Spatial_",1:2)
  assay.name=DefaultAssay(seurat.object)
  rotated.object<-CreateSeuratObject(counts=GetAssayData(seurat.object,assay = assay.name,layer = 'counts'),assay = 'Spatial')
  rotated.object[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(coordinates_new), key = "Spatial",assay = "Spatial")
  rotated.object<-AddMetaData(rotated.object,seurat.object@meta.data)
  return(rotated.object)
}


###2.transition from pixel to physical location, with built-in coord_rotate function
coord_trans<-function(seurat.object,conversion,theta,clockwise=TRUE){
  seurat.object$xcoord <- seurat.object@reductions$spatial@cell.embeddings[,1]
  seurat.object$ycoord <- seurat.object@reductions$spatial@cell.embeddings[,2]
  coordinates <- data.frame(xcoord=seurat.object$xcoord*conversion,ycoord=seurat.object$ycoord*conversion)
  
  # coordinates <- seurat.object@meta.data[c("xcoord","ycoord")]
  
  theta=theta*pi/180
  if(clockwise){
    rotation_matrix <- matrix(c(cos(theta), sin(theta),-sin(theta), cos(theta)),nrow = 2)
  }else{
    rotation_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)),nrow = 2)
  }
  
  
  # rotated coordinates
  rotated_coordinates <- as.matrix(coordinates) %*% rotation_matrix
  # new x and y
  seurat.object$x_rotated <- rotated_coordinates[, 1]
  seurat.object$y_rotated <- rotated_coordinates[, 2]
  #set min[x,y] as [0,0]
  if (min(seurat.object$x_rotated)<0){
    seurat.object$x_rotated <- abs(min(seurat.object$x_rotated)) + seurat.object$x_rotated
  }
  if (min(seurat.object$y_rotated)<0){
    seurat.object$y_rotated <- abs(min(seurat.object$y_rotated)) + seurat.object$y_rotated
  }
  if (min(seurat.object$x_rotated)>0){
    seurat.object$x_rotated <-  seurat.object$x_rotated-abs(min(seurat.object$x_rotated))
  }
  if (min(seurat.object$y_rotated)>0){
    seurat.object$y_rotated <-  seurat.object$y_rotated-abs(min(seurat.object$y_rotated))
  }
  
  coordinates_new <- data.frame(xcoord=seurat.object$x_rotated,ycoord=seurat.object$y_rotated)
  # coordinates_new <-data.frame(xcoord=as.numeric(coordinates_new$xcoord),ycoord=as.numeric(coordinates_new$ycoord))
  rownames(coordinates_new) <- colnames(seurat.object)
  colnames(coordinates_new) <- paste0("Spatial_",1:2)
  assay.name=DefaultAssay(seurat.object)
  rotated.object<-CreateSeuratObject(counts=GetAssayData(seurat.object,assay = assay.name,layer = 'counts'),assay = 'Spatial')
  rotated.object[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(coordinates_new), key = "Spatial",assay = "Spatial")
  rotated.object<-AddMetaData(rotated.object,seurat.object@meta.data)
  return(rotated.object)
}

###3. Scale gene expression
calc_integral <- function(x, y) {
  # 输入检查
  if (length(x) != length(y)) {
    stop("x and y must have the same length.")
  }
  
  if (length(x) < 2) {
    stop("x and y must have at least two points.")
  }
  
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors.")
  }
  
  if (any(diff(x) <= 0)) {
    stop("x must be strictly increasing.")
  }
  
  # 计算相邻点之间的差异（dx）
  diff_x <- diff(x)
  # 计算相邻y值的平均值（高度）
  avg_y <- (y[-1] + y[-length(y)]) / 2
  # 计算每个梯形的面积并求和
  integral <- sum(diff_x * avg_y)
  return(integral)
}

###4. calculate diffusion distance
calculate_hwhm <- function(distances, values) {
  # 输入检查
  if (length(distances) != length(values)) stop("x和y的长度必须相同")
  
  # 找到峰值位置和半高值
  # max_y <- max(values)
  # half_max <- max_y / 2
  
  interpolation <- approx(distances, values, n = 100)
  plot(distances,values)
  
  x_interp<-interpolation$x
  y_interp<-interpolation$y
  plot(x_interp,y_interp,cex = 0.5)
  
  
  # Find the two x values (distances) at which y equals the half maximum.
  # Note that we have to check if the values cross the half_max from below and from above
  peak_idx=which.max(y_interp)
  half_max <- max(y_interp) / 2
  cross_below <- which(diff(y_interp[1:peak_idx] >= half_max) != 0)
  cross_above <- which(diff(y_interp[peak_idx:length(y_interp)] <= half_max) != 0) + (peak_idx - 1)
  # cross_below <- which(interpolation$y[-length(interpolation$y)] < half_max & interpolation$y[-1] >= half_max)
  # cross_above <- which(interpolation$y[-length(interpolation$y)] > half_max & interpolation$y[-1] <= half_max)
  cross <- sort(c(cross_below, cross_above))
  if (length(cross) < 2) stop("半高点需对称")
  # points(interpolation$x[c(cross_below,cross_above)], interpolation$y[c(cross_below,cross_above)], col = "red", pch = 19, cex = 1)
  
  points(x_interp[peak_idx], y_interp[peak_idx], col = "black", pch = 19, cex = 1)
  
  # left_idx=max(cross[cross<peak_idx])
  left_idx=median(cross[cross<peak_idx])
  left<-x_interp[left_idx]
  
  # right_idx = min(cross[cross>peak_idx])
  right_idx = median(cross[cross>peak_idx])
  right<-x_interp[right_idx]
  
  points(c(left,right), y_interp[c(left_idx,right_idx)], col = "red", pch = 19, cex = 1)
  
  peak_x<- x_interp[peak_idx] # 线性插值
  HWHM_left = peak_x - left
  HWHM_right = right - peak_x
  
  # 返回左右HWHM
  return(c(
    HWHM_left = HWHM_left, 
    HWHM_right = HWHM_right,
    HWHM_mean = (HWHM_left+HWHM_right)/2,
    FWHM = HWHM_left+HWHM_right
  ))
}

#4.1
calculate_fwhm <- function(distances, values) {
  # Compute the maximum and half-maximum.
  max_y <- max(values)
  # half_max <- max_y / 2
  plot(distances,values)
  # Linear interpolation of your data.
  interpolation <- approx(distances, values,n = 200)
  plot(interpolation$x,interpolation$y)
  half_max <- max(interpolation$y) / 2
  # Find the two x values (distances) at which y equals the half maximum.
  # Note that we have to check if the values cross the half_max from below and from above
  cross_below <- which(interpolation$y[-length(interpolation$y)] < half_max & interpolation$y[-1] >= half_max)
  cross_above <- which(interpolation$y[-length(interpolation$y)] > half_max & interpolation$y[-1] <= half_max)
  cross <- sort(c(cross_below, cross_above))
  
  points(interpolation$x[which.max(interpolation$y)], interpolation$y[which.max(interpolation$y)], col = "black", pch = 19, cex = 1.5)
  
  if(length(cross) < 2){
    stop("Couldn't find two crossing points. Check your data or increase the number of points in interpolation.")
  }
  
  points(interpolation$x[cross_below], interpolation$y[cross_below], col = "red", pch = 19, cex = 1.5)
  # The FWHM is the difference between the second and the first crossing points.
  #interpolation$x[cross[cross>=which.max(interpolation$y)]][1] -> maxv
  #maxv[maxv<220] ->  maxv
  interpolation$x[cross[cross<=which.max(interpolation$y)]][1] -> minv
  fwhm <-  interpolation$x[which.max(interpolation$y)]-minv
  return(fwhm)
}

