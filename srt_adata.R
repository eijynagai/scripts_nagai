srt2adata <- function(so){
    require(SCP)
    if (class(so) != 'Seurat'){
        print("Error: Not a Seurat object")
    } else {
        adata <- srt_to_adata(so)
        return(adata)
    }
}
