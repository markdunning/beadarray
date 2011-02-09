checkSwathStatus <- function(swathFileNames, txtSuff = ".txt", grnSuff = "_Grn.tif") {
    
    ## check wether we have created a text file for each swath
    ## gets passed a list of files with "Swath" in them.
    ## Sorts out any green channel images and tries to match them with text files

    greenImages <- swathFileNames[grep(grnSuff, swathFileNames)]
    txtFiles <- swathFileNames[grep(txtSuff, swathFileNames)]

    ## if there aren't any text files just return false
    if(length(txtFiles) == 0) {
        return(FALSE)
    }
    else {
        ## strip the suffixes, leaving only the array address id's and swaths
        greenImages <- unlist(strsplit(greenImages, grnSuff));
        txtFiles <- unlist(strsplit(txtFiles, txtSuff));
        
        ## if they match perfectly, return TRUE
        if(length(match(txtFiles, greenImages)) == length(greenImages)) {
            return(TRUE);
        }
        else {
            return(FALSE);
        }
    }
}


