### from https://github.com/GreenleafLab/ArchR/blob/master/R/ValidationUtils.R

validInput <- function(input = NULL, name = NULL, valid = NULL){

  valid <- unique(valid)

  if(is.character(valid)){
    valid <- tolower(valid)
  }else{
    stop("Validator must be a character!")
  }

  if(!is.character(name)){
    stop("name must be a character!")
  }

  if("null" %in% tolower(valid)){
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }

  av <- FALSE

  for(i in seq_along(valid)){

    vi <- valid[i]

    if(vi == "integer" | vi == "wholenumber"){

      if(all(is.numeric(input))){
        #https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
        cv <- min(abs(c(input%%1, input%%1-1)), na.rm = TRUE) < .Machine$double.eps^0.5
      }else{
        cv <- FALSE
      }

    }else if(vi == "null"){

      cv <- is.null(input)

    }else if(vi == "bool" | vi == "boolean" | vi == "logical"){

      cv <- is.logical(input)

    }else if(vi == "numeric"){

      cv <- is.numeric(input)

    }else if(vi == "vector"){

      cv <- is.vector(input)

    }else if(vi == "matrix"){

      cv <- is.matrix(input)

    }else if(vi == "sparsematrix"){

      cv <- is(input, "dgCMatrix")

    }else if(vi == "character"){

      cv <- is.character(input)

    }else if(vi == "factor"){

      cv <- is.factor(input)

    }else if(vi == "rlecharacter"){

      cv1 <- is(input, "Rle")
      if(cv1){
        cv <- is(input@values, "factor") || is(input@values, "character")
      }else{
        cv <- FALSE
      }

    }else if(vi == "timestamp"){

      cv <- is(input, "POSIXct")

    }else if(vi == "dataframe" | vi == "data.frame" | vi == "df"){

      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)

    }else if(vi == "fileexists"){

      cv <- all(file.exists(input))

    }else if(vi == "direxists"){

      cv <- all(dir.exists(input))

    }else if(vi == "granges" | vi == "gr"){

      cv <- is(input, "GRanges")

    }else if(vi == "list" | vi == "simplelist"){

      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)

    }else if(vi == "bsgenome"){

      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text=input))
      }, error = function(e){
        FALSE
      })
      cv <- any(cv1, cv2)

    }else if(vi == "se" | vi == "summarizedexperiment"){

      cv <- is(input, "SummarizedExperiment")

    }else if(vi == "seurat" | vi == "seuratobject"){

      cv <- is(input, "Seurat")

    }else if(vi == "txdb"){

      cv <- is(input, "TxDb")

    }else if(vi == "orgdb"){

      cv <- is(input, "OrgDb")

    }else if(vi == "bsgenome"){

      cv <- is(input, "BSgenome")

    }else if(vi == "parallelparam"){

      cv <- is(input, "BatchtoolsParam")

    }else if(vi == "archrproj" | vi == "archrproject"){

      cv <- is(input, "ArchRProject")
      ###validObject(input) check this doesnt break anything if we
      ###add it. Useful to make sure all ArrowFiles exist! QQQ

    }else{

      stop("Validator is not currently supported by ArchR!")

    }

    if(cv){
      av <- TRUE
      break
    }

  }

  if(av){

    return(invisible(TRUE))

  }else{

    stop("Input value for '", name,"' is not a ", paste(valid, collapse="," ), ", (",name," = ",class(input),") please supply valid input!")

  }

}
