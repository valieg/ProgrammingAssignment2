##
## Here we try to build a kind of matrix and matrix inverse caching object.
## The base matrix have to be numeric, squared and inversable, in order to talk
## about an inverse(so det(A)!=0 and the inverse does exist).
##
## Parts:
##       makeCacheMatrix() - the object of type "function" that "encapsulate"
##                           the matrix, his inverse and the "classic setters and getters".
##                         - SET the base matrix(if a proper one is provided!)
##                         - WARNING - this function(more exactly his inner function 'setMatrix()')
##                                     will NOT generate automatically the inverse of a NEW matrix!
##                                     So NO synchronization.
##
##       cacheSolve()      - provide the inverse for a proper and new(!) matrix
##                         - the function is called by 'setMatrixInverse()' inner function of
##                           makeCacheMatrix()
##
##       isNumericSquareMatrix() - kind of helper, validator,...
##                               - return TRUE if the provided parameter is a numeric and
##                                 squared matrix; FALSE otherwise.
##
##       equalMatrices()         - kind of helper, validator,...
##                               - return TRUE if both provided parameters(matrices) are equals;
##                                 FALSE otherwise.
##
## Usage:
##
##       Make the object of matrix cache:
##
##       >oMatrixCache <- makeCacheMatrix(A)
##                        - if A is not a numeric squared matrix you get an 'warning message'
##
##       Ask for the inverse(for a STORED matrix!):
##
##       >cacheSolve(oMatrixCache)
##
##       Ask for the inverse of a NEW matrix:
##
##       >cacheSolve(oMatrixCache, newA)
##
##       If you want to check that the stored pair(matrix, inverse) is right, then use:
##
##       >round(oMatrixCache$getMatrixInverse() %*% oMatrixCache$getMatrix(), 3)
##                        - the response have to be an identity matrix(unit matrix)
##
##       You may also solve systems of linear equations. For a system Ax=b, the solution is x=inverseA %*% b
##

#' This function creates a special "matrix" object that can cache its inverse.
#'
#' @param x matrix
#'
#' @return list of available functions
#'
#' @examples
makeCacheMatrix <- function(x = matrix()) {
    
    # Initialize some local variables.
    oMatrix<-NULL
    oMatrixInverse<-NULL
    
    #
    # Setters & Getters
    #
    setMatrix <- function(o){
        if(TRUE == isNumericSquareMatrix(o)){
            oMatrix <<- o
        } else{
            message("The matrix is not numeric and/or squared!")
        }
    }
    
    getMatrix <- function(){
        return(oMatrix)
    }
    
    setMatrixInverse <- function(o){
        if(TRUE == isNumericSquareMatrix(o)){
            oMatrixInverse <<- o
        } else{
            message("Wrong inverse matrix!")
        }
    }
    
    getMatrixInverse <- function(){
        return(oMatrixInverse)
    }
    #
    # Set the base matrix.
    #
    setMatrix(x)
    
    return(list(setMatrix = setMatrix,
                getMatrix = getMatrix,
                setMatrixInverse = setMatrixInverse,
                getMatrixInverse = getMatrixInverse))
    
}

#' This function computes the inverse of the special "matrix" returned by 'makeCacheMatrix'.
#' If the inverse has already been calculated (and the matrix has not changed), then
#' 'cacheSolve' should retrieve the inverse from the cache.
#'
#' @param x object of type 'makeCacheMatrix'
#' @param newMatrix matrix
#' @param ...
#'
#' @return matrix
#'
#' @examples
cacheSolve <- function(x, newMatrix = NULL,...) {
    
    result = FALSE
    
    if(is.null(newMatrix)){
        #
        # Okay, we'll return the cached inverse(if does exist).
        #
        ##tmpMatrix <- x$getMatrix()
        if(is.null(x$getMatrix())){
            # No matrix has been stored, so no inverse.
            message("No matrix has been stored!")
            return(result)
        }
        #
        # So we've a proper matrix stored('setMatrix()' of 'makeCacheMatrix'
        # makes sure that we don't store garbage).
        #
        #----------------
        tmpMatrixInverse <- x$getMatrixInverse()
        if(is.null(x$getMatrixInverse())){
            # No inverse, but we have the matrix! Let's get the inverse.
            tmpMatrixInverse <- solve(x$getMatrix())
            if(TRUE == isNumericSquareMatrix(tmpMatrixInverse)){
                # Yesss, we got a right inverse. Let's store and return it!
                x$setMatrixInverse(tmpMatrixInverse)
                return(x$getMatrixInverse())
            } else {
                message("Something went wrong during matrix inversion process!")
                return(result)
            }
        }
        #
        # If we get here we've a proper matrix inverse stored('setMatrixInverse()'
        # of 'makeCacheMatrix' makes sure that we don't store garbage).
        #
        return(tmpMatrixInverse)
        #----------------
    } else {
        #
        # Okay, we've a NEW matrix; let's do the job!
        #
        if(FALSE == isNumericSquareMatrix(newMatrix)){
            message("Wrong matrix. Inversion process impossible!")
            return(result)
        }
        
        if(FALSE == equalMatrices(x$getMatrix(), newMatrix)){
            #
            # Store the NEW matrix!
            # ('setMatrix' will do the matrix check for us)
            #
            x$setMatrix(newMatrix)
        }
        # Get the inverse.
        tmpMatrixInverse <- solve(x$getMatrix())
        if(TRUE == isNumericSquareMatrix(tmpMatrixInverse)){
            # Got a right inverse.
            x$setMatrixInverse(tmpMatrixInverse)
            return(x$getMatrixInverse())
        } else {
            message("Something went wrong during the new matrix inversion process!")
            return(result)
        }
        return(tmpMatrixInverse)
    }
    
    return(result)
    
}

#' Check if the input object is a numeric squared matrix.
#'
#' @param o matrix
#'
#' @return boolean
#'
#' @examples
isNumericSquareMatrix <- function(o = matrix()){
    
    result = FALSE
    
    # Check if the input have the proper class(matrix) and the proper type(numeric).
    if(FALSE == is.matrix(o) || FALSE == is.numeric(o)){
        return(result)
    }
    # Check if the input(here is already a numeric matrix) is a squared matrix.
    if(nrow(o) == ncol(o)){
        result <- TRUE
    }
    
    return(result)
}

#' Make sure that the input are proper matrices and
#' check if they are identical(equal).
#'
#' @param o1 matrix
#' @param o1 matrix
#'
#' @return boolean
#'
#' @examples
equalMatrices <- function(o1 = matrix(), o2 = matrix()){
    
    result = FALSE
    
    # Check if both input parameters are proper matrices.
    if(FALSE == isNumericSquareMatrix(o1) || FALSE == isNumericSquareMatrix(o2)){
        return(result)
    }
    #
    # Here the "double check" is on purpose!
    # If we pass the first condition(identical) it's ideal, perfect..., but
    # if we pass the second one it's not bad too. In this case we have an "epsilon"
    # difference on some elements(.Machine$double.eps)...
    if(identical(o1, o2) || isTRUE(all.equal(o1, o2))){
        result <- TRUE
    }
    
    return(result)
}
