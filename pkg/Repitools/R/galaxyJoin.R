
.writeFileIfNeeded <- function(tabOrFile,label) {
  if( is.matrix(tabOrFile) | is.data.frame(tabOrFile) ) {
    fn<-paste(label,"data",sep=".")
    cat("Writing file:",fn,"\n")
    write.table(tabOrFile,fn,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
    return(list(file=fn,written=TRUE))
  } else if (file.exists(tabOrFile)) {
    return(list(file=tabOrFile,written=FALSE))
  } else {
    stop("Either 'tabOrFile' is not a matrix/data.frame or the file specified does not exist")
  }
}


galaxyJoin <- function(tabOrFile1=NULL,tabOrFile2=NULL,cols1=1:3,cols2=1:3,mincols=1,fill=c("none","right","left","both"),delete=TRUE,
                       script="gops_join.py") {

  options(digits.secs=6)
  
  label1<-gsub("[-: .]","",Sys.time())
  fn1<-.writeFileIfNeeded(tabOrFile1,label1)
  label2<-gsub("[-: .]","",Sys.time())
  fn2<-.writeFileIfNeeded(tabOrFile2,label2)
  
  fill<-match.arg(fill)
  
  if( is.null(options()$galaxyPath) )
    stop("Need to set options(galaxyPath=\"path/to/galaxy/dist\")")
  
  systemCall<-paste( options()$galaxyPath, script, sep="/" )
  systemCall<-paste(systemCall,"--cols1",paste(paste(cols1,collapse=","),",",sep=""))
  systemCall<-paste(systemCall,"--cols2",paste(paste(cols2,collapse=","),",",sep=""))
  systemCall<-paste(systemCall,"--mincols",mincols)
  systemCall<-paste(systemCall,"--fill",fill)
  systemCall<-paste(systemCall,fn1$file)
  systemCall<-paste(systemCall,fn2$file)
  systemCall<-paste(systemCall,label1)
  systemCall<-paste("python",systemCall)
 
  # set environment variable
  #do.call(Sys.setenv,args=list(PYTHONPATH=paste(options()$galaxyPath,"lib",sep="/")))

  cat("Calling:",systemCall,"\n")
  system(systemCall) 
  
  if(!file.exists(label1))
    stop("python call did not produce output")
	
  join<-read.table(label1,sep="\t",comment.char="",header=FALSE,stringsAsFactors=FALSE)
  
  if (delete) {
    if(fn1$written)
      unlink(fn1$file)
    if(fn2$written)
      unlink(fn2$file)
	unlink(label1)
  }
  
  join

}

