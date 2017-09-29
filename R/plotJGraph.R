plotJGraph <- function(A,fn="Net",th = 1e-06, mylayout = NULL){
M <- dim(A)[1]
p <- dim(A)[2]

#### plot case ####
for(iter in 1:M)
{
  plotGraph(A[iter,,],th = th, fn=paste(fn,iter,sep=""),mylayout = mylayout)
}
}