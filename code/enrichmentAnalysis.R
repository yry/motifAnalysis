##enrichment analysis using a Fisher exact test
##the output consits of entries with p.val < 0.05

args = commandArgs(trailingOnly = TRUE)
if(length(args) < 3)
{
  stop("three params are required, min and max win length and the pval threshold")
  }else
  {
    win.min   = args[1]
    win.max   = args[2]
    p.val.thr = args[3]
  }

#win.min <- readline(prompt="enter min win length")
#win.max <- readline(prompt="enter max win length")
#win.min = 5
#win.max = 30

print("min window length is:")
print(win.min)
print("max window length is:")
print(win.max)

for(w in (win.min:win.max))
{
  inp <- read.csv(paste("motifs",w,".csv",sep=''),header=TRUE)
  fisher.table <- matrix(rep(0,4),ncol=2)
  fisher.pval <- rep(0,nrow(inp))
  
  for(i in 1:nrow(inp))
  {
    fisher.table[1,1] = inp[i,2]
    fisher.table[2,1] = inp[i,3]
    fisher.table[1,2] = inp[i,4]
    fisher.table[2,2] = inp[i,5]
    
    fish.test <- fisher.test(fisher.table)
    fisher.pval[i] <- fish.test$p.value
    
    #bar.test <- barnard.test(inp[i,2],inp[i,3],inp[i,4],inp[i,5])
    #bar.pval[i] <- bar.test$p.value
  }
  inpFish <- cbind(inp, fisher.pval)
  inpFish.signficant = inpFish[inpFish[,6]< p.val.thr,]
  write.csv(inpFish.signficant,file=paste("motifs",w,"_significant.csv",sep=''))
}