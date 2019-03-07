pcs = as.matrix(read.table('PC_ests.txt'))

pcs = cbind(pcs,1-pchisq(pcs[,3]^2,1))
dimnames(pcs) [[2]] = c('est','SE','t','p')

c = 100

BF = pcs[1:c,]
dimnames(BF)[[1]] = 1:c

WF = pcs[(c+1):(2*c),]
dimnames(WF)[[1]] = 1:c



BF[BF[,4]<0.05,]

WF[WF[,4]<0.05,]

plot(BF[BF[,4]<0.1,1],WF[BF[,4]<0.1,1])

1-pchisq(sum(WF[,3]^2),40)