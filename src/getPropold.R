getProp <-
function(Y, G, p0=NULL){

# Y: L x 2 allele specific count matrix for L biallelic SNPs
# G: L x N genotype dosage matrix for N individuals mixed in a pool
# p0: N-dimentional vector for the initial mixture propotions. Note that sum(p0)=1.

G=G[apply(Y,1,sum)>0,]
Y=Y[apply(Y,1,sum)>0,]


N=ncol(G)
y=apply(Y,1,sum)
y1=Y[,2]
y0=Y[,1]
if(is.null(p0)){
	p=rep(1/N,N)
	b=rep(0,N-1)
}else{
	p=p0
	b=log(p[1:(N-1)]/p[N])
}
pj = c(G%*%p/2)
qj = 1-pj
lkhd0=sum(y1*log(pj) + y0*log(qj))
d=rep(0,N-1)

hoge=lkhd0
col=1

convstat=0
for(itr in 1:1000){
	g = (apply((G-2*pj)*(qj*y1-pj*y0)/(2*pj*qj),2,sum)*p)[1:(N-1)]
	H = - diag(g) - ((sum(y) - t((2-G)*y0/qj^2/4)%*%(2-G) - t(G*y1/pj^2/4)%*%G) * outer(p,p,"*"))[1:(N-1),1:(N-1)]
	d = solve(H,g)
	r=1
	#if(itr<5){r=10}else{r = 1}
	flag=0
	for(jtr in 1:20){
		b1 = b+d/r
		p1=c(exp(b1),1)/sum(c(exp(b1),1))
		pj = c(G%*%p1/2)
		qj = 1-pj
		lkhd1 = sum((y1*log(pj) + y0*log(qj))[pj>0&pj<1]); hoge=c(hoge,lkhd1); 
		if(!is.na(lkhd1) && lkhd1+0.01>lkhd0){
			#print("lkhd update")
			col=c(col,1)
			b=b1
			lkhd=lkhd1
			if(lkhd1>lkhd0)flag=1
			break
		}else{
			#print("lkhd decreased")
			col=c(col,2)
			lkhd=lkhd1
			r=r*(-2)
		}
	}
	p  = c(exp(b),1)/sum(c(exp(b),1))
	pj = c(G%*%p/2)
	qj = 1-pj
	if(itr>7 && abs(lkhd0-lkhd)<1e-7 && flag==1){
		print(c(lkhd0,lkhd))
		print("converged")
		convstat=1
		break
	}else{
		print(c(lkhd0,lkhd))
		lkhd0=lkhd
	}
}
#plot(hoge,col=col)

print(solve(H))
B=cbind(b,b-sqrt(diag(solve(H)))*2,b+sqrt(diag(solve(H)))*2)
P=t(t(exp(rbind(B,0)))/apply(exp(rbind(B,0)),2,sum))
#plot(1:N,P[,1])
#segments(1:N,P[,2],1:N,P[,3])
#abline(0,1)

list(lkhd=lkhd,P=P,convergence=c(convstat==1))

}
GetProp <-
function(Y, G){ getProp(Y, G, rep(1,ncol(G))/ncol(G)) }
