getProp <-
function(Y, G, p0=NULL){

    # Y: L x 2 allele specific count matrix for L biallelic SNPs
    # G: L x N genotype dosage matrix for N individuals mixed in a pool
    # p0: N-dimentional vector for the initial mixture propotions. Note that sum(p0)=1.

    G=G[apply(Y,1,sum)>0,]
    Y=Y[apply(Y,1,sum)>0,]


    N=ncol(G) # number of cell lines to be mixed
    y=apply(Y,1,sum) # total coverage depth at each SNP
    y1=Y[,2] # alt counts
    y0=Y[,1] # ref counts
    if(is.null(p0)){
        p=rep(1/N,N) # mixture proportion: default = 1/N
        b=rep(0,N-1) # softmax parametrisation
    }else{
        p=p0
        b=log(p[1:(N-1)]/p[N])
    }
    pj = c(G%*%p/2) # alt allele frequency at each SNP after mixing cell lines
    qj = 1-pj # ref allele frequency at each SNP after mixing cell lines
    lkhd0=sum(y1*log(pj) + y0*log(qj)) - sum(b^2)/2/10000 # likelihood
    d=rep(0,N-1) # optimum direction of each Newton-Raphson step

    lkhd_all=lkhd0 # record of all likelihood values
    col=1 # record of convergence status

    convstat=0
    for(itr in 1:100){
        g = (apply((G-2*pj)*(qj*y1-pj*y0)/(2*pj*qj),2,sum)*p)[1:(N-1)] # gradient
        H = - diag(g) - ((sum(y) - t((2-G)*y0/qj^2/4)%*%(2-G) - t(G*y1/pj^2/4)%*%G) * outer(p,p,"*"))[1:(N-1),1:(N-1)] + diag(N-1)/10000 # hessian
        d = solve(H,g) # optimal direction
        r=1 # step size
        flag=0 # likelihood is improved if flag=1
        for(jtr in 1:20){ # line search
            b1 = b+d/r # updated softmax parameters
            p1=c(exp(b1),1)/sum(c(exp(b1),1)) # updated mixture proportions
            pj = c(G%*%p1/2) # updated alt allele frequency
            qj = 1-pj
            lkhd1 = sum((y1*log(pj) + y0*log(qj))[pj>0&pj<1])-sum(b1^2)/2/10000;
            lkhd_all = c(lkhd_all,lkhd1); # updated likelihood
            #plot(lkhd_all)
            if(!is.na(lkhd1) && lkhd1+0.01>lkhd0){ # if likelihood is imporoved
                #print("lkhd update")
                col=c(col,1)
                b=b1
                lkhd=lkhd1
                if(lkhd1+1e-7>lkhd0 & jtr==1)flag=1
                break
            }else{ # likelihood is not improved
                #print("lkhd decreased")
                col=c(col,2)
                lkhd=lkhd1
                r=r*(-2)
            }
        }
        p  = c(exp(b),1)/sum(c(exp(b),1)) # set new mixture proportions
        pj = c(G%*%p/2) # set new ref allele freq
        qj = 1-pj
        if(itr>7 && abs(lkhd0-lkhd)<1e-7 && flag==1){ # conditions of convergence
            print(c(lkhd0,lkhd))
            print("converged")
            convstat=1
            break
        }else{ # not converged
            print(c(lkhd0,lkhd))
            lkhd0=lkhd
        }
    }

    #print(solve(H))
    B=cbind(b,b-sqrt(diag(solve(H)))*2,b+sqrt(diag(solve(H)))*2) # softmax parameters and their confidence intervals
    P=t(t(exp(rbind(B,0)))/apply(exp(rbind(B,0)),2,sum)) # mixture proportions of and their conf intervals


    list(lkhd=lkhd,P=P,convergence=c(convstat==1))

}
GetProp <-
function(Y, G){ getProp(Y, G, rep(1,ncol(G))/ncol(G)) }
