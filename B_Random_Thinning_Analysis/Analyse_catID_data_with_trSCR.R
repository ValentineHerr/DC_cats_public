### Analyse identified cats with random thining SCR###
# clear environment ####
rm(list = ls())

# increase memory limit to 7Tb####
memory.limit(7000000)

# load libraries ####
library(nimble)
library(coda)
library(snow)
library(doSNOW)
library(foreach)

# load data ####
all_blocks <- sapply(paste0("B_Random_Thinning_Analysis/data_block", 1:4, ".Rdata"), function(x) mget(load(x)), simplify = F)



for(b in seq_along(all_blocks)) {

  list2env(all_blocks[[b]], envir = .GlobalEnv)

  cores=chains=5
  seeds=round(runif(chains,1,10000000))


  #make cluster of "cores" cores to run "chains" chains (1 per core)
  cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
  registerDoSNOW(cl.tmp)

  clusterExport(cl=cl.tmp, list("data", "I.ID", "lam0", "sigma", "gamma", "n.cat", "n.levels", "M1_factor", "M2_factor"),
                envir=environment())


  start.time<-Sys.time()

  #run in parallel for loop
  out=foreach(c=1:chains) %dopar% {
    library(nimble)
    source("RandomThinning_Scripts/NimbleModelRT Cats NB Jcov.R")
    source("RandomThinning_Scripts/NimbleFunctionsRT Cats NB Jcov.R")
    source("RandomThinning_Scripts/init.catRT.df.R")
    source('RandomThinning_Scripts/Restart MCMC Functions.R')

    set.seed(seeds[c]) #set seed for reproducibility
    M1=ceiling(I.ID*16*M1_factor) #Augmentation level for marked.
    M2=ceiling(M1/2*M2_factor) #Unmarked individual data augmentation level
    M=M1+M2
    n.cat=4

    inits=list(lam0=lam0,sigma=sigma,gamma=gamma)

    #This function structures the simulated data to fit the model in Nimble (some more restructing below)
    #Also checks some inits
    nimbuild=init.catRT(data,inits,M=M)

    gammaMat=nimbuild$gammaMat
    G.true.init=nimbuild$G.true
    G.true.data=G.true.init*NA

    #Trap covariates
    covs.trap=model.matrix(~placement_type+distance_to_obstruction_scaled+
                             distance_to_obstruction_scaled_squared,data=data$J.cov) # no bait
    covs.trap=covs.trap[,-c(1)]#remove intercept

    #inits for nimble
    #full inits. Nimble can initialize psi1 and psi2, but if sigma and lam0 initialized too far away
    #from truth, it can stop adapting before convergence and mix very poorly.
    Niminits <-
      list(
        z = rep(1, M),
        s = nimbuild$s,
        G.true = G.true.init,
        ID = nimbuild$ID,
        capcounts = rowSums(nimbuild$y.true),
        y.true = nimbuild$y.true,
        gammaMat = gammaMat,
        G.latent = nimbuild$G.latent,
        lam0.beta0 = -1.5,
        lam0.beta.sex = 0.25,
        sigma.beta0 = 4.2,
        sigma.beta.sex = 0.2,
        thin.beta0 = -1.5,
        thin.beta.coat = c(0.25, 3.5, 2, 3.5),
        thin.beta.bi = 3.5,
        thin.beta.hair = 1.5,
        lam0.beta.trap =  c("placement_typeprivate" = 0,
                            "placement_typepublic" = -3,
                            "distance_to_obstruction_scaled" = -0.35,
                            "distance_to_obstruction_scaled_squared" = 0.125)[colnames(covs.trap)], # originally c(0, -3.5, -2, -0.35, 0.125), # private, public forest, public open, dist obstr, dist obstr squared
        r = 0.005,
        psi = 0.9
      )


    #constants for Nimble
    J=nrow(data$X)
    constants<-list(M=M,J=J,K=data$K,K1D=data$K1D,n.samples=nimbuild$n.samples,n.cat=n.cat,n.levels=n.levels,
                    xlim=data$xlim,ylim=data$ylim,covs.trap=covs.trap, n.trap.cov = ncol(covs.trap))

    # Supply data to Nimble. Note, y.true and y.true.event are completely latent.
    z.data=c(rep(1,data$n.ID),rep(NA,M-data$n.ID))

    Nimdata<-list(y.true=matrix(NA,nrow=M,ncol=J),y.ID=nimbuild$y.ID,
                  G.true=G.true.data,ID=rep(NA,nimbuild$n.samples),z=z.data,X=as.matrix(data$X),capcounts=rep(NA,M))

    # set parameters to monitor
    parameters<-c('psi','lam0.beta0','lam0.beta.sex','sigma.beta0','sigma.beta.sex','gammaMat','N','n',
                  'thin.beta0','thin.beta.coat','thin.beta.bi','thin.beta.hair','lam0.beta.trap','r')

    nt=1 #thinning rate

    #can also monitor a different set of parameters with a different thinning rate
    parameters2 <- c("s","z")
    nt2=nt*10 #thinning 10x normal here to reduce posterior file size

    # Build the model, configure the mcmc, and compile
    start.time<-Sys.time()
    Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                          inits=Niminits)

    conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt, useConjugacy = TRUE,
                          nodes=c("z","gammaMat","s","psi",'lam0.beta0','lam0.beta.sex',
                                  'sigma.beta0','sigma.beta.sex','thin.beta0','thin.beta.coat',
                                  'thin.beta.bi','thin.beta.hair','lam0.beta.trap','r'),
                          monitors2=parameters2,thin2=nt2)

    ###Two *required* sampler replacements

    ##Here, we remove the default samplers for y.true and y.event, which are not correct
    #and replace it with the custom "IDSampler"
    conf$removeSampler("y.true")
    conf$addSampler(target = paste0("y.true[1:",M,",1:",J,"]"),
                    type = 'IDSampler',control = list(M=M,J=J,K1D=data$K1D,n.cat=n.cat,
                                                      n.samples=nimbuild$n.samples,n.ID=data$n.ID,
                                                      this.j=nimbuild$this.j,G.noID=nimbuild$G.noID,
                                                      G.latent.ID=nimbuild$G.latent.ID),
                    silent = TRUE)

    #replace default G.true sampler, which is not correct, with custom sampler for G.true, "GSampler"
    conf$removeSampler("G.true")
    for(i in 1:M){
      for(m in 1:n.cat){ #don't need to update first cat bc it is mark status
        conf$addSampler(target = paste("G.true[",i,",",m,"]", sep=""),
                        type = 'GSampler',
                        control = list(i = i,m=m,M=M,n.cat=n.cat,n.samples=nimbuild$n.samples,
                                       n.levels=n.levels,G.noID=nimbuild$G.noID), silent = TRUE)
      }
    }

    ###Two *optional* sampler replacements:

    #replace default activity center sampler that updates x and y locations separately with a joint update
    #should be a little more efficient. Slice seems better than block random walk
    conf$removeSampler(paste("s[1:",M,", 1:2]", sep=""))
    for(i in 1:M){
      conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                      type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
    }

    conf$removeSampler(c("lam0.beta0","sigma.beta0","lam0.beta.sex","sigma.beta.sex"))
    conf$addSampler(target = c("lam0.beta0","sigma.beta0","lam0.beta.sex","sigma.beta.sex"),
                    type = 'AF_slice',
                    control = list(adaptive=TRUE),
                    silent = TRUE)

    conf$removeSampler(c("lam0.beta.trap"))
    conf$addSampler(target = c("lam0.beta.trap"),
                    type = 'AF_slice',
                    control = list(adaptive=TRUE),
                    silent = TRUE)

    # Build and compile
    Rmcmc <- buildMCMC(conf)
    # runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


    # Run the model.
    start.time2<-Sys.time()
    Cmcmc$run(30000, reset = TRUE)
    end.time<-Sys.time()
    end.time-start.time  # total time for compilation, replacing samplers, and fitting
    end.time-start.time2 # post-compilation run time


    # get the first set of monitored parameters
    mvSamples = as.matrix(Cmcmc$mvSamples)

    # get the second set of monitored parameters (s and z)
    mvSamples2 = as.matrix(Cmcmc$mvSamples2)

    #save MCMC and sampler states
    stateList <- list(modelState = getModelState(Cmodel),
                      mcmcState = getMCMCstate(conf, Cmcmc))

    # get some summary info

    area <- diff(range(data$X[,1])) * diff(range(data$X[,2]))
    summary_results <- data.frame(blockID = b,
                                  total_time = paste(round(end.time-start.time, 2), units(end.time-start.time)) ,
                                  run_tim = paste(round(end.time-start.time2), units(end.time-start.time2)) ,
                                  area = paste(round(area/1000000), "sqkm"),
                                  n_traps = J,
                                  M1,
                                  M2,
                                  M,
                                  n_ID = I.ID,
                                  n_noID = dim(data$y.noID)[1],
                                  N_total = mean(mvSamples[-c(1:1000),"N"]),
                                  sigma_F  = mean(mvSamples[-c(1:1000),"sigma.beta0"]),
                                  sigma_M  = mean(mvSamples[-c(1:1000),"sigma.beta0"]) + mean(mvSamples[-c(1:1000),"sigma.beta.sex"]),
                                  n_spatial_recaptures = sum(apply(apply(data$y.ID, c(1,2), sum), 1, function(x) sum(x>0)) > 1),
                                  n_spatial_recaptures_M = sum(apply(apply(data$y.ID, c(1,2), sum)[data$G.ID[,3] %in% 2, , drop = FALSE], 1, function(x) sum(x>0)) > 1),
                                  n_spatial_recaptures_F = sum(apply(apply(data$y.ID, c(1,2), sum)[data$G.ID[,3] %in% 1, , drop = FALSE], 1, function(x) sum(x>0)) > 1),
                                  n_spatial_recaptures_UNK = sum(apply(apply(data$y.ID, c(1,2), sum)[data$G.ID[,3] %in% 0, , drop = FALSE], 1, function(x) sum(x>0)) > 1),

                                  placement_type  = deparse(c(colSums(covs.trap[,grepl("place", colnames(covs.trap))]), intercept = nrow(covs.trap) - sum(covs.trap[,grepl("place", colnames(covs.trap))])))[1],
                                  # rotation_theta,
                                  D = round(mean(mvSamples[-c(1:1000),"N"])/(area/1000000),2)
    )

    return(
      list(
        mvSamples = mvSamples,
        mvSamples2 = mvSamples2,
        time = end.time - start.time2,
        run.seed = seeds[c],
        stateList = stateList,
        summary_results = summary_results,
        gammaMat = gammaMat,
        covs.trap = covs.trap
      )
    )

  }

  end.time<-Sys.time()

  # save output
  save(list = c("out", "blockID", "data", "category_key"), file= paste0("B_Random_Thinning_Analysis/trSCR_block", b, ".Rdata"))

  stopCluster(cl.tmp)

  # plot traces in pdf ####
  n.chain=length(out)

  ## gamma key
  gammaMat_translation <- rbind(("irrelevant" = rep("irrelevant", ncol(out[[1]]$gammaMat))),
                                t(sapply(category_key, function(x) paste(c(as.character(x$value), rep("irrelevant", ncol(out[[1]]$gammaMat)-length(x$value)))))))

  gammaMat_translation

  ## lambda key
  lam0.beta.trap <- colnames(out[[1]]$covs.trap)

  #extract chains ####


  stateList=vector("list",n.chain)


  mvSamplesList=mvSamples2List=vector("list",n.chain)


  file_to_append_to <- "B_Random_Thinning_Analysis/all_summary_results_trSCR_NB.csv"

  for(i in 1:n.chain){

    mvSamplesList[[i]]=mcmc(out[[i]]$mvSamples[-1,]) #remove first iteration bc screws up plotting
    mvSamples2List[[i]]=mcmc(out[[i]]$mvSamples2)
    stateList[[i]]=out[[i]]$stateList


    colnames(mvSamplesList[[i]]) <- gsub("gammaMat", "gammaMat_translation",   colnames(mvSamplesList[[i]]))

    colnames(mvSamplesList[[i]])[grepl("gammaMat_translation",   colnames(mvSamplesList[[i]]))] <- sapply(  colnames(mvSamplesList[[i]])[grepl("gammaMat_translation",   colnames(mvSamplesList[[i]]))], function(y) eval(parse(text = y)))

    colnames(mvSamplesList[[i]]) <- gsub("^F$", "Female ",   colnames(mvSamplesList[[i]]))
    colnames(mvSamplesList[[i]]) <- gsub("^M$", "Male ",   colnames(mvSamplesList[[i]]))
    colnames(mvSamplesList[[i]]) <- gsub("^TRUE$", "Bicolor TRUE ",   colnames(mvSamplesList[[i]]))
    colnames(mvSamplesList[[i]]) <- gsub("^FALSE$", "Bicolor FALSE ",   colnames(mvSamplesList[[i]]))


    colnames(mvSamplesList[[i]])[grepl("lam0.beta.trap",   colnames(mvSamplesList[[i]]))] <- sapply(  colnames(mvSamplesList[[i]])[grepl("lam0.beta.trap",   colnames(mvSamplesList[[i]]))], function(x) eval(parse(text = x)))


    mvSamplesList[[i]] <- as.mcmc(mvSamplesList[[i]][complete.cases(mvSamplesList[[i]]),])

    # append summary results together
    if(i ==1)   write.csv(cbind(out[[i]]$summary_results, chain = i), file =file_to_append_to, row.names = F) else  write.csv(unique(rbind(read.csv(file_to_append_to), cbind(out[[i]]$summary_results, chain = i))), file =file_to_append_to, row.names = F)


  }


  # plot trace ####

  pdf(paste0("B_Random_Thinning_Analysis/trSCR_block", b, ".pdf"))
  gelman.plot(mcmc.list(mvSamplesList)[,1], main = "Gelman convergence diagnostic - N")
  abline(h = 1.1, lty = 2)
  plot(mcmc.list(mvSamplesList))
  dev.off()

}


