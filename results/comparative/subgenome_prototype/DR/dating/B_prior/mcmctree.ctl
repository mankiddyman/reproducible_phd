          seed = 104729
       seqfile = ../../tree/B_full.phy
      treefile = tree.trees
       outfile = out.txt
      mcmcfile = mcmc.txt

         ndata = 1
       seqtype = 0
       usedata = 0
         clock = 2    * independent rates -- lineage rates vary ~3x
       RootAge = 'B(0.672,1.003)'
finetune = 1: .1 .1 .1 .1 .1 .1

         model = 7    * GTR
         alpha = 0.5
         ncatG = 5

     cleandata = 0
       BDparas = 1 1 0.1
   kappa_gamma = 6 2
   alpha_gamma = 1 1
   rgene_gamma = 2 20 1
  sigma2_gamma = 1 10 1

         print = 1
        burnin = 500000
      sampfreq = 20
       nsample = 200000
