import sys
from Bio.Phylo.PAML import codeml
from Bio.Phylo.PAML.chi2 import cdf_chi2

inphy =  sys.argv[1] + ".alltaxa.paml.phy"
dir = sys.argv[1] + "_mafft"
intree = sys.argv[2]
treename = intree.split(".")[1]
altdir = treename + "_alt"
nulldir = "./" + treename + "_null"
outlikelihood = treename + "_paml_delta_likelihood.txt"
resout =  sys.argv[1] +"_" + treename + ".paml.out"
resoutnull = sys.argv[1] + "_" + treename + ".paml.null.out"
# branches of interest in newick are specified by "#1"
cml = codeml.Codeml(alignment = inphy, tree =  intree,
                    out_file = resout, working_dir = altdir)

# Set options for branch site model analysis:w
#The branch-site model M8
#Model A: model = M8, NSsites = 8, fix_omega = 0 
cml.set_options(cleandata=0)
cml.set_options(NSsites=[2])
cml.set_options(fix_omega=0)
cml.set_options(clock=0)
cml.set_options(runmode=0)
cml.set_options(icode=0)
cml.set_options(seqtype=1)
cml.set_options(omega=1)
cml.set_options(Mgene=0)
cml.set_options(model=2)
cml.set_options(CodonFreq=2)


cmlnull = codeml.Codeml(alignment = inphy, tree =  intree,
                    out_file = resoutnull, working_dir = nulldir)
#Null Model: model = M8, NSsites = 8, fix_omega = 1
cmlnull.set_options(cleandata=0)
cmlnull.set_options(NSsites=[2])
cmlnull.set_options(fix_omega=1)
cmlnull.set_options(clock=0)
cmlnull.set_options(runmode=0)
cmlnull.set_options(icode=0)
cmlnull.set_options(seqtype=1)
cmlnull.set_options(omega=1)
cmlnull.set_options(Mgene=0)
cmlnull.set_options(model=2)
cml.set_options(CodonFreq=2)

# Run codeml
alt_resultdict = cml.run(verbose=False)
null_resultdict = cmlnull.run(verbose=False)
# Get model likelihoods
alt_ln = alt_resultdict.get('NSsites').get(2).get('lnL')
null_ln = null_resultdict.get('NSsites').get(2).get('lnL')
# Get model param number
alt_df  = len(alt_resultdict.get('NSsites').get(2).get('parameters').get('parameter list').split(" "))
null_df  = len(null_resultdict.get('NSsites').get(2).get('parameters').get('parameter list').split(" "))
# Calculate degrees freedom
# in our case this is always 1
df = alt_df - alt_df
# Calculate delta likelihood
delta_lrt = 2 * (alt_ln - null_ln)
# Use chi square to get pvalue
# this function has issues converging when df=1 and takes very long
#pvalue = cdf_chi2(df, delta_lrt)
# chi2 cutoff for 0.05 significance
#pval_cutoff = 3.8415
try:
    print("Null_dict",null_resultdict)
    print("Alt_dict",alt_resultdict)
except:
    pass
with open(outlikelihood,'w') as out:
    out.write("Delta likelihood = " + str(delta_lrt))
print("Delta likelihood = " + str(delta_lrt))
