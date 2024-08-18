import sys
import tree_reader
import numpy as np
import numpy.ma as ma
from scipy.stats import multivariate_normal as mvnorm
#from numpy.random import multivariate_normal as np_mvnorm
#from numpy.random import default_rng
from sklearn.covariance import GraphicalLassoCV, ledoit_wolf, LedoitWolf
from sklearn.covariance import ShrunkCovariance
import matplotlib.pyplot as plt
import warnings
#warnings.filterwarnings("ignore")
import seaborn as sns

#rng = default_rng(12345)
#np_mvnorm = rng.multivariate_normal


def DEPcov_plot(cov1,cov2):
    plt.figure(figsize=(10, 6))
    plt.subplots_adjust(left=0.02, right=0.98)
    # plot the covariances
    covs = [
        ("base model", cov1),
        ("shift node", cov2),
    ]
    vmax = cov1.max()
    for i, (name, this_cov) in enumerate(covs):
        plt.subplot(2, 2, i + 1)
        plt.imshow(
            this_cov, interpolation="nearest", vmin=-vmax, vmax=vmax, cmap=plt.cm.RdBu_r
        )
        plt.xticks(())
        plt.yticks(())
        plt.title("%s covariance" % name)
    plt.show()

def cov_plot(cov1,cov2):
    plt.figure(figsize=(10, 6))
    plt.subplots_adjust(left=0.02, right=0.98)
    # plot the covariances
    covs = [
        ("base model", cov1),
        ("shift node", cov2),
    ]
    vmax = cov1.max()
    mask = np.zeros_like(cov1)
    mask[np.triu_indices_from(mask)] = True
    for i, (name, co) in enumerate(covs):
        plt.subplot(2, 2, i + 1)
        with sns.axes_style("white"):
            ax=sns.heatmap(co,mask=mask,square=True, vmin=-vmax,vmax=vmax, cmap="PiYG")
            #ax=sns.heatmap(co,vmax=1,square=True,  cmap="YlGnBu")
            ax.set_title(name+" evolutionary correlation")
            #plt.savefig(+".shared.svg",format="svg")
    plt.show()



def plot_all_models(covd):
    plt.figure(figsize=(10, 6))
    plt.subplots_adjust(left=0.02, right=0.98)
    # plot the covariances
    covs = []
    for node in covd:
        covs.append((node,covd[node]))
    #vmax = cov1.max()
    vmax = 1.0
    #vmax = 0.1
    mask = np.zeros_like(covs[0][1])
    mask[np.triu_indices_from(mask)] = True
    for i, (name, co) in enumerate(covs):
        plt.subplot(1, len(covs), i + 1)
        with sns.axes_style("white"):
            if i+1 == len(covs):
                ax=sns.heatmap(co,mask=mask,square=True, vmin=-vmax,vmax=vmax, cmap="PiYG")
            else:
                ax=sns.heatmap(co,mask=mask,square=True, vmin=-vmax,vmax=vmax, cmap="PiYG",cbar=False)

            #ax=sns.heatmap(co,vmax=1,square=True,  cmap="YlGnBu")
            ax.set_title(name+" evolutionary correlation")
            #plt.savefig(+".shared.svg",format="svg")
    plt.show()



def test_plot(emp_cov,lw_cov_,cov_,cov):
    plt.figure(figsize=(10, 6))
    plt.subplots_adjust(left=0.02, right=0.98)

    # plot the covariances
    covs = [
        ("Empirical", emp_cov),
        ("Ledoit-Wolf", lw_cov_),
        ("GraphicalLassoCV", cov_),
        ("True", cov),
    ]
    vmax = cov_.max()
    for i, (name, this_cov) in enumerate(covs):
        plt.subplot(2, 4, i + 1)
        plt.imshow(
            this_cov, interpolation="nearest", vmin=-vmax, vmax=vmax, cmap=plt.cm.RdBu_r
        )
        plt.xticks(())
        plt.yticks(())
        plt.title("%s covariance" % name)
    plt.show()


def read_traits(tree,traitflnm,h=False):
    traitfl = open(traitflnm,"r")
    traitdic = {}
    header = ""
    if h == True:
        header=traitfl.readline()
    for line in traitfl:
        spls=line.strip().split(",")
        tax = spls[0]
        traits = []
        for i in spls[1:]:
            if i != "?":
                traits.append(float(i))
            else:
                traits.append(i)
        traitdic[tax] = traits
    traitfl.close()
    for n in tree.iternodes():
        if n.istip:
            try:
                n.data["cont"]=traitdic[n.label]
            except:
                print(n.label," represented in the tree has no data in the matrix!")
                sys.exit()

def postorder(tree):#,traitdic):
    nodenum=0
    for n in tree.iternodes("postorder"):
        n.data["prunelens"] = []
        if len(n.children) != 0:
            nodeid="node"+str(nodenum)
            n.label=nodeid
            c0 = n.children[0]
            c1 = n.children[1]
            c0traits = c0.data["cont"]
            c1traits = c1.data["cont"]
            tmpvals = []
            for i in range(len(c0traits)):
                """for ch in n.children:
                    if ch.istip:
                        ch.data["prunelens"].append(ch.length)"""
                t0=c0traits[i]
                t1=c1traits[i]
                if t0 == "?" and t1 != "?":
                    n.data["prunelens"].append(n.length+c1.length)
                    tmpvals.append(t1)
                elif t1 == "?" and t0 != "?":
                    n.data["prunelens"].append(n.length+c0.length)
                    tmpvals.append(t0)
                elif t1 == "?" and t0 == "?":
                    tmpvals.append("?")
                    n.data["prunelens"].append( n.length + ((len0 * len1) / (len0 + len1)) )
                else:
                    #print(n.label,c0.label,c0.data["prunelens"])
                    len0 = c0.data["prunelens"][i]
                    len1 = c1.data["prunelens"][i]
                    curVar = len0+len1
                    tempChar = ((len0 * t1) + (len1 * t0)) / curVar
                    tmpvals.append(tempChar)
                    n.data["prunelens"].append( n.length + ((len0 * len1) / (len0 + len1)) )
            n.data["cont"] = tmpvals
            nodenum+=1
        else:
            for i in range(len(n.data["cont"])):
                n.data["prunelens"].append(n.length)
        #print(n.label,n.data["cont"])

def preorder(tree):
    for n in tree.iternodes("preorder"):
        n.data["ancstate"]=[]
        n.data["revprnln"]=[]
        n.data["origpic"]=[]
        n.data["change"]=[]
        for i in range(len(n.data["cont"])):
            n.data["revprnln"].append(n.length)
        if n == tree:
            n.data["ancstate"] = n.data["cont"]
            n.data["cont"] = []
            for i in range(len(n.data["ancstate"])):
                n.data["cont"].append("?")
            continue
        elif n.istip:
            for i,tr in enumerate(n.data["cont"]):
                if tr == "?":
                    n.data["change"].append( np.nan )
                else:
                    n.data["change"].append(tr-n.parent.data["ancstate"][i])
            continue
        c0 = n.sib
        c1 = n.parent
        c0traits = c0.data["cont"]
        if n.parent == tree:
            try:
                c0traits = c0.data["origpics"]
            except:
                pass
        c1traits = c1.data["cont"]
        revpics = []
        for i in range(len(c0traits)):
            t0=c0traits[i]
            t1=c1traits[i]
            len0 = c0.data["prunelens"][i]
            len1 = c1.data["revprnln"][i]
            if t0 == "?" and t1 == "?":
                revprnln = n.length + ((len0 * len1) / (len0 + len1))
                n.data["revprnln"].append(revprnln) ## IS THIS OKAY ??
                n.data["ancstate"].append(np.nan)
                revpics.append("?")
            elif t0 == "?" and t1 != "?":
                n.data["revprnln"][i]=n.length+c1.data["revprnln"][i]
                revpics.append(t1)
                #n.data["ancstate"].append(np.nan)
                n.data["ancstate"].append(c1.data["ancstate"][i])
                n.data["change"].append( np.nan )
            elif t1 == "?" and t0 != "?": # parent == root case covered here too
                n.data["revprnln"][i]=n.length+c0.data["prunelens"][i]
                revpics.append(t0)
                #if n.parent == tree:
                revprnln = n.data["revprnln"][i]
                forprnln = n.data["prunelens"][i] - n.length
                n.data["ancstate"].append( ((revprnln * n.data["cont"][i]) + (forprnln * revpics[i])) / ( revprnln + forprnln ) )
                if n.parent == tree:
                    #print("CHILD",len(n.data["ancstate"]),n.data["ancstate"][i])
                    #print("PARENT",len(n.parent.data["ancstate"]),n.parent.data["ancstate"][i])
                    n.data["change"].append(n.data["ancstate"][i]-n.parent.data["ancstate"][i])
                else:
                    n.data["change"].append( np.nan )
                #else:
                #    n.data["ancstate"].append("?")
            else:
                if len0<0.00000001:
                    len0 = 0.000001
                if len1<0.00000001:
                    len1 = 0.000001
                curVar = len0+len1
                tempChar = ((len0 * t1) + (len1 * t0)) / curVar
                revpics.append(tempChar)
                revprnln = n.length + ((len0 * len1) / (len0 + len1))
                forprnln = n.data["prunelens"][i] - n.length
                n.data["revprnln"][i] =  revprnln
                #print(revprnln,forprnln,n.data["cont"][i],revpics[i])
                if n.data["cont"][i] != "?":
                    n.data["ancstate"].append( ((revprnln * n.data["cont"][i]) + (forprnln * revpics[i])) / ( revprnln + forprnln ) )
                else:
                    n.data["ancstate"].append( np.nan )
                n.data["change"].append(n.data["ancstate"][i]-n.parent.data["ancstate"][i])
        n.data["origpics"] = n.data["cont"]
        n.data["cont"] = revpics
    for n in tree.iternodes():
        if n.istip == False:
            n.data["cont"] = n.data["ancstate"]


def changes_array(tree):
    trw = {}
    cols = []
    for n in tree.iternodes():
        if n.parent == None or n.data["reg"]!=tree.data["reg"]:#or n.note==False:
            continue
        cols.append(n.label)
        for i,tr in enumerate(n.data["change"]):
            try:
                trw[i].append(tr)
            except:
                trw[i] = []
                trw[i].append(tr)
    #print(cols)
    changes = np.array(list(trw.values()))
    return changes

def node_trait_corrs(subtree,type="cov"):
    if type != "cov" and type != "corr":
        print("need to specify either `cov` or `corr` for covariance/correlation matrix between trait changes")
        sys.exit()
    chg = changes_array(subtree)
    chg = [ma.masked_invalid(i) for i in chg]
    if type == "cov":
        corrs = ma.cov(chg)
    elif type == "corr":
        corrs = ma.corrcoef(chg)
    mean = [ma.mean(i) for i in chg]
    return mean,corrs

def mean_impute(means,x):
    for i,v in enumerate(x):
        if np.isnan(v):
            if np.isnan(means[i]):
                print("ERROR: check mean_impute()", means)
                sys.exit()
            x[i] = means[i]
    return x

#def ledoit_wolf(empcov,shrink):
#    n = float(len(empcov))
#    mu = np.trace(empcov) / n
#    lw = (1.-shrink) * emp

def lassoCV(empcov):
    model = GraphicalLassoCV()
    model.fit(empcov)
    cov_ = model.covariance_
    return cov_

def node_cov_ll(subtree):
    empmeans,cov = node_trait_corrs(subtree,type="corr")
    #print(ledoit_wolf(cov))
    #lw_cov = ledoit_wolf(cov)[0]
    shrink = True
    if shrink == True:
        #print(ledoit_wolf(cov)[1])
        #lw=ledoit_wolf(cov)
        #lw_cov = lw[0]
        #lw_shrink = lw[1]
        intens=0.10
        lw_cov = ShrunkCovariance(shrinkage=intens).fit(cov)
        lw_cov =lw_cov.covariance_
        diag = np.sqrt(np.diag(np.diag(lw_cov)))
        gaid = np.linalg.inv(diag)
        lw_cov = gaid @ lw_cov @ gaid
        cov = lw_cov

    #print(cov)
    #sys.exit()
    mean = np.zeros(len(cov))
    try:
        d = mvnorm(mean,cov,allow_singular=False)
    except:
        print(subtree.label,"WARNING: you have too many traits and too few tips to get a reliable estimate of the covariance matrix. proceeding anyway, but be skeptical of your results")
        d = mvnorm(mean,cov,allow_singular=True)
    ll = 0.0
    #tc=0
    for n in subtree.iternodes():
        if n.parent == None or n.data["reg"] != subtree.data["reg"]:#or n.note == False:
            continue
        #tc+=1
        nchanges = np.array(n.data["change"]) # ma.masked_invalid(n.data["change"])
        if np.isnan(np.sum(nchanges)):
            nchanges = mean_impute(empmeans,nchanges)
        curll = d.logpdf(nchanges)
        ll += curll
    #print(tc)
    subtree.data["cov"] = cov
    return ll

def anc_changes(tree):
    postorder(tree)
    preorder(tree)
    outfl=open("NODELABELS.tre","w")
    outfl.write(tree.get_newick_repr(True))
    outfl.close() 

def recon_subtree_models(tree,min):
    for n in tree.iternodes():
        n.data["reg"] = 0
    basell = node_cov_ll(tree)
    ntraits = len(tree.data["cov"])
    ntip = len(tree.leaves())
    k = (ntraits*(ntraits-1))/2.
    nnode = float(len([n for n in tree.iternodes()])-1)
    #baseBIC = (-2. * basell) + (np.log(nnode)) * k
    crit = "aic"
    if crit == "aic":
        baseAIC = ( 2. * k ) - ( 2. * basell )
    elif crit == "bic":
        baseAIC = (-2 * basell) + np.log(nnode*ntraits) * k  #( 2. * k ) - ( 2. * basell )
    curAIC = baseAIC
    #baseAICc = baseAIC + ( ( ( ( 2. *(k**2) ) + (2. * k) ) ) / ( nnode - k - 1) )
    nodescores={}
    for n in tree.iternodes("postorder"):
        cladesize = len(n.leaves())
        if cladesize < min or ntip-cladesize < min or n == tree:
            continue
        curll = node_cov_ll(n)
        for sn in n.iternodes():
            #if sn != n:
                #sn.note = False
            sn.data["reg"] = 1

        parll = node_cov_ll(tree)
        #sys.exit()
        combll = curll + parll
        #curBIC = (-2. * combll) + (np.log(nnode)) * (k*2)
        if crit == "aic":
            curAIC = (2. * (k * 2.)) - (2. * combll)
        elif crit == "bic":
            curAIC = (-2 * basell) + np.log(nnode*ntraits)
        #curAICc = curAIC + (( ( ( 2*((k*2.)**2) ) + (2 * (k * 2)) ) ) / ( nnode - (k*2.) - 1))
        for sn in n.iternodes():
            #if sn != n:
            sn.data["reg"] = 0
                #sn.note = True

        #print(n.label)
        #print("BIC",baseBIC,curBIC)
        #print("AIC",baseAIC,curAIC)
        #print("AICc",baseAICc,curAICc)
        #cov_plot(tree.data["cov"],n.data["cov"])
        if curAIC < baseAIC:
            nodescores[n.label] = curAIC
    nodescores = sorted(nodescores,key=nodescores.get)
    return nodescores

def assemble_model(tree,nodescores):
    for n in tree.iternodes():
        n.data["reg"] = 0
    ntraits = len(tree.data["cov"])
    k = (ntraits*(ntraits-1))/2.
    basell = node_cov_ll(tree)
    lastAIC = ( 2. * k ) - ( 2. * basell )
    firstAIC=lastAIC
    nreg=1
    covd={}
    for nlab in nodescores:
        for n in tree.iternodes():
            if n.istip:
                continue
            if n.label == nlab:
                regcount = 0
                for sn in n.iternodes():
                    #if sn != n:
                    if sn.data["reg"] != 0:
                        continue
                    regcount+=1
                    sn.data["oldreg"] = sn.data["reg"]
                    sn.data["reg"] = nreg
                if regcount < 10:
                    continue # this is to deal with nested shifts. the ordering guarantees that we will pick the best shift location and we won't consider nested shifts if breaking up results in one subregime too small to estimate a covariance matrix
                curll = node_cov_ll(n)
                parll = node_cov_ll(tree)
                combll = curll + parll
                curAIC = (2. * (k * nreg)) - (2. * combll)                
                if curAIC < lastAIC:
                    lastAIC = curAIC
                    nreg+=1
                else:
                    for sn in n.iternodes():
                        sn.data["reg"] = sn.data["oldreg"]
                lab = "reg"+str(n.data["reg"])
                mat = n.data["cov"]
                covd[lab]=mat
    print("FINAL MODEL AIC....\n\n")
    print("SINGLE_REG:",firstAIC)
    print("MULTI_REG: ",curAIC)
    covd["reg"+str(tree.data["reg"])]=tree.data["cov"]
    for i in covd:
        print(i,np.mean([abs(j) for j in covd[i]]))
    plot_all_models(covd)

    for n in tree.iternodes():
        #print(n.data["reg"])
        if n.istip==False:
            n.label = str(n.data["reg"])
    outfl=open("regimes.tre","w")
    outfl.write(tree.get_newick_repr(True)+";")
    outfl.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: "+sys.argv[0]+"<tree> <measurements_table> <min_subtree_size>")
        sys.exit()

    treefl=sys.argv[1]
    nwk=open(treefl,"r").readline()
    tree=tree_reader.read_tree_string(nwk)

    read_traits(tree,sys.argv[2],True)
    anc_changes(tree)
    """for n in tree.iternodes("postorder"):
        print(n.label,n.data["change"])
        #if n != tree:
        #    n.length = sum([abs(i) for i in n.data["change"]])/float(len(n.data["change"]))
        #if n.istip == False:
        #    print(n.label,n.data["ancstate"])
        #else:
        #    print(n.label,n.data["cont"])"""
    subtr_cutoff = int(sys.argv[3])
    nodescores = recon_subtree_models(tree,subtr_cutoff)
    assemble_model(tree,nodescores)
