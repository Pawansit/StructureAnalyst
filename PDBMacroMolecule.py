from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt


def BfactorPerChain(PDBID, Chain):
    BfactorMean = PDBID.df['ATOM'][PDBID.df['ATOM']['chain_id']== Chain]['b_factor'].mean()
    BfactorMin  = PDBID.df['ATOM'][PDBID.df['ATOM']['chain_id']== Chain]['b_factor'].min()
    BfactorMax  = PDBID.df['ATOM'][PDBID.df['ATOM']['chain_id']== Chain]['b_factor'].max()
    return(BfactorMean,BfactorMin, BfactorMax)



def BfactorPerChainLinePlot(PDBID, Chain):
    PDBID.df['ATOM'][PDBID.df['ATOM']['chain_id']== Chain]['b_factor'].plot(kind="line")
    plt.title('B-Factors Along the Amino Acid Chain')
    plt.xlabel('Residue Number')
    plt.ylabel('B-factor in $A^2$')
    plt.savefig("3eiy_BfatorLine.png")



ppdb = PandasPdb().read_pdb('../../ComPath-DataDownload/PDBAPI/Compath_db/PDB_downlaod_rsync/pdb/pdb7l1z.ent')
BfactorPerChain(PDBID= ppdb, Chain="A")
BfactorPerChainLinePlot(PDBID= ppdb, Chain="A")




