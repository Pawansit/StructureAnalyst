from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter



def BfactorPerChain(PDBID, Chain):
    BfactorMean = PDBID.df['ATOM'][PDBID.df['ATOM']['chain_id']== Chain]['b_factor'].mean()
    BfactorMin  = PDBID.df['ATOM'][PDBID.df['ATOM']['chain_id']== Chain]['b_factor'].min()
    BfactorMax  = PDBID.df['ATOM'][PDBID.df['ATOM']['chain_id']== Chain]['b_factor'].max()
    return(BfactorMean,BfactorMin, BfactorMax)



def BfactorPerChainLinePlot(PDBID, Chain):
    PDB_Chain = PDBID.df['ATOM'][(PDBID.df['ATOM']['chain_id']== Chain) & (PDBID.df['ATOM']['atom_name']== 'CA')].reset_index()
    PDB_Chain['b_factor'].plot(kind="line")
    plt.title('B-Factors Along the Amino Acid '+ Chain + ' Chain and CA atoms')
    plt.xlabel('Residue Number')
    plt.ylabel('B-factor in $A^2$')
    plt.savefig("3eiy_BfatorLine.png")



if __name__ == '__main__':
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_pdb",   type=str, default="No", help="")
    parser.add_argument("--Chain", type=str, default="No", help="")

    commandline_args = vars(parser.parse_args())
    BfactorPerChain(commandline_args['input_pdb'], commandline_args['Chain'] )
    BfactorPerChainLinePlot(commandline_args['input_pdb'], commandline_args['Chain'] )



#ppdb = PandasPdb().read_pdb('../../ComPath-DataDownload/PDBAPI/Compath_db/PDB_downlaod_rsync/pdb/pdb7l1z.ent')
#BfactorPerChain(PDBID= ppdb, Chain="A")
#BfactorPerChainLinePlot(PDBID= ppdb, Chain="A")




