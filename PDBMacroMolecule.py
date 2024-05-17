import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pdbecif.mmcif_io import CifFileReader
import pandas as pd


def orgfromSRCfull(data):
	framefull = []
	resultdict = {}
	df_entity_TC = df_entitySYNrename  = df_entityNATrename  = []
	for k,v in data.items():
		if '_entity_src_gen' not in v: continue
		entity = v['_entity_src_gen']
		if 'pdbx_gene_src_scientific_name' not in entity:continue
		if isinstance(entity['entity_id'], (str, float,int)):
			df_entity = pd.DataFrame.from_dict([entity])
			df_entity_TC = df_entity[['entity_id','pdbx_gene_src_scientific_name','pdbx_gene_src_ncbi_taxonomy_id']]
			framefull.append(df_entity_TC)
		else:
			df_entity = pd.DataFrame.from_dict(entity)
			df_entity_TC = df_entity[['entity_id','pdbx_gene_src_scientific_name','pdbx_gene_src_ncbi_taxonomy_id']] 
			framefull.append(df_entity_TC)
			
	for k,v in data.items():
		if '_pdbx_entity_src_syn' not in v: continue
		entitySYN = v['_pdbx_entity_src_syn']
		if 'organism_scientific' not in entitySYN:continue
		if isinstance(entitySYN['entity_id'], (str, float,int)):
			df_entitySYN = pd.DataFrame.from_dict([entitySYN])
			df_entitySYN_TC = df_entitySYN[['entity_id','organism_scientific','ncbi_taxonomy_id']]
			df_entitySYNrename = df_entitySYN_TC.rename(columns={'organism_scientific': 'pdbx_gene_src_scientific_name', 'ncbi_taxonomy_id': 'pdbx_gene_src_ncbi_taxonomy_id'})
			framefull.append(df_entitySYNrename)
		else:
			df_entitySYN = pd.DataFrame.from_dict(entitySYN)
			df_entitySYN_TC = df_entitySYN[['entity_id','organism_scientific','ncbi_taxonomy_id']]
			df_entitySYNrename = df_entitySYN_TC.rename(columns={'organism_scientific': 'pdbx_gene_src_scientific_name', 'ncbi_taxonomy_id': 'pdbx_gene_src_ncbi_taxonomy_id'})
			framefull.append(df_entitySYNrename)
			

	for k,v in data.items():
		if '_entity_src_nat' not in v:continue
		entityNAT = v['_entity_src_nat']
		if 'pdbx_organism_scientific' not in entityNAT:continue
		if isinstance(entityNAT['entity_id'], (str, float,int)):
			df_entityNAT = pd.DataFrame.from_dict([entityNAT])
			df_entityNAT_TC = df_entityNAT[['entity_id','pdbx_organism_scientific','pdbx_ncbi_taxonomy_id']]
			df_entityNATrename = df_entityNAT_TC.rename(columns={'pdbx_organism_scientific': 'pdbx_gene_src_scientific_name', 'pdbx_ncbi_taxonomy_id': 'pdbx_gene_src_ncbi_taxonomy_id'})
			framefull.append(df_entityNATrename)
		else:
			df_entityNAT = pd.DataFrame.from_dict(entityNAT)
			df_entityNAT_TC = df_entityNAT[['entity_id','pdbx_organism_scientific','pdbx_ncbi_taxonomy_id']] 
			df_entityNATrename = df_entityNAT_TC.rename(columns={'pdbx_organism_scientific': 'pdbx_gene_src_scientific_name', 'pdbx_ncbi_taxonomy_id': 'pdbx_gene_src_ncbi_taxonomy_id'})
			framefull.append(df_entityNATrename)
			
	framefull = [df_entity_TC, df_entitySYNrename , df_entityNATrename]
	filtered_list = [ ele for ele in framefull if len(ele) != 0 ]
	if len(filtered_list) >= 1:
		result = pd.concat(filtered_list, ignore_index=True)
		resultdict["Oragnism"] = result.to_dict("records")
		return(resultdict)
	else:
		resultdict["Oragnism"] = [{'entity_id': '1','pdbx_gene_src_scientific_name': 'NA', 'pdbx_gene_src_ncbi_taxonomy_id': 'NA'}, {'entity_id': '2','pdbx_gene_src_scientific_name': 'NA', 'pdbx_gene_src_ncbi_taxonomy_id': 'NA'}]
		return(resultdict)



def Entitypoly(data):
	
	for k, v in data.items():
		
		try:
			entitypoly = v['_entity_poly']
			
			if isinstance(entitypoly['entity_id'], (str, float, int)):
				df_entitypoly = pd.DataFrame.from_dict([entitypoly])
				df_entitypoly_TC = df_entitypoly[['entity_id', 'type', 'pdbx_strand_id']]
				dict_entitypoly_TC = df_entitypoly_TC.to_dict("records")
				return(dict_entitypoly_TC)
			else:
				df_entitypoly = pd.DataFrame.from_dict(entitypoly)
				df_entitypoly_TC = df_entitypoly[['entity_id', 'type', 'pdbx_strand_id']]
				dict_entitypoly_TC = df_entitypoly_TC.to_dict("records")
				return(dict_entitypoly_TC)
		except KeyError:
			dict_entitypoly_TC = [{'entity_id': '1', 'type': 'NA', 'pdbx_strand_id': 'NA'}]
			return(dict_entitypoly_TC)
			

def StruRef(data):
	
	for k,v in data.items():
		try:
			valueslist = v['_struct_ref']
			if isinstance(valueslist['id'], (str, float,int)):
				df_strurefseq = pd.DataFrame.from_dict([valueslist])
				df_strurefseq_TC = df_strurefseq[['pdbx_db_accession','entity_id']]
				dict_strurefseq_TC = df_strurefseq_TC.to_dict('records')
				return(dict_strurefseq_TC)
			else:
				df_strurefseq = pd.DataFrame.from_dict(valueslist)
				df_strurefseq_TC = df_strurefseq[['pdbx_db_accession','entity_id']]
				#groupdf = df_strurefseq_TC.groupby('entity_id').agg(lambda x: list(x))
				#groupdfReset = groupdf.reset_index()
				dict_strurefseq_TC = df_strurefseq_TC.to_dict('records')
				return(dict_strurefseq_TC)
		except KeyError:
			dict_strurefseq_TC = [{'entity_id': '1', 'pdbx_db_accession': 'NA'}]
			return(dict_strurefseq_TC)


def EntityDetail(data):
	
	for k,v in data.items():
		try:
			valueslist = v['_entity']
			if isinstance(valueslist['id'], (str, float,int)):
				df_EntityDetail = pd.DataFrame.from_dict([valueslist])
				df_EntityDetail_TC = df_EntityDetail[['id','type','pdbx_description','pdbx_number_of_molecules','pdbx_ec','pdbx_mutation']]
				df_EntityDetail_RN = df_EntityDetail_TC.rename(columns={'id': 'entity_id'})
				dict_EntityDetail = df_EntityDetail_RN.to_dict('records')
				return(dict_EntityDetail)
			else:
				df_EntityDetail = pd.DataFrame.from_dict(valueslist)
				df_EntityDetail_TC = df_EntityDetail[['id','type','pdbx_description','pdbx_number_of_molecules','pdbx_ec','pdbx_mutation']]
				df_EntityDetail_RN = df_EntityDetail_TC.rename(columns={'id': 'entity_id'})
				dict_EntityDetail = df_EntityDetail_RN.to_dict('records')
				return(dict_EntityDetail)
		except KeyError:
			dict_EntityDetail = [{'entity_id': '1', 'type': 'NA', 'pdbx_description': 'NA', 'pdbx_number_of_molecules': 'NA', 'pdbx_ec': 'NA', 'pdbx_mutation': 'NA'}, {'entity_id': '2', 'type': 'NA', 'pdbx_description': 'NA', 'pdbx_number_of_molecules': 'NA', 'pdbx_ec': 'NA', 'pdbx_mutation': 'NA'}]
			return(dict_EntityDetail)
		

def BfactorStats(data):
	
	for k,v in data.items():
		try:
			AtomRecord = v['_atom_site']
			if isinstance(AtomRecord['id'], (str, float,int)):
				df_BF = pd.DataFrame.from_dict([AtomRecord])
				df_BF_TC = df_BF[['group_PDB','label_entity_id','B_iso_or_equiv','auth_comp_id','auth_asym_id','auth_atom_id']]
				df_BF_RN = df_BF_TC.rename(columns={'label_entity_id': 'entity_id'})
				ChainRecords = df_BF_RN['auth_asym_id'].agg(lambda x: x.unique().tolist())
				BfactorRecord = df_BF_RN[(df_BF_RN['auth_asym_id'] == i) & (df_BF_RN['group_PDB'] == "ATOM") & (df_BF_RN['auth_atom_id'] == "CA")].reset_index()
				EntityID = df_BF_RN[(df_BF_RN['auth_asym_id'] == i) & (df_BF_RN['group_PDB'] == "ATOM")]['entity_id'].agg(lambda x: x.unique().tolist())
				BfactorMean = pd.to_numeric(BfactorRecord["B_iso_or_equiv"]).mean()
				BfactorMin  = pd.to_numeric(BfactorRecord["B_iso_or_equiv"]).min()
				BfactorMax  = pd.to_numeric(BfactorRecord["B_iso_or_equiv"]).max()
				Bfactor = {"entity_id": EntityID[0], "chain": i, "Mean": round(BfactorMean,2), "Min": round(BfactorMin,2), "Max": round(BfactorMax,2)}
				return([Bfactor])
			else:
				BF = []
				df_BF = pd.DataFrame.from_dict(AtomRecord)
				df_BF_TC = df_BF[['group_PDB','label_entity_id','B_iso_or_equiv','auth_comp_id','auth_asym_id','auth_atom_id']]
				df_BF_RN = df_BF_TC.rename(columns={'label_entity_id': 'entity_id'})
				ChainRecords = df_BF_RN['auth_asym_id'].agg(lambda x: x.unique().tolist())
				for i in ChainRecords:
					BfactorRecord = df_BF_RN[(df_BF_RN['auth_asym_id'] == i) & (df_BF_RN['group_PDB'] == "ATOM") & (df_BF_RN['auth_atom_id'] == "CA")].reset_index()
					EntityID = df_BF_RN[(df_BF_RN['auth_asym_id'] == i) & (df_BF_RN['group_PDB'] == "ATOM")]['entity_id'].agg(lambda x: x.unique().tolist())
					BfactorMean = pd.to_numeric(BfactorRecord["B_iso_or_equiv"]).mean()
					BfactorMin  = pd.to_numeric(BfactorRecord["B_iso_or_equiv"]).min()
					BfactorMax  = pd.to_numeric(BfactorRecord["B_iso_or_equiv"]).max()
					Bfactor = {"entity_id": EntityID[0], "chain": i, "Mean": round(BfactorMean,2), "Min": round(BfactorMin,2), "Max": round(BfactorMax,2)}
					BF.append(Bfactor)
				
				return(BF)
		except KeyError:
			dict_BF = [{'group_PDB': '1', 'entity_id': '1', 'B_iso_or_equiv': '0', 'auth_comp_id': 'NA', 'auth_asym_id': 'NA', 'auth_atom_id': 'NA'}]
			return(dict_BF)
		

def MicroMoleculeRecord(D):
	dict_group = {}
	Chain_len  = {}
	Org_list = orgfromSRCfull(D)
	Ent_poly = Entitypoly(D)
	Stru_Ref = StruRef(D)
	Ent_Det	 = EntityDetail(D)
	BFlist	 = BfactorStats(D)
	df_EP	= pd.DataFrame.from_dict(Ent_poly)
	df_OL	= pd.DataFrame.from_dict(Org_list["Oragnism"])
	df_SR	= pd.DataFrame.from_dict(Stru_Ref)
	df_ET	= pd.DataFrame.from_dict(Ent_Det)
	df_BF   = pd.DataFrame.from_dict(BFlist)
	
	result1    = pd.merge(df_EP,df_OL , on=['entity_id'])
	result2    = pd.merge(result1,df_SR , on=['entity_id'])
	result3    = pd.merge(result2,df_ET , on=['entity_id'])
	result4    = pd.merge(result3,df_BF , on=['entity_id'])
	result5 = result4.groupby('pdbx_strand_id').agg(lambda x: x.unique().tolist())
	result6 = result5.reset_index()
	
	for k, v in D.items():
		try:
			ChainAA = v['_entity_poly']['pdbx_seq_one_letter_code']
			if isinstance(ChainAA, (str, float,int)):
				ChainLen = len(ChainAA.replace("\n",""))
				Chain_len["ChainLen"] = [ChainLen]

			else:
				ChainList = []
				for pp in range(0, len(ChainAA)):
					clen = len(ChainAA[pp].replace("\n",""))
					ChainList.append(clen)
				
				Chain_len["ChainLen"] = ChainList
		except KeyError:
			Chain_len["ChainLen"] = ["NA"]

	result6["ChainLength"] = pd.DataFrame.from_dict(Chain_len)
	dict_group["MicroMolecule"] = result6.to_dict("records")
	return(result6)


def BfactorPerChainLinePlotCif(PDBData, Chain, PDBID):
	

	for k, v in PDBData.items():
		try:
			AtomRecord = v['_atom_site']
			if isinstance(AtomRecord, (str, float,int)):
				df_ATOM = pd.DataFrame.from_dict([AtomRecord])
				df_ATOM_Modified = df_ATOM[(df_ATOM['auth_asym_id'] == Chain) & (df_ATOM['group_PDB'] == "ATOM") & (df_ATOM['auth_atom_id'] == "CA")].reset_index()
				Bfactor = pd.to_numeric(df_ATOM_Modified["B_iso_or_equiv"])
				Bfactor.plot(kind="line")
				plt.title('B-Factors Along the Amino Acid '+ Chain + ' Chain and CA atoms')
				plt.xlabel('Residue Number')
				plt.ylabel('B-factor in $A^2$')
				plt.savefig(PDBID[0]+"_BfatorLine.png")
				
			else:
				df_ATOM = pd.DataFrame.from_dict(AtomRecord)
				df_ATOM_Modified = df_ATOM[(df_ATOM['auth_asym_id'] == Chain) & (df_ATOM['group_PDB'] == "ATOM") & (df_ATOM['auth_atom_id'] == "CA")].reset_index()
				Bfactor = pd.to_numeric(df_ATOM_Modified["B_iso_or_equiv"])
				Bfactor.plot(kind="line")
				plt.title('B-Factors Along the Amino Acid '+ Chain + ' Chain and CA atoms')
				plt.xlabel('Residue Number')
				plt.ylabel('B-factor in $A^2$')
				plt.savefig(PDBID[0]+"_BfatorLine.png")
				
		except KeyError:
			df_ATOM = [{'entity_id': '1', 'pdbx_db_accession': 'NA'}]
			return(df_ATOM)



if __name__ == '__main__':
	parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("--input_pdb",   type=str, default="No", help="")
	parser.add_argument("--Chain", type=str, default="No", help="")
	commandline_args = vars(parser.parse_args())
	datadict = CifFileReader().read(commandline_args['input_pdb'])
	PDBID = list(datadict.keys())
	BfactorPerChainLinePlotCif(datadict, commandline_args['Chain'],PDBID)
	MicroMoleculeRecord(datadict).to_csv("MM.csv")