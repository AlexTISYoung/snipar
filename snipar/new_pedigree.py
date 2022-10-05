from distutils.log import error
from genericpath import exists
import itertools
import pickle
from tkinter.messagebox import NO
from unittest import result
import pandas as pd
class Individual():
    newid = itertools.count()
    def __init__(self, is_male, id=None, gen_index=None, father=None, mother=None) -> None:
        self.is_male = is_male
        self.father = father
        self.mother = mother
        if id is None:
            self.id = f"NoneExistantPerson_{next(Individual.newid)}"
        else:
            self.id = id
        self.gen_index = gen_index
    
    def is_real(self):
        return not (self.gen_index is None)
    
    def __repr__(self) -> str:
        text = f"{self.id}: male:{self.is_male}, gen_index:{self.gen_index}"
        if self.father:
            text += f', father_id:{self.father.id}'
        if self.mother:
            text += f', mother_id:{self.mother.id}'
        return text
    
    def __str__(self) -> str:
        return self.__repr__()

class PedigreeComponent():
    newid = itertools.count()
    def __init__(self, base_individual_is_male, base_individual_id, base_individual_gen_index) -> None:
        ind = Individual(base_individual_is_male, id=base_individual_id, gen_index=base_individual_gen_index)
        self.individuals = {ind.id: ind}
        self.id = f'comp_{next(PedigreeComponent.newid)}'
    
    def get_individual(self, id, default=None):
        return self.individuals.get(id, default)

    def add_sibling(self, ind, is_male, id=None, gen_index=None):
        if id in self:
            ind2 = self.individuals(id)
            return self.add_existing_sibling(ind, ind2)
        else:
            return self.add_new_sibling(ind, is_male, id, gen_index)


    def add_new_sibling(self, ind, is_male, id=None, gen_index=None):
        if ind.father is None:
            self.add_father(ind)
        
        if ind.mother is None:
            self.add_mother(ind)            
        
        sib = Individual(is_male, id, gen_index, ind.father, ind.mother)
        self.individuals[sib.id] = sib
        return sib
    
    def add_existing_sibling(self, ind1, ind2):
        if (ind1.mother is None) and (ind2.mother is None):
            ind2.mother = self.add_mother(ind1)            
        elif not (ind1.mother is None) and (ind2.mother is None):
            ind2.mother=ind1.mother
        elif (ind1.mother is None) and not (ind2.mother is None):
            ind1.mother=ind2.mother
        elif not (ind1.mother is None) and not (ind2.mother is None):
            self.merge_nodes(ind1.mother, ind2.mother)

        if (ind1.father is None) and (ind2.father is None):
            ind2.father = self.add_father(ind1)            
        elif not (ind1.father is None) and (ind2.father is None):
            ind2.father=ind1.father
        elif (ind1.father is None) and not (ind2.father is None):
            ind1.father=ind2.father
        elif not (ind1.father is None) and not (ind2.father is None):
            self.merge_nodes(ind1.father, ind2.father)
            
        return ind2
    
    def add_mother(self, ind, id=None, gen_index=None):
        if id in self:
            given_mother = self.individuals[id]
        else:
            given_mother = Individual(False, id, gen_index=gen_index)
            self.individuals[given_mother.id] = given_mother

        if ind.mother is None:
            if given_mother.is_male:
                raise Exception("mother should be female")
            ind.mother = given_mother
        else:            
            self.merge_nodes(ind.mother, given_mother)
            
        return ind.mother

    def add_father(self, ind, id=None, gen_index=None):
        if id in self:
            given_father = self.individuals[id]
        else:
            given_father = Individual(True, id, gen_index=gen_index)
            self.individuals[given_father.id] = given_father

        if ind.father is None:
            if not given_father.is_male:
                raise Exception("father sNonehould be male")
            ind.father = given_father
        else:            
            self.merge_nodes(ind.father, given_father)
            
        return ind.father

    
    def add_offspring(self, p1, is_male, id=None, gen_index=None, p2=None):
        if not (p2 is None) and (p1.is_male == p2.is_male):
            raise Exception(f"masculinity of parents {p1.is_male} and {p2.is_male} are the same")
        
        if p1.is_male:
            given_father = p1
            given_mother = p2
        else:
            given_father = p2
            given_mother = p1

        if id in self:
            offspring = self.individuals[id]
            if offspring.father is None:
                offspring.father = given_father                
            else:
                self.merge_nodes(offspring.father, given_father)
            
            if offspring.mother is None:
                offspring.mother = given_mother                
            else:
                self.merge_nodes(offspring.mother, given_mother)
        else:        
            if p1.is_male:
                offspring = Individual(is_male, id=id, gen_index=gen_index, father=p1, mother=p2)
            else:
                offspring = Individual(is_male, id=id, gen_index=gen_index, father=p2, mother=p1)
        
        self.individuals[offspring.id] = offspring
        return offspring
    
    def merge_nodes(self, ind1, ind2):
        #only functions on none nans
        if ind1 is None or ind2 is None:
            return
        if ind1.is_male != ind2.is_male:
            raise Exception("sex of merging nodes are not equal")
        is_male = ind1.is_male

        if ind1.gen_index is None:
            dest = ind1
            source = ind2
        elif ind2.gen_index is None:
            dest = ind2
            source = ind1
        else:
            raise Exception("can not merge two real nodes")
        
        for o in self.get_offsprings(dest):
            if is_male:
                o.father = source
            else:
                o.mother = source
        
        if not (source.mother is None) and not (dest.mother is None):
            self.merge_nodes(source.mother, dest.mother)
        elif (source.mother is None) and not (dest.mother is None):
            source.mother = dest.mother
        elif not (source.mother is None) and (dest.mother is None):
            dest.mother = source.mother

        if not (source.father is None) and not (dest.father is None):
            self.merge_nodes(source.father, dest.father)
        elif (source.father is None) and not (dest.father is None):
            source.father = dest.father
        elif not (source.father is None) and (dest.father is None):
            dest.father = source.father
        
        self.individuals.pop(dest.id)
            
    def merge(self, other):
        for id, ind in other.individuals.items():
            self.individuals[id] = ind

    def get_offsprings(self, parent1, parent2=None):
        if parent2 is None:
            return [i for i in self.individuals.values() if (i.father==parent1) or (i.mother==parent1)]
        else:
            p1 = (i.father==parent1) and (i.mother==parent2)
            p2 = (i.father==parent2) and (i.mother==parent1)
            return [i for i in self.individuals if p1 or p2]
    
    def get_siblings(self, offspring):
        if offspring.father is None or offspring.mother is None:
            return []
        return self.get_offsprings(offspring.father, offspring.mother)

    def get_possible_imputations(self):
        info_to_imputation_ids = {}
        for ind in self.individuals:
            if ind.gen_index is None:
                offsprings = [o for o in self.get_offsprings(ind) if not (o.gen_index is None)]
                siblings = [s for s in self.get_siblings(ind) if not (s.gen_index is None)]
                parents =  [p for p in [ind.father, ind.mother] if not (p is None) and not (p.gen_index is None)]
                information = (offsprings, siblings, parents)
                ids = info_to_imputation_ids.get(information, [])
                ids.append(ind.id)
                info_to_imputation_ids[information] = ids
        return info_to_imputation_ids
    
    def __contains__(self, key):
        return key in self.individuals
    
    def __str__(self) -> str:
        return self.id
    
    def __repr__(self) -> str:
        return self.id

    def get_desc(self):
        result = f"component {self.id}:\n"
        for key, val in self.individuals.items():
            result += f"{val}\n"
        return result
    
        

class Pedigree():
    def __init__(self) -> None:
        self.components = []
        self.ind_to_components = {}
    
    def get_comp_individual(self, id):
        comp = self.ind_to_components[id]
        ind = comp.get_individual(id)
        return comp, ind

    def add_individual(self, id, is_male, gen_index=None):
        if id in self.ind_to_components:
            comp, ind = self.get_comp_individual(id)
            if ind.is_male != is_male or ind.gen_index != gen_index:
                breakpoint()
                raise Exception("An individual with the same ID and different information exists")
        else:
            comp = PedigreeComponent(is_male, id, gen_index)
            ind = list(comp.individuals.values())[0]
            self.components.append(comp)
            self.ind_to_components[ind.id] = comp
        return ind


    def add_siblings(self, id1, is_male1, gen_index1, id2,  is_male2,  gen_index2):
        if (id1 in self.ind_to_components) and not (id2 in self.ind_to_components):
            comp, ind = self.get_comp_individual(id1)
            sib = comp.add_sibling(ind, is_male2, id2, gen_index2)
            self.ind_to_components[sib.id] = comp
        elif not (id1 in self.ind_to_components) and (id2 in self.ind_to_components):
            comp, ind = self.get_comp_individual(id2)
            sib = comp.add_sibling(ind, is_male1, id1, gen_index1)
            self.ind_to_components[sib.id] = comp
        elif not (id1 in self.ind_to_components) and not (id2 in self.ind_to_components):
            comp = PedigreeComponent(is_male1, id1, gen_index1)
            ind = list(comp.individuals.values())[0]
            sib = comp.add_sibling(ind, is_male2, id2, gen_index2)
            self.ind_to_components[ind.id] = comp
            self.ind_to_components[sib.id] = comp
            self.components.append(comp)
        else:
            comp1 = self.ind_to_components[id1]
            comp2 = self.ind_to_components[id2]
            ind1 = comp1.get_individual(id1)
            ind2 = comp2.get_individual(id2)
            if comp1 == comp2:
                comp1.add_existing_sibling(ind1, ind2)
            else:
                comp1.merge(comp2)
                for id, ind in comp2.individuals.items():
                    self.ind_to_components[id] = comp1
                self.components.remove(comp2)
                comp1.add_existing_sibling(ind1, ind2)
    
    def add_parent_offspring(self, id_p, is_male_p, gen_index_p, id_o,  is_male_o,  gen_index_o):
        c_o = self.ind_to_components.get(id_o, None)
        c_p = self.ind_to_components.get(id_p, None)
        if c_o == c_p == None:            
            comp = PedigreeComponent(is_male_p, id_p, gen_index_p)
            par = list(comp.individuals.values())[0]
            offspring = comp.add_offspring(par, is_male_o, id_o, gen_index_o)
            self.ind_to_components[par.id] = comp
            self.ind_to_components[offspring.id] = comp
            self.components.append(comp)
        elif c_o and (c_p ==None or c_p == c_o):
            comp, offspring = self.get_comp_individual(id_o)
            if is_male_p:
                par = comp.add_father(offspring, id_p, gen_index_p)
            else:
                par = comp.add_mother(offspring, id_p, gen_index_p)
            self.ind_to_components[par.id] = comp
        elif c_p and (c_o ==None or c_o == c_p):
            comp, par = self.get_comp_individual(id_p)
            offspring = comp.add_offspring(par, is_male_o, id_o, gen_index_o)
            self.ind_to_components[offspring.id] = comp
        else:
            _, par = self.get_comp_individual(id_p)
            _, offspring = self.get_comp_individual(id_o)                
            c_p.merge(c_o)
            for id, ind in c_o.individuals.items():
                self.ind_to_components[id] = c_p
            self.components.remove(c_o)
            offspring = c_p.add_offspring(par, is_male_o, id_o, gen_index_o)
    
    def get_desc(self):
        result = ""
        for c in set(self.ind_to_components.values()):
            result += c.get_desc()+"\n***************************\n"
        result+="------------------------------"
        return result

    def save(self, address):
        with open(address, "wb") as f:        
            pickle.dump(self, f)
    
    def load(address):
        result = None
        with open(address, "rb") as f:        
            result = pickle.load(f)
        return result

    def from_classic_pedigree(classic_ped_address, fam_address=None, pedigree_nan = '0'):
        ped = Pedigree()
        classic_ped = pd.read_csv(classic_ped_address, delim_whitespace=True)
        if fam_address:
            fam = pd.read_csv(fam_address, delim_whitespace=True, header=None)
            fam = {i:index for index,i in enumerate(fam[1].values)}
        else:
            fam = {}
        for index, row in classic_ped.iterrows():
            # FID IID FATHER_ID MOTHER_ID
            #TODO get an agesex file maybe? right now we are assuming unkown sex is female
            id = row['IID']
            father_id = row['FATHER_ID']
            mother_id = row['MOTHER_ID']
            is_male = id in classic_ped["FATHER_ID"].values
            gen_index = fam.get(id, None)
            father_gen_index = fam.get(father_id, None)
            mother_gen_index = fam.get(mother_id, None)
            ped.add_individual(id, is_male, gen_index)
            if row["FATHER_ID"] != pedigree_nan:
                ped.add_parent_offspring(father_id, True, father_gen_index,
                                        id, is_male, gen_index)
            if row["MOTHER_ID"] != pedigree_nan:
                ped.add_parent_offspring(mother_id, False, mother_gen_index,
                                        id, is_male, gen_index)
        return ped

    def from_king(king_address, agesex_address, fam_address=None):
        ped = Pedigree()
        king = pd.read_csv(king_address, delim_whitespace=True)[["FID1", "ID1",	"FID2",	"ID2", "InfType"]]
        agesex_df = pd.read_csv(agesex_address, delim_whitespace=True)[["FID", "IID", "age", "sex"]]
        age = {}
        is_male = {}
        for index, row in agesex_df.iterrows():
            id = row["IID"]
            age[id] = row['age']
            is_male[id] = (row['sex']=='M')

        if fam_address:
            fam = pd.read_csv(fam_address, delim_whitespace=True, header=None)
            fam = {i:index for index,i in enumerate(fam[1].values)}
        else:
            fam = {}
        
        for index, row in king.iterrows():
            id1 = row['ID1']
            id2 = row['ID2']
            if row['InfType'] == 'FS':
                ped.add_siblings(id1, is_male[id1], fam.get(id1,None),
                                id2, is_male[id2], fam.get(id2,None),
                )
            elif row['InfType'] == 'PO':
                if age[id2]>age[id1]:
                    id1, id2 = id2, id1
                ped.add_parent_offspring(id1, is_male[id1], fam.get(id1, None),
                                        id2, is_male[id2], fam.get(id2, None)
                )
            else:
                #half sibs, cousins, MZ, DZ
                pass
        return ped
    
    def get_components(self):
        return list(self.ind_to_components.values())


            
        
        
p = Pedigree()
p.add_siblings("ali", True, 12, "moeen", True, 22)
p.components
p.components[0]
c1 = p.components[0]
p.add_parent_offspring("forough", False, 13, "ali", True, 22)

p.add_parent_offspring("zahra", False, 13, "mohsen", True, 22)
p.add_parent_offspring("mohammad", True, 13, "hossein", True, 22)
p.add_siblings("mohsen", True, 12, "hossein", True, 22)

p.add_parent_offspring("kamal", True, 13, "jamal", True, 22)
p.add_parent_offspring("katal", True, 13, "kamal", True, 22)

p.add_siblings("zahra", False, 12, "forough", False, 22)

p.add_siblings("c3_2", False, 12, "c3_3", False, 22)
p.add_siblings("kamal", False, 12, "c3_3", False, 22)

print(p.get_desc())


ped = Pedigree.from_classic_pedigree("sample.ped", "sample_sib1.fam")
ped.save('tmp')
a = Pedigree.load('tmp')
print(a.ind_to_components)


ped = Pedigree.from_king("sample.king", 'sample.agesex', "sample_sib1.fam")
print(a.ind_to_components['1738_0'].get_desc())
print(a.ind_to_components['12_0'].get_desc())
print(a.ind_to_components['2278_0'].get_desc())
breakpoint()
