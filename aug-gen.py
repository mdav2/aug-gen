#!/usr/bin/python3
# coding: utf-8

# In[1]:


import copy
import sedre
import sys




# In[2]:




# In[7]:


class basis:
    #Object to contain information about an atomic orbital basis set
    #The I/O functionality is limited to G94 format
    #
    def __init__(self):
        self.atomics = {}
        self.misc = {}
    
    def add_atom(self,text):
        #takes as input a block of text w/ basis set
        #for a single nucleus (e.g. H, or C, etc)
        #
        stext=text.strip().split('\n') #give text as long string
        nucleus=stext[0].split()[0] #Grab the nucleus from the first line (e.g. H  0)
        primitives=[] #local variable for holding collections of primitives to be
                      #condensed into contracted fns vvv
        contracted=[] #local var, for holding contracted functions before modifying
                      #self.atomics w/ the new info
        temp_contracted={} #?
        ang_mom_hi=0 #highest ang momentum, used for convenience elsewhere
        for oline in stext[1:]:
            
             
            line=oline
            sline=oline.strip().split()
            if len(sline) == 2: #then this defining a primitive w/in a contracted fn
                line=oline.replace('D','E') #make them readable as exponential format by Python
                #print(line)
                sline=line.strip().split()
                primitives.append((float(sline[0]),float(sline[1]))) #add the exp and coeff to primitives
                
            elif len(sline) == 3: #then this is defining a new contracted fn
                if len(temp_contracted) != 0: #then this is NOT the first contracted fn for this atom
                    temp_contracted['primitives'] = primitives #can I just get rid of the primitives varaible?
                    contracted.append(temp_contracted) #save the previous contracted fn
                    primitives=[] #clear primitives
                    temp_contracted={} #clear temp_contracted
                sline=line.strip().split()
                ang_mom_qnumber = ang_mom_symbols.index(str(sline[0]).upper()) #convert S, P, D ... to l=0,1,2 etc
                temp_contracted['angular_momentum']   = {'symbol':str(sline[0]).upper(),'qnumber':ang_mom_qnumber}
                if ang_mom_qnumber > ang_mom_hi: #then this is the new highest ang mom
                    ang_mom_hi = ang_mom_qnumber 
                temp_contracted['num_primitives']     = int(sline[1])
                temp_contracted['gen_contract_coeff'] = float(sline[2])
                
            elif (len(sline) > 0) and (sline[0] == '****'): #delimiter for end of atom in basis set
                temp_contracted['primitives'] = primitives
                contracted.append(temp_contracted) #save the previous contracted fn
                break
                
        self.atomics[nucleus] = contracted #save basis fns to self.atomics
        self.misc[nucleus] = [nucleus,ang_mom_hi]
        
    def print_out(self,atom=None):
        #prints out the appropriate atomic basis fns in G94 format
        ostr = ''
        if not atom: #default is to print all of the atoms in self.atomics
            for key in self.atomics:
                olines = [f'{key}  0'] #create e.g. H   0 line
                for basis_fn in self.atomics[key]:
                    symbol=basis_fn['angular_momentum'] #<v add contracted basis fn line
                                                        #order of contracted fns does not seem to matter.
                    olines.append(f'{basis_fn["angular_momentum"]["symbol"]}   '+                                    f'{basis_fn["num_primitives"]}   '+                                    f'{basis_fn["gen_contract_coeff"]}')
                    
                    for i in range(basis_fn["num_primitives"]): #add all the primitive lines
                            olines.append(' '*6 + '{:.6E}      {:.6E}'.format(*basis_fn['primitives'][i]).replace('E','D'))
                olines.append('****')
                ostr=ostr+'\n'.join(olines)+'\n'
        return ostr.replace('E','D') #this is hacky, but I don't think E ever needs to be in a gbs



# In[111]:


def naug(basis,naug=2):
    #takes in a basis and a desired n-aug number (usually d-aug-, so default is 2)
    #returns a new basis with added exponents for next augmentation level
    #"the usual way" 😭😭😭
    newbasis=copy.deepcopy(basis)
    for atom in basis.atomics:
        amhi=basis.misc[atom][1] #highest angmom
        print(f'ATOM: {atom}')
        print(f'AMHI: {amhi}')
        for idx,l in enumerate(ang_mom_symbols[:amhi+1]):
            print(f'->AM:{l}')
            expm1 = None
            exp0 = None
            
            for idx2,basis_fn in enumerate(basis.atomics[atom]):
                if (idx2+1) == len(basis.atomics[atom]): #or if it's at the end, marks the last fn as appropriate
                    last = idx2
                    break
                elif idx2 == 0:
                    last=None
                    continue
                elif (basis_fn['angular_momentum']['qnumber'] < basis.atomics[atom][idx2+1]['angular_momentum']['qnumber'])                      and (basis_fn['angular_momentum']['qnumber'] == idx):
                    last=idx2
                    break

               # if (basis_fn['angular_momentum']['qnumber'] > idx):
               #     print('->route1')
               #     last = idx2 -1 #finds the split between S and P, P and D, etc
               #     break
                
            print(f'LAST:{last}')
            if last is None:
                continue
            expm1 = basis.atomics[atom][last-1]['primitives'][-1][0]
            exp0 = basis.atomics[atom][last]['primitives'][-1][0]
            newexp=exp0*(exp0/expm1)**(naug-1) #alpha*(beta^(x-1)) generator, alpha is smallest exp, beta is ratio of last two exp
            print(f'->NEW EXPONENT: {newexp}')
            basis_fn=copy.deepcopy(basis.atomics[atom][last])
            basis_fn['primitives'][-1] = (newexp,1.0)
            newbasis.atomics[atom].insert(last+1+idx, basis_fn)
        print('\n')
    return newbasis


# In[112]:


ang_mom_symbols = ['S','P','D','F','G','H','I','K','L','M','N'] #better safe than sorry!
delim_line = "sed -En \'/[A-Z][ ]+0/=\' {}"
if __name__ == "__main__":
    aug=basis()
    inputfile=sys.argv[1]
    start = sedre.common.line_num('[A-Z][ ]+0',inputfile)[0] - 1
    with open(inputfile,'r') as f: text = f.readlines()
    text = ''.join(text[start:]).strip().split('****')
    atoms = []
    for atom in text[:-1]:
        #print(atom)
        atoms.append(atom+'****')
    for atom in atoms:
        aug.add_atom(atom)
    
    
    # In[113]:
    
    
    daug=naug(aug,naug=2)
    
    
    # In[115]:
    
    
    taug=naug(daug,naug=3)
    
    
    # In[118]:
    
    
    qaug=naug(taug,naug=4)
    
    
    # In[121]:
    
    
    for idx,naug in enumerate([daug,taug,qaug]):
        with open('my-'+str(idx+2)+'-'+inputfile,'w') as f: f.write(naug.print_out())
    
