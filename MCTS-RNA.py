from subprocess import Popen, PIPE
from math import *
import random
import numpy as np
import RNA
from copy import deepcopy
from types import IntType, ListType, TupleType, StringTypes
import itertools
import time
import math
import argparse
import subprocess

class RNAstructure:

    def __init__(self):
        self.search = [0,0,0,0,0,0,0,0,0,0]
        self.basepairs=["AU", "CG", "GC", "UA","GU","UG"]
        self.bases=["A","C","G","U","AU", "CG", "GC", "UA","GU","UG"]
        self.base=["A","C","G","U"]
        self.position=str_uindex+str_index
        self.n= len(str_uindex+str_index)

    def count(self):
        self.number+=1

    def Clone(self):

        st = RNAstructure()
        st.search = self.search[:]
        st.basepairs = self.basepairs[:]
        st.position= self.position[:]
        return st


    def Rewards(self,k):
        #copy_unpairedposition=list(unpairedposition)
        #copy_bppused=list(bppused)
        if k > len(str_uindex)-1:
            posbasep=self.position[len(str_uindex):self.n]
            posbase=self.position[0:len(str_uindex)]
            e=list(itertools.chain(*posbasep))
            for i in range(len(a)):
                posbase.insert(b[i],e[c[i]])
            mutated_s= ''.join(map(str, posbase))
            mutated_str1=RNA.fold(mutated_s)
            mutated_str=mutated_str1[0]

            d=0.0
            g=0.0
            n=len(s)
            for i in range(len(s)):
                if mutated_str[i]!=s[i]:
                    d=d+1
            g=(n-d)/n
            if g==1.0:
                solution.append(mutated_s)
                return g
            else:
                return g




        if k <= len(str_uindex)-1:
            posbasep=self.position[len(str_uindex):self.n]
            posbase=self.position[0:len(str_uindex)]
            e=list(itertools.chain(*posbasep))
            for i in range(len(a)):
                posbase.insert(b[i],e[c[i]])
            mutated_s= ''.join(map(str, posbase))
            mutated_str1=RNA.fold(mutated_s)
            mutated_str=mutated_str1[0]
            d=0.0
            g=0.0
            n=len(s)
            for i in range(len(s)):
                if mutated_str[i]!=s[i]:
                    d=d+1
            g=(n-d)/n
            if g==1.0:
                solution.append(mutated_s)
                return g
            else:
                return g



    def SelectBasePairs(self,move):
        self.search[move]=self.basepairs[move]

    def Selectstateposition(self,l):
        if l > len(str_uindex)-1:
            self.position[l]=midea[l-len(str_uindex)]
        else:
            self.position[l]=copy_str_uindex[l]

    def Simulation(self,l):
        if l>len(str_uindex)-1:
            self.position[l]=random.choice(BASEPAIRS)
        else:

            self.position[l]=random.choice(bases)
    def simulation1(self, l,needgc,count_number):
        if l>len(str_uindex)-1 and count_number<=needgc:
            self.position[l]=random.choice(["CG","GC"])
            count_number=count_number+1
        if l>len(str_uindex)-1 and count_number>needgc:
            self.position[l]=random.choice(["AU","UA"])
        if l<=len(str_uindex)-1:
            self.position[l]=random.choice(["A","U"])
    def simulationGC(self,l):
        self.position[l]=random.choice(["CG","GC"])

    def simulationAU(self,l):
         self.position[l]=random.choice(["AU","UA"])
    def simulationunpairedAU(self,l):
        self.position[l]=random.choice(["A","U"])

    def simulationunpairedGC(self,l):
        self.position[l]=random.choice(["G","C"])

    def SelectPosition(self,m,k):
        self.position[k]=self.bases[m]


    def Getubpp(self):
        return [i for i in np.arange(0,4) if self.search[i] not in ["A","U","C","G"]]
    def Getbpp(self):
        return [i for i in np.arange(4,10) if self.search[i] not in ["AU", "CG", "GC", "UA","GU","UG"]]
    def GetSearch(self):
        return [i for i in range(len(self.search)) if self.search[i] not in ["A","U","C","G","AU", "CG", "GC", "UA","GU","UG"]]
    def GetPositions(self):
        return[i for i in range(len(self.position)) if self.position[i] not in ["A","U","C","G","AU", "CG", "GC", "UA","GU","UG"]]


class Node:

    def __init__(self, position = None, pt = None , parent = None, state = None):
        self.position = position
        self.pt = pt
        self.parentNode = parent
        self.childNodes = []
        self.child=None
        self.wins = 0
        self.visits = 0
        self.untriedSearches = state.GetSearch()
        self.untriedubpp=state.Getubpp()
        self.untriedbpp=state.Getbpp()
        self.untriedPositions=state.GetPositions()


    def Selectnode(self):

        s = sorted(self.childNodes, key = lambda c: c.wins/c.visits + 0.1*sqrt(2*log(self.visits)/c.visits))[-1]
        return s

    def Addnode(self, m, k, s):

        n = Node(position = m, pt=k, parent = self, state = s)
        if k in self.untriedPositions:
            self.untriedPositions.remove(k)
        self.childNodes.append(n)
        self.child=n
        return n

    def Update(self, result):

        self.visits += 1
        self.wins += result


def MCTS(root, k, verbose = False):


    running_time=time.time()
    out_time=running_time+60*10
    rootnode = Node(state = root)
    state = root.Clone() # but this state is the state of the initialization .  too important !!!


    while time.time()<=out_time:

        node = rootnode # important !    this node is different with state / node is the tree node
        state = root.Clone() # but this state is the state of the initialization .  too important !!!
        posi=[]
        posl=[]
        poslalpha=[]
        need=[]
        count_number=0
        count_number1=0
        pa=[]
        upa=[]
        while node.untriedubpp == [] or node.untriedbpp==[]:

            node = node.Selectnode()

            state.SelectPosition(node.position,node.pt)

        if node.untriedPositions != []:
            if k > len(str_uindex)-1:
                if len(node.untriedbpp)==6:
                    k = random.choice(node.untriedPositions)
            else:
                if len(node.untriedubpp)==4:
                    k = random.choice(node.untriedPositions)
            if node.untriedbpp != 6 or node.untriedubpp!=4:
                if node.child!=None:
                    k=node.child.pt

        if k > len(str_uindex)-1:
            if node.untriedbpp!=[]:
                #print node.untriedbpp
                m=random.choice(node.untriedbpp)
                #print m
                node.untriedbpp.remove(m)
                state.SelectPosition(m,k)
                node=node.Addnode(m,k,state)
        else:
            if node.untriedubpp!=[]:
                m=random.choice(node.untriedubpp)
                node.untriedubpp.remove(m)
                state.SelectPosition(m,k)
                node=node.Addnode(m,k,state)

        posi=state.position
        goal=str_index+str_uindex
        for i in range(len(state.position)):
            if goal[i] not in posi:
                posl.append(goal[i])
            else:
                need.append(goal[i])
            if posi[i] not in goal:
                    poslalpha.append(posi[i])
        if len(poslalpha)<=len(str_index):
            eposl=list(itertools.chain(*poslalpha))
        else:
            eposl111=poslalpha[len(poslalpha)-len(str_index):len(poslalpha)]
            eposl1=list(itertools.chain(*eposl111))
            eposl=poslalpha[0:len(poslalpha)-len(str_index)]+eposl1
        need_GC=calculate_GC_numbers(eposl,defined_GC,need,poslalpha)

        y=state.GetPositions()


        for i in range(len(y)):
            if y[i]>len(str_uindex)-1:
                pa.append(y[i])
            else:
                upa.append(y[i])

        while pa !=[]:
            cpa=random.choice(pa)

            if count_number<=need_GC:

                state.simulationGC(cpa)
                count_number=count_number+1
                pa.remove(cpa)
            else:
                state.simulationAU(cpa)
                pa.remove(cpa)
        while upa!=[]:
            ucpa=random.choice(upa)

            if count_number<=need_GC:
                new_count=need_GC-count_number
                if count_number1<=new_count*2:
                    state.simulationunpairedGC(ucpa)
                    upa.remove(ucpa)
                    count_number1=count_number1+1
                else:
                    state.simulationunpairedAU(ucpa)
                    upa.remove(ucpa)
            else:
                state.simulationunpairedAU(ucpa)
                upa.remove(ucpa)


        posbasep=state.position[len(str_uindex):state.n]
        posbase=state.position[0:len(str_uindex)]
        e=list(itertools.chain(*posbasep))
        for i in range(len(a)):
            posbase.insert(b[i],e[c[i]])
        mutated_s= ''.join(map(str, posbase))
        ini_seq_pool=[]
        ini_str_pool=[]
        GC_pool=[]
        index_seq=0
        if defined_pseudo==1:
            some_str_mfe,some_str_value=calculate__pseudo_mfe_and_str_pkiss(mutated_s)
            some_str_distance=calculate_structure_distance_pKiss(str_index,len(str_index),some_str_value)


        else:
            some_str_mfe,some_str_value=calculate_mfe_and_str(mutated_s)#this is the nest structures
            some_str_distance=calculate_structure_distance(s,len(s),some_str_value)

        ini_seq_pool.append(mutated_s)
        ini_str_pool.append(some_str_distance)
        GCnum=measureGC(mutated_s)

        GC_pool.append(GCnum)
        for i in range(50):
            paired_pos,dif_ini=dif_str(some_str_value)
            mutated_seq=GCcontent(defined_GC,GCnum,paired_pos,posbase,posl)
            mutated_seq1=check_GC_base3(dif_ini,mutated_seq,posl,defined_GC)
            mutated_seq2=''.join(map(str, mutated_seq1))
            GCnum=measureGC(mutated_seq2)
            GC_pool.append(GCnum)
            if defined_pseudo==1:
                mfe,kkk=pseudoknot_pkiss(mutated_seq2)
                new_str_distance=calculate_structure_distance_pKiss(str_index+str_uindex,len(str_index+str_uindex),kkk)
            else:
                kkk=RNA.fold(mutated_seq2)[0]
                new_str_distance=calculate_structure_distance(s,len(s),kkk)
            some_str_value=kkk
            some_ini_seq=mutated_seq2
            ini_seq_pool.append(mutated_seq2)
            ini_str_pool.append(new_str_distance)
            index_seq=index_seq+1
            ggg=abs(defined_GC-GCnum)

            if ini_str_pool[index_seq]==1.0:
                break
        max_idx = np.argmax(ini_str_pool)
        GCnew=GC_pool[max_idx]


        max_val = ini_str_pool[max_idx]
        seq=ini_seq_pool[index_seq]
        ggg=abs(defined_GC-GCnum)
        gggg=abs(defined_GC-GCnew)
        if ini_str_pool[index_seq]==1.0 and ggg<=defined_gd:
            break

        if ini_str_pool[index_seq]==1.0:
            if ggg<=0.01:
                re=1.0+1.0

            else:
                re=1.0+0.0

        if max_val<1.0:
            if gggg<=0.01:
                re=1.0+max_val
            else:
                re=0.0+max_val


        while node != None:
            node.Update(re)
            node = node.parentNode

    return seq, max_val, GCnum



def UCTRNA():
    one_search_start_time=time.time()
    time_out=one_search_start_time+60*10
    state = RNAstructure()
    print "search length:" + str(state.n) + "\n"
    k=random.choice(state.GetPositions())
    m,goal,GCC = MCTS(root = state, k=k, verbose = False)
    print "Solution:" + str(m)
    if goal==1.0:
        finish_time=time.time()-one_search_start_time
    else:
        finish_time=0.0


    return goal,GCC,finish_time


def MCTSnoGC(root, itermax, k, verbose = False):


    running_time=time.time()
    out_time=running_time+60*10
    rootnode = Node(state = root)

    for i in range(itermax):
        if time.time() >= out_time:
            break

        node = rootnode # important !    this node is different with state / node is the tree node
        state = root.Clone() # but this state is the state of the initialization .  too important !!!
        posi=[]
        posl=[]

        while node.untriedubpp == [] or node.untriedbpp==[]:

            node = node.Selectnode()
            state.SelectPosition(node.position,node.pt)
        if node.untriedPositions != []:
            if k > len(str_uindex)-1:
                if len(node.untriedbpp)==6:
                    k = random.choice(node.untriedPositions)
            else:
                if len(node.untriedubpp)==4:
                    k = random.choice(node.untriedPositions)

        if k > len(str_uindex)-1:
            if node.untriedbpp!=[]:
                m=random.choice(node.untriedbpp)
                node.untriedbpp.remove(m)
                state.SelectPosition(m,k)
                node=node.Addnode(m,k,state)
        else:
            if node.untriedubpp!=[]:
                m=random.choice(node.untriedubpp)
                node.untriedubpp.remove(m)
                state.SelectPosition(m,k)
                node=node.Addnode(m,k,state)


        posi=state.position
        goal=str_uindex+str_index

        for i in range(len(state.position)):
            if goal[i] not in posi:
                posl.append(goal[i])
        while state.GetPositions() != []:
            state.Simulation(random.choice(state.GetPositions()))


        posbasep=state.position[len(str_uindex):state.n]
        posbase=state.position[0:len(str_uindex)]
        e=list(itertools.chain(*posbasep))
        for i in range(len(a)):
            posbase.insert(b[i],e[c[i]])
        mutated_s= ''.join(map(str, posbase))
        ini_seq_pool=[]
        ini_str_pool=[]
        GC_pool=[]
        index_seq=0
        if defined_pseudo==1:
            some_str_mfe,some_str_value=calculate__pseudo_mfe_and_str(mutated_s)# this is the pseudoknot structure
        else:
            some_str_mfe,some_str_value=calculate_mfe_and_str(mutated_s)#this is the nest structures

        some_str_distance=calculate_structure_distance(s,len(s),some_str_value)
        ini_seq_pool.append(mutated_s)
        ini_str_pool.append(some_str_distance)
        GCnum=measureGC(mutated_s)
        GC_pool.append(GCnum)
        for i in range(50):


            paired_pos,dif_ini=dif_str(some_str_value)
            mutated_seq=check_seq_base(paired_pos,posbase,posl)
            mutated_seq1=check_GC_base(dif_ini,mutated_seq,posl)
            mutated_seq2=''.join(map(str, mutated_seq1))
            GCnum=measureGC(mutated_seq2)
            GC_pool.append(GCnum)
            if defined_pseudo==1:
                kkk=pseudoknot(mutated_seq2)[0]
            else:
                kkk=RNA.fold(mutated_seq2)[0]
            new_str_distance=calculate_structure_distance(s,len(s),kkk)
            some_str_value=kkk
            some_ini_seq=mutated_seq2
            ini_seq_pool.append(mutated_seq2)
            ini_str_pool.append(new_str_distance)
            index_seq=index_seq+1
            if ini_str_pool[index_seq]==1.0:
                break
        max_idx = np.argmax(ini_str_pool)
        GCnew=GC_pool[max_idx]
        max_val = ini_str_pool[max_idx]
        seq=ini_seq_pool[index_seq]
        if ini_str_pool[index_seq]==1.0:
            break
        if max_val<1.0:
            re=max_val
        while node != None:
            node.Update(re)
            node = node.parentNode

    return seq,ini_str_pool[index_seq], GCnew



def UCTRNAnoGC():
    one_search_start_time=time.time()
    time_out=one_search_start_time+60*10
    state = RNAstructure()
    print "search length:" + str(state.n) + "\n"

    k=state.GetPositions()
    m,goal,GC= MCTSnoGC(root = state, itermax = 100000, k=k, verbose = False)

    if goal==1.0:
        finish_time=time.time()-one_search_start_time
    else:
        finish_time=0.0

    print "solution:"+ str(m)
    print "running time:" + str(finish_time)
    print "GC-content:"+str(GC)
    print "structure distance:" + str(goal)


def calculate_structure_distance(structure_s, str_length ,some_str_value):
    sdt=0.0
    sd=0.0
    for i in range(len(structure_s)):
        if some_str_value[i]!=s[i]:
            sd=sd+1
        sdt=(str_length-sd)/str_length
    return sdt


def calculate_structure_distance_pKiss(structure_s, str_length ,some_str_value):
    paired_str=str_index
    unpaired_str=str_uindex
    struc,ustruc=calculate__pseudo_sequence_position_pKiss(some_str_value)
    structure_s_new=struc+ustruc
    sdt=0.0
    sd=0.0
    for i in range(len(str_index)):
        if paired_str[i] not in struc:
            sd=sd+1
    for i in range(len(str_uindex)):
        if unpaired_str[i] not in ustruc:
            sd=sd+1
    sdt=(len(structure_s)-sd)/len(structure_s)
    return sdt


def intialization(structure_s):

    BASEPAIRS = ["AU", "CG", "GC", "UA", "GU", "UG"]
    basepro=[0.2,0.3,0.3,0.2,0.1,0.1]
    CGbases=["CG","GC"]
    AUbases=["AU","UA"]
    GUbases=["GU","UG"]
    CGbases=["CG","GC"]
    CGbases1=random.choice(CGbases)
    AUbases1=random.choice(AUbases)
    GUbases1=random.choice(GUbases)
    j=[CGbases1,AUbases1,GUbases1]
    bases="AGCU"

    return



def pick_with_probility(some_list, probabilities):
    x = random.uniform(0, 1)
    cumulative_probability = 0.0
    for item, item_probability in zip(some_list, probabilities):
        cumulative_probability += item_probability
        if x < cumulative_probability: break
    return item


def calculate_sequence_position(seq):
    stack = []
    struc = []
    ustruc=[]
    for i in xrange(len(seq)):
        if seq[i] == '(':
            stack.append(i)
        if seq[i] == ')':
            struc.append((stack.pop(), i))
        elif seq[i]=='.':
            ustruc.append(i)
    return struc,ustruc


def getinput():
    return input ("percentage of GC : ").lower


def calculate_a(some_str_index):
    a = list(itertools.chain(*some_str_index))
    return a


def calculate_b(some_a):
    b=sorted(some_a)
    return b


def calculate_c(some_a):
    c=sorted(range(len(some_a)),key=lambda x: a[x])
    return c

def getbasepairs(some_str_index):
    midea=[]
    for i in range(len(some_str_index)):
        midea.append(random.choice(BASEPAIRS))
    return midea


def getunbases(some_str_uindex):
    some_copy_str_uindex=list(some_str_uindex)
    for i in range(len(some_str_uindex)):
        some_copy_str_uindex[i]=random.choice(bases)
    return some_copy_str_uindex


def getwholesequence(some_b,some_c,some_d,some_copy_str_uindex):
    for i in range(len(some_c)):
        some_copy_str_uindex.insert(some_b[i],some_d[some_c[i]])
    wholesequence = ''.join(map(str, some_copy_str_uindex))
    return wholesequence,some_copy_str_uindex


def calculate_d(some_midea):
    d = list(itertools.chain(*some_midea))
    return d


def calculate_mfe_and_str(sequence):
    rnafold= RNA.fold(sequence)
    mfe=rnafold[1]
    str_v=rnafold[0]
    return mfe,str_v


def error_diagnosis():
    for i in range(len(paired)):
        if paired[i]:
            pass
    return


def identical_position(input_str,predicted_str):## calculate the some position between initial and mutated
    modified_seq=[]
    modified_pos=[]
    for i in range(len(input_str)):
        if predicted_str[i]==input_str[i]:
            modified_pos.append(i)
            #modified_seq[i]=initial_seq[i]

    return modified_pos


def find_dif_pos_ini():
    dif_pos_ini=[]
    for i in range(len(str_index)):
        if str_index[i] not in paired[i]:
            dif_pos_ini.append(str_index[i])
    return

def find_dif_str_position_between_target_and_predicted(target_seq,predicted_seq):
    break_pairs=["AA","CC","GG","AG","CU","UC","UU","GA"]
    #break_pairs=["UU"]
    comp=["GC","CG"]
    save_paired_pos=[]
    save_unpaired_pos=[]
    ori_paired_pos=[]
    ori_unpaired_pos=[]
    paired_dif_pos=[]

    paired,unpaired=calculate_sequence_position(predicted_seq)

    for i in range(len(paired)):
        if paired[i] not in str_index:
            save_paired_pos.append(paired[i])
            paired_dif_pos.append(random.choice(break_pairs))
        else:
            ori_paired_pos.append(paired[i])


    for i in range(len(unpaired)):
        if unpaired[i] not in str_uindex:
            save_unpaired_pos.append(unpaired[i])
        else:
            ori_unpaired_pos.append(unpaired[i])

    dif_pos_ini=[]
    dif_pos_base=[]

    for i in range(len(str_index)):
        if str_index[i] not in paired:
            dif_pos_ini.append(str_index[i])
            dif_pos_base.append(random.choice(comp))



    return save_paired_pos,dif_pos_ini, paired_dif_pos,dif_pos_base


def dif_str(predicted_seq):

    save_paired_pos=[]
    save_unpaired_pos=[]
    ori_paired_pos=[]
    ori_unpaired_pos=[]
    paired_dif_pos=[]
    paired,unpaired=calculate_sequence_position(predicted_seq)
    for i in range(len(paired)):
        if paired[i] not in str_index:
            save_paired_pos.append(paired[i])

    dif_pos_ini=[]
    dif_pos_base=[]

    for i in range(len(str_index)):
        if str_index[i] not in paired:
            dif_pos_ini.append(str_index[i])

    return save_paired_pos,dif_pos_ini



def assign_to_paired(some_save_paired_pos, predicted_seq):

    break_pairs=["AA","CC","GG","AG","CU","UC","UU","GA"]
    break_pairs=["AA","UU","AC"]

    paired,unpaired=calculate_sequence_position(predicted_seq) #get the position of the new structure


    some_midea=getbasepairs(some_save_paired_pos)
    some_copy_str_uindex=getunbases(unpaired)

    a1=calculate_a(some_save_paired_pos)
    b1=calculate_a(a1)
    c1=calculate_c(a1)
    d1=calculate_(some_midea)
    muta_seq=getwholesequence(a1,c1,d1,some_copy_str_uindex)


    return muta_seq

def assign_to_unpaired():

    return

def CGmonitor():

    return

def pair_replace(initial_seq,save_paired_pos,some_paired):# this function used to replace basepairs

    to_modify = initial_seq
    a1=calculate_a(save_paired_pos) #connect the paired position into one sequence
    b1=calculate_b(a1)              #sorted the index of the position descend order
    c1=calculate_c(a1)              #calculate the original index of the paired positon
    d1=calculate_d(some_paired)     #connect the base paires into one sequence

    replacements=[]
    indexes=b1


    for i in range(len(c1)):
        replacements.append(d1[c1[i]])
        to_modify[indexes[i]] = replacements[i]

    return to_modify


def GC_pairreplace(ini_seq, dif_ini_GC,some_paired):
    GC_to_modify = initial_seq
    a1=calculate_a(dif_ini_GC) #connect the paired position into one sequence
    b1=calculate_b(a1)              #sorted the index of the position descend order
    c1=calculate_c(a1)              #calculate the original index of the paired positon
    d1=calculate_d(some_paired)


    return

def check_seq_base(some_paired_pos, predicted_seq,posl):
    check_even=[]
    check_odd=[]
    A_change=["G","C"]
    C_change=["A","U"]
    U_change=["U","C"]
    G_change=["A","G"]
    a1=calculate_a(some_paired_pos)
    even=a1[::2]
    odd=a1[1::2]
    new_even=[]
    new_odd=[]

    for i in range(len(even)):
        if even[i] not in posl:
            new_even.append(even[i])
            new_odd.append(odd[i])
        if odd[i] not in posl:
            new_even.append(even[i])
            new_odd.append(odd[i])

    for i in range(len(new_odd)):
        check_even.append(predicted_seq[new_even[i]])
        check_odd.append(predicted_seq[new_odd[i]])

    for i in range(len(new_odd)):
        if predicted_seq[new_odd[i]]=="A":
            predicted_seq[new_even[i]]=random.choice(A_change)

        elif predicted_seq[new_odd[i]]=="U":
            predicted_seq[new_even[i]]=random.choice(U_change)

        elif predicted_seq[new_odd[i]]=="C":
            predicted_seq[new_even[i]]=random.choice(C_change)

        elif predicted_seq[new_odd[i]]=="G":
            predicted_seq[new_even[i]]=random.choice(G_change)


    return predicted_seq

def check_GC_base(some_dif_ini,predicted_seq,posl):## assign GC or CG to predicted sequence
    GC=["G","C"]
    new_dif_ini=[]
    for i in range(len(some_dif_ini)):
        if some_dif_ini[i] not in posl:
            new_dif_ini.append(some_dif_ini[i])
    a1=calculate_a(new_dif_ini)




    even=a1[::2]
    odd=a1[1::2]
    for i in range(len(odd)):
        predicted_seq[odd[i]]=random.choice(GC)
        if predicted_seq[odd[i]]=="G":
            predicted_seq[even[i]]="C"
        if predicted_seq[odd[i]]=="C":
            predicted_seq[even[i]]="G"



    return predicted_seq


def update(paired_pos,dif_ini,posl,predicted_seq):
    a1=calculate_a(paired_pos)
    a2=calculate_a(str_index)
    #a3=calculate_a(posl)
    even=a1[::2]
    odd=a1[1::2]
    new_paired_pos=[]
    updated_weaken_pairs_position=[]
    check_odd=[]
    check_even=[]
    #new_even=[]
    #new_odd=[]
    for i in range(len(even)):
        if even[i] and odd[i] not in posl:
            new_paired_pos.append(paired_pos[i])
    #print new_paired_pos

    a_new=calculate_a(new_paired_pos)
    #print a_new
    new_even=a_new[::2]
    new_odd=a_new[1::2]
    #print new_odd
    #print new_even
    for i in range(len(new_even)):
        if new_even[i] and new_odd[i] not in a2:
            updated_weaken_pairs_position.append(new_paired_pos[i])

    a_final=calculate_a(updated_weaken_pairs_position)
    final_even=a_final[::2]
    final_odd=a_final[1::2]
    A_change=["G","C"]
    C_change=["A","U"]
    U_change=["U","C"]
    G_change=["A","G"]

    for i in range(len(final_even)):
        #print predicted_seq
        if predicted_seq[final_even[i]]=="A":
            #print predicted_seq
            predicted_seq[final_odd[i]]=random.choice(A_change)
        if predicted_seq[final_even[i]]=="U":
            predicted_seq[final_odd[i]]=random.choice(U_change)
        if predicted_seq[final_even[i]]=="C":
            predicted_seq[final_odd[i]]=random.choice(C_change)
        if predicted_seq[final_even[i]]=="G":
            predicted_seq[final_odd[i]]=random.choice(G_change)


    GC=["G","C"]
    AU=["A","U"]
    new_dif_ini=[]
    for i in range(len(dif_ini)):
        if dif_ini[i] not in posl:
            new_dif_ini.append(dif_ini[i])
    a_ini=calculate_a(new_dif_ini)
    #for i in range(len(some_dif_ini)):




    even_ini=a_ini[::2]
    odd_ini=a_ini[1::2]

    for i in range(len(even_ini)):
        if predicted_seq[odd_ini[i]]=="G":
            predicted_seq[even_ini[i]]="C"
        if predicted_seq[odd_ini[i]]=="C":
            predicted_seq[even_ini[i]]="G"

        if predicted_seq[odd_ini[i]]=="A":
            predicted_seq[even_ini[i]]="U"
        if predicted_seq[odd_ini[i]]=="U":
            predicted_seq[even_ini[i]]="A"






    return predicted_seq









def check_GC_base3(some_dif_ini,predicted_seq,posl,defined_GC):## assign GC or CG to predicted sequence
    GC=["G","C"]
    AU=["A","U"]
    newgc=measureGC(predicted_seq)
    new_dif_ini=[]
    for i in range(len(some_dif_ini)):
        if some_dif_ini[i] not in posl:
            new_dif_ini.append(some_dif_ini[i])
    a1=calculate_a(new_dif_ini)




    even=a1[::2]
    odd=a1[1::2]

    for i in range(len(odd)):
        if defined_GC>=newgc:
            predicted_seq[odd[i]]=random.choice(GC)
            if predicted_seq[odd[i]]=="G":
                predicted_seq[even[i]]="C"
            if predicted_seq[odd[i]]=="C":
                predicted_seq[even[i]]="G"
        else:
            predicted_seq[odd[i]]=random.choice(AU)
            if predicted_seq[odd[i]]=="A":
                predicted_seq[even[i]]="U"
            if predicted_seq[odd[i]]=="U":
                predicted_seq[even[i]]="A"




    return predicted_seq


def obtain_initial_sequence(input_structure_s):##obtain some good initial sequence over 0.8
    ini_seq_pool=[]
    ini_str_pool=[]

    some_str_index,some_str_uindex=calculate_sequence_position(input_structure_s)
    some_midea=getbasepairs(some_str_index)#### this is global varable
    some_copy_str_uindex=getunbases(some_str_uindex)# unpaired bases ## this is global varable
    some_a=calculate_a(some_str_index)
    some_b=calculate_b(some_a)
    some_c=calculate_c(some_a)
    some_d=calculate_d(some_midea)
    some_ini_seq,some_ini_str_seq=getwholesequence(some_b,some_c ,some_d , some_copy_str_uindex)
    some_str_mfe,some_str_value=calculate_mfe_and_str(some_ini_seq)
    some_str_distance=calculate_structure_distance(input_structure_s,len(input_structure_s),some_str_value)
    ini_seq_pool.append(some_ini_seq)
    ini_str_pool.append(some_str_distance)
    print some_str_value



    for i in range(10):


        paired_pos,dif_ini=dif_str(some_str_value)
        mutated_seq=check_seq_base(paired_pos,some_ini_str_seq)
        mutated_seq1=check_GC_base(dif_ini,mutated_seq)
        mutated_seq2=''.join(map(str, mutated_seq1))

        kkk=RNA.fold(mutated_seq2)[0]
        some_str_value=kkk
        some_ini_seq=mutated_seq2
        new_str_distance=calculate_structure_distance(s,len(s),kkk)
        ini_seq_pool.append(mutated_seq2)
        ini_str_pool.append(new_str_distance)


    max_idx = np.argmax(ini_str_pool)
    max_val = ini_str_pool[max_idx]


    seq=ini_seq_pool[max_idx]


    return seq,max_val,mutated_seq1


def GCcontent(defined_GC,getGC,some_paired_pos, predicted_seq,posl):


    check_even=[]
    check_odd=[]
    #print predicted_seq
    A_change=["G","C"]
    C_change=["A","U"]
    U_change=["U","C"]
    G_change=["A","G"]
    a1=calculate_a(some_paired_pos)
    even=a1[::2]
    odd=a1[1::2]
    new_even=[]
    new_odd=[]

    for i in range(len(even)):
        if even[i] not in posl:
            new_even.append(even[i])
            new_odd.append(odd[i])
        if odd[i] not in posl:
            new_even.append(even[i])
            new_odd.append(odd[i])

    for i in range(len(new_odd)):
        check_even.append(predicted_seq[new_even[i]])
        check_odd.append(predicted_seq[new_odd[i]])

    for i in range(len(new_odd)):
        if getGC<defined_GC:
            if predicted_seq[new_odd[i]]=="A":
                predicted_seq[new_even[i]]=random.choice(A_change)
            elif predicted_seq[new_odd[i]]=="U":
                predicted_seq[new_even[i]]="C"
            elif predicted_seq[new_odd[i]]=="C":
                predicted_seq[new_even[i]]="C"
            elif predicted_seq[new_odd[i]]=="G":
                predicted_seq[new_even[i]]="G"
        else:
            if predicted_seq[new_odd[i]]=="A":
                predicted_seq[new_even[i]]="A"
            elif predicted_seq[new_odd[i]]=="U":
                predicted_seq[new_even[i]]="U"
            elif predicted_seq[new_odd[i]]=="C":
                predicted_seq[new_even[i]]=random.choice(C_change)
            elif predicted_seq[new_odd[i]]=="G":
                predicted_seq[new_even[i]]="A"

    return predicted_seq


def GCcontent1(defined_GC,getGC,some_paired_pos, predicted_seq,posl):


    check_even=[]
    check_odd=[]
    #print predicted_seq
    A_change=["G","C"]
    C_change=["A","U"]
    U_change=["U","C"]
    G_change=["A","G"]
    a1=calculate_a(some_paired_pos)
    even=a1[::2]
    odd=a1[1::2]
    new_even=[]
    new_odd=[]

    for i in range(len(even)):
        if even[i] not in posl:
            new_even.append(even[i])
            new_odd.append(odd[i])
        if odd[i] not in posl:
            new_even.append(even[i])
            new_odd.append(odd[i])

    for i in range(len(new_odd)):
        check_even.append(predicted_seq[new_even[i]])
        check_odd.append(predicted_seq[new_odd[i]])

    for i in range(len(new_odd)):

        if predicted_seq[new_odd[i]]=="A":
            predicted_seq[new_even[i]]="A"
        if predicted_seq[new_odd[i]]=="U":
            if predicted_seq[new_even[i]]=="G":
                predicted_seq[new_even[i]]="C"
            else:
                predicted_seq[new_even[i]]="U"
        if predicted_seq[new_odd[i]]=="C":
            predicted_seq[new_even[i]]="C"
        if predicted_seq[new_odd[i]]=="G":
            if predicted_seq[new_even[i]]=="U":
                predicted_seq[new_even[i]]="A"
            else:
                predicted_seq[new_even[i]]="G"

    return predicted_seq

def measureGC(generate_seq):
    n=len(generate_seq)
    cont=0.0
    indexGC=[]
    indexnotGC=[]
    getGC=0.0
    for i in range(len(generate_seq)):
        if generate_seq[i]=="C":
            cont=cont+1
            indexGC.append(i)
        if generate_seq[i]=="G":
            cont=cont+1
            indexGC.append(i)

        else:
            indexnotGC.append(i)

    getGC=cont/n
    return getGC

def  calculate_GC_numbers(eposl,defined_GC,need,poslalpha):
    cunnt=0.0
    for i in range(len(eposl)):
        if eposl[i]=='G' or 'C':
            cunnt=cunnt+1

    needGC=len(s)*defined_GC-cunnt
    if needGC>0:
        real=round(needGC/2,0)
    else:
        real=0
    return real



def error_check(defined_GC1,defined_gd1,s1,d1):
    if defined_GC1>1.0 or defined_GC1<0.0:
        print "Error,please input a right range in [0, 1.0]"

def calculate__pseudo_sequence_position(seq):
    stack = []
    struc = []
    pseu=[]
    pseu1=[]
    ustruc=[]
    for i in xrange(len(seq)):
        if seq[i] == '(':
            stack.append(i)
        if seq[i] == ')':
            struc.append((stack.pop(), i))
        if seq[i]=='.':
            ustruc.append(i)
        if seq[i]=='[':
            pseu.append(i)
        if seq[i]==']':
            pseu1.append((pseu.pop(),i))

    return struc,ustruc,pseu1

def calculate__pseudo_mfe_and_str(sequence):
    rnafold= pseudoknot(sequence)
    mfe=rnafold[1]
    str_v=rnafold[0]
    return mfe,str_v

def pseudoknot(se):

    cmd = ["RNAPKplex","-e","-8.10"]
    #tmpdir = mkdtemp()

    p = Popen(cmd, stdin = PIPE, stdout = PIPE)
    print >> p.stdin, se
    p.stdin.close()
    t = p.stdout.readlines()[-1].strip().split(None, 1)
    p.stdout.close()
    return t


def calculate__pseudo_sequence_position_pKiss(seq):
    stack = []
    struc = []
    pseu=[]
    pseu1=[]
    ustruc=[]
    pseularg=[]
    pseularg1=[]
    pseupk=[]
    pseupk1=[]
    for i in xrange(len(seq)):
        if seq[i] == '(':
            stack.append(i)
        if seq[i] == ')':
            struc.append((stack.pop(), i))
        if seq[i]=='.':
            ustruc.append(i)
        if seq[i]=='[':
            pseu.append(i)
        if seq[i]==']':
            pseu1.append((pseu.pop(),i))
        if seq[i]=='{':
            pseularg.append(i)
        if seq[i]=='}':
            pseularg1.append((pseularg.pop(),i))
        if seq[i]=='<':
            pseupk.append(i)
        if seq[i]=='>':
            pseupk1.append((pseupk.pop(),i))

    return struc+pseu1+pseularg1+pseupk1,ustruc

def calculate__pseudo_mfe_and_str_RNApKplex(sequence):
    rnafold= pseudoknot_RNApKplex(sequence)
    mfe=rnafold[1]
    str_v=rnafold[0]
    return mfe,str_v

def calculate__pseudo_mfe_and_str_pkiss(sequence):
    mfe,str_v= pseudoknot_pkiss(sequence)
    return mfe,str_v

def checkpKiss():

  	#pKiss_output = subprocess.Popen(["which", "pKiss_mfe"], stdout=subprocess.PIPE).communicate()[0].strip()
  	pKiss_output = subprocess.Popen(["which", "pKiss_mfe"], stdout=subprocess.PIPE, shell=True).communicate()[0].strip()
	if len(pKiss_output) > 0 and pKiss_output.find("found") == -1 and pKiss_output.find(" no ") == -1:
		return True
	else:
		print "Please install pKiss"
		print "Download from http://bibiserv2.cebitec.uni-bielefeld.de/pkiss"
		exit(0)


def pseudoknot_RNApKplex(se):

    cmd = ["RNAPKplex","-e","-8.10"]
    p = Popen(cmd, stdin = PIPE, stdout = PIPE)
    print >> p.stdin, se
    p.stdin.close()
    t = p.stdout.readlines()[-1].strip().split(None, 1)
    p.stdout.close()
    return t



def checkRNAfold():
	
	RNAfold_output = subprocess.Popen(["which", "RNAfold"], stdout=subprocess.PIPE).communicate()[0].strip()
	if len(RNAfold_output) > 0 and RNAfold_output.find("found") == -1 and RNAfold_output.find(" no ") == -1:
		return True
	else:
		print "Please install RNAfold"
		print "Download from http://www.tbi.univie.ac.at/"
		exit(0)

def pseudoknot_pkiss(se):

    cmd = ["pKiss_mfe",se]
    p = Popen(cmd, stdin = PIPE, stdout = PIPE)
    t = p.stdout.read().split("\n");
    if (len(t) > 1):
        mfe = "".join(t[1].split(" ")[1])
        structure= "".join(t[1].split(" ")[3])
    p.stdout.close()
    return mfe, structure







if __name__ == "__main__":

    BASEPAIRS = ["AU", "CG", "GC", "UA"]
    bases="AGCU"
    parser = argparse.ArgumentParser()
    #parser.add_argument('-f', dest='action', action='store_const',const=None,help="monte carlo tree search for RNA inverse folding")
    parser.add_argument('-s',help="input the dot-branket representation of the RNA secondary structure")
    parser.add_argument('-GC',default=0.5,help="input the target GC content,the default GC-content is 0.5")
    parser.add_argument('-d',default=0.01,help="input the GC content error range [0, 0.02],the default GC-content error is 0.01 for nested structures, 0.02 for pseudoknot structures")
    parser.add_argument('-pk',default=0, help="this is for handling pseduoknot structures, you can use different pseduoknot prediction software, value=1 means choose pKiss, value=0 means choose RNAfold")

    parsed_args = parser.parse_args()
    s=getattr(parsed_args, 's')
    defined_GC=float(getattr(parsed_args, 'GC'))
    defined_gd=float(getattr(parsed_args,'d'))
    defined_pseudo=float(getattr(parsed_args,'pk'))
    if defined_pseudo==1:
        checkpKiss()
        str_index1,str_uindex1=calculate__pseudo_sequence_position_pKiss(s)
        str_index=str_index1
        str_uindex=str_uindex1
    else:
        checkRNAfold()
        str_index,str_uindex=calculate_sequence_position(s)
    midea=getbasepairs(str_index)#### this is global varable
    copy_str_uindex=getunbases(str_uindex)# unpaired bases ## this is global varable
    a= calculate_a(str_index)
    b= calculate_b(a)
    c=calculate_c(a)
    d=calculate_d(midea)
    ini_seq,ini_str_seq=getwholesequence(b,c ,d , copy_str_uindex)
    if defined_GC<=1.0 and defined_GC>=0.0:
        best_str,GC,run_time=UCTRNA()
        if best_str==1.0:
            print "running time:"+str(run_time)
            print "GC content:"+str(GC)
            print "GC distance:"+str(abs(GC-defined_GC))
            print "structure distance:" +str(best_str)

        else:
            print "running time:"+str(run_time)
            print "GC content:"+str(GC)
            print "GC distance:"+str(abs(GC-defined_GC))
            print "structure distance:" +str(best_str)

    else:
        UCTRNAnoGC()
