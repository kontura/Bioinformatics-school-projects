#!/usr/bin/env python2.7
import code

import sys
import math

class Select:
    start, end, score, source, phase, participating, previous_result = range(7)

class GRangeList:
    def __init__(self):
        self.list = []
        self.sorted_indexes_by_end_time = []
        self.sorted_indexes_by_start_time = []
        self.previous_non_overlap = {}

    def amount_of_bp(self):
        amount = 0
        for rang in self.list:
            amount = amount + (rang[Select.end] - rang[Select.start])
        return amount

    def print_result(self, key):
        values = key.split('/')  #list_name = values[0]+"/"+values[2]+"/"+values[6]
        current = self.list[-1]
        while(current[Select.previous_result] != -1):
            if current[Select.participating]:
                sys.stdout.write(values[0] + "\t" + current[Select.source] + "\t" + values[1] + "\t" + str(current[Select.start]) + "\t" + str(current[Select.end]) + "\t" + str(current[Select.score]) + "\t" + values[2] + "\t" + current[Select.phase])
            current = self.list[current[Select.previous_result]]
        if current[Select.participating]:
            sys.stdout.write(values[0] + "\t" + current[Select.source] + "\t" + values[1] + "\t" + str(current[Select.start]) + "\t" + str(current[Select.end]) + "\t" + str(current[Select.score]) + "\t" + values[2] + "\t" + current[Select.phase])
        sys.stdout.flush()

    def count_suma(self, i, list_of_sumas):
        if i < 0:
            return 0
        else:
            return list_of_sumas[i]

    def get_result(self, i, list_of_results):
        if i < 0:
            return []
        else:
            return list_of_results[i]

    def resolve_overlaps(self):
        #results = {}
        #results[0] = [self.sorted_indexes_by_end_time[0]]
        suma = {}
        suma[0] = self.list[self.sorted_indexes_by_end_time[0]][Select.score]
        last_i = -1
        first_non_deleted = 0

        for i in self.sorted_indexes_by_end_time:
            #print(str(i) + " :i")
            prev = self.previous_non_overlap[i]
            if (self.count_suma(last_i, suma) > (self.list[i][Select.score] + self.count_suma(prev, suma))):
                #results[i] = self.get_result(last_i, results)
                suma[i] = self.count_suma(last_i, suma)
                self.list[i][Select.participating] = None
                self.list[i][Select.previous_result] = last_i
            else:
                #code.interact(local=locals())
                #results[i] = [i] + self.get_result(prev, results)
                suma[i] = self.list[i][Select.score] + self.count_suma(prev, suma)
                self.list[i][Select.participating] = True
                self.list[i][Select.previous_result] = prev
            last_i = i

        return suma[self.sorted_indexes_by_end_time[-1]]

    def bin_search_insert_start(self, value):
        index = len(self.list)-1
        left = 0
        right = index-1
        middle = index;
        while left <= right:
            #code.interact(local=locals())
            middle = int(math.floor((left+right)/2))
            if self.list[self.sorted_indexes_by_start_time[middle]][Select.start] == value:
                break
            else:
                if self.list[self.sorted_indexes_by_start_time[middle]][Select.start] < value:
                    middle = middle + 1
                    left = middle
                else:
                    if middle == left: break
                    middle = middle - 1
                    right = middle

         
        if (middle < 0): middle = 0
        self.sorted_indexes_by_start_time.insert(middle, index)
        return middle

    def bin_search_insert_end(self, value):
        index = len(self.list)-1
        left = 0
        right = index-1
        middle = index;
        while left <= right:
            middle = int(math.floor((left+right)/2))
            if self.list[self.sorted_indexes_by_end_time[middle]][Select.end] == value:
                break
            else:
                if self.list[self.sorted_indexes_by_end_time[middle]][Select.end] < value:
                    middle = middle + 1
                    left = middle
                else:
                    if middle == left: break
                    middle = middle - 1
                    right = middle

         
        if (middle < 0): middle = 0
        self.sorted_indexes_by_end_time.insert(middle, index)
        return middle

    def find_all_previous_non_overlaps(self):
        i = 0
        closest_end = self.sorted_indexes_by_end_time[i]

        if (s_switch):
            start_range = range(0, len(self.list))
        else:
            start_range = self.sorted_indexes_by_start_time

        for start in start_range:
           # print("find non overlap for: " + str(start))
           # print(str(self.list[closest_end][Select.end]) + " compating " + str(self.list[start][Select.start]))
           # print(str(i) + " :is i")
            while(self.list[closest_end][Select.end] < self.list[start][Select.start]):
                i = i+1
                closest_end = self.sorted_indexes_by_end_time[i]
            if (self.sorted_indexes_by_end_time[i-1] != start):
                if (i-1 < 0):
                    self.previous_non_overlap[start] = -1
                else:
                    self.previous_non_overlap[start] = self.sorted_indexes_by_end_time[(i-1)]
            else:
                if (i-2 < 0):
                    self.previous_non_overlap[start] = -1
                else:
                    self.previous_non_overlap[start] = self.sorted_indexes_by_end_time[(i-2)]
        self.sorted_indexes_by_start_time = [] #to free memory

    def add_range(self, gr):
        self.list.append(gr)
        if (not s_switch):
            start_index = self.bin_search_insert_start(gr[Select.start])
        end_index = self.bin_search_insert_end(gr[Select.end])

    def count(self):
        return len(self.list)

    def print_starts(self):
        print("starts:")
        print(self.sorted_indexes_by_start_time)

    def print_ends(self):
        print("ends:")
        print(self.sorted_indexes_by_end_time)
       # for index in self.sorted_indexes_by_end_time:
       #     print(str(index) + ": " + str(self.list[index][Select.end]))

    def print_previous_overlaps(self):
        print("overlaps:")
        for overlap in self.previous_non_overlap:
            print(str(overlap) + " :is not first overlaped by: " + str(self.previous_non_overlap[overlap]));

    def print_list(self):
        i = 0
        for item in self.list:
            print(str(i)+". : " + str(item[Select.start]) + ":" + str(item[Select.end]) + " = " + str(item[Select.previous_result]) + " , participating: " + str(item[Select.participating]))
            i += 1

def argument():
    amount_of_arguments = len(sys.argv)
    if (amount_of_arguments != 2) and (amount_of_arguments != 3):
        print("Usage: resolve_overlaps.py <file.gff> [-s]")
        quit()
    else:
        if (amount_of_arguments == 2):
            return sys.argv[1], None
        else:
            return sys.argv[1], True


def parse_line(line, lists, sources_list):
   values = line.split('\t') 
   if (len(values) == 8):
       list_name = values[0]+"/"+values[2]+"/"+values[6]

       #prefering list over object
       #list of objects takes too much memory, compared to list of lists
       interval = []
       interval.append(int(values[3]))   #start
       interval.append(int(values[4]))   #end
       interval.append(float(values[5])) #score
       interval.append(values[1])        #sources
       interval.append(values[7])        #phase
       interval.append(None)             #participating(if is in the end result)
       interval.append(int(-1))          #index to next in result

       if list_name in lists:
           lists[list_name].add_range(interval)
       else:
           lists[list_name] = GRangeList()
           lists[list_name].add_range(interval)

def get_source_index(source_str, list_of_sources):
    for i in range(0,len(list_of_sources)):
        if (list_of_sources[i] == source_str):
            return i
    list_of_sources.append(source_str)
    return len(list_of_sources)

def parse_by_line(file_name, dict_of_GRangeLists, list_of_sources):
    with open(file_name) as file_gff:
        for line in file_gff:
            parse_line(line, dict_of_GRangeLists, list_of_sources)

#main
file_name, s_switch = argument()
sources_list = []
dict_of_GRangeLists = {}

#
# Pokud predpokladame ze vstup je serazeny podle pozice zacatku intervalu,
# muzeme si usetrit cely jeden seznam(sorted_indexes_by_start_time) a jedno
# binarni vyhledavani(razeni) pri jeho vytvareni. 
# Jeho roli prevezme hlavni list itervalu(reprezentovanych listy, tedy list listu),
# ktery uz bude automaticky takto serazen.
#
# Tento predpoklad na mem notebooku zrychlil vypocet z 4m45s na 2m44s (orientacne)
# Uspory pameti jsem s -s prepinacem dosahl jen zanedbatelne
#

parse_by_line(file_name, dict_of_GRangeLists, sources_list)

for key in dict_of_GRangeLists:
    dict_of_GRangeLists[key].find_all_previous_non_overlaps()
    val = dict_of_GRangeLists[key].resolve_overlaps()
    dict_of_GRangeLists[key].print_result(key)

#    sys.stdout.write(str(val))
#    sys.stdout.write("\n")

    sys.stdout.flush()

